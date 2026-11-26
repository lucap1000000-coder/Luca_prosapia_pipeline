"""
Set up a new project directory for a linker design campaign.

Reads the input PDB, extracts protein info (subunit length, chain count,
residue range), and scaffolds the project directory with:
  - protein_info.json
  - .env  
  - activation/ scripts
  - specs/ with suggested rfd3 yamls
  - tools/rfdiffusion3/ sbatch

Usage:
    sapia setup_project --pdb SlyB_trim.pdb --run-dir /path/to/project
    sapia setup_project --pdb SlyB_trim.pdb  # uses current directory
"""

import argparse
import json
import shutil
from pathlib import Path

import gemmi


ACTIVATION_TEMPLATES = {
    "rfdiffusion3.sh": """# Activation script for RFdiffusion3
module load Foundry/0.2.0-foss-2025a-CUDA-12.8.0
export FOUNDRY_CHECKPOINT_DIRS="$VSC_DATA_VO/foundry/foundry_checkpoints"
""",
    "python.sh": """# Activation script for CPU Python tools
source /data/brussel/118/vsc11875/prosapia-env/bin/activate
""",
    "boltz.sh": """# Activation script for Boltz on VUB HPC
module load Boltz-2/2.2.1-foss-2024a-CUDA-12.6.0
export BOLTZ_CACHE="/data/brussel/118/vsc11875/.boltz"
""",
    "proteinmpnn.sh": """# Activation script for ProteinMPNN on VUB HPC
module load ProteinMPNN/1.0.1-20230627-foss-2025a-CUDA-12.8.0
export PROTEIN_MPNN="/apps/brussel/RL9/zen2-ib/software/ProteinMPNN/1.0.1-20230627-foss-2025a-CUDA-12.8.0"
""",
    "usalign.sh": """# Activation script for USalign
export USALIGN_BIN="/data/brussel/118/vsc11875/bin/USalign"
""",
}

RFD3_SBATCH_SRC = Path("/data/brussel/118/vsc11875/prosapia/src/prosapia/tools/rfdiffusion3/rfdiffusion3.sbatch")

SPEC_TEMPLATE_NOHOT = """contig: "{chain_a}{start}-{end},{linker_min}-{linker_max},{chain_b}{start}-{end}"
is_non_loopy: true
"""

RFD4_SPEC_TEMPLATE = """slyb_linker:
  input: {input_cif}
  contig: "{chain_a}{start}-{end}, LINK:{linker_min}-{linker_max}, /0, {chain_b}{start}-{end}"
  conditions:
    C_CRD: {{true: [{{select: "sel(\'{chain_a}\') and is_protein_backbone"}}, {{select: "sel(\'{chain_b}\') and is_protein_backbone"}}]}}
    C_SEQ: {{true: [{{select: "sel(\'{chain_a}\')"}} , {{select: "sel(\'{chain_b}\')"}}]}}
"""

SPEC_TEMPLATE_HOT = """contig: "{chain_a}{start}-{end},{linker_min}-{linker_max},{chain_b}{start}-{end}"
is_non_loopy: true
select_hotspots:
{hotspots}
infer_ori_strategy: hotspots
"""


def get_protein_info(pdb_path: Path, tm_start: int | None = None, tm_end: int | None = None) -> dict:
    st = gemmi.read_structure(str(pdb_path))
    model = st[0]
    chains = list(model)
    
    chain_a = chains[0]
    residues = [r for r in chain_a if r.find_atom("CA", "*")]
    
    native_sequence = gemmi.one_letter_code([r.name for r in residues])

    return {
        "subunit_len": len(residues),
        "n_chains": len(chains),
        "contig_start": residues[0].seqid.num,
        "contig_end": residues[-1].seqid.num,
        "chain_a": chain_a.name,
        "chain_b": chains[1].name if len(chains) > 1 else "B",
        "native_sequence": native_sequence,
        "tm_start": tm_start,
        "tm_end": tm_end,
        "input_pdb": str(pdb_path.resolve()),
    }



def build_c12_templates(full_c12_path: Path, info: dict, run_dir: Path) -> dict:
    """Build fused dimer, context, and trimer CIFs from a full C12 structure."""
    import gemmi

    st = gemmi.read_structure(str(full_c12_path))
    model = st[0]
    chains = [c.name for c in model]

    s = info["contig_start"]
    e = info["contig_end"]
    a = info["chain_a"]
    b = info["chain_b"]

    # Neighbors: L (before A) and C (after B) in C12 ring A-B-C-D-E-F-G-H-I-J-K-L
    idx_a = chains.index(a)
    chain_l = chains[(idx_a - 1) % len(chains)]
    chain_c = chains[(chains.index(b) + 1) % len(chains)]

    inputs_dir = run_dir / "inputs"
    inputs_dir.mkdir(exist_ok=True)

    def trim_chain(src_chain, name, offset=0):
        new_chain = gemmi.Chain(name)
        for res in model[src_chain]:
            if s <= res.seqid.num <= e:
                r = res.clone()
                r.seqid = gemmi.SeqId(res.seqid.num + offset, " ")
                r.label_seq = res.seqid.num + offset
                r.subchain = name
                new_chain.add_residue(r)
        return new_chain

    def fix_entity_ids(st):
        ent_map = {}
        for i, ent in enumerate(st.entities, start=1):
            ent.name = str(i)
            for sub in ent.subchains:
                ent_map[sub] = str(i)
        for model in st:
            for chain in model:
                for res in chain:
                    res.entity_id = ent_map.get(res.subchain, "1")

    def write_cif(chains_list, path):
        new_st = gemmi.Structure()
        new_model = gemmi.Model("1")
        for ch in chains_list:
            new_model.add_chain(ch)
        new_st.add_model(new_model)
        new_st.setup_entities()
        fix_entity_ids(new_st)
        doc = new_st.make_mmcif_document()
        block = doc.sole_block()
        # Add _entity_poly_seq (required by boltz template parser)
        entity_seqs = {}
        for chain in new_model:
            subchain = list(chain)[0].subchain if list(chain) else chain.name
            entity_id = None
            for i, ent in enumerate(new_st.entities, start=1):
                if subchain in ent.subchains:
                    entity_id = str(i)
                    break
            if entity_id and entity_id not in entity_seqs:
                entity_seqs[entity_id] = [r.name for r in chain if r.find_atom("CA", "*")]
        loop = block.init_loop("_entity_poly_seq.", ["entity_id", "num", "mon_id", "hetero"])
        for eid, seq in entity_seqs.items():
            for i, res_name in enumerate(seq, start=1):
                loop.add_row([eid, str(i), res_name, "n"])
        doc.write_file(str(path))

    subunit_len = info["subunit_len"]
    offset = subunit_len + 1  # 1-residue gap

    # 1. Fused dimer CIF (A + B with 1-residue gap, single chain A)
    chain_a_part = trim_chain(a, "A", offset=0)
    chain_b_part = trim_chain(b, "A", offset=offset)
    # Combine into one chain
    fused_chain = gemmi.Chain("A")
    for res in chain_a_part:
        fused_chain.add_residue(res.clone())
    for res in chain_b_part:
        fused_chain.add_residue(res.clone())
    dimer_fused_path = inputs_dir / "C12_dimer_fused_rfd4.cif"
    write_cif([fused_chain], dimer_fused_path)

    # 2. Context CIF (fused dimer + L + C)
    chain_l_part = trim_chain(chain_l, chain_l)
    chain_c_part = trim_chain(chain_c, chain_c)
    context_path = inputs_dir / "C12_dimer_context_rfd4.cif"
    fused_chain2 = gemmi.Chain("A")
    for res in chain_a_part:
        fused_chain2.add_residue(res.clone())
    for res in chain_b_part:
        fused_chain2.add_residue(res.clone())
    write_cif([fused_chain2, chain_l_part, chain_c_part], context_path)

    # 3. Trimer fused CIF (A + B + C concatenated, continuous numbering)
    trimer_chain = gemmi.Chain("A")
    new_resnum = 1
    for src_chain in [a, b, chains[(chains.index(b) + 1) % len(chains)]]:
        for res in model[src_chain]:
            if s <= res.seqid.num <= e:
                r = res.clone()
                r.seqid = gemmi.SeqId(new_resnum, " ")
                r.label_seq = new_resnum
                r.subchain = "A"
                trimer_chain.add_residue(r)
                new_resnum += 1
    trimer_path = inputs_dir / "C12_trimer_fused.cif"
    write_cif([trimer_chain], trimer_path)

    print(f"Written {dimer_fused_path}")
    print(f"Written {context_path}")
    print(f"Written {trimer_path}")

    return {
        "dimer_fused_cif": str(dimer_fused_path.resolve()),
        "context_cif": str(context_path.resolve()),
        "trimer_fused_cif": str(trimer_path.resolve()),
        "context_chain_l": chain_l,
        "context_chain_c": chain_c,
    }

def setup_project(pdb_path: Path, run_dir: Path, tm_start: int | None = None, tm_end: int | None = None, full_c12: Path | None = None) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)

    # 1. protein_info.json
    info = get_protein_info(pdb_path, tm_start=tm_start, tm_end=tm_end)
    info_path = run_dir / "protein_info.json"
    info_path.write_text(json.dumps(info, indent=2))
    print(f"Written {info_path}")

    # 1b. build C12 templates if full_c12 provided
    if full_c12 is not None:
        template_info = build_c12_templates(full_c12, info, run_dir)
        info.update(template_info)
        info_path.write_text(json.dumps(info, indent=2))
        print(f"Updated {info_path} with template paths")

    # 2. activation scripts
    act_dir = run_dir / "activation"
    act_dir.mkdir(exist_ok=True)
    for name, content in ACTIVATION_TEMPLATES.items():
        (act_dir / name).write_text(content)
    print(f"Written activation scripts to {act_dir}")

    # 3. .env
    scratch = str(run_dir.resolve()).replace("/vscmnt/brussel_pixiu_data/_data_brussel", "/data/brussel")
    rhea = str(run_dir.resolve())
    env_content = f"""# prosapia config for VUB HPC

SAPIA_ACTIVATE_RFDIFFUSION3="{rhea}/activation/rfdiffusion3.sh"
SAPIA_ACTIVATE_BBCHECK="{rhea}/activation/python.sh"
SAPIA_ACTIVATE_CHAINLEN="{rhea}/activation/python.sh"
SAPIA_ACTIVATE_LINKER_CLASH="{rhea}/activation/python.sh"
SAPIA_ACTIVATE_BOLTZ="{rhea}/activation/boltz.sh"
SAPIA_ACTIVATE_PROTEINMPNN="{rhea}/activation/proteinmpnn.sh"
SAPIA_ACTIVATE_USALIGN="{rhea}/activation/usalign.sh"
SAPIA_ACTIVATE_C12_CLASH="{rhea}/activation/python.sh"
PROSAPIA_TOOLS_DIR=/data/brussel/118/vsc11875/prosapia/src/prosapia/tools
"""
    (run_dir / ".env").write_text(env_content)
    print(f"Written .env")

    # 4. tools/rfdiffusion3 sbatch
    tools_rfd3 = run_dir / "tools" / "rfdiffusion3"
    tools_rfd3.mkdir(parents=True, exist_ok=True)
    if RFD3_SBATCH_SRC.exists():
        shutil.copy2(RFD3_SBATCH_SRC, tools_rfd3 / "rfdiffusion3.sbatch")
        print(f"Copied rfdiffusion3.sbatch to {tools_rfd3}")

    # 5. specs
    specs_dir = run_dir / "specs"
    specs_dir.mkdir(exist_ok=True)
    
    s = info["contig_start"]
    e = info["contig_end"]
    a = info["chain_a"]
    b = info["chain_b"]
    
    for name, lmin, lmax in [("l20_65", 20, 65), ("l65_104", 65, 104)]:
        spec = SPEC_TEMPLATE_NOHOT.format(
            chain_a=a, chain_b=b, start=s, end=e,
            linker_min=lmin, linker_max=lmax
        )
        (specs_dir / f"{name}_nohot.yaml").write_text(spec)
    
    # Write rfd4 specs (two-chain format, no hotspots)
    input_cif = str(pdb_path.with_suffix(".cif").resolve())
    for name, lmin, lmax in [("l20_65", 20, 65), ("l65_104", 65, 104)]:
        spec = RFD4_SPEC_TEMPLATE.format(
            input_cif=input_cif, chain_a=a, chain_b=b,
            start=s, end=e, linker_min=lmin, linker_max=lmax
        )
        (specs_dir / f"rfd4_{name}.yaml").write_text(spec)

    print(f"Written specs to {specs_dir}")
    print(f"\nProtein info:")
    for k, v in info.items():
        print(f"  {k}: {v}")
    print(f"\nProject ready at {run_dir}")


def build_setup_project_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--pdb", type=Path, required=True, help="Input PDB file.")
    parser.add_argument("--run-dir", type=Path, default=Path("."), help="Project directory to set up. Defaults to current directory.")
    parser.add_argument("--tm-start", type=int, default=None, help="First residue of TM helix (original numbering).")
    parser.add_argument("--tm-end", type=int, default=None, help="Last residue of TM helix (original numbering).")
    parser.add_argument("--full-c12", type=Path, default=None, help="Full C12 PDB/CIF (all 12 chains) to build rfd4 templates automatically.")
    return parser


def setup_project_from_args(args: argparse.Namespace) -> None:
    setup_project(pdb_path=args.pdb, run_dir=args.run_dir, tm_start=args.tm_start, tm_end=args.tm_end, full_c12=args.full_c12)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Set up a new linker design project directory.",
        parents=[build_setup_project_parser()],
    )
    setup_project_from_args(parser.parse_args())


if __name__ == "__main__":
    main()
