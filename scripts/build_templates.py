"""
Build the per-protein structure files for a design campaign:

  inputs/C12_dimer_context_rfd4.cif   RFD4 input (fused dimer + neighbours L, C)
  templates/dimer.cif + dimer.yaml    Boltz template, fused dimer, residues 1..2n
  templates/trimer.cif + trimer.yaml  Boltz template, fused trimer, residues 1..3n
  protein_info.json                   with the paths added

Usage:
  python build_templates.py --pdb trimmed.pdb --full-c12 full_c12.pdb \
      --tm-start 63 --tm-end 100 --out PROJECT_DIR
RFD4 specs are written by hand; this only builds structures and checks them.
"""
import argparse
import json
import shutil
import sys
from pathlib import Path

import gemmi

from setup_project import build_c12_templates, get_protein_info


def renumber_1_to_n(src: Path, dst: Path) -> None:
    st = gemmi.read_structure(str(src))
    for i, res in enumerate(st[0][0], start=1):
        res.seqid = gemmi.SeqId(i, " ")
        res.label_seq = i
    st.make_mmcif_document().write_file(str(dst))


def check_boltz_cif(path: Path, expected_nres: int) -> list[str]:
    """The three defects from check_templates.py, plus residue count."""
    problems = []
    st = gemmi.read_structure(str(path))
    if len(st) == 0:
        return ["zero models"]
    chains = [c.name for c in st[0]]
    asym = [c[0].subchain for c in st[0]]
    ents = [len(e.full_sequence) for e in st.entities]
    if not ents or min(ents) == 0:
        problems.append("empty full_sequence")
    if "_entity_poly_seq" not in path.read_text(errors="ignore"):
        problems.append("no _entity_poly_seq")
    if chains != asym:
        problems.append(f"label_asym {asym} != chains {chains}")
    if chains != ["A"]:
        problems.append(f"expected a single chain A, got {chains}")
    elif len(st[0][0]) != expected_nres:
        problems.append(f"{len(st[0][0])} residues, expected {expected_nres}")
    nums = [r.seqid.num for r in st[0][0]]
    if nums != list(range(1, len(nums) + 1)):
        problems.append("residue numbering is not 1..N")
    return problems


def arc_angle(path: Path) -> float:
    """Rotation between the two halves of the fused dimer (should be 360/ring_size)."""
    import numpy as np
    ch = gemmi.read_structure(str(path))[0][0]
    half = len(ch) // 2
    pos = lambda rs: [gemmi.Position(*r["CA"][0].pos.tolist()) for r in rs]
    sp = gemmi.superpose_positions(pos(list(ch)[:half]), pos(list(ch)[half:]))
    m = np.array(sp.transform.mat.tolist())
    return float(np.degrees(np.arccos((np.trace(m) - 1) / 2)))


def write_yaml(path: Path, cif: Path) -> None:
    path.write_text(
        f"templates:\n  - cif: {cif.resolve()}\n    chain_id: [A]\n    template_id: [A]\n"
    )


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pdb", type=Path, required=True)
    p.add_argument("--full-ring", "--full-c12", dest="full_c12", type=Path, required=True)
    p.add_argument("--tm-start", type=int)
    p.add_argument("--tm-end", type=int)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--ring-size", type=int, required=True,
                   help="subunits in the full ring (12 for C12, 11 for C11); aborts if the file differs")
    a = p.parse_args()

    ring_chains = [c.name for c in gemmi.read_structure(str(a.full_c12))[0]]
    print(f"full ring: {len(ring_chains)} chains ({' '.join(ring_chains)})")
    if len(ring_chains) != a.ring_size:
        sys.exit(f"ERROR: --ring-size {a.ring_size} but the file has {len(ring_chains)} chains")

    out = a.out
    out.mkdir(parents=True, exist_ok=True)
    tdir = out / "templates"
    tdir.mkdir(exist_ok=True)

    info = get_protein_info(a.pdb, tm_start=a.tm_start, tm_end=a.tm_end)
    ca, cb = info["chain_a"], info["chain_b"]
    if ca not in ring_chains or cb not in ring_chains:
        sys.exit(f"ERROR: trimmed PDB chains {ca},{cb} not found in the ring file ({' '.join(ring_chains)})")
    expected_b = ring_chains[(ring_chains.index(ca) + 1) % len(ring_chains)]
    if cb != expected_b:
        sys.exit(f"ERROR: chain after {ca} in the ring is {expected_b}, but the trimmed PDB's second chain is {cb}")
    info.update(build_c12_templates(a.full_c12, info, out))
    n = info["subunit_len"]

    renumber_1_to_n(Path(info["dimer_fused_cif"]), tdir / "dimer.cif")
    shutil.copy2(info["trimer_fused_cif"], tdir / "trimer.cif")
    write_yaml(tdir / "dimer.yaml", tdir / "dimer.cif")
    write_yaml(tdir / "trimer.yaml", tdir / "trimer.cif")

    failed = False
    for name, nres in [("dimer", 2 * n), ("trimer", 3 * n)]:
        probs = check_boltz_cif(tdir / f"{name}.cif", nres)
        print(f"{'FAIL' if probs else 'PASS'} templates/{name}.cif ({nres} residues)")
        for pr in probs:
            print("   ->", pr)
        failed = failed or bool(probs)

    ang, want = arc_angle(tdir / "dimer.cif"), 360.0 / a.ring_size
    ok = abs(ang - want) < 2.0
    print(f"{'PASS' if ok else 'FAIL'} dimer arc {ang:.3f} deg (expected {want:.3f} for a {a.ring_size}-fold ring)")
    failed = failed or not ok

    info["boltz_dimer_template"] = str((tdir / "dimer.yaml").resolve())
    info["boltz_trimer_template"] = str((tdir / "trimer.yaml").resolve())
    (out / "protein_info.json").write_text(json.dumps(info, indent=2))
    print(f"protein_info.json written to {out}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
