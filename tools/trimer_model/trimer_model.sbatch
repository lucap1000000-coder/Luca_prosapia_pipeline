#!/bin/bash
#SBATCH --job-name=trimer_model
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -uo pipefail
source "${SAPIA_PRELUDE:?}"

set +u
source "${SAPIA_ACTIVATE_TRIMER_MODEL:?set SAPIA_ACTIVATE_TRIMER_MODEL_BUILDER in your .env}"
set -u

name=$(echo "$SAPIA_LINE" | cut -f1)
src=$(echo "$SAPIA_LINE" | cut -f2)
subunit_len=$(echo "$SAPIA_LINE" | cut -f3)
mkdir -p "$OUT_DIR"

python - "$src" "$OUT_DIR/$name.json" "$OUT_DIR/$name.pdb" "$subunit_len" <<'PYEOF'
import json, sys
import numpy as np
import gemmi

src = sys.argv[1]
out_json = sys.argv[2]
out_pdb = sys.argv[3]
subunit_len = int(sys.argv[4])

st = gemmi.read_structure(src)
chain = list(st[0])[0]
residues = [r for r in chain if r.find_atom("CA", "*")]
total = len(residues)
linker_len = total - 2 * subunit_len

if linker_len <= 0:
    json.dump({"error": f"linker_len={linker_len} <= 0"}, open(out_json, "w"))
    sys.exit(0)

# Extract CA coords for subunit A and subunit B
def get_ca(res_list):
    coords = []
    for r in res_list:
        ca = r.find_atom("CA", "*")
        if ca:
            coords.append(gemmi.Position(ca.pos.x, ca.pos.y, ca.pos.z))
    return coords

sub_a = residues[:subunit_len]
sub_b = residues[subunit_len + linker_len:]

pos_a = get_ca(sub_a)
pos_b = get_ca(sub_b)

# Compute A->B transform
result = gemmi.superpose_positions(pos_b, pos_a)
T = result.transform

# Apply T to entire chain to get chain B (subunit_B + linker2 + subunit_C)
def apply_transform(residues, T, start_resnum):
    new_res_list = []
    for i, res in enumerate(residues):
        new_res = gemmi.Residue()
        new_res.name = res.name
        new_res.seqid = gemmi.SeqId(start_resnum + i, ' ')
        new_res.subchain = 'A'
        for atom in res:
            new_atom = atom.clone()
            p = T.apply(atom.pos)
            new_atom.pos = gemmi.Position(p.x, p.y, p.z)
            new_res.add_atom(new_atom)
        new_res_list.append(new_res)
    return new_res_list

# Build trimer: original chain + transformed(original) = A+link+B+link+C
transformed = apply_transform(residues, T, total + 1)

new_st = gemmi.Structure()
new_model = gemmi.Model('1')
new_chain = gemmi.Chain('A')

# Add original residues (A+linker+B)
for res in residues:
    new_chain.add_residue(res.clone())

# Add transformed residues (B_copy+linker2+C)
# We only want linker2+C: skip the first subunit_len (which is B_copy)
for res in transformed[subunit_len:]:
    new_chain.add_residue(res)

new_model.add_chain(new_chain)
new_st.add_model(new_model)
new_st.setup_entities()
new_st.write_pdb(out_pdb)

json.dump({
    "trimer_model_path": out_pdb,
    "linker_len": linker_len,
    "subunit_len": subunit_len,
    "superpose_rmsd": round(result.rmsd, 4),
}, open(out_json, "w"))
PYEOF
