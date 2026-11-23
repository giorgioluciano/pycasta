import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--poc", type=Path, required=True)
parser.add_argument("--poc-info", type=Path, required=True)
args = parser.parse_args()

run = Path(
    "output_isolated_tests/castp_reference_1tcd/"
    "cast_pdb_1tcd_20261118_214940_651772"
)
pdb_path = Path("data/cast_historical/1tcd.pdb")
report = json.loads(
    (run / "alpha_results/cast_alpha_20261118_214942_316566/results.json")
    .read_bytes()
)
raw = (run / "cast_complex_all_dims.json").read_bytes()
prov = report["provenance"]
if hashlib.sha256(raw).hexdigest() not in prov["input_sha256"].values():
    raise SystemExit("STOP: il complesso non corrisponde al report.")

pdb = pdb_path.read_bytes()
if hashlib.sha256(pdb).hexdigest() != prov["manifest"]["input_sha256"]:
    raise SystemExit("STOP: il PDB non corrisponde al manifest.")

atoms = [l for l in pdb.decode("ascii").splitlines() if l.startswith("ATOM  ")]
coords = [[float(l[i:i + 8]) for i in (30, 38, 46)] for l in atoms]
atom_index = {l[:54]: i for i, l in enumerate(atoms)}
if (
    len(atoms) != 3818
    or len(atom_index) != len(atoms)
    or coords != prov["manifest"]["coordinates"]
):
    raise SystemExit("STOP: identita o ordine degli atomi non verificati.")

pockets = defaultdict(set)
for line in args.poc.read_text().splitlines():
    if line.startswith(("ATOM  ", "HETATM")):
        if line[:54] not in atom_index:
            raise SystemExit(f"STOP: atomo CASTp non trovato: {line[:54]}")
        pockets[int(line.split()[-2])].add(atom_index[line[:54]])

mouths = {}
for line in args.poc_info.read_text().splitlines()[1:]:
    fields = line.split()
    if fields and fields[0] == "POC:":
        mouths[int(fields[2])] = int(fields[3])
if len(mouths) != 78 or set(pockets) != set(mouths):
    raise SystemExit("STOP: gli ID di .poc e .pocInfo non coincidono con 1TCD.")

g = json.loads(raw)
outside = set(g["outside_ancestors"])
region_of = {
    i: r["region_id"]
    for r in report["regions"]
    for i in r["tetra_indices"]
}
cells_by_atom = defaultdict(list)
labels = []
for i, (cell, in_alpha) in enumerate(zip(g["cells_idx"], g["tetra_in_alpha"])):
    for atom in cell:
        cells_by_atom[atom].append(i)
    if in_alpha:
        label = "alpha"
    elif i in outside:
        label = "flow_exterior"
    elif i in region_of:
        label = "region_0" if region_of[i] == 0 else "other_region"
    else:
        label = "unclassified"
    labels.append(label)

categories = ("alpha", "flow_exterior", "region_0", "other_region", "unclassified")
out = Path("1tcd_audit_locale.csv")
with out.open("x", encoding="utf-8", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(
        ["castp_id", "castp_n_mouths", "castp_n_atoms"]
        + [f"ge{k}_{cat}" for k in (1, 2, 3) for cat in categories]
    )
    for pid, vertices in sorted(pockets.items()):
        hits = Counter(i for atom in vertices for i in cells_by_atom[atom])
        counts = {
            k: Counter(labels[i] for i, n in hits.items() if n >= k)
            for k in (1, 2, 3)
        }
        writer.writerow(
            [pid, mouths[pid], len(vertices)]
            + [counts[k][cat] for k in (1, 2, 3) for cat in categories]
        )

print(f"{len(pockets)} regioni CASTp; CSV: {out.resolve()}")