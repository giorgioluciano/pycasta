import hashlib
import json
from pathlib import Path

data = Path("data/cast_historical")
run = Path(
    "output_isolated_tests/castp_reference_1tcd/"
    "cast_pdb_1tcd_20261118_214940_651772"
)
report_path = run / "alpha_results/cast_alpha_20261118_214942_316566/results.json"
report = json.loads(report_path.read_bytes())
prov = report["provenance"]
manifest = prov["manifest"]

raw = (run / "cast_complex_all_dims.json").read_bytes()
digest = hashlib.sha256(raw).hexdigest()
if digest not in prov["input_sha256"].values():
    raise SystemExit("STOP: il complesso non corrisponde al report.")

pdb = (data / "1tcd.pdb").read_bytes()
if hashlib.sha256(pdb).hexdigest() != manifest["input_sha256"]:
    raise SystemExit("STOP: il PDB non corrisponde al manifest.")

atoms = [line for line in pdb.decode("ascii").splitlines()
         if line.startswith("ATOM  ")]
coords = [[float(line[k:k + 8]) for k in (30, 38, 46)] for line in atoms]
index = {line[:54]: i for i, line in enumerate(atoms)}
if coords != manifest["coordinates"] or len(index) != len(atoms):
    raise SystemExit("STOP: identita o ordine degli atomi non verificati.")

ref_ids = sorted({
    index[line[:54]]
    for line in (data / "1tcd.poc").read_text().splitlines()
    if line.startswith(("ATOM  ", "HETATM")) and line.split()[-2] == "6"
})
ref_set = set(ref_ids)

g = json.loads(raw)
cells = g["cells_idx"]
region_0 = set(next(
    r["tetra_indices"] for r in report["regions"] if r["region_id"] == 0
))
local = {
    i for i, cell in enumerate(cells)
    if sum(atom in ref_set for atom in cell) >= 3
}
if len(local) != 107 or not local <= region_0:
    raise SystemExit("STOP: l'estratto non riproduce il risultato dell'audit.")

links = []
for i in sorted(local):
    for k, j in enumerate(g["neighbors"][i]):
        if j in region_0 and j not in local and not g["face_in_alpha"][i][k]:
            face = cells[i][:k] + cells[i][k + 1:]
            links.append({"local_cell": i, "other_cell": j, "face": face})

vertex_ids = set(ref_ids)
for i in local:
    vertex_ids.update(cells[i])
for link in links:
    vertex_ids.update(cells[link["other_cell"]])

out = Path("1tcd_cavita6_spaziale.json")
with out.open("x", encoding="utf-8") as f:
    json.dump({
        "castp_id": 6,
        "complex_sha256": digest,
        "pdb_sha256": manifest["input_sha256"],
        "reference_atom_indices": ref_ids,
        "vertices": {
            i: {"xyz": coords[i], "radius_A": manifest["atomic_radii"][i]}
            for i in sorted(vertex_ids)
        },
        "local_tetrahedra": {i: cells[i] for i in sorted(local)},
        "adjacent_tetrahedra": {
            i: cells[i] for i in sorted({link["other_cell"] for link in links})
        },
        "connections_to_rest_of_region_0": links,
    }, f, separators=(",", ":"))

print(f"Atomi CASTp: {len(ref_ids)}; tetraedri locali: {len(local)}")
print(f"Facce verso il resto della regione 0: {len(links)}")
print(f"Estratto: {out.resolve()} ({out.stat().st_size} byte)")