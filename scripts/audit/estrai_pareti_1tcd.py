import hashlib
import json
from pathlib import Path

run = Path(
    "output_isolated_tests/castp_reference_1tcd/"
    "cast_pdb_1tcd_20261118_214940_651772"
)
complex_path = run / "cast_complex_all_dims.json"
results_path = run / "alpha_results/cast_alpha_20261118_214942_316566/results.json"

raw = complex_path.read_bytes()
digest = hashlib.sha256(raw).hexdigest()
report = json.loads(results_path.read_bytes())
provenance = report["provenance"]

if digest not in provenance["input_sha256"].values():
    raise SystemExit("STOP: il complesso non corrisponde al report.")
if len(report["regions"]) != 14:
    raise SystemExit("STOP: attese 14 regioni nel run 1TCD.")

geometry = json.loads(raw)
alpha_faces = set()
for cell, flags in zip(
    geometry["cells_idx"], geometry["face_in_alpha"], strict=True
):
    for k, included in enumerate(flags):
        if included:
            alpha_faces.add(tuple(sorted(cell[:k] + cell[k + 1:])))

regions = []
for region in report["regions"]:
    boundary = {tuple(sorted(f)) for f in region["closure_boundary"]["triangles"]}
    wall_faces = boundary & alpha_faces
    wall_atoms = {i for f in wall_faces for i in f}
    mouth_atoms = {i for f in boundary - alpha_faces for i in f}
    regions.append({
        "region_id": region["region_id"],
        "n_tetrahedra": len(region["tetra_indices"]),
        "n_mouths": region["n_mouths"],
        "n_alpha_boundary_faces": len(wall_faces),
        "wall_candidate_atom_indices": sorted(wall_atoms),
        "mouth_only_atom_indices": sorted(mouth_atoms - wall_atoms),
    })

out = Path("1tcd_pycasta_wall_candidates.json")
with out.open("x", encoding="utf-8") as f:
    json.dump({
        "complex_sha256": digest,
        "pdb_sha256": provenance["manifest"]["input_sha256"],
        "regions": regions,
    }, f, separators=(",", ":"))

print(f"{len(regions)} regioni; estratto: {out.resolve()} ({out.stat().st_size} byte)")