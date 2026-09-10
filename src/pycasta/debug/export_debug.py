import json
from pathlib import Path


def export_pocket_debug_json(output_path, protein_coords, candidates):
    payload = {
        "protein_coords": protein_coords.tolist(),
        "pockets": [
            {
                "pocket_id": c.pocket_id,
                "representative_point": c.representative_point,
                "tetra_indices": list(c.tetra_indices),
                "features": dict(c.features),
            }
            for c in candidates
        ],
    }

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)