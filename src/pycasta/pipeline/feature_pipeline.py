import numpy as np

from .models import PocketCandidate
from ..features.geometry_features import compute_geometry_features
from ..features.mouth_features import compute_mouth_features_for_all_pockets


def build_pocket_candidates(ranked_pockets, tetra_positions):
    candidates = []

    for idx, pocket in enumerate(ranked_pockets, start=1):
        pocket_points = tetra_positions[pocket].reshape(-1, 3)
        representative_point = np.mean(pocket_points, axis=0).tolist()

        candidates.append(
            PocketCandidate(
                pocket_id=idx,
                tetra_indices=list(pocket),
                representative_point=representative_point,
            )
        )

    return candidates


def enrich_with_geometry(candidates, tetra_positions):
    for candidate in candidates:
        feats = compute_geometry_features(tetra_positions, candidate.tetra_indices)
        candidate.volume = feats["volume"]
        candidate.depth = feats["depth"]
        candidate.features.update(feats)
    return candidates


def enrich_with_mouth(candidates, ranked_pockets, simplices, protein_coords, ligand_coords=None):
    mouth_data = compute_mouth_features_for_all_pockets(
        ranked_pockets=ranked_pockets,
        simplices=simplices,
        protein_coords=protein_coords,
        ligand_coords=ligand_coords,
    )

    for candidate, mouth in zip(candidates, mouth_data):
        candidate.mouth_area = mouth.get("mouth_area")
        candidate.mouth_perimeter = mouth.get("mouth_circumference")
        candidate.mouth_center = mouth.get("mouth_center")
        candidate.rim_atoms = mouth.get("rim_atoms", [])
        candidate.ligand_distance_to_mouth = mouth.get("ligand_distance")

        candidate.features.update(
            {
                "mouth_area": candidate.mouth_area,
                "mouth_perimeter": candidate.mouth_perimeter,
                "ligand_distance_to_mouth": candidate.ligand_distance_to_mouth,
                "rim_atom_count": len(candidate.rim_atoms),
            }
        )

    return candidates