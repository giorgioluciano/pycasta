import numpy as np

from ..pocket_detection import compute_analytic_pocket_volume


def compute_pocket_depth_from_tetra(tetra_positions, pocket):
    if not pocket:
        return 0.0

    pocket_points = tetra_positions[pocket].reshape(-1, 3)
    centroid = np.mean(pocket_points, axis=0)
    distances = np.linalg.norm(pocket_points - centroid, axis=1)
    return float(np.max(distances))


def compute_geometry_features(tetra_positions, pocket):
    volume = compute_analytic_pocket_volume(tetra_positions, pocket)
    depth = compute_pocket_depth_from_tetra(tetra_positions, pocket)
    n_tetra = len(pocket)

    return {
        "volume": float(volume),
        "depth": float(depth),
        "n_tetra": int(n_tetra),
    }