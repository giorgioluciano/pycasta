# features/occlusion_features.py
import numpy as np
from scipy.spatial import KDTree


def fibonacci_sphere(n_points):
    points = []
    phi = np.pi * (3.0 - np.sqrt(5.0))
    for i in range(n_points):
        y = 1.0 - (i / float(max(n_points - 1, 1))) * 2.0
        radius = np.sqrt(max(0.0, 1.0 - y * y))
        theta = phi * i
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        points.append([x, y, z])
    return np.asarray(points, dtype=float)


def ray_sphere_intersection(origin, direction, center, radius):
    oc = origin - center
    b = 2.0 * np.dot(oc, direction)
    c = np.dot(oc, oc) - radius * radius
    disc = b * b - 4.0 * c

    if disc < 0:
        return None

    sqrt_disc = np.sqrt(disc)
    t1 = (-b - sqrt_disc) / 2.0
    t2 = (-b + sqrt_disc) / 2.0
    valid_ts = [t for t in (t1, t2) if t > 1e-6]
    return min(valid_ts) if valid_ts else None


def _normalize_01(value, low, high):
    if high <= low:
        return 0.0
    return float(np.clip((value - low) / (high - low), 0.0, 1.0))


def compute_hybrid_occlusion_features(
    candidate,
    protein_coords,
    protein_atoms,
    atom_radii,
    sasa_by_residue=None,
    num_rays=96,
    max_distance=12.0,
    neighbor_radius=16.0,
    low_sasa_threshold=5.0,
    kd_tree=None,
):
    center = np.asarray(candidate.representative_point, dtype=float)
    protein_coords = np.asarray(protein_coords, dtype=float)
    atom_radii = np.asarray(atom_radii, dtype=float)

    if kd_tree is None:
        kd_tree = KDTree(protein_coords)

    neighbor_indices = kd_tree.query_ball_point(center, r=neighbor_radius)

    if not neighbor_indices:
        return {
            "ray_occlusion_fraction": 0.0,
            "ray_escape_fraction": 1.0,
            "ray_mean_hit_distance": float(max_distance),
            "ray_min_hit_distance": float(max_distance),
            "hit_residue_sasa_mean": 0.0,
            "low_sasa_hit_fraction": 0.0,
            "hybrid_occlusion_score": 0.0,
            "occlusion_debug": {
                "num_hits": 0,
                "num_rays": int(num_rays),
                "hit_atom_indices": [],
            },
        }

    local_coords = protein_coords[neighbor_indices]
    local_radii = atom_radii[neighbor_indices]
    directions = fibonacci_sphere(num_rays)

    hit_distances = []
    hit_atom_indices = []
    hit_residues = []

    for direction in directions:
        direction = direction / np.linalg.norm(direction)
        best_t = None
        best_atom_index = None

        for local_i, (atom_center, atom_radius) in enumerate(zip(local_coords, local_radii)):
            t_hit = ray_sphere_intersection(center, direction, atom_center, atom_radius)
            if t_hit is None or t_hit > max_distance:
                continue
            if best_t is None or t_hit < best_t:
                best_t = t_hit
                best_atom_index = neighbor_indices[local_i]

        if best_atom_index is not None:
            hit_distances.append(best_t)
            hit_atom_indices.append(int(best_atom_index))
            residue_id = str(protein_atoms[best_atom_index].get("resid", ""))
            hit_residues.append(residue_id)

    num_hits = len(hit_distances)
    ray_occlusion_fraction = num_hits / float(num_rays)
    ray_escape_fraction = 1.0 - ray_occlusion_fraction
    ray_mean_hit_distance = float(np.mean(hit_distances)) if hit_distances else float(max_distance)
    ray_min_hit_distance = float(np.min(hit_distances)) if hit_distances else float(max_distance)

    sasa_values = []
    if sasa_by_residue:
        for residue_id in set(hit_residues):
            sasa_values.append(float(sasa_by_residue.get(str(residue_id), 0.0)))

    hit_residue_sasa_mean = float(np.mean(sasa_values)) if sasa_values else 0.0
    low_sasa_hit_fraction = (
        float(np.mean([value <= low_sasa_threshold for value in sasa_values]))
        if sasa_values else 0.0
    )

    mean_hit_distance_norm = _normalize_01(ray_mean_hit_distance, 0.0, max_distance)
    sasa_norm = _normalize_01(hit_residue_sasa_mean, 0.0, 30.0)

    hybrid_occlusion_score = (
        0.55 * ray_occlusion_fraction +
        0.25 * (1.0 - mean_hit_distance_norm) +
        0.20 * (1.0 - sasa_norm)
    )

    return {
        "ray_occlusion_fraction": float(ray_occlusion_fraction),
        "ray_escape_fraction": float(ray_escape_fraction),
        "ray_mean_hit_distance": float(ray_mean_hit_distance),
        "ray_min_hit_distance": float(ray_min_hit_distance),
        "hit_residue_sasa_mean": float(hit_residue_sasa_mean),
        "low_sasa_hit_fraction": float(low_sasa_hit_fraction),
        "hybrid_occlusion_score": float(hybrid_occlusion_score),
        "occlusion_debug": {
            "num_hits": int(num_hits),
            "num_rays": int(num_rays),
            "hit_atom_indices": hit_atom_indices[:128],
        },
    }
    
def compute_hybrid_occlusion_features_v2(**kwargs):
    return compute_hybrid_occlusion_features(**kwargs)