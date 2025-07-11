#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pocket_utils.py

Module containing utility functions for pockets:
- Compute distances and mouth parameters.
- Merge nearby pocket clusters.
- Ranking functions and saving pocket properties.
"""

import numpy as np
from scipy.spatial.distance import cdist
import logging
from config import (
    SIMPLE_SPLIT,
    ATOM_SPLITTING,
    ATOM_SPLITTING_MODE,
    SAVE_ALPHA,
    VERSION_TAG,
    OUTPUT_DIR,
    LIGAND_THRESHOLD,
    BETA_RANKING,
    ALPHA_RANKING,
)
import os
import pandas as pd


def compute_ligand_contact_distance(
    protein_coords, ligand_coords, threshold=LIGAND_THRESHOLD
):
    if protein_coords.size == 0 or ligand_coords.size == 0:
        return {}
    distances = cdist(protein_coords, ligand_coords)
    contacting = {}
    for i in range(protein_coords.shape[0]):
        min_dist = distances[i].min()
        if min_dist < threshold:
            contacting[i] = min_dist
    return contacting


def compute_mouth_parameters(
    ranked_pockets, simplices, protein_coords, ligand_coords=None
):
    ligand_center = None
    if ligand_coords is not None and len(ligand_coords) > 0:
        ligand_center = np.mean(ligand_coords, axis=0)
    mouth_params_list = []
    for pocket in ranked_pockets:
        face_count = {}
        for tetra_idx in pocket:
            atoms = simplices[tetra_idx]
            faces = [
                tuple(sorted((atoms[0], atoms[1], atoms[2]))),
                tuple(sorted((atoms[0], atoms[1], atoms[3]))),
                tuple(sorted((atoms[0], atoms[2], atoms[3]))),
                tuple(sorted((atoms[1], atoms[2], atoms[3]))),
            ]
            for face in faces:
                face_count[face] = face_count.get(face, 0) + 1
        boundary_faces = [face for face, count in face_count.items() if count == 1]
        total_area, total_circumference, rim_atoms = 0, 0, set()
        face_centroids = []
        from utils.geometry_utils import compute_triangle_area

        for face in boundary_faces:
            pts = protein_coords[list(face)]
            area = compute_triangle_area(pts)
            total_area += area
            centroid = np.mean(pts, axis=0)
            face_centroids.append(centroid)
            edges = [(face[i], face[j]) for i, j in [(0, 1), (1, 2), (2, 0)]]
            total_circumference += sum(
                np.linalg.norm(protein_coords[e[0]] - protein_coords[e[1]])
                for e in edges
            )
            rim_atoms.update(face)
        mouth_center = (
            np.mean(face_centroids, axis=0) if face_centroids else np.zeros(3)
        )
        if ligand_coords is not None and len(ligand_coords) > 0:
            dists = np.linalg.norm(ligand_coords - mouth_center[None, :], axis=-1)
            ligand_distance = float(np.min(dists))
        else:
            ligand_distance = None
        mouth_params_list.append(
            {
                "mouth_area": round(total_area, 2),
                "mouth_circumference": round(total_circumference, 2),
                "rim_atoms": list(rim_atoms),
                "mouth_center": mouth_center.tolist(),
                "ligand_distance": (
                    round(ligand_distance, 2) if ligand_distance is not None else None
                ),
            }
        )
    return mouth_params_list


def compute_ranking_scores(
    pocket_volumes, pockets, alpha_ranking=ALPHA_RANKING, beta_ranking=BETA_RANKING
):
    ranking_scores = []
    for volume, pocket in zip(pocket_volumes, pockets):
        rank_score = volume + 0.1 * len(pocket)
        ranking_scores.append(rank_score)
    return ranking_scores


def merge_nearby_clusters(pockets, tetra_positions, threshold):
    if not pockets:
        return []
    clusters = list(pockets)
    changed = True
    while changed:
        changed = False
        new_clusters = []
        merged_flags = [False] * len(clusters)
        for i in range(len(clusters)):
            if merged_flags[i]:
                continue
            cluster_i = clusters[i]
            if not cluster_i:
                continue
            tetra_coords_i = np.concatenate(
                [tetra_positions[idx] for idx in cluster_i], axis=0
            )
            centroid_i = np.mean(tetra_coords_i, axis=0)
            for j in range(i + 1, len(clusters)):
                if merged_flags[j]:
                    continue
                cluster_j = clusters[j]
                if not cluster_j:
                    continue
                tetra_coords_j = np.concatenate(
                    [tetra_positions[idx] for idx in cluster_j], axis=0
                )
                centroid_j = np.mean(tetra_coords_j, axis=0)
                if np.linalg.norm(centroid_i - centroid_j) < threshold:
                    cluster_i = list(set(cluster_i) | set(cluster_j))
                    merged_flags[j] = True
                    changed = True
                    tetra_coords_i = np.concatenate(
                        [tetra_positions[idx] for idx in cluster_i], axis=0
                    )
                    centroid_i = np.mean(tetra_coords_i, axis=0)
            new_clusters.append(cluster_i)
        clusters = new_clusters
    return clusters


def save_pocket_properties(result: dict, pdb_id: str, output_dir: str):

    n = len(result.get("volumes", []))
    data = []
    for i in range(n):
        data.append(
            {
                "pocket_index": i + 1,
                "volume": result.get("volumes", [None])[i],
                "depth": result.get("depths", [None])[i],
                "mouth_area": result.get("mouth_area", [None])[i],
                "mouth_perimeter": result.get("mouth_perimeter", [None])[i],
                "min_ligand_mouth_dist": result.get(
                    "ligand_mouth_min_distance", [None]
                )[i],
                "ligand_inside_mesh": result.get("ligand_containment_mesh", [None])[i],
            }
        )

    df = pd.DataFrame(data)
    # Updated output directory to include VERSION_TAG:
    outdir = os.path.join(output_dir, VERSION_TAG, "pocket_properties")
    os.makedirs(outdir, exist_ok=True)
    df.to_csv(os.path.join(outdir, f"{pdb_id}_pocket_properties.csv"), index=False)


def compute_volume_for_tetra_group(tetra_positions, group_indices):
    total_volume = 0.0
    for idx in group_indices:
        tetra = tetra_positions[idx]
        a, b, c, d = tetra[0], tetra[1], tetra[2], tetra[3]
        vol = np.abs(np.dot(np.cross(b - a, c - a), d - a)) / 6.0
        total_volume += vol
    return total_volume


def apply_splitting_to_shared_atoms(pockets, protein_coords, atom_radii):
    if not ATOM_SPLITTING:
        return pockets
    for i in range(len(pockets)):
        for j in range(i + 1, len(pockets)):
            pocket_i = set(pockets[i])
            pocket_j = set(pockets[j])
            common_atoms = pocket_i.intersection(pocket_j)
            if common_atoms:
                pocket_center1 = np.mean(protein_coords[list(pocket_i)], axis=0)
                pocket_center2 = np.mean(protein_coords[list(pocket_j)], axis=0)
                for atom in common_atoms:
                    atom_coord = protein_coords[atom]
                    atom_radius = atom_radii[atom] if atom_radii is not None else 1.0
                    if ATOM_SPLITTING_MODE.upper() == "CAST":
                        frac1, frac2 = split_atom_contribution(
                            atom_coord, atom_radius, pocket_center1, pocket_center2
                        )
                    else:
                        frac1, frac2 = 0.5, 0.5
                    logging.info(
                        f"Atom {atom} shared between pockets {i} and {j}: split fractions {frac1:.2f}, {frac2:.2f}"
                    )
    return pockets


def split_atom_contribution(atom_coord, atom_radius, pocket_center1, pocket_center2):
    midpoint = (pocket_center1 + pocket_center2) / 2
    plane_normal = pocket_center2 - pocket_center1
    norm = np.linalg.norm(plane_normal)
    if norm == 0:
        return (0.5, 0.5)
    plane_normal = plane_normal / norm
    d = np.abs(np.dot(atom_coord - midpoint, plane_normal))
    frac = compute_spherical_cap_fraction(atom_radius, d)
    signed_d = np.dot(atom_coord - midpoint, plane_normal)
    if signed_d >= 0:
        return (1 - frac, frac)
    else:
        return (frac, 1 - frac)


def compute_spherical_cap_fraction(R, d):
    if d >= R:
        return 0
    V_cap = (np.pi * (R - d) ** 2 * (2 * R + d)) / 3.0
    V_sphere = (4 / 3) * np.pi * R**3
    return V_cap / V_sphere
