#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
alpha_shape.py

Module for computing the alpha shape and for the flow_detection function.
"""

import os
import numpy as np
import logging
from config import (
    SIMPLE_SPLIT,
    ATOM_SPLITTING,
    ATOM_SPLITTING_MODE,
    SAVE_ALPHA,
    VERSION_TAG,
    OUTPUT_DIR,
    FILTER_AFTER_MERGE,
    SASA_CONTACT_THRESHOLD,
)
from utils.geometry_utils import compute_all_circumsphere_radii, log_radii_statistics
from scipy.spatial import KDTree
from utils.pocket_utils import apply_splitting_to_shared_atoms


def filter_alpha_complex_by_sasa(
    tetrahedra, alpha_values, protein_coords, pdb_file, sasa_threshold
):
    from utils.sasa_utils import compute_sasa
    from biopandas.pdb import PandasPdb

    logging.info("Computing SASA for filtering the alpha complex...")
    sasa_dict = compute_sasa(pdb_file, method="pymol", remove_ligand=True)
    if sasa_dict is None:
        logging.error("SASA computation failed; no filtering applied.")
        return tetrahedra, alpha_values, np.ones_like(alpha_values, dtype=bool)
    ppdb = PandasPdb().read_pdb(pdb_file)
    atoms_df = ppdb.df["ATOM"]
    atom_residues = atoms_df["residue_number"].values
    from scipy.spatial import KDTree

    tree = KDTree(protein_coords)
    filtered_indices = []
    total_tetra = len(alpha_values)
    for i, tet in enumerate(tetrahedra):
        centroid = np.mean(tet, axis=0)
        dist, idx = tree.query(centroid)
        residue = atom_residues[idx]
        sasa = sasa_dict.get(str(residue), 0)
        if sasa < sasa_threshold:
            filtered_indices.append(i)
    logging.info(
        f"Tetrahedra: {total_tetra} total, {len(filtered_indices)} pass the SASA filter (threshold={sasa_threshold})."
    )
    filtered_tetrahedra = tetrahedra[filtered_indices]
    filtered_alpha_values = alpha_values[filtered_indices]
    filtered_mask = np.zeros_like(alpha_values, dtype=bool)
    filtered_mask[filtered_indices] = True
    return filtered_tetrahedra, filtered_alpha_values, filtered_mask


def compute_alpha_complex_from_tetrahedra(
    simplices,
    tetra_positions,
    alpha_value,
    molecule_name=None,
    protein_coords=None,
    filter_by_sasa=False,
    sasa_threshold=SASA_CONTACT_THRESHOLD,
    pdb_file=None,
):
    """
    Compute the alpha complex for a set of tetrahedra.

    Returns:
        alpha_mask: Boolean mask of tetrahedra with circumsphere radii below alpha_value.
        radii: The circumsphere radii used.
        tetra_positions: Possibly updated tetra_positions after SASA filtering.
        simplices: The filtered simplices array, so that indices align with tetra_positions.
    """
    assert tetra_positions.ndim == 3 and tetra_positions.shape[1:] == (
        4,
        3,
    ), f"Invalid tetra_positions shape: {tetra_positions.shape}"

    # Compute circumsphere radii and initial alpha mask.
    radii = compute_all_circumsphere_radii(tetra_positions)
    log_radii_statistics(radii)
    alpha_mask = (radii < alpha_value).astype(bool)

    # Optionally, save the unfiltered alpha shape data.
    if SAVE_ALPHA and molecule_name and protein_coords is not None:
        version_dir = os.path.join(OUTPUT_DIR, VERSION_TAG)
        os.makedirs(version_dir, exist_ok=True)
        out_file = os.path.join(version_dir, f"{molecule_name}.alpha.npz")
        np.savez(
            out_file,
            alpha_values=radii,
            alpha_mask=alpha_mask,
            protein_coords=protein_coords,
            tetrahedra=tetra_positions,
        )
        logging.info(f"Alpha shape data saved to {out_file}")

    # If filtering by SASA is enabled, update tetra_positions and simplices.
    if filter_by_sasa:
        if pdb_file is None:
            logging.error("pdb_file must be provided when filtering by SASA.")
        else:
            filtered_tetra, filtered_alpha_values, filtered_mask = (
                filter_alpha_complex_by_sasa(
                    tetra_positions, radii, protein_coords, pdb_file, sasa_threshold
                )
            )
            # Update tetra_positions and radii based on the SASA filter.
            tetra_positions = filtered_tetra
            radii = filtered_alpha_values
            # Recalculate alpha_mask from the filtered radii.
            alpha_mask = radii < alpha_value
            # Also update simplices using the same boolean mask.
            simplices = simplices[filtered_mask]
            version_dir = os.path.join(OUTPUT_DIR, VERSION_TAG)
            os.makedirs(version_dir, exist_ok=True)
            out_file_filtered = os.path.join(
                version_dir, f"{molecule_name}.alpha_filtered.npz"
            )
            np.savez(
                out_file_filtered,
                alpha_values=radii,
                alpha_mask=alpha_mask,
                protein_coords=protein_coords,
                tetrahedra=tetra_positions,
            )
            logging.info(f"Filtered alpha shape data saved to {out_file_filtered}")

    return alpha_mask, radii, tetra_positions, simplices


def flow_detection(
    simplices, alpha_mask, tetra_positions, protein_coords, flow_params=None
):
    """
    Apply discrete flow on the tetrahedra to grow pockets.
    This version uses a relative tolerance (tol_fraction) to decide if a neighbor is lower.
    """
    if flow_params is None:
        flow_params = {}
    sigma_p = flow_params.get("sigma_p", 1.0)
    # Use a relative tolerance fraction (e.g., 0.01 means a 1% drop)
    tol_fraction = flow_params.get("tol_fraction", 0.01)
    max_steps = flow_params.get("max_steps", 100)
    adaptive = flow_params.get("adaptive", False)
    adaptive_factor = flow_params.get("adaptive_factor", 0.9)
    min_steps = flow_params.get("min_steps", 0)

    # Compute proxy values for all tetrahedra (e.g., based on circumsphere radii)
    proxies = compute_all_circumsphere_radii(tetra_positions)
    proxy_values = sigma_p * proxies

    # Get indices of tetrahedra that are outside the alpha complex
    empty_indices = np.where(~alpha_mask)[0]

    # Build a face-to-tetrahedron mapping (neighbors sharing a face)
    face_to_tetra = {}
    for i, simplex in enumerate(simplices):
        faces = [
            tuple(sorted(simplex[[a, b, c]]))
            for (a, b, c) in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
        ]
        for face in faces:
            face_to_tetra.setdefault(face, []).append(i)

    # Build neighbor dictionary for each empty tetrahedron
    neighbor_dict = {}
    for i in empty_indices:
        faces = [
            tuple(sorted(simplices[i, [a, b, c]]))
            for (a, b, c) in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
        ]
        neighbors = set()
        for face in faces:
            neighbors.update(face_to_tetra.get(face, []))
        # Exclude the tetrahedron itself and those that are not in the empty region
        neighbor_dict[i] = [j for j in neighbors if j != i and not alpha_mask[j]]

    # Flow: for each tetrahedron in empty_indices, try to follow a gradient descent
    flow_target = {}
    flow_steps = {}
    total_steps = 0

    for i in empty_indices:
        current = i
        steps = 0
        current_proxy = proxy_values[current]
        while steps < max_steps:
            # Use relative tolerance: neighbor qualifies if its proxy is at least tol_fraction lower than current
            lower_neighbors = [
                j
                for j in neighbor_dict[current]
                if (current_proxy - proxy_values[j]) > tol_fraction * current_proxy
            ]
            if not lower_neighbors:
                if adaptive:
                    # Optionally, if no neighbor is found, increase tolerance (or adapt in some way)
                    tol_fraction *= adaptive_factor
                    lower_neighbors = [
                        j
                        for j in neighbor_dict[current]
                        if (current_proxy - proxy_values[j])
                        > tol_fraction * current_proxy
                    ]
                    if not lower_neighbors:
                        break
                else:
                    break
            next_current = min(lower_neighbors, key=lambda j: proxy_values[j])
            current = next_current
            current_proxy = proxy_values[current]
            steps += 1
        if steps >= min_steps:
            flow_target[i] = current
            flow_steps[i] = steps
            total_steps += steps

    average_steps = total_steps / max(len(flow_target), 1)
    logging.info(f"Average number of flow steps per tetrahedron: {average_steps:.2f}")

    # Group tetrahedra by their final sink (flow target)
    sink_groups = {}
    for i, sink in flow_target.items():
        sink_groups.setdefault(sink, []).append(i)

    final_pockets = []
    connectivity_graph = {}
    # For each sink group, build a connectivity graph among tetrahedra
    for sink, tetra_indices in sink_groups.items():
        graph = {i: [] for i in tetra_indices}
        for idx, i in enumerate(tetra_indices):
            for j in tetra_indices[idx + 1 :]:
                if len(set(simplices[i]) & set(simplices[j])) == 3:
                    graph[i].append(j)
                    graph[j].append(i)
        visited = set()
        for node in tetra_indices:
            if node not in visited:
                component = []
                stack = [node]
                while stack:
                    cur = stack.pop()
                    if cur not in visited:
                        visited.add(cur)
                        component.append(cur)
                        stack.extend(graph[cur])
                final_pockets.append(component)
                connectivity_graph[len(final_pockets) - 1] = graph

    logging.info(f"Flow detection found {len(final_pockets)} preliminary pockets.")
    final_pockets = apply_splitting_to_shared_atoms(
        final_pockets, protein_coords, flow_params.get("atom_radii")
    )
    return final_pockets, flow_steps, connectivity_graph
