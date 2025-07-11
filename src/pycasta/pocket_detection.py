#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pocket_detection.py

Module for detecting and ranking pockets in a protein structure.
Uses alpha shape calculation, discrete flow, and cluster merging.
"""

import logging
import numpy as np
from alpha_shape import compute_alpha_complex_from_tetrahedra, flow_detection
from ranking import compute_ranking_scores
from config import MIN_POCKET_VOLUME, FLOW_KWARGS, MERGE_CLUSTERS, MERGE_THRESHOLD
from utils.pocket_utils import merge_nearby_clusters
from dual_sets import compute_dual_set_for_pocket, compute_dual_mouths
from config import get_alpha, get_merge_threshold


# Analytic pocket volume computation defined here
def compute_analytic_pocket_volume(tetra_positions, pocket_indices):
    pocket_tetra = tetra_positions[pocket_indices]
    vec1 = pocket_tetra[:, 1] - pocket_tetra[:, 0]
    vec2 = pocket_tetra[:, 2] - pocket_tetra[:, 0]
    vec3 = pocket_tetra[:, 3] - pocket_tetra[:, 0]
    volumes = np.abs(np.einsum("ij,ij->i", np.cross(vec1, vec2), vec3)) / 6.0
    return float(np.sum(volumes))


def detect_pockets(
    protein_coords,
    simplices=None,
    tetra_positions=None,
    alpha_mask=None,
    min_volume_threshold=None,
    flow_params=None,
    merge_clusters=None,
    merge_threshold=None,
    molecule_name=None,
    radii=None,
):
    if min_volume_threshold is None:
        min_volume_threshold = MIN_POCKET_VOLUME
    if flow_params is None:
        flow_params = FLOW_KWARGS
    if merge_clusters is None:
        merge_clusters = MERGE_CLUSTERS
    if merge_threshold is None:
        merge_threshold = get_merge_threshold()
    logging.info(f"[CONFIG] Merge threshold currently set to: {merge_threshold}")
    if tetra_positions is None:
        tetra_positions = protein_coords[simplices]
    logging.info(f"Delaunay triangulation: {len(simplices)} tetrahedra found.")
    max_index = len(protein_coords) - 1
    if np.max(simplices) > max_index:
        raise ValueError(
            f"ERROR: Some simplex indices exceed protein_coords size ({max_index})."
        )

    from config import get_alpha

    alpha_value = get_alpha()

    if alpha_mask is None or radii is None:
        alpha_mask, radii, tetra_positions, simplices = (
            compute_alpha_complex_from_tetrahedra(
                simplices, tetra_positions, alpha_value, molecule_name, protein_coords
            )
        )

    logging.info(f"Using alpha shape with alpha = {alpha_value}")
    if alpha_mask is None:
        alpha_mask, _ = compute_alpha_complex_from_tetrahedra(
            tetra_positions, alpha_value, molecule_name, protein_coords
        )
    logging.info(
        f"Alpha shape retained {int(np.sum(alpha_mask))} tetrahedra out of {len(alpha_mask)}."
    )
    if np.sum(alpha_mask) == 0:
        logging.warning("Alpha shape filtered out all tetrahedra. No pockets detected.")
        return [], [], []
    pockets, flow_steps, connectivity_graph = flow_detection(
        simplices, alpha_mask, tetra_positions, protein_coords, flow_params
    )
    logging.info(f"Flow detection found {len(pockets)} preliminary pockets.")
    valid_pockets = []
    for pocket in pockets:
        if all(0 <= idx <= max_index for idx in pocket):
            valid_pockets.append(pocket)
    if merge_clusters:
        pockets = merge_nearby_clusters(valid_pockets, tetra_positions, merge_threshold)
        logging.info(
            f"After merging (threshold={merge_threshold}), {len(pockets)} pockets remain."
        )
    else:
        pockets = valid_pockets
    pocket_volumes = [
        compute_analytic_pocket_volume(tetra_positions, pocket) for pocket in pockets
    ]
    filtered = [
        (p, v) for p, v in zip(pockets, pocket_volumes) if v >= min_volume_threshold
    ]
    if not filtered:
        logging.warning(
            f"No pockets exceed the minimum volume threshold of {min_volume_threshold} Å³."
        )
        return {
            "ranked_pockets": [],
            "representative_points": [],
            "ranking_scores": [],
            "dual_sets_info": [],
            "mouths_info": [],
        }
    pockets, pocket_volumes = zip(*filtered)
    ranking_scores = compute_ranking_scores(
        pocket_volumes, pockets, flow_steps, connectivity_graph, tetra_positions
    )
    ranked_data = sorted(
        zip(pockets, pocket_volumes, ranking_scores), key=lambda x: x[2], reverse=True
    )
    ranked_pockets, pocket_volumes, ranking_scores = zip(*ranked_data)

    for pocket in ranked_pockets:
        dual_info = compute_dual_set_for_pocket(
            pocket, tetra_positions, radii, alpha_value=alpha_value, decimals=3
        )
    representative_points = [
        np.mean(tetra_positions[p].reshape(-1, 3), axis=0) for p in ranked_pockets
    ]
    logging.info(f"Final number of pockets in {molecule_name}: {len(ranked_pockets)}")
    # Compute dual set information for each pocket:
    dual_sets_info = (
        []
    )  # Will store the dual set info (closure, dual set, boundary, etc.)
    mouths_info = []  # Will store the computed mouth dual components for each pocket

    for pocket in ranked_pockets:
        # Compute the dual set info for this pocket.
        # Here, we use the current MANUAL_ALPHA as a proxy threshold.
        dual_info = compute_dual_set_for_pocket(
            pocket, tetra_positions, radii, alpha_value=alpha_value, decimals=3
        )
        dual_sets_info.append(dual_info)
        # Compute the dual mouth components from the boundary of the closure.
        mouth_duals = compute_dual_mouths(dual_info["boundary"], global_dual=None)
        mouths_info.append(mouth_duals)

    # Optionally, store these dual representations in the result dictionary:
    # Return all relevant information in a dictionary:
    return {
        "ranked_pockets": ranked_pockets,
        "representative_points": representative_points,
        "ranking_scores": ranking_scores,
        "dual_sets_info": dual_sets_info,
        "mouths_info": mouths_info,
    }
