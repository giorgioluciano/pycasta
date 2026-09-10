#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.spatial import KDTree


def compute_tetra_max_edge_lengths(tetra_positions):
    tetra_positions = np.asarray(tetra_positions, dtype=float)

    if tetra_positions.ndim != 3 or tetra_positions.shape[1:] != (4, 3):
        raise ValueError(
            f"Expected tetra_positions shape (n, 4, 3), got {tetra_positions.shape}"
        )

    edge_pairs = [
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 2),
        (1, 3),
        (2, 3),
    ]

    edge_lengths = []
    for i, j in edge_pairs:
        edge_lengths.append(
            np.linalg.norm(tetra_positions[:, i, :] - tetra_positions[:, j, :], axis=1)
        )

    edge_lengths = np.stack(edge_lengths, axis=1)
    max_edge_lengths = np.max(edge_lengths, axis=1)
    return max_edge_lengths


def compute_local_atom_scale(protein_coords, k=6):
    protein_coords = np.asarray(protein_coords, dtype=float)

    if protein_coords.ndim != 2 or protein_coords.shape[1] != 3:
        raise ValueError(
            f"Expected protein_coords shape (n, 3), got {protein_coords.shape}"
        )

    if len(protein_coords) < 2:
        return np.ones(len(protein_coords), dtype=float)

    k_eff = min(max(int(k), 2), len(protein_coords))
    tree = KDTree(protein_coords)
    dists, _ = tree.query(protein_coords, k=k_eff)

    if dists.ndim == 1:
        dists = dists[:, None]

    # skip self-distance in column 0
    if dists.shape[1] > 1:
        local_scale = np.mean(dists[:, 1:], axis=1)
    else:
        local_scale = np.ones(len(protein_coords), dtype=float)

    local_scale = np.where(local_scale > 1e-8, local_scale, 1.0)
    return local_scale


def compute_tetra_local_scale(tetra_positions, protein_coords, k=6):
    tetra_positions = np.asarray(tetra_positions, dtype=float)
    protein_coords = np.asarray(protein_coords, dtype=float)

    atom_local_scale = compute_local_atom_scale(protein_coords, k=k)

    tree = KDTree(protein_coords)
    flat_vertices = tetra_positions.reshape(-1, 3)
    _, nearest_idx = tree.query(flat_vertices, k=1)

    vertex_scale = atom_local_scale[np.asarray(nearest_idx, dtype=int)]
    vertex_scale = vertex_scale.reshape(-1, 4)

    tetra_local_scale = np.mean(vertex_scale, axis=1)
    tetra_local_scale = np.where(tetra_local_scale > 1e-8, tetra_local_scale, 1.0)
    return tetra_local_scale


def compute_tetra_edge_ratio(
    tetra_positions,
    protein_coords,
    local_scale_k=6,
):
    max_edge_lengths = compute_tetra_max_edge_lengths(tetra_positions)
    tetra_local_scale = compute_tetra_local_scale(
        tetra_positions=tetra_positions,
        protein_coords=protein_coords,
        k=local_scale_k,
    )

    edge_ratio = max_edge_lengths / tetra_local_scale
    return edge_ratio, max_edge_lengths, tetra_local_scale


def build_graph_constraint_mask(
    tetra_positions,
    protein_coords,
    max_edge_factor=2.5,
    local_scale_k=6,
):
    edge_ratio, max_edge_lengths, tetra_local_scale = compute_tetra_edge_ratio(
        tetra_positions=tetra_positions,
        protein_coords=protein_coords,
        local_scale_k=local_scale_k,
    )

    constraint_mask = edge_ratio <= float(max_edge_factor)

    stats = {
        "num_total_tetrahedra": int(len(edge_ratio)),
        "num_kept_by_graph_filter": int(np.sum(constraint_mask)),
        "num_removed_by_graph_filter": int(np.sum(~constraint_mask)),
        "kept_fraction_by_graph_filter": (
            float(np.mean(constraint_mask)) if len(constraint_mask) > 0 else 0.0
        ),
        "edge_ratio_mean": float(np.mean(edge_ratio)) if len(edge_ratio) > 0 else 0.0,
        "edge_ratio_median": float(np.median(edge_ratio)) if len(edge_ratio) > 0 else 0.0,
        "edge_ratio_max": float(np.max(edge_ratio)) if len(edge_ratio) > 0 else 0.0,
        "max_edge_length_mean": (
            float(np.mean(max_edge_lengths)) if len(max_edge_lengths) > 0 else 0.0
        ),
        "local_scale_mean": (
            float(np.mean(tetra_local_scale)) if len(tetra_local_scale) > 0 else 0.0
        ),
        "max_edge_factor": float(max_edge_factor),
        "local_scale_k": int(local_scale_k),
    }

    return constraint_mask, stats


def apply_alpha_constraints(
    alpha_mask,
    tetra_positions,
    protein_coords,
    use_graph_filter=True,
    max_edge_factor=2.5,
    local_scale_k=6,
):
    alpha_mask = np.asarray(alpha_mask, dtype=bool)

    if not use_graph_filter:
        stats = {
            "graph_filter_enabled": False,
            "num_alpha_before_graph_filter": int(np.sum(alpha_mask)),
            "num_alpha_after_graph_filter": int(np.sum(alpha_mask)),
            "num_removed_from_alpha_by_graph_filter": 0,
        }
        return alpha_mask, stats

    constraint_mask, constraint_stats = build_graph_constraint_mask(
        tetra_positions=tetra_positions,
        protein_coords=protein_coords,
        max_edge_factor=max_edge_factor,
        local_scale_k=local_scale_k,
    )

    constrained_alpha_mask = alpha_mask & constraint_mask

    stats = {
        "graph_filter_enabled": True,
        "num_alpha_before_graph_filter": int(np.sum(alpha_mask)),
        "num_alpha_after_graph_filter": int(np.sum(constrained_alpha_mask)),
        "num_removed_from_alpha_by_graph_filter": int(
            np.sum(alpha_mask & (~constraint_mask))
        ),
    }
    stats.update(constraint_stats)

    return constrained_alpha_mask, stats