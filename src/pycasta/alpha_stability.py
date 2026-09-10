#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import numpy as np

from .alpha_shape import compute_alpha_complex_from_tetrahedra, flow_detection
from .pocket_detection import compute_analytic_pocket_volume
from .utils.pocket_utils import merge_nearby_clusters
from . import config as cfg


def _round_float_list(values, ndigits=4):
    return [round(float(v), ndigits) for v in values]


def _top_volume_signature(volumes, top_k=3, decimals=1):
    if not volumes:
        return tuple()
    top = sorted((float(v) for v in volumes), reverse=True)[:top_k]
    return tuple(round(v, decimals) for v in top)


def _signature_similarity(sig_a, sig_b):
    a = set(sig_a)
    b = set(sig_b)
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def _single_alpha_run(
    protein_coords,
    simplices,
    tetra_positions,
    alpha_value,
    molecule_name=None,
    min_volume_threshold=None,
    flow_params=None,
    merge_clusters=None,
    merge_threshold=None,
):
    if min_volume_threshold is None:
        min_volume_threshold = cfg.MIN_POCKET_VOLUME
    if flow_params is None:
        flow_params = cfg.FLOW_KWARGS
    if merge_clusters is None:
        merge_clusters = cfg.MERGE_CLUSTERS
    if merge_threshold is None:
        merge_threshold = cfg.get_merge_threshold()

    alpha_mask, local_radii, local_tetra_positions, local_simplices = (
        compute_alpha_complex_from_tetrahedra(
            simplices=simplices.copy(),
            tetra_positions=tetra_positions.copy(),
            alpha_value=alpha_value,
            molecule_name=molecule_name,
            protein_coords=protein_coords,
            save_alpha_data=False,
        )
    )

    if np.sum(alpha_mask) == 0:
        return {
            "alpha": float(alpha_value),
            "num_alpha_tetrahedra": 0,
            "num_preliminary_pockets": 0,
            "num_final_pockets": 0,
            "volumes": [],
            "top_signature": tuple(),
            "mean_volume": 0.0,
            "max_volume": 0.0,
        }

    # Use the same radii that were used for alpha filtering
    # (circumradii for standard Delaunay, power radii for weighted Delaunay)
    pockets, flow_steps, connectivity_graph = flow_detection(
        simplices,
        alpha_mask,
        tetra_positions,
        protein_coords,
        flow_params,
        proxy_values=radii,  # radii is already computed earlier in this function
    )

    max_index = len(protein_coords) - 1
    valid_pockets = [p for p in pockets if all(0 <= idx <= max_index for idx in p)]

    if merge_clusters:
        valid_pockets = merge_nearby_clusters(
            valid_pockets,
            local_tetra_positions,
            merge_threshold,
        )

    pocket_volumes = [
        compute_analytic_pocket_volume(local_tetra_positions, pocket)
        for pocket in valid_pockets
    ]

    filtered = [
        (pocket, volume)
        for pocket, volume in zip(valid_pockets, pocket_volumes)
        if volume >= min_volume_threshold
    ]

    if not filtered:
        return {
            "alpha": float(alpha_value),
            "num_alpha_tetrahedra": int(np.sum(alpha_mask)),
            "num_preliminary_pockets": len(valid_pockets),
            "num_final_pockets": 0,
            "volumes": [],
            "top_signature": tuple(),
            "mean_volume": 0.0,
            "max_volume": 0.0,
        }

    _, filtered_volumes = zip(*filtered)
    filtered_volumes = list(filtered_volumes)

    return {
        "alpha": float(alpha_value),
        "num_alpha_tetrahedra": int(np.sum(alpha_mask)),
        "num_preliminary_pockets": len(valid_pockets),
        "num_final_pockets": len(filtered_volumes),
        "volumes": _round_float_list(filtered_volumes, 4),
        "top_signature": _top_volume_signature(filtered_volumes, top_k=3, decimals=1),
        "mean_volume": float(np.mean(filtered_volumes)),
        "max_volume": float(np.max(filtered_volumes)),
    }


def alpha_stability(
    protein_coords,
    simplices,
    tetra_positions,
    molecule_name=None,
    base_alpha=None,
    alpha_offsets=None,
    min_volume_threshold=None,
    flow_params=None,
    merge_clusters=None,
    merge_threshold=None,
):
    if base_alpha is None:
        base_alpha = cfg.get_alpha()

    if alpha_offsets is None:
        alpha_offsets = [-0.2, -0.1, 0.0, 0.1, 0.2]

    alpha_values = sorted(
        {
            round(float(base_alpha + offset), 4)
            for offset in alpha_offsets
            if (base_alpha + offset) > 0
        }
    )

    runs = []
    for alpha_value in alpha_values:
        logging.info(f"[alpha_stability] Testing alpha={alpha_value:.4f}")
        run = _single_alpha_run(
            protein_coords=protein_coords,
            simplices=simplices.copy(),
            tetra_positions=tetra_positions.copy(),
            alpha_value=alpha_value,
            molecule_name=molecule_name,
            min_volume_threshold=min_volume_threshold,
            flow_params=flow_params,
            merge_clusters=merge_clusters,
            merge_threshold=merge_threshold,
        )
        runs.append(run)

    base_run = min(runs, key=lambda r: abs(r["alpha"] - float(base_alpha)))

    pocket_counts = [r["num_final_pockets"] for r in runs]
    max_volumes = [r["max_volume"] for r in runs]
    mean_volumes = [r["mean_volume"] for r in runs]

    count_range = int(max(pocket_counts) - min(pocket_counts)) if pocket_counts else 0
    max_volume_cv = (
        float(np.std(max_volumes) / np.mean(max_volumes))
        if max_volumes and np.mean(max_volumes) > 0
        else 0.0
    )
    mean_volume_cv = (
        float(np.std(mean_volumes) / np.mean(mean_volumes))
        if mean_volumes and np.mean(mean_volumes) > 0
        else 0.0
    )

    signature_similarities = [
        _signature_similarity(base_run["top_signature"], r["top_signature"])
        for r in runs
    ]
    avg_signature_similarity = (
        float(np.mean(signature_similarities)) if signature_similarities else 0.0
    )

    stable = (
        count_range <= 1
        and max_volume_cv <= 0.25
        and avg_signature_similarity >= 0.6
    )

    summary = {
        "base_alpha": float(base_alpha),
        "tested_alphas": [float(r["alpha"]) for r in runs],
        "runs": [
            {
                "alpha": r["alpha"],
                "num_alpha_tetrahedra": r["num_alpha_tetrahedra"],
                "num_preliminary_pockets": r["num_preliminary_pockets"],
                "num_final_pockets": r["num_final_pockets"],
                "volumes": r["volumes"],
                "top_signature": list(r["top_signature"]),
                "mean_volume": round(float(r["mean_volume"]), 4),
                "max_volume": round(float(r["max_volume"]), 4),
            }
            for r in runs
        ],
        "metrics": {
            "count_range": count_range,
            "max_volume_cv": round(max_volume_cv, 4),
            "mean_volume_cv": round(mean_volume_cv, 4),
            "avg_top_signature_similarity": round(avg_signature_similarity, 4),
        },
        "stable": bool(stable),
    }

    logging.info(
        "[alpha_stability] stable=%s | count_range=%s | max_volume_cv=%.4f | avg_top_signature_similarity=%.4f",
        summary["stable"],
        summary["metrics"]["count_range"],
        summary["metrics"]["max_volume_cv"],
        summary["metrics"]["avg_top_signature_similarity"],
    )

    return summary