#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


from . import config as cfg
from .pocket_detection import compute_analytic_pocket_volume
from .utils.pocket_utils import merge_nearby_clusters
from .alpha_connectivity import compute_alpha_connectivity_metrics
from .alpha_shape import compute_alpha_complex_from_tetrahedra, flow_detection
from .alpha_constraints import apply_alpha_constraints
from .config import (
    FLOW_KWARGS,
    MERGE_CLUSTERS,
    MIN_POCKET_VOLUME,
    OUTPUT_DIR,
    VERSION_TAG,
    USE_GRAPH_FILTERED_ALPHA,
    ALPHA_MAX_EDGE_FACTOR,
    ALPHA_LOCAL_SCALE_K,
    get_alpha,
    get_merge_threshold,
)

from .config import (
    FLOW_KWARGS,
    MERGE_CLUSTERS,
    MIN_POCKET_VOLUME,
    OUTPUT_DIR,
    VERSION_TAG,
    USE_GRAPH_FILTERED_ALPHA,
    ALPHA_MAX_EDGE_FACTOR,
    ALPHA_LOCAL_SCALE_K,
    get_alpha,
    get_merge_threshold,
)



def _ensure_output_dir(base_name):
    out_dir = os.path.join(str(cfg.OUTPUT_DIR), cfg.VERSION_TAG, "alpha_diagnostics", base_name)
    os.makedirs(out_dir, exist_ok=True)
    return out_dir


def _sanitize_flow_params(flow_params=None):
    params = dict(cfg.FLOW_KWARGS if flow_params is None else flow_params)
    if "tol_fraction" not in params and "tol" in params:
        params["tol_fraction"] = params["tol"]
    return params


def _compute_histogram_stats(radii):
    radii = np.asarray(radii, dtype=float)
    if radii.size == 0:
        return {
            "count": 0,
            "min": 0.0,
            "max": 0.0,
            "mean": 0.0,
            "median": 0.0,
            "std": 0.0,
            "q05": 0.0,
            "q25": 0.0,
            "q50": 0.0,
            "q75": 0.0,
            "q95": 0.0,
        }

    return {
        "count": int(radii.size),
        "min": float(np.min(radii)),
        "max": float(np.max(radii)),
        "mean": float(np.mean(radii)),
        "median": float(np.median(radii)),
        "std": float(np.std(radii)),
        "q05": float(np.quantile(radii, 0.05)),
        "q25": float(np.quantile(radii, 0.25)),
        "q50": float(np.quantile(radii, 0.50)),
        "q75": float(np.quantile(radii, 0.75)),
        "q95": float(np.quantile(radii, 0.95)),
    }


def _alpha_grid_from_base(base_alpha, radii, alpha_values=None, alpha_offsets=None):
    radii = np.asarray(radii, dtype=float)

    if alpha_values is not None:
        vals = sorted({round(float(a), 4) for a in alpha_values if float(a) > 0})
        return vals

    if alpha_offsets is None:
        alpha_offsets = [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]

    vals = sorted(
        {
            round(float(base_alpha + off), 4)
            for off in alpha_offsets
            if float(base_alpha + off) > 0
        }
    )

    if radii.size > 0:
        q25 = float(np.quantile(radii, 0.25))
        q50 = float(np.quantile(radii, 0.50))
        q75 = float(np.quantile(radii, 0.75))
        vals = sorted(set(vals + [round(q25, 4), round(q50, 4), round(q75, 4)]))

    return vals


def _run_single_alpha_diagnostics(
    protein_coords,
    simplices,
    tetra_positions,
    alpha_value,
    molecule_name=None,
    radii=None,
    min_volume_threshold=None,
    flow_params=None,
    merge_clusters=None,
    merge_threshold=None,
    filter_by_sasa=False,
    sasa_threshold=None,
    pdb_file=None,
):
    if min_volume_threshold is None:
        min_volume_threshold = cfg.MIN_POCKET_VOLUME
    if merge_clusters is None:
        merge_clusters = cfg.MERGE_CLUSTERS
    if merge_threshold is None:
        merge_threshold = cfg.get_merge_threshold()

    flow_params = _sanitize_flow_params(flow_params)

    alpha_mask, local_radii, local_tetra_positions, local_simplices = (
        compute_alpha_complex_from_tetrahedra(
            simplices=simplices.copy(),
            tetra_positions=tetra_positions.copy(),
            alpha_value=alpha_value,
            molecule_name=molecule_name,
            protein_coords=protein_coords,
            radii=radii,
            filter_by_sasa=False,
            pdb_file=None,
            save_alpha_data=False,
        )
    )

    total_tetra = int(len(local_radii))
    alpha_tetra = int(np.sum(alpha_mask))
    alpha_fraction = float(alpha_tetra / total_tetra) if total_tetra > 0 else 0.0

    graph_constraint_stats = {
        "graph_filter_enabled": False,
        "num_alpha_before_graph_filter": alpha_tetra,
        "num_alpha_after_graph_filter": alpha_tetra,
        "num_removed_from_alpha_by_graph_filter": 0,
        "kept_fraction_by_graph_filter": 1.0 if alpha_tetra > 0 else 0.0,
        "edge_ratio_mean": 0.0,
        "edge_ratio_median": 0.0,
        "edge_ratio_max": 0.0,
        "max_edge_length_mean": 0.0,
        "local_scale_mean": 0.0,
    }

    if (
        total_tetra > 0
        and protein_coords is not None
        and len(local_tetra_positions) > 0
    ):
        base_alpha_mask = (local_radii < alpha_value).astype(bool)

        _, graph_constraint_stats = apply_alpha_constraints(
            alpha_mask=base_alpha_mask,
            tetra_positions=local_tetra_positions,
            protein_coords=protein_coords,
            use_graph_filter=cfg.USE_GRAPH_FILTERED_ALPHA,
            max_edge_factor=cfg.ALPHA_MAX_EDGE_FACTOR,
            local_scale_k=cfg.ALPHA_LOCAL_SCALE_K,
        )

        connectivity_metrics = {
            "num_alpha_components": 0,
            "largest_component_size": 0,
            "largest_component_fraction": 0.0,
        }

    if total_tetra > 0 and len(local_simplices) > 0:
        connectivity_metrics = compute_alpha_connectivity_metrics(
            simplices=local_simplices,
            alpha_mask=alpha_mask,
        )

    if total_tetra == 0:
        return {
            "alpha": float(alpha_value),
            "total_tetrahedra": 0,
            "num_alpha_tetrahedra": 0,
            "alpha_fraction": 0.0,
            "graph_filter_enabled": False,
            "num_alpha_before_graph_filter": 0,
            "num_alpha_after_graph_filter": 0,
            "num_removed_from_alpha_by_graph_filter": 0,
            "kept_fraction_by_graph_filter": 0.0,
            "edge_ratio_mean": 0.0,
            "edge_ratio_median": 0.0,
            "edge_ratio_max": 0.0,
            "max_edge_length_mean": 0.0,
            "local_scale_mean": 0.0,
            "num_alpha_components": 0,
            "largest_component_size": 0,
            "largest_component_fraction": 0.0,
            "num_preliminary_pockets": 0,
            "num_merged_pockets": 0,
            "num_final_pockets": 0,
            "mean_flow_steps": 0.0,
            "max_flow_steps": 0,
            "mean_volume": 0.0,
            "median_volume": 0.0,
            "max_volume": 0.0,
            "top3_volumes": [],
            "all_volumes": [],
        }

    if alpha_tetra == 0:
        return {
            "alpha": float(alpha_value),
            "total_tetrahedra": total_tetra,
            "num_alpha_tetrahedra": 0,
            "alpha_fraction": 0.0,
            "graph_filter_enabled": bool(
                graph_constraint_stats.get("graph_filter_enabled", False)
            ),
            "num_alpha_before_graph_filter": int(
                graph_constraint_stats.get("num_alpha_before_graph_filter", 0)
            ),
            "num_alpha_after_graph_filter": int(
                graph_constraint_stats.get("num_alpha_after_graph_filter", 0)
            ),
            "num_removed_from_alpha_by_graph_filter": int(
                graph_constraint_stats.get(
                    "num_removed_from_alpha_by_graph_filter",
                    0,
                )
            ),
            "kept_fraction_by_graph_filter": float(
                graph_constraint_stats.get("kept_fraction_by_graph_filter", 0.0)
            ),
            "edge_ratio_mean": float(
                graph_constraint_stats.get("edge_ratio_mean", 0.0)
            ),
            "edge_ratio_median": float(
                graph_constraint_stats.get("edge_ratio_median", 0.0)
            ),
            "edge_ratio_max": float(
                graph_constraint_stats.get("edge_ratio_max", 0.0)
            ),
            "max_edge_length_mean": float(
                graph_constraint_stats.get("max_edge_length_mean", 0.0)
            ),
            "local_scale_mean": float(
                graph_constraint_stats.get("local_scale_mean", 0.0)
            ),
            "num_alpha_components": int(
                connectivity_metrics.get("num_alpha_components", 0)
            ),
            "largest_component_size": int(
                connectivity_metrics.get("largest_component_size", 0)
            ),
            "largest_component_fraction": float(
                connectivity_metrics.get("largest_component_fraction", 0.0)
            ),
            "num_preliminary_pockets": 0,
            "num_merged_pockets": 0,
            "num_final_pockets": 0,
            "mean_flow_steps": 0.0,
            "max_flow_steps": 0,
            "mean_volume": 0.0,
            "median_volume": 0.0,
            "max_volume": 0.0,
            "top3_volumes": [],
            "all_volumes": [],
        }

    # Use the same radii that were used for alpha filtering
    # (circumradii for standard Delaunay, power radii for weighted Delaunay)
    pockets, flow_steps, connectivity_graph = flow_detection(
        local_simplices,
        alpha_mask,
        local_tetra_positions,
        protein_coords,
        flow_params,
        proxy_values=local_radii,
    )

    max_index = len(protein_coords) - 1
    valid_pockets = [p for p in pockets if all(0 <= idx <= max_index for idx in p)]
    num_preliminary_pockets = len(valid_pockets)
    max_index = len(protein_coords) - 1
    valid_pockets = [p for p in pockets if all(0 <= idx <= max_index for idx in p)]
    num_preliminary_pockets = len(valid_pockets)

    if merge_clusters:
        merged_pockets = merge_nearby_clusters(
            valid_pockets,
            local_tetra_positions,
            merge_threshold,
        )
    else:
        merged_pockets = valid_pockets

    num_merged_pockets = len(merged_pockets)

    pocket_volumes = [
        compute_analytic_pocket_volume(local_tetra_positions, pocket)
        for pocket in merged_pockets
    ]

    filtered_volumes = [
        float(v) for v in pocket_volumes if float(v) >= float(min_volume_threshold)
    ]
    filtered_volumes = sorted(filtered_volumes, reverse=True)

    flow_values = list(flow_steps.values())
    mean_flow_steps = float(np.mean(flow_values)) if flow_values else 0.0
    max_flow_steps = int(np.max(flow_values)) if flow_values else 0

    return {
        "alpha": float(alpha_value),
        "total_tetrahedra": total_tetra,
        "num_alpha_tetrahedra": alpha_tetra,
        "alpha_fraction": alpha_fraction,
        "graph_filter_enabled": bool(
            graph_constraint_stats.get("graph_filter_enabled", False)
        ),
        "num_alpha_before_graph_filter": int(
            graph_constraint_stats.get("num_alpha_before_graph_filter", alpha_tetra)
        ),
        "num_alpha_after_graph_filter": int(
            graph_constraint_stats.get("num_alpha_after_graph_filter", alpha_tetra)
        ),
        "num_removed_from_alpha_by_graph_filter": int(
            graph_constraint_stats.get("num_removed_from_alpha_by_graph_filter", 0)
        ),
        "kept_fraction_by_graph_filter": float(
            graph_constraint_stats.get("kept_fraction_by_graph_filter", 0.0)
        ),
        "edge_ratio_mean": float(graph_constraint_stats.get("edge_ratio_mean", 0.0)),
        "edge_ratio_median": float(
            graph_constraint_stats.get("edge_ratio_median", 0.0)
        ),
        "edge_ratio_max": float(graph_constraint_stats.get("edge_ratio_max", 0.0)),
        "max_edge_length_mean": float(
            graph_constraint_stats.get("max_edge_length_mean", 0.0)
        ),
        "local_scale_mean": float(
            graph_constraint_stats.get("local_scale_mean", 0.0)
        ),
        "num_preliminary_pockets": num_preliminary_pockets,
        "num_merged_pockets": num_merged_pockets,
        "num_final_pockets": len(filtered_volumes),
        "mean_flow_steps": mean_flow_steps,
        "max_flow_steps": max_flow_steps,
        "mean_volume": float(np.mean(filtered_volumes)) if filtered_volumes else 0.0,
        "median_volume": float(np.median(filtered_volumes)) if filtered_volumes else 0.0,
        "max_volume": float(np.max(filtered_volumes)) if filtered_volumes else 0.0,
        "top3_volumes": [round(v, 4) for v in filtered_volumes[:3]],
        "all_volumes": [round(v, 4) for v in filtered_volumes],
    }


def _save_radii_histogram(radii, alpha_values, out_png, title):
    radii = np.asarray(radii, dtype=float)

    plt.figure(figsize=(10, 6))
    if radii.size > 0:
        plt.hist(radii, bins=50, color="#4C78A8", alpha=0.8, edgecolor="white")
        for alpha in alpha_values:
            plt.axvline(alpha, linestyle="--", linewidth=1.5, label=f"alpha={alpha:.2f}")
    plt.xlabel("Circumsphere radius")
    plt.ylabel("Count")
    plt.title(title)
    if radii.size > 0 and len(alpha_values) <= 10:  # <-- FIX: aggiungi radii.size > 0
        plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    plt.close()

def _save_alpha_sweep_plot(df, out_png, title):
    if df.empty:
        return

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    axes[0, 0].plot(df["alpha"], df["num_alpha_tetrahedra"], marker="o", color="#4C78A8")
    axes[0, 0].set_title("Alpha tetrahedra")
    axes[0, 0].set_xlabel("alpha")
    axes[0, 0].set_ylabel("count")

    axes[0, 1].plot(df["alpha"], df["num_preliminary_pockets"], marker="o", label="pre", color="#F58518")
    axes[0, 1].plot(df["alpha"], df["num_merged_pockets"], marker="o", label="merged", color="#54A24B")
    axes[0, 1].plot(df["alpha"], df["num_final_pockets"], marker="o", label="final", color="#E45756")
    axes[0, 1].set_title("Pocket counts")
    axes[0, 1].set_xlabel("alpha")
    axes[0, 1].set_ylabel("count")
    axes[0, 1].legend()

    axes[1, 0].plot(df["alpha"], df["max_volume"], marker="o", color="#72B7B2")
    axes[1, 0].set_title("Top pocket volume")
    axes[1, 0].set_xlabel("alpha")
    axes[1, 0].set_ylabel("volume")

    axes[1, 1].plot(df["alpha"], df["mean_flow_steps"], marker="o", color="#B279A2")
    axes[1, 1].set_title("Mean flow steps")
    axes[1, 1].set_xlabel("alpha")
    axes[1, 1].set_ylabel("steps")

    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    plt.close()


def _save_top3_volume_plot(df, out_png, title):
    if df.empty:
        return

    plot_df = df.copy()
    for i in range(3):
        col = f"top{i+1}_volume"
        plot_df[col] = plot_df["top3_volumes"].apply(lambda x: x[i] if len(x) > i else 0.0)

    plt.figure(figsize=(10, 6))
    plt.plot(plot_df["alpha"], plot_df["top1_volume"], marker="o", label="top1")
    plt.plot(plot_df["alpha"], plot_df["top2_volume"], marker="o", label="top2")
    plt.plot(plot_df["alpha"], plot_df["top3_volume"], marker="o", label="top3")
    plt.xlabel("alpha")
    plt.ylabel("volume")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    plt.close()


def run_alpha_diagnostics(
    protein_coords,
    simplices,
    tetra_positions,
    molecule_name,
    pdb_file=None,
    radii=None,
    base_alpha=None,
    alpha_values=None,
    alpha_offsets=None,
    min_volume_threshold=None,
    flow_params=None,
    merge_clusters=None,
    merge_threshold=None,
    filter_by_sasa=False,
    sasa_threshold=None,
):
    if base_alpha is None:
        base_alpha = cfg.get_alpha()
    if min_volume_threshold is None:
        min_volume_threshold = cfg.MIN_POCKET_VOLUME
    if merge_clusters is None:
        merge_clusters = cfg.MERGE_CLUSTERS
    if merge_threshold is None:
        merge_threshold = cfg.get_merge_threshold()

    # For weighted-only, skip raw_radii computation (not needed)
    raw_radii = np.array([], dtype=float)
    stats = _compute_histogram_stats(raw_radii)
    alpha_grid = _alpha_grid_from_base(
        base_alpha=base_alpha,
        radii=raw_radii,
        alpha_values=alpha_values,
        alpha_offsets=alpha_offsets,
    )

    out_dir = _ensure_output_dir(molecule_name)

    runs = []
    for alpha in alpha_grid:
        logging.info(f"[alpha_diagnostics] {molecule_name}: alpha={alpha:.4f}")
        run = _run_single_alpha_diagnostics(
            protein_coords=protein_coords,
            simplices=simplices,
            tetra_positions=tetra_positions,
            alpha_value=alpha,
            molecule_name=molecule_name,
            radii=radii, 
            min_volume_threshold=min_volume_threshold,
            flow_params=flow_params,
            merge_clusters=merge_clusters,
            merge_threshold=merge_threshold,
            filter_by_sasa=filter_by_sasa,
            sasa_threshold=sasa_threshold,
            pdb_file=pdb_file,
        )
        runs.append(run)

    df = pd.DataFrame(runs)

    csv_path = os.path.join(out_dir, f"{molecule_name}_alpha_sweep_metrics.csv")
    json_path = os.path.join(out_dir, f"{molecule_name}_alpha_diagnostics.json")
    hist_png = os.path.join(out_dir, f"{molecule_name}_alpha_radii_hist.png")
    sweep_png = os.path.join(out_dir, f"{molecule_name}_alpha_sweep.png")
    top3_png = os.path.join(out_dir, f"{molecule_name}_alpha_top3_volumes.png")

    df.to_csv(csv_path, index=False)

    _save_radii_histogram(
        raw_radii,
        alpha_grid,
        hist_png,
        title=f"{molecule_name} — circumsphere radii and tested alphas",
    )
    _save_alpha_sweep_plot(
        df,
        sweep_png,
        title=f"{molecule_name} — alpha sweep diagnostics",
    )
    _save_top3_volume_plot(
        df,
        top3_png,
        title=f"{molecule_name} — top pocket volumes across alpha",
    )

    summary = {
        "molecule_name": molecule_name,
        "pdb_file": pdb_file,
        "base_alpha": float(base_alpha),
        "tested_alphas": [float(a) for a in alpha_grid],
        "radii_stats": stats,
        "config": {
            "min_volume_threshold": float(min_volume_threshold),
            "merge_clusters": bool(merge_clusters),
            "merge_threshold": float(merge_threshold),
            "filter_by_sasa": bool(filter_by_sasa),
            "sasa_threshold": None if sasa_threshold is None else float(sasa_threshold),
            "flow_params": _sanitize_flow_params(flow_params),
        },
        "runs": runs,
        "artifacts": {
            "csv": csv_path,
            "json": json_path,
            "hist_png": hist_png,
            "sweep_png": sweep_png,
            "top3_png": top3_png,
        },
    }

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    return summary