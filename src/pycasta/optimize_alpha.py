# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 21:14:08 2025

@author: giorgio
"""

import numpy as np
import pandas as pd
import os
from config import set_alpha, set_merge_threshold, get_merge_threshold
from run_entrypoint import run_all_analysis

# Optional: dynamic setter for flow kwargs
from config import FLOW_KWARGS


def set_flow_params(min_steps=None, min_volume=None, adaptive_factor=None, tol=None):
    if min_steps is not None:
        FLOW_KWARGS["min_steps"] = min_steps
    if min_volume is not None:
        FLOW_KWARGS["min_volume"] = min_volume
    if adaptive_factor is not None:
        FLOW_KWARGS["adaptive_factor"] = adaptive_factor
    if tol is not None:
        FLOW_KWARGS["tol"] = tol


def optimize_alpha_and_flow(
    alpha_values,
    min_steps_values,
    min_volume_values,
    adaptive_factor_values,
    tol_values,
    merge_thresholds,
    output_csv="alpha_flow_optimization_results.csv",
):
    """
    Runs analysis for each combination of alpha and flow parameters,
    collects performance metrics.
    """
    records = []

    for alpha in alpha_values:
        set_alpha(alpha)
        for min_steps in min_steps_values:
            for min_volume in min_volume_values:
                for adaptive_factor in adaptive_factor_values:
                    for tol in tol_values:
                        for merge_threshold in merge_thresholds:
                            set_flow_params(
                                min_steps=min_steps,
                                min_volume=min_volume,
                                adaptive_factor=adaptive_factor,
                                tol=tol,
                            )
                            set_merge_threshold(merge_threshold)
                            print(
                                f"\n🚀 Testing alpha={alpha}, min_steps={min_steps}, min_volume={min_volume}, adaptive_factor={adaptive_factor}, tol={tol}, merge_threshold={merge_threshold}"
                            )
                            try:
                                summary = run_all_analysis()
                                record = {
                                    "alpha": alpha,
                                    "min_steps": min_steps,
                                    "min_volume": min_volume,
                                    "adaptive_factor": adaptive_factor,
                                    "tol": tol,
                                    "merge_threshold": merge_threshold,
                                    "top1_percentage": summary.get(
                                        "top1_percentage", 0
                                    ),
                                    "top3_percentage": summary.get(
                                        "top3_percentage", 0
                                    ),
                                    "top5_percentage": summary.get(
                                        "top5_percentage", 0
                                    ),
                                    "top_all_percentage": summary.get(
                                        "top_all_percentage", 0
                                    ),
                                    "average_step_to_ligand": summary.get(
                                        "average_step_to_ligand", "N/A"
                                    ),
                                    "average_ligand_size": summary.get(
                                        "average_ligand_size", "N/A"
                                    ),
                                    "average_pocket_volume": summary.get(
                                        "average_pocket_volume", "N/A"
                                    ),
                                    "average_pocket_depth": summary.get(
                                        "average_pocket_depth", "N/A"
                                    ),
                                    "average_mouth_area": summary.get(
                                        "average_mouth_area", "N/A"
                                    ),
                                    "average_mouth_perimeter": summary.get(
                                        "average_mouth_perimeter", "N/A"
                                    ),
                                    "total_molecules": summary.get(
                                        "total_molecules", 0
                                    ),
                                }
                                records.append(record)
                            except Exception as e:
                                print(f"❌ Error: {e}")
                                continue

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    print(f"\n✅ Optimization completed. Results saved to: {output_csv}")
    return df


if __name__ == "__main__":
    test_alphas = [2.0, 2.2]
    test_min_steps = [3]
    test_min_volumes = [50]
    test_adaptive_factors = [0.5]
    test_tols = [0.01]
    test_merge_thresholds = [16, 18]

    optimize_alpha_and_flow(
        test_alphas,
        test_min_steps,
        test_min_volumes,
        test_adaptive_factors,
        test_tols,
        test_merge_thresholds,
    )
