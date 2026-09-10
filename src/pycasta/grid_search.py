#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
grid_search.py

Script to perform a grid search over model parameters.
"""

import itertools
import json
import logging
import os
import numpy as np

from .run_analysis import process_pdb
from .config import (
    OUTPUT_DIR,
    VERSION_TAG,
    FLOW_KWARGS,
    set_alpha,
    reset_alpha,
    set_merge_threshold,
    reset_merge_threshold,
    set_min_pocket_volume,
    reset_min_pocket_volume,
)
alpha_values = [5.0, 5.4, 6.0]
volume_thresholds = [40, 50, 60]
merge_thresholds = [3, 4, 5]
flow_tols = [0.5, 0.6]
flow_max_steps = [30, 50]
flow_adaptive_factors = [0.9, 0.8]

parameter_grid = list(
    itertools.product(
        alpha_values,
        volume_thresholds,
        merge_thresholds,
        flow_tols,
        flow_max_steps,
        flow_adaptive_factors,
    )
)

logging.info(f"Total parameter combinations to test: {len(parameter_grid)}")


def evaluate_result(result, target_min=1, target_max=3):
    num_pockets = len(result.get("ranked_pockets", []))
    if target_min <= num_pockets <= target_max:
        return 0
    else:
        penalty = abs(
            num_pockets - (target_min if num_pockets < target_min else target_max)
        )
        return penalty


test_pdb = "data/bounded_test/1acj.pdb"
best_score = float("inf")
best_params = None
results_summary = []

for combo in parameter_grid:
    alpha_val, vol_thresh, merge_thresh, flow_tol, flow_max, flow_adapt = combo

    set_alpha(alpha_val)
    set_min_pocket_volume(vol_thresh)
    set_merge_threshold(merge_thresh)

    old_tol = FLOW_KWARGS.get("tol_fraction")
    old_max_steps = FLOW_KWARGS.get("max_steps")
    old_adaptive_factor = FLOW_KWARGS.get("adaptive_factor")

    FLOW_KWARGS["tol_fraction"] = flow_tol
    FLOW_KWARGS["max_steps"] = flow_max
    FLOW_KWARGS["adaptive_factor"] = flow_adapt

    logging.info(
        f"Testing parameters: alpha={alpha_val}, volume_threshold={vol_thresh}, "
        f"merge_threshold={merge_thresh}, flow_tol={flow_tol}, "
        f"max_steps={flow_max}, adaptive_factor={flow_adapt}"
    )

    result = process_pdb(test_pdb)
    score = evaluate_result(result)

    logging.info(
        f"Result: {len(result.get('ranked_pockets', []))} pockets, score={score}"
    )

    results_summary.append(
        {
            "params": {
                "alpha": alpha_val,
                "volume_threshold": vol_thresh,
                "merge_threshold": merge_thresh,
                "flow_tol": flow_tol,
                "flow_max_steps": flow_max,
                "flow_adaptive_factor": flow_adapt,
            },
            "num_pockets": len(result.get("ranked_pockets", [])),
            "score": score,
        }
    )

    if score < best_score:
        best_score = score
        best_params = combo

    FLOW_KWARGS["tol_fraction"] = old_tol
    FLOW_KWARGS["max_steps"] = old_max_steps
    FLOW_KWARGS["adaptive_factor"] = old_adaptive_factor
    
reset_alpha()
reset_min_pocket_volume()
reset_merge_threshold()


logging.info(f"Best score: {best_score} with parameters: {best_params}")
output_file = os.path.join(
    OUTPUT_DIR, f"parameter_grid_search_summary_{VERSION_TAG}.json"
)
with open(output_file, "w") as f:
    json.dump(results_summary, f, indent=4)
logging.info(f"Grid search summary saved to: {output_file}")
