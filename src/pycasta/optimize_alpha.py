import pandas as pd

from .config import (
    set_alpha,
    set_merge_threshold,
    set_min_pocket_volume,
    reset_alpha,
    reset_merge_threshold,
    reset_min_pocket_volume,
    FLOW_KWARGS,
)
from .run_analysis import run_all_analysis


def set_flow_params(min_steps=None, adaptive_factor=None, tol_fraction=None):
    if min_steps is not None:
        FLOW_KWARGS["min_steps"] = min_steps
    if adaptive_factor is not None:
        FLOW_KWARGS["adaptive_factor"] = adaptive_factor
    if tol_fraction is not None:
        FLOW_KWARGS["tol_fraction"] = tol_fraction


def optimize_alpha_and_flow(
    alpha_values,
    min_steps_values,
    min_volume_values,
    adaptive_factor_values,
    tol_values,
    merge_thresholds,
    output_csv="alpha_flow_optimization_results.csv",
):
    records = []

    for alpha in alpha_values:
        set_alpha(alpha)
        for min_steps in min_steps_values:
            for min_volume in min_volume_values:
                for adaptive_factor in adaptive_factor_values:
                    for tol in tol_values:
                        for merge_threshold in merge_thresholds:
                            set_min_pocket_volume(min_volume)
                            set_flow_params(
                                min_steps=min_steps,
                                adaptive_factor=adaptive_factor,
                                tol_fraction=tol,
                            )
                            set_merge_threshold(merge_threshold)

                            print(
                                f"\n🚀 Testing alpha={alpha}, min_steps={min_steps}, "
                                f"min_volume={min_volume}, adaptive_factor={adaptive_factor}, "
                                f"tol={tol}, merge_threshold={merge_threshold}"
                            )

                            try:
                                summary = run_all_analysis()
                                records.append({
                                    "alpha": alpha,
                                    "min_steps": min_steps,
                                    "min_volume": min_volume,
                                    "adaptive_factor": adaptive_factor,
                                    "tol": tol,
                                    "merge_threshold": merge_threshold,
                                    "top1_percentage": summary.get("top1_percentage", 0),
                                    "top3_percentage": summary.get("top3_percentage", 0),
                                    "top5_percentage": summary.get("top5_percentage", 0),
                                    "top_all_percentage": summary.get("top_all_percentage", 0),
                                    "average_step_to_ligand": summary.get("average_step_to_ligand", "N/A"),
                                    "average_ligand_size": summary.get("average_ligand_size", "N/A"),
                                    "average_pocket_volume": summary.get("average_pocket_volume", "N/A"),
                                    "average_pocket_depth": summary.get("average_pocket_depth", "N/A"),
                                    "average_mouth_area": summary.get("average_mouth_area", "N/A"),
                                    "average_mouth_perimeter": summary.get("average_mouth_perimeter", "N/A"),
                                    "total_molecules": summary.get("total_molecules", 0),
                                })
                            except Exception as e:
                                print(f"❌ Error: {e}")
                                continue

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)

    reset_alpha()
    reset_merge_threshold()
    reset_min_pocket_volume()

    print(f"\n✅ Optimization completed. Results saved to: {output_csv}")
    return df

if __name__ == "__main__":
    test_alphas = [2.0, 2.2]
    test_min_steps = [3]
    test_min_volumes = [30,80]
    test_adaptive_factors = [0.5]
    test_tols = [0.01]
    test_merge_thresholds = [18,24]

    optimize_alpha_and_flow(
        test_alphas,
        test_min_steps,
        test_min_volumes,
        test_adaptive_factors,
        test_tols,
        test_merge_thresholds,
    )
