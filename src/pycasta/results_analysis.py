import json
import logging
import os
from collections import Counter

import numpy as np
import pandas as pd
from tabulate import tabulate

from config import VALIDATION_METHOD


def safe_sum_bool(lst):
    return sum(bool(x) for x in lst) if isinstance(lst, list) else "N/A"


def safe_any(lst):
    return any(bool(x) for x in lst) if isinstance(lst, list) else False


def safe_distances(dist_list):
    return [d for d in dist_list if isinstance(d, (int, float)) and not np.isnan(d)]


def _get_validation_containment_and_distances(result):
    if VALIDATION_METHOD == "mesh_extrusion":
        containment = result.get("ligand_containment_mesh", [])
        distances = safe_distances(result.get("ligand_mesh_distances", []))
    elif VALIDATION_METHOD == "sasa":
        containment = result.get("ligand_containment_strict", [])
        distances = safe_distances(result.get("ligand_to_pocket_distances", []))
    else:
        containment = []
        distances = []
    return containment, distances


def _is_paired(results):
    return bool(results) and "bounded_file" in results[0]


def summarize_results(results, analysis_type="single"):
    total = len(results)
    if total == 0:
        logging.warning("No valid results found to summarize.")
        print("No results to summarize.")
        return {}, pd.DataFrame()

    valid_steps = [
        res.get("step_to_ligand")
        for res in results
        if isinstance(res.get("step_to_ligand"), int)
    ]

    top1 = sum(1 for s in valid_steps if s == 1)
    top3 = sum(1 for s in valid_steps if s <= 3)
    top5 = sum(1 for s in valid_steps if s <= 5)
    top_all = len(valid_steps)

    ligand_sizes = [len(res.get("ligand_coords", [])) for res in results]
    pocket_volumes = [
        vol
        for res in results
        for vol in res.get("pocket_volumes", [])
        if isinstance(vol, (int, float))
    ]
    pocket_depths = [
        depth
        for res in results
        for depth in res.get("pocket_depths", [])
        if isinstance(depth, (int, float))
    ]
    mouth_areas = [
        area
        for res in results
        for area in res.get("mouth_area", [])
        if isinstance(area, (int, float))
    ]
    mouth_perimeters = [
        perim
        for res in results
        for perim in res.get("mouth_perimeter", [])
        if isinstance(perim, (int, float))
    ]

    summary = {
        "total_molecules": total,
        "top1_percentage": round(top1 / total * 100, 2),
        "top3_percentage": round(top3 / total * 100, 2),
        "top5_percentage": round(top5 / total * 100, 2),
        "top_all_percentage": round(top_all / total * 100, 2),
        "average_step_to_ligand": round(np.mean(valid_steps), 2) if valid_steps else "N/A",
        "average_ligand_size": round(np.mean(ligand_sizes), 2) if ligand_sizes else "N/A",
        "average_pocket_volume": round(np.mean(pocket_volumes), 2) if pocket_volumes else "N/A",
        "average_pocket_depth": round(np.mean(pocket_depths), 2) if pocket_depths else "N/A",
        "average_mouth_area": round(np.mean(mouth_areas), 2) if mouth_areas else "N/A",
        "average_mouth_perimeter": round(np.mean(mouth_perimeters), 2) if mouth_perimeters else "N/A",
    }

    table = []

    if _is_paired(results):
        header = [
            "File Pair",
            "Top1",
            "Top3",
            "Top5",
            "Rank (B/U)",
            "#Pockets",
            "#Validated",
            "Validation",
            "Contact",
            "Min. Dist. (Å)",
            "RMSD (Å)",
        ]
        for r in results:
            file_pair = f"{r.get('bounded_file', 'N/A')} / {r.get('unbounded_file', 'N/A')}"
            top1_val = "Yes" if r.get("is_top1_bounded") and r.get("is_top1_unbounded") else "No"
            top3_val = "Yes" if r.get("is_top3_bounded") and r.get("is_top3_unbounded") else "No"
            top5_val = "Yes" if r.get("is_top5_bounded") and r.get("is_top5_unbounded") else "No"
            validation_used = VALIDATION_METHOD.capitalize()

            containment, distances = _get_validation_containment_and_distances(r)
            num_validated = safe_sum_bool(containment)
            contact = "Yes" if safe_any(containment) else "No"
            rank_info = f"{r.get('bounded_pocket_rank', 'None')}/{r.get('unbounded_pocket_rank', 'None')}"
            num_pockets = len(r.get("pocket_volumes", []))
            min_distance = round(min(distances), 2) if distances else "N/A"
            rmsd = round(r.get("alignment_rmsd", float("nan")), 2)

            table.append([
                file_pair,
                top1_val,
                top3_val,
                top5_val,
                rank_info,
                num_pockets,
                num_validated,
                validation_used,
                contact,
                min_distance,
                rmsd,
            ])

        print("\nPaired Analysis Summary:")
        print(tabulate(table, headers=header, tablefmt="github"))

    else:
        header = [
            "File",
            "Top1",
            "Top3",
            "Top5",
            "Top All",
            "Rank",
            "Validation",
            "#Pockets",
            "#Validated",
            "Contact",
            "Min. Dist. (Å)",
            "Avg. Volume",
            "Avg. Depth",
        ]
        for res in results:
            filename = os.path.basename(res.get("pdb_path", res.get("file", "N/A")))
            step = res.get("step_to_ligand")
            top1_val = "Yes" if step == 1 else "No"
            top3_val = "Yes" if isinstance(step, int) and step <= 3 else "No"
            top5_val = "Yes" if isinstance(step, int) and step <= 5 else "No"
            top_all_val = "Yes" if isinstance(step, int) else "No"
            rank = step if isinstance(step, int) else "N/A"
            validation_used = VALIDATION_METHOD.capitalize()

            containment, distances = _get_validation_containment_and_distances(res)
            num_validated = safe_sum_bool(containment)
            contact = "Yes" if safe_any(containment) else "No"
            num_pockets = res.get("num_pockets", len(res.get("pocket_volumes", [])))
            min_distance = round(min(distances), 2) if distances else "N/A"
            avg_volume = round(np.mean(res.get("pocket_volumes", [])), 2) if res.get("pocket_volumes") else "N/A"
            avg_depth = round(np.mean(res.get("pocket_depths", [])), 2) if res.get("pocket_depths") else "N/A"

            table.append([
                filename,
                top1_val,
                top3_val,
                top5_val,
                top_all_val,
                rank,
                validation_used,
                num_pockets,
                num_validated,
                contact,
                min_distance,
                avg_volume,
                avg_depth,
            ])

        print("\nSummary Table:")
        print(tabulate(table, headers=header, tablefmt="github"))

    print("\nPerformance Summary:")
    perf = [
        ["Top1 %", summary["top1_percentage"]],
        ["Top3 %", summary["top3_percentage"]],
        ["Top5 %", summary["top5_percentage"]],
        ["Top All %", summary["top_all_percentage"]],
    ]
    print(tabulate(perf, tablefmt="github"))

    df = pd.DataFrame(table, columns=header)
    return summary, df


def save_summary_json(results, output_json, analysis_type="single"):
    summary, df = summarize_results(results, analysis_type=analysis_type)
    payload = {
        "analysis_type": analysis_type,
        "summary": summary,
        "table": df.to_dict(orient="records"),
    }
    os.makedirs(os.path.dirname(output_json), exist_ok=True)
    with open(output_json, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
    logging.info(f"Summary saved to JSON: {output_json}")
    return payload


def save_summary_csv(results, output_csv, analysis_type="single"):
    _, df = summarize_results(results, analysis_type=analysis_type)
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df.to_csv(output_csv, index=False)
    logging.info(f"Summary table saved to CSV: {output_csv}")
    return df