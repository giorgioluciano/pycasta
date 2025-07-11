#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
data_utils.py

Module for saving and exporting results (JSON, CSV) and converting NumPy types.
"""

import os
import json
import logging
import numpy as np
import csv
from config import VALIDATION_METHOD


def np_converter(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.integer, np.floating)):
        return obj.item()
    elif isinstance(obj, list):
        return [np_converter(item) for item in obj]
    elif isinstance(obj, dict):
        return {key: np_converter(value) for key, value in obj.items()}
    return obj


def convert_numpy_types(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, (np.bool_,)):
        return bool(obj)
    elif isinstance(obj, tuple):
        return [convert_numpy_types(x) for x in obj]
    elif isinstance(obj, set):
        return [convert_numpy_types(x) for x in obj]
    elif isinstance(obj, list):
        return [convert_numpy_types(x) for x in obj]
    elif isinstance(obj, dict):
        new_dict = {}
        for k, v in obj.items():
            # Ensure keys are acceptable
            if not isinstance(k, (str, int, float, bool)) and k is not None:
                k = str(k)
            new_dict[k] = convert_numpy_types(v)
        return new_dict
    return obj


def save_json(data, filename):
    with open(filename, "w") as f:
        json.dump(convert_numpy_types(data), f, indent=4)


def load_json(in_path):
    if not os.path.exists(in_path):
        logging.warning(f"JSON file not found: {in_path}")
        return {}
    with open(in_path, "r") as f:
        content = f.read().strip()
        if not content:
            return {}
        try:
            return json.loads(content)
        except json.decoder.JSONDecodeError as e:
            logging.warning(f"JSON decode error in {in_path}: {e}")
            return {}


def summarize_and_print_results(
    results, analysis_type="single", output_file=None, csv_output_file=None
):
    from tabulate import tabulate
    import numpy as np
    import pandas as pd
    import os
    import logging
    from config import VALIDATION_METHOD

    def safe_sum_bool(lst):
        return sum(bool(x) for x in lst) if isinstance(lst, list) else "N/A"

    def safe_any(lst):
        return any(bool(x) for x in lst) if isinstance(lst, list) else False

    def safe_distances(dist_list):
        return [d for d in dist_list if isinstance(d, (int, float)) and not np.isnan(d)]

    total = len(results)
    if total == 0:
        logging.warning("No valid results found to summarize.")
        print("No results to summarize.")
        return

    is_paired = "bounded_file" in results[0]
    valid_steps = [
        res["step_to_ligand"]
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
        "average_step_to_ligand": (
            round(np.mean(valid_steps), 2) if valid_steps else "N/A"
        ),
        "average_ligand_size": (
            round(np.mean(ligand_sizes), 2) if ligand_sizes else "N/A"
        ),
        "average_pocket_volume": (
            round(np.mean(pocket_volumes), 2) if pocket_volumes else "N/A"
        ),
        "average_pocket_depth": (
            round(np.mean(pocket_depths), 2) if pocket_depths else "N/A"
        ),
        "average_mouth_area": round(np.mean(mouth_areas), 2) if mouth_areas else "N/A",
        "average_mouth_perimeter": (
            round(np.mean(mouth_perimeters), 2) if mouth_perimeters else "N/A"
        ),
    }

    if output_file:

        save_json(summary, output_file)
        logging.info(f"Summary saved to: {output_file}")

    table = []

    if is_paired:
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
            file_pair = (
                f"{r.get('bounded_file','N/A')} / {r.get('unbounded_file','N/A')}"
            )
            top1_val = (
                "Yes"
                if r.get("is_top1_bounded") and r.get("is_top1_unbounded")
                else "No"
            )
            top3_val = (
                "Yes"
                if r.get("is_top3_bounded") and r.get("is_top3_unbounded")
                else "No"
            )
            top5_val = (
                "Yes"
                if r.get("is_top5_bounded") and r.get("is_top5_unbounded")
                else "No"
            )

            validation_used = VALIDATION_METHOD.capitalize()

            if VALIDATION_METHOD == "mesh_extrusion":
                containment = r.get("ligand_containment_mesh", [])
                distances = safe_distances(r.get("ligand_mesh_distances", []))
            elif VALIDATION_METHOD == "sasa":
                containment = r.get("ligand_containment_strict", [])
                distances = safe_distances(r.get("ligand_to_pocket_distances", []))
            else:
                containment = []
                distances = []

            num_validated = safe_sum_bool(containment)
            contact = "Yes" if safe_any(containment) else "No"
            rank_info = f"{r.get('bounded_pocket_rank', 'None')}/{r.get('unbounded_pocket_rank', 'None')}"
            num_pockets = len(r.get("pocket_volumes", []))
            min_distance = round(min(distances), 2) if distances else "N/A"
            rmsd = round(r.get("alignment_rmsd", float("nan")), 2)

            table.append(
                [
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
                ]
            )

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
            top_all = "Yes" if isinstance(step, int) else "No"
            rank = step if isinstance(step, int) else "N/A"
            validation_used = VALIDATION_METHOD.capitalize()

            if VALIDATION_METHOD == "mesh_extrusion":
                containment = res.get("ligand_containment_mesh", [])
                distances = safe_distances(res.get("ligand_mesh_distances", []))
            elif VALIDATION_METHOD == "sasa":
                containment = res.get("ligand_containment_strict", [])
                distances = safe_distances(res.get("ligand_to_pocket_distances", []))
            else:
                containment = []
                distances = []

            num_validated = safe_sum_bool(containment)
            contact = "Yes" if safe_any(containment) else "No"
            num_pockets = res.get("num_pockets", len(res.get("pocket_volumes", [])))
            min_distance = round(min(distances), 2) if distances else "N/A"
            avg_volume = (
                round(np.mean(res.get("pocket_volumes", [])), 2)
                if res.get("pocket_volumes")
                else "N/A"
            )
            avg_depth = (
                round(np.mean(res.get("pocket_depths", [])), 2)
                if res.get("pocket_depths")
                else "N/A"
            )

            table.append(
                [
                    filename,
                    top1_val,
                    top3_val,
                    top5_val,
                    top_all,
                    rank,
                    validation_used,
                    num_pockets,
                    num_validated,
                    contact,
                    min_distance,
                    avg_volume,
                    avg_depth,
                ]
            )

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

    if csv_output_file:
        df = pd.DataFrame(table, columns=header)
        df.to_csv(csv_output_file, index=False)
        logging.info(f"Summary table saved to CSV: {csv_output_file}")

    return summary


def save_pocket_results_csv(result, pdb_id, output_dir="results"):
    output_csv = os.path.join(output_dir, f"{pdb_id}_pockets.csv")
    os.makedirs(output_dir, exist_ok=True)

    headers = [
        "pocket_id",
        "volume",
        "depth",
        "mouth_area",
        "mouth_circumference",
        "ligand_mouth_distance",
        "ligand_in_mesh",
    ]

    rows = []
    num_pockets = len(result.get("volumes", []))
    for i in range(num_pockets):
        rows.append(
            {
                "pocket_id": i + 1,
                "volume": result.get("volumes", [None])[i],
                "depth": result.get("depths", [None])[i],
                "mouth_area": result.get("mouth_area", [None])[i],
                "mouth_circumference": result.get("mouth_perimeter", [None])[i],
                "ligand_mouth_distance": result.get(
                    "ligand_mouth_min_distance", [None]
                )[i],
                "ligand_in_mesh": result.get("ligand_containment_mesh", [None])[i],
            }
        )

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)
    print(f"✅ Pocket results saved to {output_csv}")
