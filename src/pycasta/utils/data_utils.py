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
