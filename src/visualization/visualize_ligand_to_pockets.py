#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
visualize_ligand_to_pockets.py

Module for visualizing pocket centroids and the ligand.
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from config import OUTPUT_DIR, VERSION_TAG

def visualize_ligand_and_pockets(pdb_id, max_pockets=5):
    json_path = os.path.join(OUTPUT_DIR, VERSION_TAG, f"{pdb_id}.json")
    if not os.path.exists(json_path):
        raise FileNotFoundError(f"Result JSON not found for {pdb_id}")
    with open(json_path, "r") as f:
        data = json.load(f)
    ligand = np.array(data.get("ligand_coords", []))
    pocket_dists = data.get("ligand_to_pocket_distances", [])
    rep_points = data.get("representative_points", [])
    if ligand.size == 0 or not rep_points:
        print("❌ No ligand or pocket data found.")
        return
    rep_points = np.array(rep_points)
    distances = np.array(pocket_dists)
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(ligand[:, 0], ligand[:, 1], ligand[:, 2], c='red', s=20, label='Ligand')
    n = min(max_pockets, len(rep_points))
    ax.scatter(rep_points[:n, 0], rep_points[:n, 1], rep_points[:n, 2], c='blue', s=40, label='Pocket centroids')
    ligand_com = np.mean(ligand, axis=0)
    for i in range(n):
        dists = np.linalg.norm(ligand - rep_points[i], axis=1)
        closest_ligand_point = ligand[np.argmin(dists)]
        ax.plot([closest_ligand_point[0], rep_points[i, 0]],
                [closest_ligand_point[1], rep_points[i, 1]],
                [closest_ligand_point[2], rep_points[i, 2]], 'gray', linestyle='--', linewidth=1)
        dist = distances[i] if i < len(distances) else "N/A"
        ax.text(rep_points[i, 0], rep_points[i, 1], rep_points[i, 2], f"{i+1} ({dist:.1f}Å)", fontsize=9)
    ax.set_title(f"Pocket Centroids vs Ligand for {pdb_id}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend()
    plt.tight_layout()
    plt.show()
    print(f"\n🔎 Ligand-to-pocket distances for top {n} pockets:")
    for i in range(n):
        d = distances[i] if i < len(distances) else "N/A"
        print(f"Pocket {i+1}: {d:.2f} Å" if isinstance(d, float) else f"Pocket {i+1}: N/A")

if __name__ == "__main__":
    visualize_ligand_and_pockets("1fbp")
