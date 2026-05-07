# -*- coding: utf-8 -*-
"""
adaptive_alpha.py
Funzioni per il calcolo e la diagnostica dell'alpha adattivo.
"""
import numpy as np
import os
import csv


def compute_alpha_gap(tetra_positions, simplices, radii, alpha_values=None):
    if alpha_values is None:
        alpha_values = np.linspace(0.5, 5.0, 50)
    if radii is None or len(radii) == 0:
        return float(np.median(alpha_values))
    radii_arr = np.array(radii)
    gaps = []
    for alpha in alpha_values:
        mask = radii_arr <= alpha
        gap = float(np.sum(mask)) / len(radii_arr)
        gaps.append(gap)
    gaps = np.array(gaps)
    diffs = np.diff(gaps)
    best_idx = int(np.argmax(diffs)) + 1
    return float(alpha_values[best_idx])


def get_adaptive_alpha(tetra_positions, simplices, radii, alpha_values=None):
    return compute_alpha_gap(tetra_positions, simplices, radii, alpha_values=alpha_values)


def save_alpha_diagnostics(alpha_values, gap_scores, best_alpha,
                            output_dir="results", molecule_name="molecule"):
    os.makedirs(output_dir, exist_ok=True)
    filename = os.path.join(output_dir, f"{molecule_name}_alpha_diagnostics.csv")
    with open(filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["alpha", "gap_score", "best_alpha"])
        for alpha, score in zip(alpha_values, gap_scores):
            writer.writerow([alpha, score, best_alpha])
    return filename