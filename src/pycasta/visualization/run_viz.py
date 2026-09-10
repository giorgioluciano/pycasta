#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_viz.py - Entry point PyCASTA Visualization

Uso:
    python run_viz.py <pdb> <json> [--mode alpha|ao|pocket|all]
                                   [--alpha FLOAT]   valore alpha custom
                                   [--rays  INT]     raggi AO (default 64)
                                   [--res   FLOAT]   risoluzione Ang (default 0.5)
                                   [--results-dir PATH]

Esempi:
    python run_viz.py 1a26A.pdb results/coach420_test/1a26A_pockets.json
    python run_viz.py 1a26A.pdb 1a26A_pockets.json --mode alpha
    python run_viz.py 1a26A.pdb 1a26A_pockets.json --alpha 1.8
    python run_viz.py 1a26A.pdb 1a26A_pockets.json --mode ao --rays 128
    python run_viz.py 1a26A.pdb 1a26A_pockets.json --results-dir results/coach420_test
"""
import sys, argparse, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from core import load_mol, guess_alpha_values
from triangulation_viewer import show_triangulation
from alpha_viewer import show_alpha
from surface_viewer import show_surface
from pocket_viewer import show_pockets


def main():
    p = argparse.ArgumentParser(description="PyCASTA Visualization")
    p.add_argument("pdb")
    p.add_argument("json")
    p.add_argument("--mode", default="all", choices=["tri", "alpha", "ao", "pocket", "all"])
    p.add_argument("--alpha", type=float, default=None)
    p.add_argument("--rays", type=int, default=64)
    p.add_argument("--res", type=float, default=0.5)
    p.add_argument("--results-dir", type=str, default=None)
    args = p.parse_args()

    mol = load_mol(args.pdb, args.json)
    av = guess_alpha_values(mol)
    if args.alpha:
        av["custom"] = args.alpha

    n_raw = (
        len(mol.triangulation_tetrahedra)
        if getattr(mol, "triangulation_tetrahedra", None) is not None
        else 0
    )

    print(f"\nMolecola   : {mol.name}")
    print(f"Atomi      : {len(mol.coords)}")
    print(f"Tetra raw  : {n_raw}")
    print(f"Tetra alpha: {len(mol.tetrahedra)}  (inside={mol.alpha_mask.sum()}, outside={(~mol.alpha_mask).sum()})")
    print(f"Pockets    : {len(mol.ranked_pockets)}  (validati={sum(mol.ligand_validation)})")
    print(f"Alpha      : tight={av['tight']:.2f}  medium={av['medium']:.2f}  loose={av['loose']:.2f}  Rg={av['_rg']:.1f}A")

    if args.mode in ("tri", "all"):
        show_triangulation(mol)
    if args.mode in ("alpha", "all"):
        show_alpha(mol, av)
    if args.mode in ("ao", "all"):
        show_surface(mol, N_rays=args.rays, resolution=args.res)
    if args.mode in ("pocket", "all"):
        show_pockets(mol, results_dir=args.results_dir)


if __name__ == "__main__":
    main()
