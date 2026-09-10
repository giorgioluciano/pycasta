#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_viz_alpha_only_no_json.py - PyCASTA alpha-only visualization without JSON.
"""
import sys, argparse, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import MolData, get_vdw, guess_alpha_values
from alpha_viewer import show_alpha


def load_protein_atoms_from_pdb(pdb_path):
    coords = []
    atom_names = []
    with open(pdb_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])
                atom_names.append(atom_name)
    coords = np.array(coords, dtype=float)
    radii = np.array([get_vdw(name) for name in atom_names], dtype=float)
    ligand_coords = np.empty((0, 3), dtype=float)
    representative_points = np.empty((0, 3), dtype=float)
    return coords, atom_names, radii, ligand_coords, representative_points


def load_alpha_only(pdb_path, alpha_npz_path, alpha_value=None):
    coords, atom_names, radii, ligand_coords, representative_points = load_protein_atoms_from_pdb(pdb_path)
    with np.load(alpha_npz_path, allow_pickle=True) as alpha_data:
        tetrahedra = np.array(alpha_data["tetrahedra"], dtype=float)
        alpha_mask = np.array(alpha_data["alpha_mask"], dtype=bool)

    mol = MolData(
        name=os.path.splitext(os.path.basename(pdb_path))[0],
        pdb_path=str(pdb_path),
        coords=coords,
        radii=radii,
        atom_names=atom_names,
        ligand_coords=ligand_coords,
        tetrahedra=tetrahedra,
        alpha_mask=alpha_mask,
        ranked_pockets=[],
        representative_points=representative_points,
        pocket_volumes=[],
        ligand_validation=[],
        ranking_scores=[],
        alpha_used=float(alpha_value) if alpha_value is not None else None,
        ao_scores=None,
        alpha_npz_path=str(alpha_npz_path),
    )
    return mol


def main():
    p = argparse.ArgumentParser(description="PyCASTA alpha-only visualization without JSON")
    p.add_argument("pdb")
    p.add_argument("alpha_npz")
    p.add_argument("--alpha", type=float, default=None)
    args = p.parse_args()

    mol = load_alpha_only(args.pdb, args.alpha_npz, args.alpha)
    av = guess_alpha_values(mol)
    if args.alpha is not None:
        av["custom"] = args.alpha

    print(f"\nMolecola  : {mol.name}")
    print(f"Atomi     : {len(mol.coords)}")
    print(f"Tetraedri : {len(mol.tetrahedra)}  (inside={mol.alpha_mask.sum()}, outside={(~mol.alpha_mask).sum()})")
    print(f"Alpha     : tight={av['tight']:.2f}  medium={av['medium']:.2f}  loose={av['loose']:.2f}")

    show_alpha(mol, av)


if __name__ == "__main__":
    main()
