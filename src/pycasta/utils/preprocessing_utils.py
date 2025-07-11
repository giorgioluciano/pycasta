#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
preprocessing_utils.py

Module for loading and preprocessing PDB files:
- Separating protein and ligand coordinates.
- Calculating atomic radii.
- Extracting connectivity (using CONECT records or a distance threshold).
"""

import logging
import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb
from config import INCLUDE_WATER, CONNECTIVITY_THRESHOLD
from atomic_radii import ATOMIC_RADII_DICT, DEFAULT_RADIUS
from scipy.spatial import KDTree


def load_and_separate_pdb(pdb_path, include_water=INCLUDE_WATER):
    ppdb = PandasPdb().read_pdb(pdb_path)

    # --- Protein coordinates ---
    protein_df = ppdb.df["ATOM"]
    protein_coords = protein_df[["x_coord", "y_coord", "z_coord"]].values
    logging.info(f"Protein coordinates loaded: {protein_coords.shape[0]} atoms.")

    # --- Ligand coordinates ---
    ligand_df = ppdb.df.get("HETATM")
    if ligand_df is not None and not ligand_df.empty:
        col_name = "residue_name" if "residue_name" in ligand_df.columns else "blank_1"
        if not include_water:
            ligand_df = ligand_df[ligand_df[col_name] != "HOH"]
        ligand_coords = (
            ligand_df[["x_coord", "y_coord", "z_coord"]].values
            if not ligand_df.empty
            else np.array([])
        )
    else:
        ligand_coords = np.array([])
    logging.info(f"Ligand coordinates loaded: {ligand_coords.shape[0]} atoms.")

    # --- Build list of protein atoms with metadata ---
    protein_atoms = []
    for _, row in protein_df.iterrows():
        protein_atoms.append(
            {
                "resid": int(row["residue_number"]),
                "atom_name": row["atom_name"].strip(),
                "coord": [row["x_coord"], row["y_coord"], row["z_coord"]],
            }
        )

    return protein_coords, ligand_coords, protein_atoms, ppdb


def calculate_atomic_radii(ppdb):
    try:
        protein_df = ppdb.df.get("ATOM")
        if protein_df is None or protein_df.empty:
            raise ValueError("No ATOM data available to compute radii.")
        atom_elements = protein_df["element_symbol"].values
        radii = np.array(
            [ATOMIC_RADII_DICT.get(elem, DEFAULT_RADIUS) for elem in atom_elements]
        )
        logging.info(f"Calculated atomic radii for {len(protein_df)} atoms.")
    except Exception as e:
        logging.error(f"Error calculating atomic radii: {e}")
        raise
    return radii


def extract_connectivity(ppdb, distance_threshold=CONNECTIVITY_THRESHOLD):
    connectivity = {}
    if "CONECT" in ppdb.df and not ppdb.df["CONECT"].empty:
        logging.info("Using CONECT records to extract connectivity.")
        conect_df = ppdb.df["CONECT"]
        for idx, row in conect_df.iterrows():
            source = int(row[0]) - 1  # 0-based index
            neighbors = [int(x) - 1 for x in row[1:] if not pd.isna(x)]
            connectivity[source] = neighbors
        num_atoms = len(ppdb.df["ATOM"])
        for i in range(num_atoms):
            if i not in connectivity:
                connectivity[i] = []
    else:
        logging.info(
            "No CONECT records found; computing connectivity based on distance threshold."
        )
        protein_df = ppdb.df.get("ATOM")
        if protein_df is None or protein_df.empty:
            raise ValueError("No ATOM data available to compute connectivity.")
        coords = protein_df[["x_coord", "y_coord", "z_coord"]].to_numpy()
        tree = KDTree(coords)
        num_atoms = coords.shape[0]
        for i in range(num_atoms):
            indices = tree.query_ball_point(coords[i], r=distance_threshold)
            indices = [j for j in indices if j != i]
            connectivity[i] = indices
    return connectivity
