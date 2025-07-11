#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
sasa_utils.py

Module for calculating SASA and evaluating ligand-pocket contact.
"""

import os
import logging
import numpy as np
from scipy.spatial.distance import cdist
from config import SASA_METHOD, USE_PYMOL2, SASA_CONTACT_THRESHOLD

try:
    import pymol

    pymol.finish_launching(["-qc"])
    PYMOL_AVAILABLE = True
except ImportError:
    PYMOL_AVAILABLE = False
    logging.warning("PyMOL not found. SASA will not be computed using PyMOL.")

try:
    import pymol2

    PYMOL2_AVAILABLE = True
except ImportError:
    PYMOL2_AVAILABLE = False
    logging.warning("pymol2 not found.")

try:
    import freesasa

    RDKit_AVAILABLE = True
except ImportError:
    RDKit_AVAILABLE = False
    logging.warning("FreeSASA not found.")


def compute_sasa(pdb_path, method=None, remove_ligand=False):
    if method is None:
        method = SASA_METHOD
    if method.lower() == "pymol2":
        return compute_sasa_pymol2(pdb_path, remove_ligand=remove_ligand)
    else:
        return compute_sasa_pymol(pdb_path, remove_ligand=remove_ligand)


def compute_sasa_pymol2(pdb_path, remove_ligand=False):
    if not PYMOL2_AVAILABLE:
        logging.error("pymol2 is not available for SASA computation.")
        return None
    with pymol2.PyMOL() as pymol_instance:
        cmd = pymol_instance.cmd
        cmd.reinitialize()
        logging.info(f"Loading {pdb_path} in pymol2 for SASA computation...")
        cmd.load(pdb_path, "protein")
        if remove_ligand:
            cmd.remove("hetatm")
            logging.info("Removed hetatm (ligands) for SASA computation.")
        cmd.set("dot_solvent", 1)
        cmd.get_area("protein", load_b=1)
        model = cmd.get_model("protein")
        sasa_values = {str(atom.resi): atom.b for atom in model.atom}
        logging.info("SASA computation using pymol2 completed.")
        return sasa_values


def compute_sasa_pymol(pdb_path, remove_ligand=False):
    from utils.common import Fore, Style

    pdb_file = os.path.abspath(pdb_path)
    if not os.path.exists(pdb_file):
        logging.error(f"File not found: {pdb_file}")
        return None
    if USE_PYMOL2:
        try:
            import pymol2

            p = pymol2.PyMOL()
            p.start()
            cmd = p.cmd
        except Exception as e:
            logging.error(
                Fore.RED + Style.BRIGHT + f"[ERROR] Failed to start PyMOL2: {e}"
            )
            return None
    else:
        try:
            from pymol import cmd

            cmd.reinitialize()
        except Exception as e:
            logging.error(
                Fore.RED + Style.BRIGHT + f"[ERROR] Failed to initialize PyMOL: {e}"
            )
            return None
    logging.info("Loading " + pdb_file + " in pymol for SASA computation...")
    try:
        cmd.load(pdb_file, "protein")
    except Exception as e:
        logging.error(
            Fore.RED
            + Style.BRIGHT
            + f"[ERROR] Failed to load PDB file '{pdb_file}': {e}"
        )
        return None
    if remove_ligand:
        cmd.remove("hetatm")
        logging.info("Ligands removed for SASA computation.")
    cmd.set("dot_solvent", 1)
    cmd.get_area("protein", load_b=1)
    sasa_values = {resi.resi: resi.b for resi in cmd.get_model("protein").atom}
    cmd.delete("protein")
    logging.info("SASA computation completed.")
    return sasa_values


def compute_ligand_contact_sasa(sasa_before, sasa_after):
    if not sasa_before or not sasa_after:
        logging.warning("Incomplete SASA data; cannot compute ligand contact.")
        return {}
    contacting = {}
    for atom_idx, area_before in sasa_before.items():
        area_after = sasa_after.get(atom_idx, 0)
        if area_before > area_after:
            contacting[atom_idx] = area_before
    return contacting


def evaluate_sasa_contact(
    pocket_atom_indices,
    sasa_before,
    sasa_after,
    ligand_coords,
    protein_coords,
    threshold=SASA_CONTACT_THRESHOLD,
):
    contact_atoms = []
    for idx in pocket_atom_indices:
        res_id = str(idx)
        sasa_b = sasa_before.get(res_id, 0)
        sasa_a = sasa_after.get(res_id, 0)
        if sasa_b > sasa_a:
            contact_atoms.append(protein_coords[idx])
    if not contact_atoms or ligand_coords.size == 0:
        return False
    distances = cdist(np.array(contact_atoms), ligand_coords)
    min_dist = np.min(distances)
    return min_dist <= threshold
