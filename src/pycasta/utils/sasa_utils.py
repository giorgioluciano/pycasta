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
from ..config import SASA_METHOD


try:
    import freesasa

    RDKit_AVAILABLE = True
except ImportError:
    RDKit_AVAILABLE = False
    logging.warning("FreeSASA not found.")


def compute_sasa(pdbpath, method=None, remove_ligand=False):
    structure = freesasa.Structure(pdbpath)
    result = freesasa.calc(structure)
    sasa_values = {}
    for i in range(structure.nAtoms()):
        resi = structure.residueNumber(i)
        area = result.atomArea(i)
        sasa_values[resi] = sasa_values.get(resi, 0) + area
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
