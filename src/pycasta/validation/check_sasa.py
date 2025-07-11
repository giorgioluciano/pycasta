# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 11:57:11 2025

@author: giorgio
"""

import logging
from utils.sasa_utils import compute_sasa
import pymol

pymol.finish_launching(["-qc"])
logging.basicConfig(level=logging.DEBUG, format="%(levelname)s: %(message)s")

sasa_vals = compute_sasa(
    "data/bounded_test/1acj.pdb", method="pymol", remove_ligand=False
)
logging.info(f"SASA values computed: {sasa_vals}")
