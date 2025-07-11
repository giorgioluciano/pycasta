# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 15:12:38 2025

@author: giorg
"""

import os
import subprocess


pdb_file = "data/bounded/1acj.pdb"
ply_file = "results/bov1/pocket_pdbs/1acj/1acj_ranked_1.pdb"
script = f"""
load {pdb_file}, protein
load {ply_file}, pocket
show surface, protein
color slate, protein
show mesh, pocket
color orange, pocket
orient
"""

with open("tmp_load_script.pml", "w") as f:
    f.write(script)

# Percorso completo al batch PyMOL trovato
PYMOL_PATH = r"C:\Users\giorg\miniconda3\envs\debug\Scripts\pymol.BAT"

# Lancia PyMOL interattivo con lo script
subprocess.run([PYMOL_PATH, "-r", "tmp_load_script.pml"])
