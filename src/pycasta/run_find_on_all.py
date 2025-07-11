#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_find_on_all.py

Script to run a batch analysis on all saved pocket PDB files.
It calls find_closest_pocket.py for each subdirectory in the pocket_pdbs folder.
"""

import os
import subprocess

base_dir = (
    r"C:\Users\giorgio\OneDrive - CNR\AAAProj\Molly\pycast\results\single\pocket_pdbs"
)

for pdb_id in os.listdir(base_dir):
    subdir = os.path.join(base_dir, pdb_id)
    if os.path.isdir(subdir):
        print(f"\n🔎 Running find_closest_pocket.py for {pdb_id}...")
        subprocess.run(["python", "find_closest_pocket.py", "--pdb", pdb_id])
