# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 15:17:35 2025

@author: giorg
"""

import os
import shutil
import subprocess

# Possibili nomi ed estensioni dell'eseguibile PyMOL
pymol_names = [
    "pymol.exe",  # Installazione classica
    "PyMOLWin.exe",  # Vecchie versioni di Windows
    "pymol-open-source.exe",  # Alcuni bundle conda
]

# Possibili directory in cui PyMOL potrebbe trovarsi
possible_dirs = [
    r"C:\Program Files\Schrodinger\PyMOL",
    r"C:\Program Files\PyMOL",
    r"C:\Program Files (x86)\PyMOL",
    r"C:\Program Files (x86)\Schrodinger\PyMOL",
    r"C:\Users\%USERNAME%\AppData\Local\Programs\PyMOL",
    r"C:\Users\%USERNAME%\AppData\Local\PyMOL",
]

# Cerca nei PATH di sistema
pymol_in_path = shutil.which("pymol")
if pymol_in_path:
    print(f"Found pymol in system PATH: {pymol_in_path}")
    pymol_path = pymol_in_path
else:
    # Prova le directory comuni
    found = False
    for dir_ in possible_dirs:
        dir_ = os.path.expandvars(dir_)
        for exe in pymol_names:
            exe_path = os.path.join(dir_, exe)
            if os.path.isfile(exe_path):
                print(f"Found pymol here: {exe_path}")
                pymol_path = exe_path
                found = True
                break
        if found:
            break
    else:
        print("PyMOL not found in common locations or system PATH.")
        pymol_path = None

# Esegui un comando di test se trovato
if pymol_path:
    print(f"Trying to launch PyMOL: {pymol_path}")
    # Lancia PyMOL senza file (solo per test)
    try:
        subprocess.run([pymol_path], check=True)
    except Exception as e:
        print(f"Could not launch PyMOL: {e}")
else:
    print("Please install PyMOL or add it to your system PATH.")
