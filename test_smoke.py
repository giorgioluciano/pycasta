"""
test_smoke.py

Smoke test minimale per pycasta.
Verifica che il package sia importabile e che process_pdb giri su un file di esempio.
"""

import sys
import traceback

def test_import():
    print("=== TEST 1: import package ===")
    import pycasta
    print(f"OK: pycasta importato da {pycasta.__file__}")

def test_config():
    print("\n=== TEST 2: config ===")
    from pycasta.config import (
        BOUNDED_DIR, OUTPUT_DIR, MANUAL_ALPHA,
        VERSION_TAG, MIN_POCKET_VOLUME
    )
    print(f"  BOUNDED_DIR   = {BOUNDED_DIR}")
    print(f"  OUTPUT_DIR    = {OUTPUT_DIR}")
    print(f"  MANUAL_ALPHA  = {MANUAL_ALPHA}")
    print(f"  VERSION_TAG   = {VERSION_TAG}")
    print(f"  MIN_VOLUME    = {MIN_POCKET_VOLUME}")
    print("OK: config caricato")

def test_process_pdb():
    print("\n=== TEST 3: process_pdb su 1 file ===")
    from pycasta.config import BOUNDED_DIR
    from pycasta.run_analysis import process_pdb
    import glob
    import os

    pdbs = glob.glob(os.path.join(str(BOUNDED_DIR), "*.pdb"))
    if not pdbs:
        print(f"SKIP: nessun .pdb trovato in {BOUNDED_DIR}")
        return

    pdb = pdbs[0]
    print(f"  File: {pdb}")
    result = process_pdb(pdb)
    pockets = result.get("ranked_pockets", [])
    print(f"  Pockets trovate: {len(pockets)}")
    if len(pockets) > 0:
        print(f"  Volume pocket 1: {result.get('pocket_volumes', [None])[0]}")
    print("OK: process_pdb eseguito")

if __name__ == "__main__":
    errori = 0
    for test in [test_import, test_config, test_process_pdb]:
        try:
            test()
        except Exception as e:
            print(f"\nFAIL: {test.__name__}")
            traceback.print_exc()
            errori += 1

    print(f"\n{'='*40}")
    if errori == 0:
        print("TUTTI I TEST PASSATI")
    else:
        print(f"{errori} TEST FALLITI")
    sys.exit(errori)