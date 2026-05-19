"""
test_batch.py
Esegue un mini-batch su 3 PDB per controllare stabilità del flusso.
"""

import glob
import os
import traceback

from pycasta.run_analysis import process_pdb
from pycasta.config import BOUNDED_DIR

def main():
    pdbs = sorted(glob.glob(os.path.join(str(BOUNDED_DIR), "*.pdb")))
    if not pdbs:
        print(f"Nessun PDB trovato in {BOUNDED_DIR}")
        return 1

    sample = pdbs[:3]
    print("Batch test sui file:")
    for p in sample:
        print(" -", p)

    rows = []
    for pdb in sample:
        print(f"\n=== RUN: {os.path.basename(pdb)} ===")
        try:
            result = process_pdb(pdb)
            pockets = result.get("ranked_pockets", [])
            volumes = result.get("pocket_volumes", [])
            top_volume = volumes[0] if volumes else None
            rows.append({
                "file": os.path.basename(pdb),
                "pockets": len(pockets),
                "top_volume": top_volume,
                "status": "OK",
            })
            print(f"OK | pockets={len(pockets)} | top_volume={top_volume}")
        except Exception:
            rows.append({
                "file": os.path.basename(pdb),
                "pockets": None,
                "top_volume": None,
                "status": "FAIL",
            })
            print("FAIL")
            traceback.print_exc()

    print("\n=== SUMMARY ===")
    for r in rows:
        print(r)

    fails = sum(1 for r in rows if r["status"] != "OK")
    return 0 if fails == 0 else 1

if __name__ == "__main__":
    raise SystemExit(main())