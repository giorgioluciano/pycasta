import glob
import os
import traceback

from pycasta.config import BOUNDED_DIR, UNBOUNDED_DIR
from pycasta.run_analysis import process_pdb


def pick_one_pdb(folder):
    if folder is None:
        return None

    folder = str(folder)
    if not os.path.exists(folder):
        return None

    pdb_files = sorted(glob.glob(os.path.join(folder, "*.pdb")))
    return pdb_files[0] if pdb_files else None


def run_case(label, pdb_path):
    if not pdb_path:
        print(f"SKIP {label}: no .pdb file found")
        return True

    print(f"\n=== {label} ===")
    print(f"Input: {pdb_path}")

    try:
        result = process_pdb(pdb_path)

        ranked_pockets = result.get("ranked_pockets", []) or []
        pocket_volumes = result.get("pocket_volumes", []) or []
        step_to_ligand = result.get("step_to_ligand")
        triangulation_backend = result.get("triangulation_backend")

        top_volume = pocket_volumes[0] if pocket_volumes else None

        print("STATUS: OK")
        print(f"triangulation_backend = {triangulation_backend}")
        print(f"num_pockets = {len(ranked_pockets)}")
        print(f"top_volume = {top_volume}")
        print(f"step_to_ligand = {step_to_ligand}")

        return True

    except Exception:
        print("STATUS: FAIL")
        traceback.print_exc()
        return False


def main():
    print(f"BOUNDED_DIR   = {BOUNDED_DIR}")
    print(f"UNBOUNDED_DIR = {UNBOUNDED_DIR}")

    bounded_pdb = pick_one_pdb(BOUNDED_DIR)
    unbounded_pdb = pick_one_pdb(UNBOUNDED_DIR)

    ok_bounded = run_case("BOUNDED", bounded_pdb)
    ok_unbounded = run_case("UNBOUNDED", unbounded_pdb)

    print("\n=== SANITY TEST ===")
    if ok_bounded and ok_unbounded:
        print("PASS")
        return 0

    print("FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())