"""
stress_test.py

Stress test per pycasta su modalità, dataset e configurazioni.
Salva output/stress_test_results.csv.
"""

import glob
import itertools
import os
import time
import traceback
from contextlib import contextmanager
from copy import deepcopy

import pandas as pd

from pycasta import config as cfg
from pycasta.run_analysis import process_pdb


@contextmanager
def temp_config(**kwargs):
    old = {}
    for k, v in kwargs.items():
        old[k] = getattr(cfg, k)
        setattr(cfg, k, v)
    try:
        yield
    finally:
        for k, v in old.items():
            setattr(cfg, k, v)


def list_pdbs(folder):
    if folder is None or not os.path.exists(str(folder)):
        return []
    return sorted(glob.glob(os.path.join(str(folder), "*.pdb")))


def run_one(pdb_path):
    t0 = time.perf_counter()
    result = process_pdb(pdb_path)
    dt = time.perf_counter() - t0
    pockets = result.get("ranked_pockets", [])
    vols = result.get("pocket_volumes", [])
    return {
        "pdb": os.path.basename(pdb_path),
        "status": "OK",
        "pockets": len(pockets),
        "top_volume": vols[0] if vols else None,
        "elapsed_s": round(dt, 2),
    }


def paired_like_check():
    corr = str(cfg.CORRESPONDENCE_FILE)
    bounded = list_pdbs(cfg.BOUNDED_DIR)
    unbounded = list_pdbs(cfg.UNBOUNDED_DIR)

    rows = []
    if not os.path.exists(corr):
        return rows

    corr_df = pd.read_excel(corr)
    for _, row in corr_df.iterrows():
        bname = str(row["bounded"]).strip() if "bounded" in corr_df.columns else str(row.iloc[0]).strip()
        uname = str(row["unbounded"]).strip() if "unbounded" in corr_df.columns else str(row.iloc[1]).strip()

        bpath = os.path.join(str(cfg.BOUNDED_DIR), bname)
        upath = os.path.join(str(cfg.UNBOUNDED_DIR), uname)

        if not os.path.exists(bpath) or not os.path.exists(upath):
            rows.append({
                "mode": "paired",
                "pair": f"{bname} / {uname}",
                "status": "SKIP",
                "detail": "missing file",
                "pockets": None,
                "top_volume": None,
                "elapsed_s": None,
            })
            continue

        try:
            br = run_one(bpath)
            ur = run_one(upath)
            rows.append({
                "mode": "paired",
                "pair": f"{bname} / {uname}",
                "status": "OK",
                "detail": f"B:{br['pockets']} U:{ur['pockets']}",
                "pockets": br["pockets"],
                "top_volume": br["top_volume"],
                "elapsed_s": None,
            })
        except Exception:
            rows.append({
                "mode": "paired",
                "pair": f"{bname} / {uname}",
                "status": "FAIL",
                "detail": traceback.format_exc(limit=1),
                "pockets": None,
                "top_volume": None,
                "elapsed_s": None,
            })
    return rows


def main():
    bounded = list_pdbs(cfg.BOUNDED_DIR)
    unbounded = list_pdbs(cfg.UNBOUNDED_DIR)

    bounded = bounded[:3]
    unbounded = unbounded[:3]

    alpha_vals = [cfg.MANUAL_ALPHA, 2.2]
    merge_vals = [cfg.MERGE_THRESHOLD, 12, 18]
    sasa_vals = [False, True] if hasattr(cfg, "USE_SASA") else [False]
    filter_alpha_vals = [False, True] if hasattr(cfg, "FILTER_ALPHA_BY_SASA") else [False]

    combos = list(itertools.product(alpha_vals, merge_vals, sasa_vals, filter_alpha_vals))

    rows = []
    print("=== STRESS TEST START ===")
    print(f"bounded files: {len(bounded)}")
    print(f"unbounded files: {len(unbounded)}")
    print(f"correspondence: {cfg.CORRESPONDENCE_FILE}")

    for alpha, merge_th, use_sasa, filter_alpha in combos:
        print(f"\n--- CONFIG alpha={alpha} merge={merge_th} use_sasa={use_sasa} filter_alpha={filter_alpha} ---")
        with temp_config(
            MANUAL_ALPHA=alpha,
            MERGE_THRESHOLD=merge_th,
            USE_SASA=use_sasa,
            FILTER_ALPHA_BY_SASA=filter_alpha,
        ):
            for pdb in bounded:
                try:
                    r = run_one(pdb)
                    r.update({
                        "mode": "single-bounded",
                        "alpha": alpha,
                        "merge_threshold": merge_th,
                        "use_sasa": use_sasa,
                        "filter_alpha": filter_alpha,
                    })
                    rows.append(r)
                    print(f"OK bounded | {r['pdb']} | pockets={r['pockets']} | time={r['elapsed_s']}s")
                except Exception:
                    rows.append({
                        "mode": "single-bounded",
                        "pdb": os.path.basename(pdb),
                        "status": "FAIL",
                        "pockets": None,
                        "top_volume": None,
                        "elapsed_s": None,
                        "alpha": alpha,
                        "merge_threshold": merge_th,
                        "use_sasa": use_sasa,
                        "filter_alpha": filter_alpha,
                        "detail": traceback.format_exc(limit=1),
                    })
                    print(f"FAIL bounded | {os.path.basename(pdb)}")

            for pdb in unbounded:
                try:
                    r = run_one(pdb)
                    r.update({
                        "mode": "single-unbounded",
                        "alpha": alpha,
                        "merge_threshold": merge_th,
                        "use_sasa": use_sasa,
                        "filter_alpha": filter_alpha,
                    })
                    rows.append(r)
                    print(f"OK unbounded | {r['pdb']} | pockets={r['pockets']} | time={r['elapsed_s']}s")
                except Exception:
                    rows.append({
                        "mode": "single-unbounded",
                        "pdb": os.path.basename(pdb),
                        "status": "FAIL",
                        "pockets": None,
                        "top_volume": None,
                        "elapsed_s": None,
                        "alpha": alpha,
                        "merge_threshold": merge_th,
                        "use_sasa": use_sasa,
                        "filter_alpha": filter_alpha,
                        "detail": traceback.format_exc(limit=1),
                    })
                    print(f"FAIL unbounded | {os.path.basename(pdb)}")

    rows.extend(paired_like_check())

    df = pd.DataFrame(rows)
    out_csv = os.path.join(str(cfg.OUTPUT_DIR), "stress_test_results.csv")
    os.makedirs(str(cfg.OUTPUT_DIR), exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"\nSaved CSV: {out_csv}")
    print("\n=== SUMMARY ===")
    print(df.head(30).to_string(index=False))
    fails = int((df["status"] == "FAIL").sum()) if "status" in df.columns else 0
    print(f"\nFailures: {fails}")
    return 0 if fails == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())