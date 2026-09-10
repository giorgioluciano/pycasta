from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from radial_closure_test_k import analyze_pdb


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_files", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, default=Path("results/radial_closure_k_sensitivity.csv"))
    parser.add_argument("--k-values", nargs="+", type=int, default=[5, 10, 20])
    parser.add_argument("--rays", type=int, default=16)
    parser.add_argument("--max-distance", type=float, default=12.0)
    parser.add_argument("--neighbor-radius", type=float, default=None)
    args = parser.parse_args()

    rows = []
    for pdb_file in args.pdb_files:
        for k in args.k_values:
            for row in analyze_pdb(
                path=pdb_file,
                k_nearest=k,
                n_rays=args.rays,
                max_distance=args.max_distance,
                neighbor_radius=args.neighbor_radius,
            ):
                row["k_nearest"] = k
                rows.append(row)

    if not rows:
        raise SystemExit("Nessun gruppo HETATM trovato.")

    df = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"Salvate {len(df)} righe in {args.output}")

    print("\nMedie per k:")
    print(
        df.groupby("k_nearest")[[
            "radial_closure",
            "hit_fraction",
            "distance_score",
            "mean_hit_distance",
        ]].mean().round(3)
    )


if __name__ == "__main__":
    main()
