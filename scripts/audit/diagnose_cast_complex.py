"""Isolated exact CAST-complex diagnostic; no legacy pipeline/ranking."""
import argparse
from collections import Counter
from datetime import datetime
import hashlib
import importlib
import json
from pathlib import Path

from pycasta.utils.preprocessing_utils import (
    load_and_separate_pdb,
    calculate_atomic_radii,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", default=r"F:\datasets\protein_pockets\pocketpicker\holo\1a6w.pdb")
    parser.add_argument("--alpha-squared", type=float, default=0.0)
    parser.add_argument("--output-root", default=r"C:\github\pycasta\output_isolated_tests")
    args = parser.parse_args()
    source = Path(args.pdb).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    core = importlib.import_module("pocketlab.core")
    coords, ligand, atoms, ppdb = load_and_separate_pdb(str(source))
    radii = calculate_atomic_radii(ppdb)
    result = core.build_cast_complex_exact(coords, radii, args.alpha_squared)
    dest = Path(args.output_root).resolve() / (
        "cast_exact_" + source.stem + "_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    )
    dest.mkdir(parents=True, exist_ok=False)
    manifest = {
        "input": str(source),
        "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "backend": str(core.__file__),
        "alpha_squared": args.alpha_squared,
        "model_note": "Base balls use calculate_atomic_radii output unchanged; no additive solvent inflation.",
        "n_atoms": len(coords),
        "n_ligand_atoms_loaded_not_used_for_detection": len(ligand),
        "coordinates": coords.tolist() if hasattr(coords, "tolist") else coords,
        "atomic_radii": radii.tolist() if hasattr(radii, "tolist") else radii,
    }
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (dest / "cast_complex.json").write_text(json.dumps(result), encoding="utf-8")
    report = {
        "status": result["status"],
        "alpha_squared": args.alpha_squared,
        "n_cells": len(result["cells_idx"]),
        "n_nonvertex_inputs": len(result["nonvertex_input_ids"]),
        "zero_face_predicates": len(result["zero_face_predicates"]),
        "equal_power_flow_faces": len(result["equal_power_flow_faces"]),
    }
    if result["ok"]:
        groups = result["pocket_tetra_components"]
        flat = [i for group in groups for i in group]
        outside = set(result["outside_ancestors"])
        expected = {
            i for i in range(len(result["cells_idx"]))
            if i not in outside
            and not result["tetra_in_alpha"][i]
        }
        if len(flat) != len(set(flat)) or set(flat) != expected:
            raise RuntimeError("Components do not partition eligible tetrahedra")
        unique_faces = {}
        for i, cell in enumerate(result["cells_idx"]):
            for k in range(4):
                key = tuple(sorted(cell[j] for j in range(4) if j != k))
                state = result["face_alpha_class"][i][k]
                if key in unique_faces and unique_faces[key] != state:
                    raise RuntimeError("Inconsistent shared face classification")
                unique_faces[key] = state
        report.update({
            "tetra_in_alpha": sum(result["tetra_in_alpha"]),
            "outside_ancestors": len(result["outside_ancestors"]),
            "finite_sinks": len(result["finite_sinks"]),
            "eligible_tetrahedra": result["eligible_tetrahedra"],
            "n_tetra_components": len(groups),
            "largest_component_sizes": sorted(map(len, groups), reverse=True)[:20],
            "unique_face_classes": dict(Counter(
                ["EXTERIOR", "SINGULAR", "REGULAR", "INTERIOR"][v]
                for v in unique_faces.values()
            )),
            "partition_check": "passed",
        })
    (dest / "summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print("OUTPUT:", dest)
    if not result["ok"]:
        raise SystemExit("Degenerate case saved; no pocket classification accepted.")


if __name__ == "__main__":
    main()
