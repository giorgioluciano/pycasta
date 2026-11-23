"""Bridge PyCASTa PDB preprocessing to the exact PocketLab backend and alpha CLI."""
import argparse
import hashlib
import json
import math
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdb", type=Path)
    parser.add_argument("--out", type=Path, default=Path("output_isolated_tests"))
    parser.add_argument("--samples", type=int, default=128)
    parser.add_argument("--alpha-squared", type=float, default=0.0)
    parser.add_argument("--reference-manifest", type=Path)
    parser.add_argument("--reference-topology", type=Path)
    a = parser.parse_args()
    alpha_cli = Path(__file__).with_name("cast_alpha.py")
    if not alpha_cli.is_file():
        parser.error("Place this runner next to cast_alpha.py")
    if not a.pdb.is_file() or not math.isfinite(a.alpha_squared) or a.samples < 0 or a.samples == 1:
        parser.error("Need an existing PDB, finite alpha_squared and samples=0 or >=2")
    from pycasta.utils.preprocessing_utils import load_and_separate_pdb, calculate_atomic_radii
    import pocketlab.core as core
    import numpy as np
    coord, ligand, _, ppdb = load_and_separate_pdb(str(a.pdb.resolve()))
    coord = np.asarray(coord, dtype=float)
    radii = np.asarray(calculate_atomic_radii(ppdb), dtype=float)
    if radii.ndim != 1 or coord.shape != (len(radii), 3) or len(radii) < 4:
        parser.error("Invalid atom array dimensions")
    if not np.isfinite(coord).all() or not np.isfinite(radii).all() or (radii <= 0).any():
        parser.error("Nonfinite coordinates or invalid radii")
    digest = hashlib.sha256(a.pdb.read_bytes()).hexdigest()
    if a.reference_manifest:
        ref = json.loads(a.reference_manifest.read_text(encoding="utf-8-sig"))
        if (ref.get("input_sha256") != digest or ref.get("alpha_squared") != a.alpha_squared
                or not np.array_equal(coord, ref.get("coordinates"))
                or not np.array_equal(radii, ref.get("atomic_radii"))):
            parser.error("PDB hash, alpha, coordinates/order or radii differ from reference manifest")
        print("Reference model: exact match", flush=True)
    g = core.build_cast_complex_exact(coord.tolist(), radii.tolist(), a.alpha_squared)
    if not isinstance(g, dict) or not g.get("ok") or not g.get("alpha_all_dimensions_exported"):
        parser.error("Exact backend did not produce a successful all-dimensional complex; no fallback")
    dest = a.out / ("cast_pdb_" + a.pdb.stem + "_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f"))
    dest.mkdir(parents=True, exist_ok=False)
    manifest = dict(input=str(a.pdb.resolve()), input_sha256=digest, backend=core.__file__,
                    backend_sha256=hashlib.sha256(Path(core.__file__).read_bytes()).hexdigest(),
                    alpha_squared=a.alpha_squared, n_atoms=len(radii),
                    n_ligand_atoms_loaded_not_used_for_detection=len(ligand),
                    model_note="calculate_atomic_radii unchanged; no additive solvent inflation",
                    coordinates=coord.tolist(), atomic_radii=radii.tolist())
    gp, mp = dest / "cast_complex_all_dims.json", dest / "manifest.json"
    for path, data in ((gp, g), (mp, manifest)):
        with path.open("x", encoding="utf-8") as f:
            json.dump(data, f, indent=2, allow_nan=False)
    cmd = [sys.executable, str(alpha_cli), "--complex", str(gp), "--manifest", str(mp),
           "--samples", str(a.samples), "--out", str(dest / "alpha_results")]
    if a.reference_topology:
        cmd.extend(["--reference-topology", str(a.reference_topology)])
    print("CGAL EXPORT:", dest.resolve(), flush=True)
    return subprocess.call(cmd)


if __name__ == "__main__":
    raise SystemExit(main())
