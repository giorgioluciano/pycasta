"""Run the exact CAST backend with verified per-atom ProtOr radii."""
import argparse
import csv
import hashlib
import io
import json
import math
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdb", type=Path)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    parser.add_argument("--radii-csv", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("output_isolated_tests/protor_1fkf"))
    parser.add_argument("--alpha-squared", type=float, default=0.0)
    parser.add_argument("--samples", type=int, default=0)
    args = parser.parse_args()

    alpha_cli = Path(__file__).resolve().with_name("cast_alpha.py")
    if not alpha_cli.is_file():
        parser.error("Il runner deve stare accanto a cast_alpha.py.")
    if args.alpha_squared != 0 or args.samples != 0:
        parser.error("Per questo confronto usa --alpha-squared 0 --samples 0.")

    import numpy as np
    import pocketlab.core as core
    from pycasta.utils.preprocessing_utils import (
        calculate_atomic_radii,
        load_and_separate_pdb,
    )

    try:
        pdb_bytes = args.pdb.read_bytes()
        pdb_hash = hashlib.sha256(pdb_bytes).hexdigest()
        ref = json.loads(args.reference_manifest.read_text(encoding="utf-8-sig"))
        backend_path = Path(core.__file__)
        backend_hash = hashlib.sha256(backend_path.read_bytes()).hexdigest()

        coord, ligand, _, ppdb = load_and_separate_pdb(str(args.pdb.resolve()))
        coord = np.asarray(coord, dtype=float)
        baseline = np.asarray(calculate_atomic_radii(ppdb), dtype=float)

        if (
            coord.shape != (832, 3)
            or baseline.shape != (832,)
            or pdb_hash != ref.get("input_sha256")
            or args.alpha_squared != ref.get("alpha_squared")
            or backend_hash != ref.get("backend_sha256")
            or not np.array_equal(coord, ref.get("coordinates"))
            or not np.array_equal(baseline, ref.get("atomic_radii"))
        ):
            raise ValueError("PDB, alpha, coordinate, raggi baseline o backend non coincidono.")

        atoms = [
            line for line in pdb_bytes.decode("ascii").splitlines()
            if line.startswith("ATOM  ")
        ]
        pdb_coords = [[float(line[i:i + 8]) for i in (30, 38, 46)] for line in atoms]
        if len(atoms) != 832 or not np.array_equal(coord, pdb_coords):
            raise ValueError("I record ATOM non coincidono con l'ordine del preprocessing.")

        print("Reference model: exact match", flush=True)

        csv_bytes = args.radii_csv.read_bytes()
        rows = list(csv.DictReader(io.StringIO(csv_bytes.decode("utf-8-sig"), newline="")))
        if len(rows) != 832:
            raise ValueError(f"Attese 832 righe nel CSV, trovate {len(rows)}.")

        radii = []
        for i, (row, line, old) in enumerate(zip(rows, atoms, baseline)):
            expected = (i, line[6:11].strip(), line[17:20].strip(), line[12:16].strip())
            actual = (
                int(row["atom_index"]), row["pdb_serial"].strip(),
                row["residue"].strip(), row["atom"].strip(),
            )
            new = float(row["protor_radius_A"])
            if actual != expected or float(row["baseline_radius_A"]) != old:
                raise ValueError(f"Identita, ordine o raggio baseline errato alla riga {i + 2}.")
            if not math.isfinite(new) or new <= 0:
                raise ValueError(f"Raggio ProtOr non valido alla riga {i + 2}.")
            radii.append(new)

        changed = int(np.count_nonzero(np.asarray(radii) != baseline))
        if changed != 832:
            raise ValueError(f"Attesi 832 raggi modificati, trovati {changed}.")
    except (OSError, KeyError, TypeError, ValueError) as exc:
        parser.error(str(exc))

    print(f"ProtOr override: {changed}/832 raggi modificati", flush=True)
    geometry = core.build_cast_complex_exact(coord.tolist(), radii, args.alpha_squared)
    if (
        not isinstance(geometry, dict)
        or not geometry.get("ok")
        or not geometry.get("alpha_all_dimensions_exported")
    ):
        parser.error("Il backend non ha prodotto un complesso completo; nessun fallback.")

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    dest = args.out / f"cast_protor_{args.pdb.stem}_{stamp}"
    dest.mkdir(parents=True, exist_ok=False)
    manifest = dict(
        input=str(args.pdb.resolve()),
        input_sha256=pdb_hash,
        backend=str(backend_path),
        backend_sha256=backend_hash,
        alpha_squared=args.alpha_squared,
        n_atoms=len(radii),
        n_ligand_atoms_loaded_not_used_for_detection=len(ligand),
        model_note="ProtOr from FreeSASA 2.2.1 CSV; no additive solvent inflation",
        radii_csv=str(args.radii_csv.resolve()),
        radii_csv_sha256=hashlib.sha256(csv_bytes).hexdigest(),
        coordinates=coord.tolist(),
        atomic_radii=radii,
    )
    complex_path = dest / "cast_complex_all_dims.json"
    manifest_path = dest / "manifest.json"
    for path, data in ((complex_path, geometry), (manifest_path, manifest)):
        with path.open("x", encoding="utf-8") as handle:
            json.dump(data, handle, indent=2, allow_nan=False)

    print("CGAL EXPORT:", dest.resolve(), flush=True)
    return subprocess.call([
        sys.executable, str(alpha_cli),
        "--complex", str(complex_path),
        "--manifest", str(manifest_path),
        "--samples", str(args.samples),
        "--out", str(dest / "alpha_results"),
    ])


if __name__ == "__main__":
    raise SystemExit(main())