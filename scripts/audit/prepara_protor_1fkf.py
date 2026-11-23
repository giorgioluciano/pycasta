import argparse
import csv
import hashlib
import json
import math
from importlib.metadata import version
from pathlib import Path

import freesasa


parser = argparse.ArgumentParser()
parser.add_argument("--pdb", type=Path, required=True)
parser.add_argument("--results", type=Path, required=True)
parser.add_argument("--out", type=Path, default=Path("1fkf_protor_radii.csv"))
args = parser.parse_args()

manifest = json.loads(args.results.read_text(encoding="utf-8"))["provenance"]["manifest"]
pdb_bytes = args.pdb.read_bytes()

if hashlib.sha256(pdb_bytes).hexdigest() != manifest["input_sha256"]:
    raise SystemExit("STOP: il PDB non ha lo stesso hash della baseline.")

atoms = [line for line in pdb_bytes.decode("ascii").splitlines()
         if line.startswith("ATOM  ")]
coords = [[float(line[i:i + 8]) for i in (30, 38, 46)] for line in atoms]
baseline_radii = manifest["atomic_radii"]

if len(atoms) != 832 or len(baseline_radii) != 832 or coords != manifest["coordinates"]:
    raise SystemExit("STOP: numero di atomi o coordinate diversi dalla baseline 1FKF.")

classifier = freesasa.Classifier.getStandardClassifier("protor")
rows, missing = [], []

for index, (line, old_radius) in enumerate(zip(atoms, baseline_radii)):
    serial = line[6:11].strip()
    residue = line[17:20].strip()
    atom = line[12:16].strip()
    new_radius = classifier.radius(residue, atom)

    if not math.isfinite(new_radius) or new_radius <= 0:
        missing.append((index, serial, residue, atom))
    else:
        rows.append((
            index, serial, residue, atom, old_radius, new_radius,
            round(new_radius - old_radius, 6),
        ))

print(f"FreeSASA {version('freesasa')} - copertura ProtOr: {len(rows)}/832")

if missing:
    for row in missing:
        print("NON ASSEGNATO:", *row)
    raise SystemExit("STOP: nessun CSV scritto; nessun raggio di ripiego.")

if args.out.exists():
    raise SystemExit(f"STOP: {args.out} esiste gia. Usa un altro percorso --out.")

with args.out.open("x", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle)
    writer.writerow([
        "atom_index", "pdb_serial", "residue", "atom",
        "baseline_radius_A", "protor_radius_A", "delta_A",
    ])
    writer.writerows(rows)

print(f"Raggi modificati: {sum(row[-1] != 0 for row in rows)}/832")
print(f"CSV: {args.out.resolve()}")