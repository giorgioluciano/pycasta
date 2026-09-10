from pathlib import Path
from urllib.request import urlretrieve

pairs = [
    ("1bid", "3tms"), ("1cdo", "8adh"), ("1dwd", "1hxf"),
    ("1fbp", "2fbp"), ("1gca", "1gcg"), ("1hew", "1hel"),
    ("1hyt", "1npc"), ("1inc", "1esa"), ("1rbp", "1brq"),
    ("1rob", "8rat"), ("1stp", "1swb"), ("1ulb", "1ula"),
    ("2ifb", "1ifb"), ("3ptb", "3ptn"), ("2ypi", "1ypi"),
    ("4dfr", "5dfr"), ("4phv", "3phv"), ("5cna", "2ctv"),
    ("7cpa", "5cpa"), ("1a6w", "1a6u"), ("1acj", "1qif"),
    ("1apu", "3app"), ("1blh", "1djb"), ("1byb", "1bya"),
    ("1hfc", "1cge"), ("1ida", "1hsi"), ("1igj", "1a4j"),
    ("1imb", "1ime"), ("1ivd", "1nna"), ("1mrg", "1ahc"),
    ("1mtw", "2tga"), ("1okm", "4ca2"), ("1pdz", "1pdy"),
    ("1phd", "1phc"), ("1pso", "1psn"), ("1qpe", "3lck"),
    ("1rne", "1bbs"), ("1snc", "1stn"), ("1srf", "1pts"),
    ("2ctc", "2ctb"), ("2h4n", "2cba"), ("2pk4", "1krn"),
    ("2sim", "2sil"), ("2tmn", "1l3f"), ("3gch", "1chg"),
    ("3mth", "6ins"), ("5p2p", "3p2p"), ("6rsa", "7rat"),
]

root = Path(__file__).resolve().parent
holo_dir = root / "data" / "pocketpicker" / "holo"
apo_dir = root / "data" / "pocketpicker" / "apo"
holo_dir.mkdir(parents=True, exist_ok=True)
apo_dir.mkdir(parents=True, exist_ok=True)

failures = []

for label, entries, target_dir in (("holo", [p[0] for p in pairs], holo_dir), ("apo", [p[1] for p in pairs], apo_dir)):
    for pdb_id in entries:
        out = target_dir / f"{pdb_id}.pdb"
        if out.exists() and out.stat().st_size > 0:
            print(f"[{label}] presente: {out.name}")
            continue
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        try:
            print(f"[{label}] scarico {pdb_id} ...", end=" ")
            urlretrieve(url, out)
            print("OK")
        except Exception as exc:
            failures.append((label, pdb_id, str(exc)))
            if out.exists():
                out.unlink()
            print(f"ERRORE: {exc}")

(root / "data" / "pocketpicker" / "pocketpicker_pairs.csv").write_text(
    "holo_pdb,apo_pdb\n" + "\n".join(f"{h},{a}" for h, a in pairs) + "\n",
    encoding="utf-8",
)

print(f"\nHolo scaricati: {len(list(holo_dir.glob('*.pdb')))} / {len(pairs)}")
print(f"Apo scaricati: {len(list(apo_dir.glob('*.pdb')))} / {len(pairs)}")
if failures:
    print("\nDownload falliti:")
    for label, pdb_id, error in failures:
        print(f"{label},{pdb_id},{error}")
    raise SystemExit(1)
