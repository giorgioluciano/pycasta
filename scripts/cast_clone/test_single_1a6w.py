from datetime import datetime
from pathlib import Path
import hashlib
import logging
import subprocess

from pycasta import config as cfg


def main():
    # Seleziona un solo input prima di avviare l'analisi.
    matches = [
        p.resolve()
        for p in Path(cfg.BOUNDED_DIR).glob("*.pdb")
        if p.stem.lower() == "1a6w"
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"Atteso un solo 1a6w.pdb in {cfg.BOUNDED_DIR}; trovati: {matches}"
        )
    pdb_path = matches[0]

    # Directory nuova: se esiste già, interrompe invece di riutilizzarla.
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    root = Path("output_isolated_tests").resolve() / f"1a6w_fixed_alpha_{stamp}"
    root.mkdir(parents=True, exist_ok=False)

    # Imposta la radice PRIMA di importare i moduli dell'analisi.
    # VERSION_TAG e parametri scientifici restano invariati.
    cfg.OUTPUT_DIR = str(root)
    cfg.USE_EXISTING_RESULTS = False

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        handlers=[
            logging.FileHandler(root / "run.log", encoding="utf-8"),
            logging.StreamHandler(),
        ],
        force=True,
    )

    import pycasta.run_analysis as ra

    out_dir = root / cfg.VERSION_TAG
    out_dir.mkdir(parents=True, exist_ok=False)

    # Conserva configurazione, input e modifica al sorgente.
    settings = "\n".join(
        f"{name} = {getattr(cfg, name)!r}"
        for name in sorted(dir(cfg))
        if name.isupper()
    )
    (root / "config_snapshot.txt").write_text(settings, encoding="utf-8")

    provenance = [
        f"input = {pdb_path}",
        f"input_sha256 = {hashlib.sha256(pdb_path.read_bytes()).hexdigest()}",
        f"run_analysis = {ra.__file__}",
        f"alpha = {cfg.get_alpha()!r}",
        f"min_pocket_volume = {cfg.get_min_pocket_volume()!r}",
    ]
    for args in (
        ["git", "rev-parse", "HEAD"],
        ["git", "diff", "--", "src/pycasta/run_analysis.py"],
    ):
        completed = subprocess.run(
            args, capture_output=True, text=True, errors="replace"
        )
        provenance.extend(
            ["\nCOMMAND: " + " ".join(args), completed.stdout, completed.stderr]
        )
    (root / "provenance.txt").write_text(
        "\n".join(provenance), encoding="utf-8"
    )

    logging.info("Single input: %s", pdb_path)
    logging.info("Isolated output root: %s", root)
    logging.info("Alpha: %s", cfg.get_alpha())
    logging.info("Validation: %s", cfg.VALIDATION_METHOD)

    try:
        results = ra.run_analysis_mode(
            pdb_sources=[str(pdb_path)],
            analysis_type="single",
        )
        ra.save_summary_json(
            results,
            str(out_dir / "single_summary.json"),
            analysis_type="single",
        )
        ra.save_summary_csv(
            results,
            str(out_dir / "single_summary.csv"),
            analysis_type="single",
        )
    except Exception:
        logging.exception("Single-structure test failed")
        raise

    if not results:
        logging.warning(
            "Nessun risultato nel riepilogo: controllare il log; "
            "run_analysis_mode esclude i casi senza pocket."
        )

    for result in results:
        for key in (
            "file",
            "num_pockets",
            "pocket_rank",
            "step_to_ligand",
            "step_to_ligand_mesh",
            "pocket_volumes",
        ):
            logging.info("%s = %s", key, result.get(key))

    logging.info("Test complete. Output root: %s", root)


if __name__ == "__main__":
    main()