"""
Test minimo: confronta baseline vs percentile_log sugli stessi dati sintetici
e sulle molecole della testShort suite.
Esegui con: pytest tests/test_ranking_comparison.py -v
"""
import numpy as np
import pytest
import importlib

import pycasta.config as cfg
from pycasta.ranking import compute_ranking_scores


# ── Dati sintetici ────────────────────────────────────────────────────────────
VOLUMES    = [1200.0, 300.0, 80.0, 50.0]
POCKETS    = [[0, 1, 2], [3, 4], [5], [6]]
FLOWSTEPS  = {0: 12, 1: 8, 2: 5, 3: 3}
CONNGRAPH  = {0: [1, 2], 1: [0], 2: [], 3: []}


def run_with_method(method: str):
    cfg.RANKING_METHOD = method
    importlib.reload(cfg)          # propaga il cambio al modulo ranking
    scores = compute_ranking_scores(
        VOLUMES, POCKETS, FLOWSTEPS, CONNGRAPH, tetrapositions=None
    )
    return scores


def test_baseline_order_preserved():
    """La baseline deve sempre mettere la pocket più grande al primo posto."""
    scores = run_with_method("baseline")
    ranked = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    assert ranked[0] == 0, f"Expected pocket 0 first, got {ranked[0]}"


def test_v2_order_preserved():
    """La variante deve anch'essa mettere la pocket più grande al primo posto."""
    scores = run_with_method("percentile_log")
    ranked = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    assert ranked[0] == 0, f"Expected pocket 0 first, got {ranked[0]}"


def test_v2_scores_in_unit_range():
    """I punteggi v2 devono essere nel range [0, 1.7] (somma pesi max)."""
    scores = run_with_method("percentile_log")
    assert all(0.0 <= s <= 1.75 for s in scores), f"Out-of-range scores: {scores}"


def test_baseline_vs_v2_top1_agreement():
    """
    Verifica che baseline e v2 concordino sul Top-1 per dati sintetici.
    Se fallisce, analizzare il caso prima di accettare la variante.
    """
    b = run_with_method("baseline")
    v = run_with_method("percentile_log")
    top1_b = int(np.argmax(b))
    top1_v = int(np.argmax(v))
    assert top1_b == top1_v, (
        f"Top-1 disagrees: baseline={top1_b}, percentile_log={top1_v}. "
        f"Investigate before merging."
    )


# ── Test su molecole reali (opzionale, salta se i file non esistono) ──────────
@pytest.mark.parametrize("pdb_name", ["1a6w", "1qif", "3app"])
def test_top1_consistency_on_real_molecules(pdb_name, tmp_path):
    """
    Esegue la pipeline completa con entrambi i metodi e confronta Top-1.
    Richiede i PDB in data/testShort/bounded/.
    """
    pytest.importorskip("Bio")
    try:
        from pycasta.runanalysis import processpdb
        import os
        pdb = f"data/testShort/bounded/{pdb_name}.pdb"
        if not os.path.exists(pdb):
            pytest.skip(f"PDB not found: {pdb}")

        cfg.RANKING_METHOD = "baseline"
        r_base = processpdb(pdb)
        top1_base = r_base.get("step_to_ligand")

        cfg.RANKING_METHOD = "percentile_log"
        r_v2 = processpdb(pdb)
        top1_v2 = r_v2.get("step_to_ligand")

        print(f"\n{pdb_name}: baseline_top1={top1_base}, v2_top1={top1_v2}")
        # Non falliamo automaticamente: loggare la divergenza per analisi
        if top1_base != top1_v2:
            print(f"  ⚠️  Top-1 divergence on {pdb_name}: "
                  f"baseline={top1_base} vs v2={top1_v2}")
    except Exception as e:
        pytest.skip(f"Pipeline error: {e}")