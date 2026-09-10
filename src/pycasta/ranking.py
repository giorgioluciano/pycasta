#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ranking.py

Module for computing the ranking scores of detected pockets.
"""


# ranking.py
import numpy as np
import logging
from .config import RANKING_METHOD


def compute_ranking_scores(volumes, pockets, flowsteps, connectivitygraph, tetrapositions):
    """Dispatcher: sceglie il metodo in base a config.RANKING_METHOD."""
    method = RANKING_METHOD
    logging.info(f"[ranking] Using method: {method}")
    if method == "percentile_log":
        return _ranking_percentile_log(volumes, pockets, flowsteps, connectivitygraph)
    else:
        return _ranking_baseline(volumes, pockets, flowsteps, connectivitygraph)


# ── Baseline (metodo originale, invariato) ───────────────────────────────────
def _ranking_baseline(volumes, pockets, flowsteps, connectivitygraph):
    scores = []
    for i, (volume, pocket) in enumerate(zip(volumes, pockets)):
        nsteps       = flowsteps.get(i, 1)
        nconnections = len(connectivitygraph.get(i, []))
        norm_volume      = volume
        norm_steps       = np.log1p(nsteps)
        norm_connectivity = np.log1p(nconnections)
        score = 1.0 * norm_volume + 0.5 * norm_steps + 0.2 * norm_connectivity
        logging.debug(f"Pocket {i+1} volume={volume:.2f}, steps={nsteps}, "
                      f"connections={nconnections}, score={score:.2f}")
        scores.append(score)
    return scores


# ── Variante: percentile-log normalizzata ────────────────────────────────────
def _ranking_percentile_log(volumes, pockets, flowsteps, connectivitygraph):
    """
    Normalizza ogni componente sul range osservato (percentile 5-95),
    poi combina con pesi identici alla baseline ma su scala [0,1].
    Più stabile su distribuzioni skewed di volumi.
    """
    if not volumes:
        return []

    vols  = np.array(volumes, dtype=float)
    steps = np.array([flowsteps.get(i, 1) for i in range(len(pockets))], dtype=float)
    conns = np.array([len(connectivitygraph.get(i, [])) for i in range(len(pockets))], dtype=float)

    def _pct_norm(arr):
        lo, hi = np.percentile(arr, 5), np.percentile(arr, 95)
        if hi == lo:
            return np.ones_like(arr) * 0.5
        return np.clip((arr - lo) / (hi - lo), 0.0, 1.0)

    norm_v = _pct_norm(np.log1p(vols))
    norm_s = _pct_norm(np.log1p(steps))
    norm_c = _pct_norm(np.log1p(conns))

    scores = (1.0 * norm_v + 0.5 * norm_s + 0.2 * norm_c).tolist()

    for i, score in enumerate(scores):
        logging.debug(f"[v2] Pocket {i+1} score={score:.4f}")

    return scores