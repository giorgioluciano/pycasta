#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ranking.py

Module for computing the ranking scores of detected pockets.
"""


def compute_ranking_scores(
    volumes, pockets, flow_steps, connectivity_graph, tetra_positions
):
    import numpy as np
    import logging

    scores = []
    for i, (volume, pocket) in enumerate(zip(volumes, pockets)):
        n_steps = flow_steps.get(i, 1)
        n_connections = len(connectivity_graph.get(i, []))
        norm_volume = volume
        norm_steps = np.log1p(n_steps)
        norm_connectivity = np.log1p(n_connections)
        score = 1.0 * norm_volume + 0.5 * norm_steps + 0.2 * norm_connectivity
        logging.debug(
            f"Pocket {i+1}: volume={volume:.2f}, steps={n_steps}, connections={n_connections}, score={score:.2f}"
        )
        scores.append(score)
    return scores
