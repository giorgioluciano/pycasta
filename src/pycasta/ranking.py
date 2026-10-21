# -*- coding: utf-8 -*-
"""
ranking.py

Pocket ranking and scoring algorithms.
"""

import logging
import numpy as np

from . import config as cfg


def compute_ranking_scores(
    volumes,
    pockets,
    flow_steps,
    connectivity_graph,
    tetra_positions,
    protein_coords=None,
    compute_sasa=True,
):
    """
    Compute ranking scores for pockets.
    
    Parameters
    ----------
    volumes : list of float
        Pocket volumes.
    pockets : list of list of int
        Tetrahedron indices for each pocket.
    flow_steps : dict
        Flow steps for each tetrahedron.
    connectivity_graph : list
        Connectivity graph for each pocket.
    tetra_positions : np.ndarray
        Tetrahedron vertex coordinates.
    protein_coords : np.ndarray, optional
        Protein atom coordinates (required for SASA).
    compute_sasa : bool
        Whether to compute SASA metrics.
    
    Returns
    -------
    scores : list of float
        Ranking scores for each pocket.
    """
    # Compute SASA if requested and protein_coords available
    sasa_areas = []
    if compute_sasa and protein_coords is not None:
        from .utils.sasa_utils import compute_pocket_surface_area
        
        logging.info("Computing SASA for pocket ranking...")
        for i, pocket in enumerate(pockets):
            area = compute_pocket_surface_area(
                pocket,
                tetra_positions,
                protein_coords,
                probe_radius=cfg.SASA_PROBE_RADIUS,
                n_spheres=cfg.SASA_N_SPHERES,
            )
            sasa_areas.append(area)
            logging.debug(f"Pocket {i}: SASA = {area:.2f} Å²")
    else:
        sasa_areas = None
    
    return _ranking_percentile_log(
        volumes, pockets, flow_steps, connectivity_graph, sasa_areas
    )


def _ranking_percentile_log(volumes, pockets, flowsteps, connectivitygraph, sasa_areas=None):
    """
    Percentile log scoring with optional SASA.
    
    Score = volume_weight * vol_percentile + 
            sasa_weight * sasa_percentile + 
            conn_weight * conn_percentile
    
    Parameters
    ----------
    volumes : list of float
        Pocket volumes.
    pockets : list of list of int
        Tetrahedron indices for each pocket.
    flowsteps : dict
        Flow steps (not used in this version).
    connectivitygraph : list
        Connectivity graph for each pocket.
    sasa_areas : list of float, optional
        SASA for each pocket.
    
    Returns
    -------
    scores : list of float
        Ranking scores (0-1 range).
    """
    # Volume score (log-transformed)
    vol_log = np.log1p(np.array(volumes, dtype=float))
    vol_percentile = np.searchsorted(np.sort(vol_log), vol_log) / len(vol_log)
    
    # SASA score (if available)
    if sasa_areas is not None and any(sasa_areas):
        sasa_log = np.log1p(np.array(sasa_areas, dtype=float))
        sasa_percentile = np.searchsorted(np.sort(sasa_log), sasa_log) / len(sasa_log)
    else:
        sasa_percentile = np.zeros(len(volumes))
    
    # Connectivity (number of tetrahedra in pocket)
    conns = np.array([len(pocket) for pocket in pockets], dtype=float)
    conn_percentile = np.searchsorted(np.sort(conns), conns) / len(conns)
    
    # Combined score: weighted sum
    scores = (
        cfg.VOLUME_WEIGHT * vol_percentile +
        cfg.SASA_WEIGHT * sasa_percentile +
        cfg.CONNECTIVITY_WEIGHT * conn_percentile
    )
    
    return scores.tolist()