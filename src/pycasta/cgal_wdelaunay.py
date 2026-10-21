#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
cgal_wdelaunay.py

Weighted (regular) Delaunay triangulation using CGAL via pocketlab.
CGAL is required - no fallback.
"""

import logging
import numpy as np


def weighted_delaunay(protein_coords, radii):
    """
    Weighted (regular) Delaunay triangulation using CGAL via pocketlab.
    
    Parameters
    ----------
    protein_coords : np.ndarray
        Array of shape (N, 3) containing protein atom coordinates.
    radii : np.ndarray
        Array of shape (N,) containing atomic radii (e.g., van der Waals radii).
    
    Returns
    -------
    simplices : np.ndarray
        Array of shape (M, 4) containing tetrahedron vertex indices.
    triangulation_info : dict
        Additional information from the triangulation.
    
    Raises
    ------
    RuntimeError
        If pocketlab/CGAL is not available or triangulation fails.
    """
    try:
        from pocketlab import build_regular_triangulation

        result = build_regular_triangulation(protein_coords, radii)

        if isinstance(result, dict) and result.get('ok') and 'cells_idx' in result:
            simplices = np.asarray(result['cells_idx'], dtype=int)
            logging.info(
                f'[CGAL] Weighted Delaunay built successfully: {len(simplices)} cells.'
            )
            return simplices, result

        raise RuntimeError(
            f'[CGAL] pocketlab returned invalid result: {result}'
        )

    except Exception as e:
        logging.error(
            f'[weighted_delaunay] pocketlab failed: {e}'
        )
        raise RuntimeError(
            f'Weighted Delaunay triangulation failed: {e}'
        ) from e