#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
cgal_wdelaunay.py

Weighted Delaunay triangulation wrapper.
Uses CGAL via pocketlab when enabled, otherwise falls back to SciPy Delaunay.
"""

import logging
import numpy as np
from scipy.spatial import Delaunay

from config import USE_CGAL, WEIGHTED_DELAUNAY


def cgal_weighted_delaunay(protein_coords, radii):
    if USE_CGAL and WEIGHTED_DELAUNAY:
        try:
            from pocketlab import build_regular_triangulation

            result = build_regular_triangulation(protein_coords, radii)

            if isinstance(result, dict) and result.get("ok") and "cells_idx" in result:
                simplices = np.asarray(result["cells_idx"], dtype=int)
                logging.info(
                    f"[CGAL] Weighted Delaunay built successfully: {len(simplices)} cells."
                )
                return simplices, result

            logging.warning(
                "[CGAL] pocketlab returned an invalid result, falling back to SciPy Delaunay."
            )

        except Exception as e:
            logging.warning(
                f"[CGAL] pocketlab unavailable or failed ({e}). Falling back to SciPy Delaunay."
            )

    simplices = Delaunay(protein_coords).simplices
    logging.info(f"[SciPy] Standard Delaunay built: {len(simplices)} simplices.")
    return simplices, None