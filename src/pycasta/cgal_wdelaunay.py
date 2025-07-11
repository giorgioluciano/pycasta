#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
cgal_wdelaunay.py

Optional module for weighted Delaunay triangulation using CGAL.
If USE_CGAL is True in config, this function is used; otherwise, standard Delaunay is used.
"""

import numpy as np
from scipy.spatial import Delaunay


def cgal_weighted_delaunay(protein_coords, radii):
    # Dummy implementation: fall back to standard Delaunay triangulation
    simplices = Delaunay(protein_coords).simplices
    return simplices, None
