#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
check_geometry.py

Module for performing geometric checks on a pocket.
"""

import numpy as np
from utils.geometry_utils import compute_analytic_pocket_volume

# Example test: compute the volume of a standard tetrahedron
tetrahedron = np.array([[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]])
pocket_indices = [0]
volume = compute_analytic_pocket_volume(tetrahedron, pocket_indices)
print("Calculated volume of tetrahedron:", volume)

tetrahedron2 = tetrahedron + np.array([2, 2, 2])
tetrahedra = np.concatenate([tetrahedron, tetrahedron2], axis=0)
pocket_indices = [0, 1]
volume_two = compute_analytic_pocket_volume(tetrahedra, pocket_indices)
print("Calculated volume of two tetrahedra:", volume_two)
