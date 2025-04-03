# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 13:00:08 2025

@author: giorgio
"""

from plyfile import PlyData

def load_ply_to_pymol(name, ply_path, color=(1.0, 0.5, 0.0)):
    plydata = PlyData.read(ply_path)
    vertices = np.vstack([plydata['vertex'][axis] for axis in ['x', 'y', 'z']]).T
    faces = np.vstack([f[0] for f in plydata['face'].data])
    draw_extruded_pocket(name, vertices, faces, color=color)