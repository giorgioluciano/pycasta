#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ply_to_obj.py

Script to convert .ply files to .obj format.
"""

import trimesh
import os

ply_path_1 = "/mnt/data/1rob_pocket_1.ply"
ply_path_2 = "/mnt/data/1rob_pocket_2.ply"

mesh1 = trimesh.load(ply_path_1, file_type="ply")
mesh2 = trimesh.load(ply_path_2, file_type="ply")

obj_path_1 = ply_path_1.replace(".ply", ".obj")
obj_path_2 = ply_path_2.replace(".ply", ".obj")

mesh1.export(obj_path_1)
mesh2.export(obj_path_2)

print(obj_path_1, obj_path_2)
