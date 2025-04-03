#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
visualize_extruded_mesh.py

Module for visualizing the extruded mesh of a pocket along with the ligand.
"""

import os
import json
import numpy as np
import trimesh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from sklearn.decomposition import PCA
from config import OUTPUT_DIR, VERSION_TAG, MESH_EXTRUSION_DISTANCE

def load_pocket_mesh_data(pdb_id, pocket_idx=0):
    json_path = os.path.join(OUTPUT_DIR, VERSION_TAG, f"{pdb_id}.json")
    if not os.path.exists(json_path):
        raise FileNotFoundError(f"Result JSON not found: {json_path}")
    with open(json_path, "r") as f:
        data = json.load(f)
    ligand = np.array(data["ligand_coords"])
    protein = np.array(data["protein_coords"])
    mouth_rim_list = data.get("mouth_rim_atoms") or []
    if not mouth_rim_list or pocket_idx >= len(mouth_rim_list):
        raise ValueError("Missing rim atoms for the selected pocket.")
    rim_atoms = mouth_rim_list[pocket_idx]
    rim_coords = np.array([protein[i] for i in rim_atoms])
    return ligand, rim_coords

from scipy.spatial import ConvexHull

def build_extruded_mesh(mouth_coords, extrusion_distance):
    if len(mouth_coords) < 3:
        raise ValueError("Not enough rim atoms to build the mesh.")
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    flat_2d = pca.fit_transform(mouth_coords)
    hull_2d = ConvexHull(flat_2d)
    ordered_indices = hull_2d.vertices
    cap_3d = mouth_coords[ordered_indices]
    normal_3d = np.cross(cap_3d[1] - cap_3d[0], cap_3d[2] - cap_3d[0])
    normal_3d /= np.linalg.norm(normal_3d)
    extruded = cap_3d + extrusion_distance * normal_3d
    all_vertices = np.vstack([cap_3d, extruded])
    num = len(cap_3d)
    faces = []
    for i in range(1, num - 1):
        faces.append([0, i, i + 1])
    for i in range(1, num - 1):
        faces.append([0 + num, i + 1 + num, i + num])
    for i in range(num):
        j = (i + 1) % num
        faces.append([i, j, j + num])
        faces.append([i, j + num, i + num])
    return all_vertices, faces

def plot_extruded_mesh(vertices, faces, ligand_coords):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    mesh = Poly3DCollection(vertices[faces], alpha=0.4, facecolor='skyblue', edgecolor='k')
    ax.add_collection3d(mesh)
    ax.scatter(ligand_coords[:, 0], ligand_coords[:, 1], ligand_coords[:, 2],
               c='red', s=20, label='Ligand')
    ligand_center = np.mean(ligand_coords, axis=0)
    pocket_center = np.mean(vertices, axis=0)
    ax.plot([ligand_center[0], pocket_center[0]],
            [ligand_center[1], pocket_center[1]],
            [ligand_center[2], pocket_center[2]], 'gray', linestyle='--')
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Extruded Pocket Mesh and Ligand")
    ax.legend()
    plt.tight_layout()
    plt.show()

def run(pdb_id, pocket_idx=0):
    ligand, rim_coords = load_pocket_mesh_data(pdb_id, pocket_idx)
    verts, faces = build_extruded_mesh(rim_coords, MESH_EXTRUSION_DISTANCE)
    plot_extruded_mesh(verts, faces, ligand)

if __name__ == "__main__":
    run("1fbp", pocket_idx=0)
