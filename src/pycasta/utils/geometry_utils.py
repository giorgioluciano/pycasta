#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
geometry_utils.py

Module containing functions for:
- Computing circumsphere radii of tetrahedra.
- Estimating pocket volume using Monte Carlo integration.
- Validating pockets via mesh extrusion and fake sphere tests.
"""

import os
import numpy as np
from scipy.spatial import Delaunay, ConvexHull
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
import trimesh
from config import (
    SAVE_EXTRUDED_MESHES,
    OUTPUT_DIR,
    MESH_EXTRUSION_DISTANCE,
    LIGAND_THRESHOLD,
)


def get_unique_vertices(tetra_positions, pocket_indices):
    vertices = []
    for idx in pocket_indices:
        vertices.extend(tetra_positions[idx])
    vertices = np.array(vertices)
    return np.unique(vertices, axis=0)


def get_convex_hull_points(points):
    if len(points) < 4:
        return points
    hull = ConvexHull(points)
    return points[hull.vertices]


def validate_ligand_in_extruded_mesh(
    mouth_coords, ligand_coords, extrusion_distance=MESH_EXTRUSION_DISTANCE
):
    if len(mouth_coords) < 3 or ligand_coords.size == 0:
        return False, None, None

    from sklearn.decomposition import PCA
    from scipy.spatial import ConvexHull
    import logging
    import trimesh

    try:
        # Flatten the mouth coordinates using PCA
        pca = PCA(n_components=2)
        flat_2d = pca.fit_transform(mouth_coords)

        # Create 2D convex hull of the flattened points
        hull_2d = ConvexHull(flat_2d)
        ordered_indices = hull_2d.vertices
        cap_3d = mouth_coords[ordered_indices]

        # Compute extrusion normal
        normal_3d = np.cross(cap_3d[1] - cap_3d[0], cap_3d[2] - cap_3d[0])
        if np.linalg.norm(normal_3d) == 0:
            return False, None, None
        normal_3d /= np.linalg.norm(normal_3d)

        # Extrude cap
        extruded = cap_3d + extrusion_distance * normal_3d
        all_vertices = np.vstack([cap_3d, extruded])
        num = len(cap_3d)

        # Construct mesh faces
        faces = []
        for i in range(1, num - 1):
            faces.append([0, i, i + 1])
        for i in range(1, num - 1):
            faces.append([num, num + i + 1, num + i])
        for i in range(num):
            j = (i + 1) % num
            faces.append([i, j, j + num])
            faces.append([i, j + num, i + num])

        # Build the mesh
        mesh = trimesh.Trimesh(vertices=all_vertices, faces=faces, process=True)

        # Check if ligand atoms are inside the mesh
        inside_mask = mesh.contains(ligand_coords)
        inside = np.any(inside_mask)

        # Compute signed distances
        signed_distances = mesh.nearest.signed_distance(ligand_coords)
        min_distance = float(np.min(np.abs(signed_distances)))

        # Key fix: treat as valid contact if inside OR very close
        in_contact = inside or min_distance <= LIGAND_THRESHOLD

        return in_contact, mesh, min_distance

    except Exception as e:
        logging.warning(f"⚠️ Failed to validate ligand in extruded mesh: {e}")
        return False, None, None


def validate_ligand_in_fake_sphere(mouth_coords, ligand_coords, radius=1.4):
    if len(mouth_coords) == 0 or ligand_coords.size == 0:
        return False
    mouth_center = np.mean(mouth_coords, axis=0)
    distances = np.linalg.norm(ligand_coords - mouth_center, axis=1)
    print(f"✅ Ligand atoms: {ligand_coords.shape[0]}")
    print(f"✅ Mouth center coordinates: {mouth_center}")
    print(f"✅ Minimum ligand-mouth distance: {np.min(distances):.2f} Å")
    return np.any(distances <= radius)


def compute_circumsphere_radius(tetra):
    if tetra.shape[0] != 4:
        return np.inf
    A, B, C, D = tetra
    AB, AC, AD = B - A, C - A, D - A
    M = np.vstack([AB, AC, AD]).T
    try:
        rhs = 0.5 * np.array([np.dot(AB, AB), np.dot(AC, AC), np.dot(AD, AD)])
        sol = np.linalg.solve(M, rhs)
        radius = np.linalg.norm(sol)
    except np.linalg.LinAlgError:
        radius = np.inf
    return radius


def log_radii_statistics(radii):
    stats = {
        "min": float(np.min(radii)),
        "max": float(np.max(radii)),
        "mean": float(np.mean(radii)),
    }
    import logging

    logging.info(f"Circumsphere radii stats: {stats}")
    return stats


def compute_all_circumsphere_radii(tetra_positions):
    radii = np.array([compute_circumsphere_radius(tetra) for tetra in tetra_positions])
    log_radii_statistics(radii)
    return radii


def monte_carlo_pocket_volume(tetra_positions, pocket_indices, num_samples=10000):
    pocket_tetra = tetra_positions[pocket_indices]
    all_points = pocket_tetra.reshape(-1, 3)
    min_bounds = np.min(all_points, axis=0)
    max_bounds = np.max(all_points, axis=0)
    bbox_volume = np.prod(max_bounds - min_bounds)
    random_points = np.random.uniform(
        low=min_bounds, high=max_bounds, size=(num_samples, 3)
    )
    count_inside = 0
    for p in random_points:
        if any(point_in_tetrahedron(p, tetra) for tetra in pocket_tetra):
            count_inside += 1
    return (count_inside / num_samples) * bbox_volume


def monte_carlo_ligand_volume(ligand_coords, num_samples=10000):
    if ligand_coords.shape[0] < 4:
        return 0.0
    min_bounds = np.min(ligand_coords, axis=0)
    max_bounds = np.max(ligand_coords, axis=0)
    bbox_volume = np.prod(max_bounds - min_bounds)
    delaunay = Delaunay(ligand_coords)
    random_points = np.random.uniform(
        low=min_bounds, high=max_bounds, size=(num_samples, 3)
    )
    inside = delaunay.find_simplex(random_points) >= 0
    return float(np.sum(inside)) / num_samples * bbox_volume


def point_in_tetrahedron(p, tetra):
    A, B, C, D = tetra
    M = np.column_stack([B - A, C - A, D - A])
    try:
        u, v, w = np.linalg.solve(M, p - A)
        return (u >= 0) and (v >= 0) and (w >= 0) and (u + v + w <= 1)
    except np.linalg.LinAlgError:
        return False


def compute_triangle_area(triangle):
    a, b, c = triangle
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a))


def compute_triangle_perimeter(triangle):
    a, b, c = triangle
    return np.linalg.norm(a - b) + np.linalg.norm(b - c) + np.linalg.norm(c - a)


def compute_circumcenter(tetra):
    A, B, C, D = tetra
    M = np.vstack([B - A, C - A, D - A]).T
    rhs = 0.5 * np.array(
        [np.dot(B - A, B - A), np.dot(C - A, C - A), np.dot(D - A, D - A)]
    )
    try:
        sol = np.linalg.solve(M, rhs)
        circumcenter = A + sol
    except np.linalg.LinAlgError:
        circumcenter = np.array([np.nan, np.nan, np.nan])
    return circumcenter


def compute_volume_for_tetra_group(tetra_positions, tetra_indices):
    total_volume = 0.0
    for idx in tetra_indices:
        verts = tetra_positions[idx]
        a, b, c, d = verts
        vol = np.abs(np.dot(a - d, np.cross(b - d, c - d))) / 6.0
        total_volume += vol
    return total_volume
