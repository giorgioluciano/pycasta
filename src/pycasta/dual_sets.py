# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 10:52:03 2025

@author: giorgio
"""

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
dual_sets.py

Module for computing the dual set of a pocket and for extracting dual mouth components.

The idea is as follows:
  1. Compute the closure of the pocket, i.e. collect all faces (triangles) from every tetrahedron
     in the pocket.
  2. For each face, assign a "filtration value" (we use the minimum circumsphere radius of any
     tetrahedron that contributed that face).
  3. Compute the dual set by keeping only those faces whose filtration value is at least a given
     threshold (usually the chosen alpha value).
  4. Compute the boundary of the dual set (those faces that appear only once) and then extract the
     connected components, each of which corresponds to a pocket mouth.
     
Note: This is a simplified version based on the pseudo‐algorithm. You may need to further refine it
for your specific application.
"""

import numpy as np
import logging
from collections import defaultdict
from scipy.spatial import ConvexHull


def get_faces_from_tetrahedron(tet, decimals=3):
    """
    Given a tetrahedron (4x3 array of coordinates), return a list of its 4 triangle faces.
    Each face is represented as a tuple of 3 vertex tuples, with coordinates rounded to a
    specified number of decimals (to allow consistent comparison).
    """
    # Round vertices
    vertices = [tuple(np.round(v, decimals=decimals)) for v in tet]
    # Tetrahedron faces: there are 4 triangles
    face_indices = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    faces = [tuple(sorted([vertices[i] for i in inds])) for inds in face_indices]
    return faces


def compute_closure_of_pocket(pocket_indices, tetra_positions, decimals=3):
    """
    Compute the closure of a pocket defined by a list of tetrahedron indices.
    The closure is the union of all faces (triangles) of the tetrahedra.

    Returns:
      closure: a dictionary mapping each face (tuple) to its occurrence count.
    """
    closure = {}
    for idx in pocket_indices:
        tet = tetra_positions[idx]
        faces = get_faces_from_tetrahedron(tet, decimals=decimals)
        for face in faces:
            closure[face] = closure.get(face, 0) + 1
    return closure


def compute_boundary_from_closure(closure):
    """
    The boundary of the closure is defined as the set of faces that appear only once.

    Returns:
      boundary: a set of faces.
    """
    boundary = {face for face, count in closure.items() if count == 1}
    return boundary


def compute_face_filtration_values(pocket_indices, tetra_positions, radii, decimals=3):
    """
    For each face in the pocket, compute its filtration value as the minimum circumsphere radius
    among the tetrahedra (in the pocket) that contributed that face.

    Parameters:
      pocket_indices: list of indices of tetrahedra in the pocket.
      tetra_positions: array of shape (n_tetra, 4, 3) with tetrahedron coordinates.
      radii: array of circumsphere radii corresponding to tetra_positions.

    Returns:
      face_filtration: dict mapping face (tuple) to a filtration value.
    """
    face_filtration = {}
    for idx in pocket_indices:
        tet = tetra_positions[idx]
        face_value = radii[
            idx
        ]  # Use the tetrahedron’s circumsphere radius as the proxy filtration value
        faces = get_faces_from_tetrahedron(tet, decimals=decimals)
        for face in faces:
            if face in face_filtration:
                face_filtration[face] = min(face_filtration[face], face_value)
            else:
                face_filtration[face] = face_value
    return face_filtration


def compute_dual_set_from_closure(closure, face_filtration, alpha_value):
    """
    Compute the dual set of a pocket by removing from the closure all faces whose filtration value is
    less than the given alpha_value.

    Returns:
      dual_set: set of faces (triangles) that form the dual representation.
    """
    dual_set = {
        face
        for face in closure
        if face_filtration.get(face, float("inf")) >= alpha_value
    }
    return dual_set


def build_connectivity_graph(simplices):
    """
    Build a connectivity graph among simplices.
    Two simplices are considered connected if they share at least two vertices (i.e. an edge).

    Parameters:
      simplices: iterable of simplices (e.g., faces represented as tuples of vertex coordinates).

    Returns:
      graph: dict mapping each simplex to a list of neighboring simplices.
    """
    graph = defaultdict(list)
    simplices = list(simplices)
    for i, s in enumerate(simplices):
        for j in range(i + 1, len(simplices)):
            t = simplices[j]
            # For triangles, sharing two vertices means they share an edge.
            if len(set(s) & set(t)) >= 2:
                graph[s].append(t)
                graph[t].append(s)
    return graph


def find_connected_components(graph):
    """
    Find connected components in a graph.

    Parameters:
      graph: dict mapping node to list of neighbors.

    Returns:
      components: list of sets, each set representing a connected component.
    """
    visited = set()
    components = []
    for node in graph:
        if node not in visited:
            component = set()
            stack = [node]
            while stack:
                cur = stack.pop()
                if cur not in visited:
                    visited.add(cur)
                    component.add(cur)
                    stack.extend(graph[cur])
            components.append(component)
    return components


def compute_dual_mouths(boundary_faces, global_dual=None):
    """
    Given the boundary faces of a pocket (from the closure) and an optional global dual complex,
    compute the connected components that represent the mouth(s) of the pocket.

    Parameters:
      boundary_faces: set of faces (triangles) forming the boundary of the pocket.
      global_dual: set of faces belonging to the dual complex of the entire structure (optional).

    Returns:
      mouth_components: list of sets, each set representing one connected mouth.
    """
    # If provided, remove faces that are in the global dual complex.
    if global_dual is not None:
        filtered_boundary = {face for face in boundary_faces if face not in global_dual}
    else:
        filtered_boundary = boundary_faces

    # Build connectivity graph among the remaining boundary faces.
    graph = build_connectivity_graph(filtered_boundary)
    # Extract connected components.
    mouth_components = find_connected_components(graph)
    return mouth_components


# Example integration function
def compute_dual_set_for_pocket(
    pocket_indices, tetra_positions, radii, alpha_value, decimals=3
):
    """
    Given a pocket (indices into tetra_positions), compute:
      - the closure (all faces),
      - the filtration value for each face,
      - the dual set (faces with filtration value >= alpha_value),
      - and the boundary of the dual set.

    Returns a dictionary with keys:
      'closure', 'face_filtration', 'dual_set', and 'boundary'
    """
    closure = compute_closure_of_pocket(
        pocket_indices, tetra_positions, decimals=decimals
    )
    face_filtration = compute_face_filtration_values(
        pocket_indices, tetra_positions, radii, decimals=decimals
    )
    dual_set = compute_dual_set_from_closure(closure, face_filtration, alpha_value)
    boundary = compute_boundary_from_closure(closure)
    return {
        "closure": closure,
        "face_filtration": face_filtration,
        "dual_set": dual_set,
        "boundary": boundary,
    }


# =============================================================================
# if __name__ == "__main__":
#     # Example demonstration
#     # Suppose we have a pocket consisting of tetrahedra with indices [0, 1, 2]
#     pocket = [0, 1, 2]
#     # Example tetra_positions (each tetrahedron is 4 points in 3D)
#     tetra_positions = np.array([
#         [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
#         [[1, 0, 0], [1, 1, 0], [0, 1, 0], [1, 0, 1]],
#         [[0, 1, 0], [1, 1, 0], [0, 1, 1], [0, 0, 1]],
#     ])
#     # For demonstration, let's compute radii for these tetrahedra
#     # (In practice, you would use your actual function compute_all_circumsphere_radii)
#     def dummy_radii(tets):
#         return np.array([1.0, 1.5, 2.0])
#     radii = dummy_radii(tetra_positions)
#
#     alpha_value = 1.5  # example threshold (could be your MANUAL_ALPHA)
#
#     dual_info = compute_dual_set_for_pocket(pocket, tetra_positions, radii, alpha_value, decimals=3)
#     logging.info("Closure of the pocket:")
#     logging.info(dual_info["closure"])
#     logging.info("Dual set:")
#     logging.info(dual_info["dual_set"])
#     # For mouth duals, one would then call compute_dual_mouths with the boundary.
#     mouth_duals = compute_dual_mouths(dual_info["boundary"], global_dual=None)
#     logging.info("Mouth dual components:")
#     logging.info(mouth_duals)
#
# =============================================================================
