#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""CGAL/pocketlab regular triangulation with input-atom index mapping."""

import logging

import numpy as np


def weighted_delaunay(protein_coords, radii):
    """Return tetrahedra indexed into the original protein_coords array."""
    try:
        from pocketlab import build_regular_triangulation

        coords = np.asarray(protein_coords, dtype=float)
        atom_radii = np.asarray(radii, dtype=float)

        if coords.ndim != 2 or coords.shape[1] != 3:
            raise ValueError("protein_coords must have shape (N, 3).")
        if atom_radii.shape != (len(coords),):
            raise ValueError("radii must have shape (N,).")
        if not np.all(np.isfinite(coords)):
            raise ValueError("protein_coords contains non-finite values.")
        if not np.all(np.isfinite(atom_radii)) or np.any(atom_radii < 0):
            raise ValueError("radii must contain finite, non-negative values.")

        input_index = {}
        for atom_id, xyz in enumerate(coords):
            key = tuple(float(x) for x in xyz)
            if key in input_index:
                raise RuntimeError(
                    "Duplicate input coordinates: cannot map vertices "
                    "unambiguously with this compatibility adapter."
                )
            input_index[key] = atom_id

        result = build_regular_triangulation(coords, atom_radii)

        if not isinstance(result, dict) or not result.get("ok"):
            raise RuntimeError(f"Invalid pocketlab result: {result}")
        if "cells_idx" not in result or "vertices" not in result:
            raise RuntimeError(
                "pocketlab result must contain cells_idx and vertices."
            )

        vertex_to_atom = {}
        for vertex in result["vertices"]:
            vertex_id = int(vertex["id"])
            key = (
                float(vertex["x"]),
                float(vertex["y"]),
                float(vertex["z"]),
            )
            if key not in input_index:
                raise RuntimeError(
                    f"CGAL vertex {vertex_id} does not match any input atom."
                )
            if vertex_id in vertex_to_atom:
                raise RuntimeError(f"Duplicate CGAL vertex ID: {vertex_id}")

            vertex_to_atom[vertex_id] = input_index[key]

        raw_cells = np.asarray(result["cells_idx"], dtype=int)
        if raw_cells.size == 0:
            raise RuntimeError("CGAL returned no finite tetrahedra.")
        if raw_cells.ndim != 2 or raw_cells.shape[1] != 4:
            raise RuntimeError(
                f"Invalid CGAL cell array shape: {raw_cells.shape}"
            )

        simplices = np.asarray(
            [
                [vertex_to_atom[int(vertex_id)] for vertex_id in cell]
                for cell in raw_cells
            ],
            dtype=int,
        )

        info = dict(result)
        info["input_cells_idx"] = simplices.tolist()
        info["vertex_to_input_atom"] = vertex_to_atom
        info["returned_simplices_index_space"] = "input_atoms"

        logging.info(
            "[CGAL] Weighted Delaunay: %d cells; "
            "vertex indices mapped to original input atoms.",
            len(simplices),
        )
        return simplices, info

    except Exception as exc:
        logging.error("[weighted_delaunay] pocketlab failed: %s", exc)
        raise RuntimeError(
            f"Weighted Delaunay triangulation failed: {exc}"
        ) from exc