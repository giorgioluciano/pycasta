# -*- coding: utf-8 -*-
"""visualization/triangulation_viewer.py"""

import numpy as np
import pyvista as pv


def _tetra_edges():
    return [
        (0, 1), (0, 2), (0, 3),
        (1, 2), (1, 3), (2, 3),
    ]


def _build_wireframe_lines(tetrahedra: np.ndarray):
    points = []
    lines = []

    for tet in tetrahedra:
        base = len(points)
        points.extend(tet.tolist())

        for i, j in _tetra_edges():
            lines.append(2)
            lines.append(base + i)
            lines.append(base + j)

    if not points:
        return None

    poly = pv.PolyData(np.asarray(points, dtype=float))
    poly.lines = np.asarray(lines, dtype=np.int32)
    return poly


def show_triangulation(mol, max_tetra=None):
    tetra = getattr(mol, "triangulation_tetrahedra", None)

    if tetra is None or len(tetra) == 0:
        print("[triangulation_viewer] No triangulation debug data available.")
        return

    tetra = np.asarray(tetra, dtype=float)

    if max_tetra is not None and len(tetra) > max_tetra:
        tetra = tetra[:max_tetra]

    line_mesh = _build_wireframe_lines(tetra)
    if line_mesh is None:
        print("[triangulation_viewer] Empty triangulation.")
        return

    pl = pv.Plotter(window_size=(1400, 900))
    pl.set_background("#0b0f14")

    if getattr(mol, "coords", None) is not None and len(mol.coords) > 0:
        pl.add_points(
            np.asarray(mol.coords, dtype=float),
            color="white",
            point_size=4,
            render_points_as_spheres=True,
            opacity=0.35,
        )

    if getattr(mol, "ligand_coords", None) is not None and len(mol.ligand_coords) > 0:
        pl.add_points(
            np.asarray(mol.ligand_coords, dtype=float),
            color="#ff7f50",
            point_size=9,
            render_points_as_spheres=True,
            opacity=1.0,
        )

    pl.add_mesh(
        line_mesh,
        color="#9aa4b2",
        line_width=1,
        opacity=0.35,
    )

    pl.add_text(
        f"{mol.name} | raw triangulation | tetra={len(tetra)}",
        font_size=12,
        color="white",
    )
    pl.add_axes()
    pl.show()