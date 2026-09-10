# -*- coding: utf-8 -*-
"""visualization/alpha_viewer.py - Finestra 1/3: Alpha Shape Debug"""

import numpy as np
import pyvista as pv
from typing import Optional
from core import MolData, guess_alpha_values


ECOLS = {
    "C": (0.5, 0.5, 0.5),
    "N": (0.2, 0.4, 0.8),
    "O": (0.9, 0.2, 0.2),
    "S": (0.9, 0.8, 0.1),
    "P": (1.0, 0.6, 0.0),
    "H": (0.9, 0.9, 0.9),
}
PCOLS = [
    (0.3, 0.75, 0.3),
    (0.2, 0.5, 0.9),
    (0.95, 0.55, 0.1),
    (0.9, 0.2, 0.55),
    (0.1, 0.8, 0.8),
    (0.8, 0.8, 0.1),
]


def _tet_edges(tets):
    ed = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    pts, lines, idx = [], [], 0
    for tet in tets:
        for a, b in ed:
            pts += [tet[a], tet[b]]
            lines += [2, idx, idx + 1]
            idx += 2
    if not pts:
        return None
    mesh = pv.PolyData()
    mesh.points = np.array(pts, dtype=np.float32)
    mesh.lines = np.array(lines, dtype=np.int32)
    return mesh


def _toggle(actor, plotter):
    if actor is not None:
        actor.SetVisibility(not bool(actor.GetVisibility()))
        plotter.render()


def show_alpha(mol: MolData, alpha_values: Optional[dict] = None):
    if alpha_values is None:
        alpha_values = guess_alpha_values(mol)
    av = alpha_values

    print(f"\n=== Alpha Shape Debug: {mol.name} ===")
    print(f" inside={mol.alpha_mask.sum()} outside={(~mol.alpha_mask).sum()}")
    print(f" alpha -> tight={av['tight']:.2f} medium={av['medium']:.2f} loose={av['loose']:.2f}")

    pl = pv.Plotter(window_size=(1000, 800), title=f"Alpha Shape - {mol.name}")
    pl.set_background("#0f1117")
    actors = {}

    if len(mol.coords):
        ec = np.array([
            ECOLS.get("".join(c for c in n if c.isalpha()).upper()[:1], (0.6, 0.6, 0.6))
            for n in mol.atom_names
        ], dtype=np.float32)
        pc = pv.PolyData(mol.coords.astype(np.float32))
        pc["colors"] = (ec * 255).astype(np.uint8)
        actors["protein"] = pl.add_mesh(
            pc,
            scalars="colors",
            rgb=True,
            point_size=5,
            render_points_as_spheres=True,
            opacity=0.4,
            label="Protein",
        )

    ins = mol.alpha_mask
    out = ~mol.alpha_mask

    if ins.any():
        m = _tet_edges(mol.tetrahedra[ins])
        if m is not None:
            actors["alpha_inside"] = pl.add_mesh(
                m, color="#1e88e5", line_width=1.2, opacity=0.55, label="alpha inside"
            )

    if out.any():
        m = _tet_edges(mol.tetrahedra[out])
        if m is not None:
            actors["alpha_outside"] = pl.add_mesh(
                m, color="#fdd835", line_width=0.8, opacity=0.20, label="alpha outside"
            )

    raw_tets = getattr(mol, "triangulation_tetrahedra", None)
    if raw_tets is not None and len(raw_tets):
        m = _tet_edges(raw_tets)
        if m is not None:
            actors["triangulation_raw"] = pl.add_mesh(
                m,
                color="#9e9e9e",
                line_width=0.35,
                opacity=0.10,
                label="raw triangulation",
                render=False,
            )
            actors["triangulation_raw"].SetVisibility(False)

    if len(mol.ligand_coords):
        actors["ligand"] = pl.add_mesh(
            pv.PolyData(mol.ligand_coords.astype(np.float32)),
            color="white",
            point_size=12,
            render_points_as_spheres=True,
            opacity=0.95,
            label="Ligand",
        )

    rep_points = []
    for i, pt in enumerate(mol.representative_points):
        rep_points.append(pt)
    rep_points = np.asarray(rep_points, dtype=np.float32) if len(rep_points) else np.empty((0, 3), dtype=np.float32)

    pocket_count = {"value": min(1, len(rep_points))}
    pocket_actor = {"value": None}

    def refresh_pockets():
        if pocket_actor["value"] is not None:
            pl.remove_actor(pocket_actor["value"], render=False)
            pocket_actor["value"] = None
        n = pocket_count["value"]
        if n > 0 and len(rep_points):
            subset = pv.PolyData(rep_points[:n])
            pocket_actor["value"] = pl.add_mesh(
                subset,
                color="#ff4d6d",
                point_size=14,
                render_points_as_spheres=True,
                opacity=1.0,
                label="representative points",
                render=False,
            )
        pl.render()

    refresh_pockets()

    pl.add_key_event("1", lambda: _toggle(actors.get("protein"), pl))
    pl.add_key_event("2", lambda: _toggle(actors.get("ligand"), pl))
    pl.add_key_event("3", lambda: _toggle(actors.get("alpha_inside"), pl))
    pl.add_key_event("4", lambda: _toggle(actors.get("alpha_outside"), pl))
    pl.add_key_event("5", lambda: _toggle(pocket_actor.get("value"), pl))
    pl.add_key_event("7", lambda: _toggle(actors.get("triangulation_raw"), pl))

    pl.add_key_event(
        "Prior",
        lambda: (
            pocket_count.__setitem__("value", min(pocket_count["value"] + 1, len(rep_points))),
            refresh_pockets(),
        ),
    )
    pl.add_key_event(
        "Next",
        lambda: (
            pocket_count.__setitem__("value", max(pocket_count["value"] - 1, 0)),
            refresh_pockets(),
        ),
    )
    pl.add_key_event(
        "PageUp",
        lambda: (
            pocket_count.__setitem__("value", min(pocket_count["value"] + 1, len(rep_points))),
            refresh_pockets(),
        ),
    )
    pl.add_key_event(
        "PageDown",
        lambda: (
            pocket_count.__setitem__("value", max(pocket_count["value"] - 1, 0)),
            refresh_pockets(),
        ),
    )
    pl.add_key_event("0", pl.reset_camera)

    pl.add_text(
        "1 protein  2 ligand  3 alpha in  4 alpha out  5 pockets  7 raw triangulation  PgUp/PgDn pockets  0 reset",
        position="upper_left",
        font_size=10,
        color="white",
    )
    pl.show()