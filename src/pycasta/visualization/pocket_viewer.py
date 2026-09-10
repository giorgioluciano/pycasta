# -*- coding: utf-8 -*-
"""visualization/pocket_viewer.py - Finestra 3/3: Pocket Overview"""
import numpy as np
import pyvista as pv
from pathlib import Path
from typing import Optional
from core import MolData

PCOLS = ["#4caf50", "#2196f3", "#ff9800", "#e91e63", "#00bcd4", "#cddc39", "#9c27b0", "#ff5722"]


def _toggle(actor, plotter):
    if actor is not None:
        actor.SetVisibility(not bool(actor.GetVisibility()))
        plotter.render()


def show_pockets(mol: MolData, results_dir: Optional[str] = None):
    print(f"\n=== Pocket Overview: {mol.name} ===")
    pl = pv.Plotter(window_size=(1000, 800), title=f"Pocket Overview - {mol.name}")
    pl.set_background("#0f1117")
    actors = {}

    if len(mol.coords):
        actors["protein"] = pl.add_mesh(
            pv.PolyData(mol.coords.astype(np.float32)),
            color="#607d8b",
            point_size=4,
            render_points_as_spheres=True,
            opacity=0.3,
            label="Protein",
        )

    if len(mol.ligand_coords):
        actors["ligand"] = pl.add_mesh(
            pv.PolyData(mol.ligand_coords.astype(np.float32)),
            color="white",
            point_size=13,
            render_points_as_spheres=True,
            opacity=0.98,
            label="Ligand",
        )

    reps = []
    for i, pocket_idx in enumerate(mol.ranked_pockets[:6]):
        col = PCOLS[i % len(PCOLS)]
        if pocket_idx < len(mol.representative_points):
            reps.append(mol.representative_points[pocket_idx])
    reps = np.asarray(reps, dtype=np.float32) if len(reps) else np.empty((0, 3), dtype=np.float32)

    pocket_count = {"value": min(1, len(reps))}
    pocket_actor = {"value": None}

    def refresh_pockets():
        if pocket_actor["value"] is not None:
            pl.remove_actor(pocket_actor["value"], render=False)
            pocket_actor["value"] = None
        n = pocket_count["value"]
        if n > 0 and len(reps):
            pdata = pv.PolyData(reps[:n])
            pocket_actor["value"] = pl.add_mesh(
                pdata,
                color="#ff4d6d",
                point_size=16,
                render_points_as_spheres=True,
                opacity=1.0,
                label="Top pockets",
                render=False,
            )
        pl.render()

    refresh_pockets()

    pl.add_key_event("1", lambda: _toggle(actors.get("protein"), pl))
    pl.add_key_event("2", lambda: _toggle(actors.get("ligand"), pl))
    pl.add_key_event("3", lambda: _toggle(pocket_actor.get("value"), pl))
    pl.add_key_event("Prior", lambda: (pocket_count.__setitem__("value", min(pocket_count["value"] + 1, len(reps))), refresh_pockets()))
    pl.add_key_event("Next", lambda: (pocket_count.__setitem__("value", max(pocket_count["value"] - 1, 0)), refresh_pockets()))
    pl.add_key_event("PageUp", lambda: (pocket_count.__setitem__("value", min(pocket_count["value"] + 1, len(reps))), refresh_pockets()))
    pl.add_key_event("PageDown", lambda: (pocket_count.__setitem__("value", max(pocket_count["value"] - 1, 0)), refresh_pockets()))
    pl.add_key_event("0", pl.reset_camera)

    pl.add_text(
        "1 protein | 2 ligand | 3 pockets | PgUp/PgDn pocket count | 0 reset",
        position="upper_left",
        font_size=10,
        color="white",
    )
    pl.show()
