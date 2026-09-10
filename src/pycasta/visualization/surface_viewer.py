# -*- coding: utf-8 -*-
"""visualization/surface_viewer.py - Finestra 2/3: Superficie VdW colorata per AO"""
import numpy as np
import pyvista as pv
from core import MolData



def _build_surface(coords, radii, res=0.5):
    pad = 3.0
    lo = coords.min(0) - pad
    hi = coords.max(0) + pad
    nx, ny, nz = int((hi[0] - lo[0]) / res) + 1, int((hi[1] - lo[1]) / res) + 1, int((hi[2] - lo[2]) / res) + 1
    gx, gy, gz = np.meshgrid(
        np.linspace(lo[0], hi[0], nx),
        np.linspace(lo[1], hi[1], ny),
        np.linspace(lo[2], hi[2], nz),
        indexing="ij",
    )
    gpts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
    sdf = np.full(len(gpts), -np.inf, dtype=np.float32)
    for s in range(0, len(coords), 64):
        e = min(s + 64, len(coords))
        d = np.linalg.norm(gpts[:, None, :] - coords[None, s:e, :], axis=2)
        sdf = np.maximum(sdf, (radii[None, s:e] - d).max(axis=1))
    grid = pv.ImageData(dimensions=(nx, ny, nz), spacing=(res, res, res), origin=lo)
    grid.point_data["sdf"] = sdf.reshape(nx, ny, nz).ravel(order="F")
    return grid.contour([0.0], scalars="sdf")


def _interp_ao(surface, coords, ao, k=4):
    from scipy.spatial import KDTree
    dists, idxs = KDTree(coords).query(np.array(surface.points), k=k)
    dists = np.clip(dists, 1e-6, None)
    w = 1 / dists**2
    w /= w.sum(1, keepdims=True)
    return (w * ao[idxs]).sum(1).astype(np.float32)


def _toggle(actor, plotter):
    if actor is not None:
        actor.SetVisibility(not bool(actor.GetVisibility()))
        plotter.render()


def show_surface(mol: MolData, N_rays=64, resolution=0.5):
    print(f"\n=== Surface AO: {mol.name} ===")
    if mol.ao_scores is None:
        print(f" Calcolo AO ({N_rays} raggi)...")
        mol.ao_scores = compute_ao(mol.coords, mol.radii, N_rays=N_rays)
    ao = mol.ao_scores
    print(f" AO min={ao.min():.3f} max={ao.max():.3f} mean={ao.mean():.3f}")
    print(f" Costruzione superficie ({resolution} A)...")
    surf = _build_surface(mol.coords, mol.radii, resolution)
    surf.point_data["AO"] = _interp_ao(surf, mol.coords, ao)
    print(f" {surf.n_points} vertici, {surf.n_cells} celle")

    pl = pv.Plotter(window_size=(1000, 800), title=f"Surface AO - {mol.name}")
    pl.set_background("#0f1117")
    actors = {}

    actors["surface"] = pl.add_mesh(
        surf,
        scalars="AO",
        cmap="RdBu_r",
        clim=(0, 1),
        smooth_shading=True,
        scalar_bar_args={
            "title": "Occlusion\n(0=esposto 1=cavita)",
            "color": "white",
            "title_font_size": 12,
            "label_font_size": 10,
        },
    )

    if len(mol.ligand_coords):
        actors["ligand"] = pl.add_mesh(
            pv.PolyData(mol.ligand_coords.astype(np.float32)),
            color="yellow",
            point_size=14,
            render_points_as_spheres=True,
            opacity=0.95,
            label="Ligand",
        )

    reps = np.asarray(mol.representative_points, dtype=np.float32) if len(mol.representative_points) else np.empty((0, 3), dtype=np.float32)
    actors["rep_points"] = None
    if len(reps):
        actors["rep_points"] = pl.add_mesh(
            pv.PolyData(reps),
            color="#ff4d6d",
            point_size=13,
            render_points_as_spheres=True,
            opacity=1.0,
            label="Representative points",
        )

    pl.add_key_event("1", lambda: _toggle(actors.get("surface"), pl))
    pl.add_key_event("2", lambda: _toggle(actors.get("ligand"), pl))
    pl.add_key_event("3", lambda: _toggle(actors.get("rep_points"), pl))
    pl.add_key_event("0", pl.reset_camera)

    pl.add_text(
        "1 surface | 2 ligand | 3 rep points | 0 reset",
        position="upper_left",
        font_size=10,
        color="white",
    )
    pl.show()
