# -*- coding: utf-8 -*-
"""visualization/core.py"""
import json, os, numpy as np
from dataclasses import dataclass
from typing import List, Optional, Dict
from pathlib import Path

VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52,
    "S": 1.80, "P": 1.80, "F": 1.47, "CL": 1.75,
    "BR": 1.85, "I": 1.98, "FE": 1.40, "ZN": 1.39,
    "MG": 1.73, "CA": 1.97, "MN": 1.61, "CU": 1.40,
    "SE": 1.90, "K": 2.75, "NA": 2.27,
}
DEFAULT_VDW = 1.50


def get_vdw(atom_name: str) -> float:
    elem = "".join(c for c in atom_name if c.isalpha()).upper()
    return VDW_RADII.get(elem, VDW_RADII.get(elem[:1], DEFAULT_VDW))


@dataclass
class MolData:
    name: str
    pdb_path: str
    coords: np.ndarray
    radii: np.ndarray
    atom_names: List[str]
    ligand_coords: np.ndarray
    tetrahedra: np.ndarray
    alpha_mask: np.ndarray
    ranked_pockets: List[List[int]]
    representative_points: np.ndarray
    pocket_volumes: List[float]
    ligand_validation: List[bool]
    ranking_scores: List[float]
    alpha_used: float = 2.0
    ao_scores: Optional[np.ndarray] = None
    alpha_npz_path: Optional[str] = None
    triangulation_simplices: np.ndarray | None = None
    triangulation_tetrahedra: np.ndarray | None = None
    triangulation_radii: np.ndarray | None = None
    triangulation_npz_path: Optional[str] = None


def _flat_coords(raw) -> np.ndarray:
    if not raw:
        return np.zeros((0, 3))
    if isinstance(raw[0], (list, tuple)):
        return np.array(raw, dtype=float)
    return np.array(raw, dtype=float).reshape(-1, 3)


def _find_alpha_npz(json_path: str, pdb_path: str) -> Optional[str]:
    name = Path(pdb_path).stem
    json_dir = Path(json_path).parent
    candidates = [
        json_dir / f"{name}.alpha_filtered.npz",
        json_dir / f"{name}.alpha.npz",
    ]
    for p in candidates:
        if p.exists():
            return str(p)
    return None

def _find_triangulation_npz(json_path: str, pdb_path: str) -> Optional[str]:
    name = Path(pdb_path).stem
    json_dir = Path(json_path).parent
    candidate = json_dir / f"{name}.triangulation_debug.npz"
    if candidate.exists():
        return str(candidate)
    return None


def load_mol(pdb_path: str, json_path: str) -> MolData:
    with open(json_path, encoding="utf-8") as f:
        d = json.load(f)

    raw = d.get("protein_atoms", [])
    if raw and isinstance(raw[0], dict):
        coords = np.array([a["coord"] for a in raw], dtype=float)
        atom_names = [a.get("atom_name", "C") for a in raw]
    else:
        coords = _flat_coords(d.get("protein_coords", []))
        atom_names = ["C"] * len(coords)

    radii = np.array([get_vdw(n) for n in atom_names], dtype=float)

    ligand_coords = _flat_coords(d.get("ligand_coords", []))
    representative_points = _flat_coords(d.get("representative_points", []))

    tetrahedra = np.zeros((0, 4, 3), dtype=float)
    alpha_mask = np.zeros((0,), dtype=bool)
    alpha_npz_path = _find_alpha_npz(json_path, pdb_path)
    alpha_used: Optional[float] = d.get("alpha_used", None)

    if alpha_npz_path:
        try:
            z = np.load(alpha_npz_path, allow_pickle=True)

            if "tetrahedra" in z:
                tetrahedra = np.array(z["tetrahedra"], dtype=float)

            if "alpha_mask" in z:
                alpha_mask = np.array(z["alpha_mask"], dtype=bool)

            if "alpha_value" in z:
                alpha_used = float(z["alpha_value"])

            print(f"[core] alpha NPZ loaded: {alpha_npz_path}")
            print(f"[core] tetrahedra={len(tetrahedra)} alpha_mask={len(alpha_mask)}")
        except Exception as e:
            print(f"[core] WARNING: failed to load alpha NPZ: {e}")
    else:
        print("[core] WARNING: no .alpha.npz found next to JSON")

    triangulation_simplices = None
    triangulation_tetrahedra = None
    triangulation_radii = None
    triangulation_npz_path = _find_triangulation_npz(json_path, pdb_path)

    if triangulation_npz_path:
        try:
            ztri = np.load(triangulation_npz_path, allow_pickle=True)

            if "simplices" in ztri:
                triangulation_simplices = np.array(ztri["simplices"], dtype=np.int32)

            if "tetrahedra" in ztri:
                triangulation_tetrahedra = np.array(ztri["tetrahedra"], dtype=float)

            if "circumsphere_radii" in ztri:
                triangulation_radii = np.array(ztri["circumsphere_radii"], dtype=float)

            print(f"[core] triangulation NPZ loaded: {triangulation_npz_path}")
            print(
                "[core] triangulation tetrahedra="
                f"{0 if triangulation_tetrahedra is None else len(triangulation_tetrahedra)}"
            )
        except Exception as e:
            print(f"[core] WARNING: failed to load triangulation NPZ: {e}")
    else:
        print("[core] WARNING: no .triangulation_debug.npz found next to JSON")

    return MolData(
        name=Path(pdb_path).stem,
        pdb_path=pdb_path,
        coords=coords,
        radii=radii,
        atom_names=atom_names,
        ligand_coords=ligand_coords,
        tetrahedra=tetrahedra,
        alpha_mask=alpha_mask,
        ranked_pockets=d.get("ranked_pockets", []),
        representative_points=representative_points,
        pocket_volumes=d.get("pocket_volumes", []),
        ligand_validation=d.get("ligand_validation", []),
        ranking_scores=d.get("ranking_scores", []),
        alpha_used=alpha_used,
        alpha_npz_path=alpha_npz_path,
        triangulation_simplices=triangulation_simplices,
        triangulation_tetrahedra=triangulation_tetrahedra,
        triangulation_radii=triangulation_radii,
        triangulation_npz_path=triangulation_npz_path,
    )


def guess_alpha_values(mol: MolData) -> Dict[str, float]:
    from scipy.spatial import KDTree

    r_mean = float(np.mean(mol.radii))
    tree = KDTree(mol.coords)
    dists, _ = tree.query(mol.coords, k=5)
    d_nn = float(np.mean(dists[:, 1:5]))
    centroid = mol.coords.mean(axis=0)
    rg = float(np.sqrt(np.mean(np.sum((mol.coords - centroid) ** 2, axis=1))))

    return {
        "tight": round(r_mean * 1.00, 2),
        "medium": round(r_mean * 1.40, 2),
        "loose": round(r_mean * 1.95, 2),
        "_r_mean": r_mean,
        "_d_nn": d_nn,
        "_rg": rg,
    }