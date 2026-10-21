#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_analysis.py

This script runs the pocket detection analysis on PDB files,
either in single (bounded) mode or paired mode (if a correspondence file is provided).
Results (in JSON and CSV formats) are saved in the OUTPUT_DIR.
"""

import glob
import json
import logging
import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd
from Bio.PDB import Superimposer
from Bio.PDB.Atom import Atom
from biopandas.pdb import PandasPdb
from scipy.spatial import Delaunay, KDTree
from scipy.spatial.distance import cdist
from tabulate import tabulate

from . import config as cfg

# Core modules
from .pocket_detection import detect_pockets, compute_analytic_pocket_volume
from .alpha_shape import compute_alpha_complex_from_tetrahedra
from .cgal_wdelaunay import weighted_delaunay
from .alpha_stability import alpha_stability
from .alpha_diagnostics import run_alpha_diagnostics
from .atomic_radii import ATOMIC_RADII_DICT
from .results_analysis import save_summary_json, save_summary_csv

# Utility modules
from .utils.preprocessing_utils import load_and_separate_pdb, calculate_atomic_radii
from .utils.geometry_utils import (
    validate_ligand_in_extruded_mesh,
    validate_ligand_in_fake_sphere,
    compute_all_circumsphere_radii,
)

from .utils.geometry_utils import (
    validate_ligand_in_extruded_mesh,
    validate_ligand_in_fake_sphere,
)
# from .utils.sasa_utils import (
#     compute_ligand_contact_sasa,
#     evaluate_sasa_contact,
# )
# Import commentati - funzioni non implementate
from .utils.pocket_utils import (
    apply_splitting_to_shared_atoms,
    compute_volume_for_tetra_group,
    merge_nearby_clusters,
    compute_mouth_parameters,
    save_pocket_properties,
)

from .utils.data_utils import (
    save_pocket_results_csv,
    save_json,
    load_json,
    convert_numpy_types,
)

from .ranking import compute_ranking_scores
from .pipeline.feature_pipeline import (
    build_pocket_candidates,
    enrich_with_geometry,
    enrich_with_mouth,
  
)
from .pipeline.scoring_pipeline import score_candidates

print("Working directory:", os.getcwd())
print("DEBUG: Looking for PDBs in BOUNDED_DIR:", os.path.abspath(cfg.BOUNDED_DIR))
if os.path.isdir(cfg.BOUNDED_DIR):
    print("DEBUG: Files found in directory:", os.listdir(cfg.BOUNDED_DIR))
else:
    print("DEBUG: Directory does not exist:", cfg.BOUNDED_DIR)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logging.info(f"🔍 Using Validation Method: {cfg.VALIDATION_METHOD}")


def compute_pocket_distances(pocket_points, ligand_coords):
    if ligand_coords.size == 0 or pocket_points.size == 0:
        logging.warning("Ligand or pocket points are missing!")
        return {"min": None, "max": None, "avg": None}

    distances = cdist(pocket_points, ligand_coords)
    return {
        "min": np.min(distances),
        "max": np.max(distances),
        "avg": np.mean(distances),
    }


class DummyAtom:
    def __init__(self, coord):
        self.coord = np.array(coord, dtype=float)

    def get_coord(self):
        return self.coord


def match_atoms_by_residue_and_name(
    fixed_atoms,
    moving_atoms,
    allowed_atom_names=["CA"],
):
    fixed_filtered = [a for a in fixed_atoms if a["atom_name"] in allowed_atom_names]
    moving_filtered = [a for a in moving_atoms if a["atom_name"] in allowed_atom_names]

    fixed_dict = {(a["resid"], a["atom_name"]): a["coord"] for a in fixed_filtered}
    moving_dict = {(a["resid"], a["atom_name"]): a["coord"] for a in moving_filtered}

    common_keys = sorted(set(fixed_dict.keys()) & set(moving_dict.keys()))
    if len(common_keys) == 0:
        logging.warning("No matching atoms found between fixed and moving sets.")
        return [], []

    fixed_coords = [fixed_dict[k] for k in common_keys]
    moving_coords = [moving_dict[k] for k in common_keys]

    logging.info(
        f"Matched {len(common_keys)} atoms for alignment based on {allowed_atom_names}."
    )
    return fixed_coords, moving_coords


def make_dummy_atom(coord):
    return Atom(
        name="X",
        coord=coord,
        bfactor=0.0,
        occupancy=1.0,
        altloc=" ",
        fullname=" X ",
        serial_number=1,
        element="X",
    )


def align_using_superimposer(fixed_coords, moving_coords, return_transform=False):
    fixed_coords = np.array(fixed_coords)
    moving_coords = np.array(moving_coords)

    n = min(len(fixed_coords), len(moving_coords))
    if len(fixed_coords) != len(moving_coords):
        logging.warning(
            f"Different number of atoms: fixed {len(fixed_coords)} vs moving {len(moving_coords)}. Using the first {n} atoms."
        )
        fixed_coords = fixed_coords[:n]
        moving_coords = moving_coords[:n]

    fixed_atoms = [make_dummy_atom(coord) for coord in fixed_coords]
    moving_atoms = [make_dummy_atom(coord) for coord in moving_coords]

    sup = Superimposer()
    sup.set_atoms(fixed_atoms, moving_atoms)

    R, t = sup.rotran
    aligned = np.dot(moving_coords, R) + t
    rmsd = sup.rms

    logging.info(f"Alignment RMSD: {rmsd:.3f} Å")

    if return_transform:
        return aligned, rmsd, (R, t)
    return aligned, rmsd


def analyze_top_pockets(result):
    """
    Determines the first validated pocket based on the chosen validation method.
    Returns a tuple: (top_index, validation_tags)
    where top_index is the 1-based index of the first pocket that validated (or None),
    and validation_tags is a list (one per pocket) indicating the chosen method if validated, or "None".
    This function checks the unified field 'ligand_validation'.
    """
    method = cfg.VALIDATION_METHOD.lower()
    validations = result.get("ligand_validation", [])

    top_index = None
    validation_tags = []

    for idx, valid in enumerate(validations):
        if valid:
            validation_tags.append(method.capitalize())
            if top_index is None:
                top_index = idx + 1
        else:
            validation_tags.append("None")

    return top_index, validation_tags


def run_all_analysis():
    out_dir = os.path.join(str(cfg.OUTPUT_DIR), cfg.VERSION_TAG)
    os.makedirs(out_dir, exist_ok=True)

    bounded_files = glob.glob(os.path.join(str(cfg.BOUNDED_DIR), "*.pdb"))
    if not bounded_files:
        raise FileNotFoundError("No PDB files found in BOUNDED_DIR.")

    correspondence_exists = (
        cfg.CORRESPONDENCE_FILE is not None
        and os.path.exists(str(cfg.CORRESPONDENCE_FILE))
    )
    unbounded_exists = (
        cfg.UNBOUNDED_DIR is not None
        and os.path.exists(str(cfg.UNBOUNDED_DIR))
    )

    if unbounded_exists and correspondence_exists:
        results = run_analysis_mode(
            pdb_sources=bounded_files,
            analysis_type="paired",
            correspondence_df=pd.read_excel(cfg.CORRESPONDENCE_FILE),
        )
        summary, _ = summarize_results(results, analysis_type="paired")
        return summary

    elif unbounded_exists and not correspondence_exists:
        unbounded_files = glob.glob(os.path.join(str(cfg.UNBOUNDED_DIR), "*.pdb"))
        if not unbounded_files:
            raise FileNotFoundError("No PDB files found in UNBOUNDED_DIR.")
        results = run_analysis_mode(
            pdb_sources=unbounded_files,
            analysis_type="unbounded",
        )
        summary, _ = summarize_results(results, analysis_type="unbounded")
        return summary

    else:
        results = run_analysis_mode(
            pdb_sources=bounded_files,
            analysis_type="single",
        )
        summary, _ = summarize_results(results, analysis_type="single")
        return summary

def save_triangulation_debug_data(
    simplices,
    tetra_positions,
    protein_coords,
    molecule_name,
    circumsphere_radii=None,
):
    try:
        out_dir = os.path.join(cfg.OUTPUT_DIR, cfg.VERSION_TAG)
        os.makedirs(out_dir, exist_ok=True)

        out_file = os.path.join(out_dir, f"{molecule_name}.triangulation_debug.npz")

        if circumsphere_radii is None:
            circumsphere_radii = compute_all_circumsphere_radii(tetra_positions)

        np.savez(
            out_file,
            simplices=np.asarray(simplices, dtype=np.int32),
            tetrahedra=np.asarray(tetra_positions, dtype=float),
            protein_coords=np.asarray(protein_coords, dtype=float),
            circumsphere_radii=np.asarray(circumsphere_radii, dtype=float),
        )

        logging.info(f"Triangulation debug data saved to {out_file}")

    except Exception as exc:
        logging.warning(f"Could not save triangulation debug data for {molecule_name}: {exc}")

# ── Core PDB processing ──────────────────────────────────────────────────────
def process_pdb(pdb_path):
    logging.info(f"=== Processing: {pdb_path} ===")
    base = os.path.splitext(os.path.basename(pdb_path))[0]
    out_dir = os.path.join(cfg.OUTPUT_DIR, cfg.VERSION_TAG)
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, f"{base}_pockets.json")

    result = {
        "pdb_path": pdb_path,
        "triangulation_backend": None,
        "ranked_pockets": [],
        "ranking_scores": [],
        "ligand_coords": [],
        "protein_coords": [],
        "protein_atoms": [],
        "step_to_ligand": None,
        "step_to_ligand_mesh": None,
        "validation_methods": [],
        "ligand_validation": [],
        "pocket_volumes": [],
        "pocket_depths": [],
        "mouth_area": [],
        "mouth_perimeter": [],
        "ligand_to_pocket_distances": [],
        "ligand_containment_mesh": [],
        "ligand_mesh_distances": [],
        "ligand_containment_strict": [],
    }

    # Cache hit
    if cfg.USE_EXISTING_RESULTS and os.path.exists(out_file):
        cached = load_json(out_file)
        if cached:
            logging.info(f"Using cached result: {out_file}")
            return cached
        logging.warning("Cached result invalid — recalculating.")

    # Load PDB
    try:
        protein_coords, ligand_coords, protein_atoms, ppdb = load_and_separate_pdb(pdb_path)
        radii = calculate_atomic_radii(ppdb)
    except Exception as e:
        logging.error(f"Failed to load {pdb_path}: {e}")
        return result

    ligand_coords = (
        np.array(ligand_coords)
        if ligand_coords is not None and len(ligand_coords) > 0
        else np.array([])
    )

    from .atomic_radii import ATOMIC_RADII_DICT
    radii = np.array([
        ATOMIC_RADII_DICT.get(atom.get('element', 'C'), ATOMIC_RADII_DICT['C'])
        for atom in protein_atoms
    ], dtype=float)

    
    result["ligand_coords"] = ligand_coords
    result["protein_coords"] = protein_coords.tolist()
    result["protein_atoms"] = protein_atoms

    if np.unique(protein_coords, axis=0).shape[0] < 4:
        logging.error(f"Too few unique atoms in {pdb_path}. Skipping.")
        return result

    sasa_by_residue = None


    # Delaunay triangulation
    logging.info("Running triangulation...")

   
    simplices, triangulation_info = weighted_delaunay(protein_coords, radii)
    logging.info("Using weighted Delaunay via CGAL/pocketlab.")


    simplices = np.asarray(simplices, dtype=int)
    if simplices.ndim != 2 or simplices.shape[1] != 4:
        logging.error(f"Invalid simplices shape: {simplices.shape}")
        return result

    result["triangulation_backend"] = ("weighted_delaunay")
        
    logging.info(f"Triangulation backend: {result['triangulation_backend']}")

    tetra_positions = protein_coords[simplices]
    if getattr(cfg, "SAVE_TRIANGULATION_DEBUG", False):
        save_triangulation_debug_data(
        simplices=simplices,
        tetra_positions=tetra_positions,
        protein_coords=protein_coords,
        molecule_name=base,
        )

    # Alpha shape
    logging.info("Computing alpha shape...")
    alpha = cfg.get_alpha()
    logging.info(f"Alpha = {alpha}")
    alpha_mask, power_radii, tetra_positions, simplices = compute_alpha_complex_from_tetrahedra(
        simplices,
        tetra_positions,
        alpha,
        base,
        protein_coords,
        radii=radii,  # <-- AGGIUNGI QUESTO PARAMETRO!
        filter_by_sasa=cfg.FILTER_ALPHA_BY_SASA,
        pdb_file=pdb_path,
    )

    # Alpha shape
    logging.info("Computing alpha shape...")
    alpha = cfg.get_alpha()
    logging.info(f"Alpha = {alpha}")
    alpha_mask, radii, tetra_positions, simplices = compute_alpha_complex_from_tetrahedra(
        simplices,
        tetra_positions,
        alpha,
        base,
        protein_coords,
        radii=radii,
        filter_by_sasa=cfg.FILTER_ALPHA_BY_SASA,
        pdb_file=pdb_path,
    )

    # Pocket detection
    logging.info("Detecting pockets...")
    pocket_info = detect_pockets(
        protein_coords,
        simplices,
        tetra_positions,
        alpha_mask,
        min_volume_threshold=cfg.get_min_pocket_volume(),
        flow_params=cfg.FLOW_KWARGS,
        molecule_name=base,
        radii=radii,
    )
    raw_ranked_pockets = pocket_info.get("ranked_pockets", [])

    candidates = build_pocket_candidates(raw_ranked_pockets, tetra_positions)
    candidates = enrich_with_geometry(candidates, tetra_positions)
    candidates = enrich_with_mouth(
        candidates,
        ranked_pockets=raw_ranked_pockets,
        simplices=simplices,
        protein_coords=protein_coords,
        ligand_coords=ligand_coords,
    )
   
    candidates = score_candidates(candidates)

    ranked_pockets = [c.tetra_indices for c in candidates]

    # Compute SASA for each pocket using FreeSASA
    # Compute SASA for each pocket using FreeSASA
    from .utils.sasa_utils import compute_pocket_sasa

    logging.info("Computing SASA for pockets...")
    pocket_sasa = []

    try:
        logging.info(f"DEBUG: ranked_pockets type={type(ranked_pockets)}, len={len(ranked_pockets) if ranked_pockets else 'None'}")
        print(f"DEBUG: ranked_pockets = {ranked_pockets}")
        print("ABOUT TO START SASA LOOP")
        
        for i, pocket in enumerate(ranked_pockets):
            print(f"ENTERING POCKET {i}")
            print(f"  - pocket type: {type(pocket)}, len: {len(pocket)}")
            print(f"  - tetra_positions shape: {tetra_positions.shape}")
            print(f"  - protein_coords shape: {protein_coords.shape}")
            print(f"  - pdb_path: {pdb_path}")
            print("  - ABOUT TO CALL compute_pocket_sasa")
            
            try:
                sasa = compute_pocket_sasa(
                    pocket_tetra_indices=pocket,
                    tetra_positions=tetra_positions,
                    protein_coords=protein_coords,
                    pdb_file=pdb_path,
                    probe_radius=cfg.SASA_PROBE_RADIUS,
                )
                print(f"DEBUG: SASA computed for pocket {i+1}: {sasa}")
                pocket_sasa.append(sasa)
                print(f"Pocket {i+1}: SASA = {sasa:.2f} Å²")
            except Exception as e:
                print(f"Failed to compute SASA for pocket {i+1}: {e}")
                import traceback
                traceback.print_exc()
                pocket_sasa.append(0.0)
        
        print("SASA LOOP COMPLETED")
        
    except Exception as e:
        print(f"FATAL ERROR in SASA computation: {e}")
        import traceback
        traceback.print_exc()
        pocket_sasa = [0.0] * len(ranked_pockets)

    result["pocket_sasa"] = pocket_sasa

    try:
        stability_info = alpha_stability(
            protein_coords=protein_coords,
            simplices=simplices,
            tetra_positions=tetra_positions,
            molecule_name=base,
            radii=radii,
            base_alpha=None,
            alpha_offsets=[-0.2, -0.1, 0.0, 0.1, 0.2],
            min_volume_threshold=cfg.get_min_pocket_volume(),
            flow_params=cfg.FLOW_KWARGS,
            merge_clusters=cfg.MERGE_CLUSTERS,
            merge_threshold=cfg.get_merge_threshold(),
        )
        result["alpha_stability"] = stability_info
    except Exception as e:
        logging.warning(f"Alpha stability analysis failed for {base}: {e}")
        result["alpha_stability"] = {
            "stable": False,
            "error": str(e),
        }

    try:
        alpha_diag = run_alpha_diagnostics(
            protein_coords=protein_coords,
            simplices=simplices,
            tetra_positions=tetra_positions,
            molecule_name=base,
            radii=radii,
            pdb_file=pdb_path,
            base_alpha=alpha,
            alpha_offsets=[-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3],
            min_volume_threshold=cfg.get_min_pocket_volume(),
            flow_params=cfg.FLOW_KWARGS,
            merge_clusters=cfg.MERGE_CLUSTERS,
            merge_threshold=cfg.get_merge_threshold(),
            filter_by_sasa=cfg.FILTER_ALPHA_BY_SASA,
            sasa_threshold=cfg.SASA_THRESHOLD,
        )
        result["alpha_diagnostics"] = alpha_diag
    except Exception as e:
        logging.warning(f"Alpha diagnostics failed for {base}: {e}")
        result["alpha_diagnostics"] = {
            "error": str(e),
        }

    result.update(
        {
            "ranked_pockets": ranked_pockets,
            # "representative_points": representative_pts,  # BUG: not defined
            "ranking_scores": [c.scores for c in candidates],
            "pocket_features": [c.features for c in candidates],
            "pocket_scores": [c.scores for c in candidates],
        }
    )

    # Volumes and depths
    pocket_volumes = [c.volume if c.volume is not None else 0.0 for c in candidates]
    pocket_depths = [c.depth if c.depth is not None else 0.0 for c in candidates]

    for idx, candidate in enumerate(candidates, start=1):
        logging.info(
            f"Pocket {idx}: volume={candidate.volume:.2f}, depth={candidate.depth:.2f}, score={candidate.scores.get('global', 0.0):.2f}"
        )

    result["pocket_volumes"] = pocket_volumes
    result["pocket_depths"] = pocket_depths

    # Mouth parameters
    result["mouth_area"] = [c.mouth_area if c.mouth_area is not None else 0 for c in candidates]
    result["mouth_perimeter"] = [c.mouth_perimeter if c.mouth_perimeter is not None else 0 for c in candidates]

    mouth_params = [
        {
            "mouth_area": c.mouth_area,
            "mouth_circumference": c.mouth_perimeter,
            "rim_atoms": c.rim_atoms,
            "mouth_center": c.mouth_center,
            "ligand_distance": c.ligand_distance_to_mouth,
        }
        for c in candidates
    ]

    # Validation
    validation_methods = []
    ligand_containment_mesh = []
    ligand_mesh_distances = []
    ligand_validation = []
    method = cfg.VALIDATION_METHOD.lower()

    if method == "mesh_extrusion":
        for i, mouth in enumerate(mouth_params):
            rim_atoms = mouth.get("rim_atoms", [])
            mouth_coords = np.array(protein_coords)[rim_atoms] if rim_atoms else np.array([])
            in_contact, mesh, min_dist = validate_ligand_in_extruded_mesh(
                mouth_coords,
                ligand_coords,
                extrusion_distance=cfg.MESH_EXTRUSION_DISTANCE,
            )
            validation_methods.append("Mesh" if in_contact else "None")
            ligand_containment_mesh.append(in_contact)
            ligand_mesh_distances.append(min_dist)
            ligand_validation.append(in_contact)

            if cfg.SAVE_EXTRUDED_MESHES and mesh:
                mesh_dir = os.path.join(
                    cfg.OUTPUT_DIR,
                    cfg.VERSION_TAG,
                    "extruded_meshes",
                    base,
                )
                os.makedirs(mesh_dir, exist_ok=True)
                mesh.export(os.path.join(mesh_dir, f"{base}_pocket_{i+1}.ply"))

    elif method == "sasa":
        sasa_before = compute_sasa(pdb_path, cfg.SASA_METHOD)
        sasa_after = compute_sasa(pdb_path, cfg.SASA_METHOD, remove_ligand=True)
        tree = KDTree(protein_coords)

        for pocket in ranked_pockets:
            tetra_pts = np.concatenate([tetra_positions[i] for i in pocket], axis=0)
            _, atom_idx = tree.query(tetra_pts, k=1)
            valid = evaluate_sasa_contact(
                list(set(atom_idx)),
                sasa_before,
                sasa_after,
                ligand_coords,
                protein_coords,
                cfg.SASA_THRESHOLD,
            )
            validation_methods.append("SASA" if valid else "None")
            ligand_validation.append(valid)

    elif method == "fake_ball":
        for mouth in mouth_params:
            rim_atoms = mouth.get("rim_atoms", [])
            mouth_coords = np.array(protein_coords)[rim_atoms] if rim_atoms else np.array([])
            valid = validate_ligand_in_fake_sphere(
                mouth_coords,
                ligand_coords,
                cfg.FAKE_SPHERE_RADIUS,
            )
            validation_methods.append("FakeBall" if valid else "None")
            ligand_validation.append(valid)

    else:
        validation_methods = ["None"] * len(ranked_pockets)
        ligand_validation = [False] * len(ranked_pockets)

    result["validation_methods"] = validation_methods
    result["ligand_validation"] = ligand_validation
    result["ligand_containment_mesh"] = ligand_containment_mesh
    result["ligand_mesh_distances"] = ligand_mesh_distances
    result["step_to_ligand"] = next(
        (i + 1 for i, m in enumerate(validation_methods) if m != "None"),
        None,
    )
    result["step_to_ligand_mesh"] = (
        result["step_to_ligand"] if method == "mesh_extrusion" else None
    )

    # Save pocket PDBs
    protein_df = ppdb.df["ATOM"]
    output_pdb_dir = os.path.join(cfg.OUTPUT_DIR, cfg.VERSION_TAG, "pocket_pdbs", base)
    os.makedirs(output_pdb_dir, exist_ok=True)
    tree = KDTree(protein_coords)

    for rank, pocket in enumerate(ranked_pockets[: cfg.NUM_POCKETS_TO_SAVE], 1):
        pocket_atoms = np.unique(
            np.concatenate([tetra_positions[i] for i in pocket], axis=0),
            axis=0,
        )
        _, idx = tree.query(pocket_atoms, k=1)
        pocket_df = protein_df.iloc[idx].copy()
        pocket_df["b_factor"] = rank
        pdb_out = PandasPdb()
        pdb_out.df["ATOM"] = pocket_df
        pdb_out.to_pdb(os.path.join(output_pdb_dir, f"{base}_ranked_{rank}.pdb"))
        logging.info(f"Pocket {rank} saved.")

    save_pocket_properties(result, base, cfg.OUTPUT_DIR)
    logging.info(
        f"=== Done: {pdb_path} | pockets={len(ranked_pockets)}, top={result['step_to_ligand']} ==="
    )
    return convert_numpy_types(result)


# ── Analysis modes ────────────────────────────────────────────────────────────


def run_analysis_mode(pdb_sources, analysis_type="single", correspondence_df=None):
    results = []
    out_dir = os.path.join(cfg.OUTPUT_DIR, cfg.VERSION_TAG)
    os.makedirs(out_dir, exist_ok=True)

    if analysis_type == "paired" and correspondence_df is not None:
        for _, row in correspondence_df.iterrows():
            b_name = row["bounded"]
            u_name = row["unbounded"]
            b_path = os.path.join(cfg.BOUNDED_DIR, f"{b_name}.pdb")
            u_path = os.path.join(cfg.UNBOUNDED_DIR, f"{u_name}.pdb")

            if not os.path.exists(b_path) or not os.path.exists(u_path):
                logging.warning(f"Skipping pair {b_name}/{u_name}: missing file(s).")
                continue

            b_res = process_pdb(b_path)
            if not b_res or not b_res.get("ranked_pockets"):
                logging.warning(f"Skipping {b_name}: no pockets.")
                continue

            u_res = process_pdb(u_path)
            if not u_res or not u_res.get("ranked_pockets"):
                logging.warning(f"Skipping {u_name}: no pockets.")
                continue

            ALIGN_ATOMS = ["CA"]
            fixed_c, moving_c = match_atoms_by_residue_and_name(
                b_res["protein_atoms"],
                u_res["protein_atoms"],
                allowed_atom_names=ALIGN_ATOMS,
            )
            _, rmsd, (R, t) = align_using_superimposer(
                fixed_c,
                moving_c,
                return_transform=True,
            )
            u_res["protein_coords"] = (
                np.dot(np.array(u_res["protein_coords"]), R) + t
            ).tolist()
            logging.info(f"Pair {b_name}/{u_name} RMSD={rmsd:.3f} Å")

            b_rank, _ = analyze_top_pockets(b_res)
            u_rank, _ = analyze_top_pockets(u_res)

            results.append(
                {
                    "bounded_file": b_name,
                    "unbounded_file": u_name,
                    "bounded_pocket_rank": b_rank,
                    "unbounded_pocket_rank": u_rank,
                    "is_top1_bounded": b_rank == 1 if b_rank else False,
                    "is_top3_bounded": b_rank is not None and b_rank <= 3,
                    "is_top5_bounded": b_rank is not None and b_rank <= 5,
                    "is_top1_unbounded": u_rank == 1 if u_rank else False,
                    "is_top3_unbounded": u_rank is not None and u_rank <= 3,
                    "is_top5_unbounded": u_rank is not None and u_rank <= 5,
                    "alignment_rmsd": rmsd,
                    "step_to_ligand": b_res.get("step_to_ligand"),
                    "step_to_ligand_mesh": b_res.get("step_to_ligand_mesh"),
                    "ligand_to_pocket_distances": b_res.get("ligand_to_pocket_distances", []),
                    "ligand_mesh_distances": b_res.get("ligand_mesh_distances", []),
                    "ligand_containment_mesh": b_res.get("ligand_containment_mesh", []),
                    "ligand_containment_strict": b_res.get("ligand_containment_strict", []),
                    "pocket_volumes": b_res.get("pocket_volumes", []),
                    "pocket_depths": b_res.get("pocket_depths", []),
                    "mouth_area": b_res.get("mouth_area", []),
                    "mouth_perimeter": b_res.get("mouth_perimeter", []),
                    "ligand_coords": b_res.get("ligand_coords", []),
                }
            )

            save_json(b_res, os.path.join(out_dir, f"bounded_{b_name}.json"))
            save_json(u_res, os.path.join(out_dir, f"unbounded_{u_name}.json"))

    else:
        for pdb_file in pdb_sources:
            logging.info(f"Processing {pdb_file}")
            result = process_pdb(pdb_file)

            if not result or not result.get("ranked_pockets"):
                logging.warning(f"Skipping {pdb_file}: no pockets.")
                continue

            pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
            save_json(result, os.path.join(out_dir, f"{pdb_name}_pockets.json"))

            rank_sasa = result.get("step_to_ligand")
            rank_mesh = result.get("step_to_ligand_mesh")
            chosen = rank_mesh if cfg.VALIDATION_METHOD == "mesh_extrusion" else rank_sasa

            results.append(
                {
                    "pdb_path": pdb_file,
                    "file": os.path.basename(pdb_file),
                    "pocket_rank": chosen,
                    "step_to_ligand": rank_sasa,
                    "step_to_ligand_mesh": rank_mesh,
                    "is_top1": chosen == 1 if isinstance(chosen, int) else False,
                    "is_top3": chosen <= 3 if isinstance(chosen, int) else False,
                    "is_top5": chosen <= 5 if isinstance(chosen, int) else False,
                    "pocket_volumes": result.get("pocket_volumes", []),
                    "pocket_depths": result.get("pocket_depths", []),
                    "mouth_area": result.get("mouth_area", []),
                    "mouth_perimeter": result.get("mouth_perimeter", []),
                    "ligand_coords": result.get("ligand_coords", []),
                    "ligand_to_pocket_distances": result.get("ligand_to_pocket_distances", []),
                    "ligand_mesh_distances": result.get("ligand_mesh_distances", []),
                    "ligand_containment_mesh": result.get("ligand_containment_mesh", []),
                    "ligand_containment_strict": result.get("ligand_containment_strict", []),
                    "num_pockets": len(result.get("ranked_pockets", [])),
                }
            )

    return results


# ── Entry point ───────────────────────────────────────────────────────────────


def main():
    logging.info(f"=== PyCASTA analysis start | version={cfg.VERSION_TAG} ===")
    out_dir = os.path.join(str(cfg.OUTPUT_DIR), cfg.VERSION_TAG)
    os.makedirs(out_dir, exist_ok=True)

    bounded_files = glob.glob(os.path.join(str(cfg.BOUNDED_DIR), "*.pdb"))
    if not bounded_files:
        logging.error("No PDB files found in BOUNDED_DIR.")
        sys.exit(1)

    correspondence_exists = (
        cfg.CORRESPONDENCE_FILE is not None
        and os.path.exists(str(cfg.CORRESPONDENCE_FILE))
    )
    unbounded_exists = (
        cfg.UNBOUNDED_DIR is not None
        and os.path.exists(str(cfg.UNBOUNDED_DIR))
    )

    if unbounded_exists and correspondence_exists:
        logging.info("Mode: paired")
        results = run_analysis_mode(
            pdb_sources=bounded_files,
            analysis_type="paired",
            correspondence_df=pd.read_excel(cfg.CORRESPONDENCE_FILE),
        )
        save_summary_json(
            results,
            os.path.join(out_dir, "paired_summary.json"),
            analysis_type="paired",
        )
        save_summary_csv(
            results,
            os.path.join(out_dir, "paired_summary.csv"),
            analysis_type="paired",
        )

    elif unbounded_exists and not correspondence_exists:
        logging.info("Mode: unbounded")
        unbounded_files = glob.glob(os.path.join(str(cfg.UNBOUNDED_DIR), "*.pdb"))
        if not unbounded_files:
            logging.error("No PDB files found in UNBOUNDED_DIR.")
            sys.exit(1)

        results = run_analysis_mode(
            pdb_sources=unbounded_files,
            analysis_type="unbounded",
        )
        save_summary_json(
            results,
            os.path.join(out_dir, "unbounded_summary.json"),
            analysis_type="unbounded",
        )
        save_summary_csv(
            results,
            os.path.join(out_dir, "unbounded_summary.csv"),
            analysis_type="unbounded",
        )

    else:
        logging.info("Mode: single (bounded)")
        results = run_analysis_mode(
            pdb_sources=bounded_files,
            analysis_type="single",
        )
        save_summary_json(
            results,
            os.path.join(out_dir, "single_summary.json"),
            analysis_type="single",
        )
        save_summary_csv(
            results,
            os.path.join(out_dir, "single_summary.csv"),
            analysis_type="single",
        )

    logging.info(f"=== Analysis complete. Results in: {out_dir} ===")

if __name__ == "__main__":
    main()
