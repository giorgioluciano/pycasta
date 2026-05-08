#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_analysis.py

Pipeline principale di pocket detection su file PDB.
Supporta tre modalità: single (bounded), unbounded, paired.
Risultati salvati in OUTPUT_DIR/<VERSION_TAG>/ come JSON e CSV.
"""

import glob
import logging
import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd
from Bio.PDB import Superimposer
from Bio.PDB.Atom import Atom
from scipy.spatial import Delaunay, KDTree
from scipy.spatial.distance import cdist

from config import (
    BOUNDED_DIR,
    CORRESPONDENCE_FILE,
    DATASET_DIR,
    FAKE_SPHERE_RADIUS,
    FILTER_ALPHA_BY_SASA,
    FLOW_KWARGS,
    INCLUDE_WATER,
    LIGAND_CONTACT_METHOD,
    LIGAND_CONTACT_THRESHOLD,
    MERGE_CLUSTERS,
    MERGE_THRESHOLD,
    MESH_EXTRUSION_DISTANCE,
    MIN_POCKET_VOLUME,
    NUM_POCKETS_TO_SAVE,
    NUM_POCKETS_TO_VALIDATE,
    OUTPUT_DIR,
    SASA_CONTACT_THRESHOLD,
    SASA_METHOD,
    SASA_THRESHOLD,
    SAVE_EXTRUDED_MESHES,
    SIMPLE_SPLIT,
    STRICT_DISTANCE_THRESHOLD,
    UNBOUNDED_DIR,
    USE_CGAL,
    USE_EXISTING_RESULTS,
    USE_SASA,
    USE_SASA_CONTACT_VALIDATION,
    VALIDATION_METHOD,
    VERSION_TAG,
    VOLUME_METHOD,
    get_alpha,
)

from alpha_shape import compute_alpha_complex_from_tetrahedra
from cgal_wdelaunay import cgal_weighted_delaunay
from pocket_detection import compute_analytic_pocket_volume, detect_pockets
from ranking import compute_ranking_scores

from utils.data_utils import (
    convert_numpy_types,
    load_json,
    save_json,
    save_pocket_results_csv,
    summarize_and_print_results,
)
from utils.geometry_utils import (
    validate_ligand_in_extruded_mesh,
    validate_ligand_in_fake_sphere,
)
from utils.pocket_utils import (
    apply_splitting_to_shared_atoms,
    compute_mouth_parameters,
    compute_volume_for_tetra_group,
    merge_nearby_clusters,
    save_pocket_properties,
)
from utils.preprocessing_utils import calculate_atomic_radii, load_and_separate_pdb
from utils.sasa_utils import compute_sasa, compute_ligand_contact_sasa, evaluate_sasa_contact

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logging.info(f"Validation method: {VALIDATION_METHOD}")


# ── Alignment helpers ────────────────────────────────────────────────────────

def _make_atom(coord):
    return Atom(
        name="X", coord=coord, bfactor=0.0, occupancy=1.0,
        altloc=" ", fullname=" X ", serial_number=1, element="X",
    )


def match_atoms_by_residue_and_name(fixed_atoms, moving_atoms, allowed_atom_names=("CA",)):
    fixed_dict  = {(a["resid"], a["atom_name"]): a["coord"] for a in fixed_atoms  if a["atom_name"] in allowed_atom_names}
    moving_dict = {(a["resid"], a["atom_name"]): a["coord"] for a in moving_atoms if a["atom_name"] in allowed_atom_names}
    common = sorted(fixed_dict.keys() & moving_dict.keys())
    if not common:
        logging.warning("No matching atoms found for alignment.")
        return [], []
    logging.info(f"Matched {len(common)} atoms for alignment ({list(allowed_atom_names)}).")
    return [fixed_dict[k] for k in common], [moving_dict[k] for k in common]


def align_using_superimposer(fixed_coords, moving_coords, return_transform=False):
    fixed_coords  = np.array(fixed_coords)
    moving_coords = np.array(moving_coords)
    n = min(len(fixed_coords), len(moving_coords))
    if len(fixed_coords) != len(moving_coords):
        logging.warning(f"Atom count mismatch ({len(fixed_coords)} vs {len(moving_coords)}), using first {n}.")
        fixed_coords, moving_coords = fixed_coords[:n], moving_coords[:n]
    sup = Superimposer()
    sup.set_atoms([_make_atom(c) for c in fixed_coords], [_make_atom(c) for c in moving_coords])
    R, t = sup.rotran
    aligned = np.dot(moving_coords, R) + t
    logging.info(f"Alignment RMSD: {sup.rms:.3f} Å")
    if return_transform:
        return aligned, sup.rms, (R, t)
    return aligned, sup.rms


# ── Pocket analysis helpers ───────────────────────────────────────────────────

def compute_pocket_distances(pocket_points, ligand_coords):
    if ligand_coords.size == 0 or pocket_points.size == 0:
        logging.warning("Ligand or pocket points missing — skipping distance computation.")
        return {"min": None, "max": None, "avg": None}
    d = cdist(pocket_points, ligand_coords)
    return {"min": float(np.min(d)), "max": float(np.max(d)), "avg": float(np.mean(d))}


def analyze_top_pockets(result):
    """Return (top_1based_index, validation_tags_list) from unified 'ligand_validation' field."""
    validations = result.get("ligand_validation", [])
    top_index, tags = None, []
    for idx, valid in enumerate(validations):
        if valid:
            tags.append(VALIDATION_METHOD.capitalize())
            if top_index is None:
                top_index = idx + 1
        else:
            tags.append("None")
    return top_index, tags


# ── Core PDB processing ──────────────────────────────────────────────────────

def process_pdb(pdb_path):
    logging.info(f"=== Processing: {pdb_path} ===")
    base    = os.path.splitext(os.path.basename(pdb_path))[0]
    out_dir = os.path.join(OUTPUT_DIR, VERSION_TAG)
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, f"{base}_pockets.json")

    result = {
        "pdb_path": pdb_path,
        "ranked_pockets": [], "ranking_scores": [],
        "ligand_coords": [], "protein_coords": [], "protein_atoms": [],
        "step_to_ligand": None, "step_to_ligand_mesh": None,
        "validation_methods": [], "ligand_validation": [],
        "pocket_volumes": [], "pocket_depths": [],
        "mouth_area": [], "mouth_perimeter": [],
        "ligand_to_pocket_distances": [],
        "ligand_containment_mesh": [], "ligand_mesh_distances": [],
        "ligand_containment_strict": [],
    }

    # Cache hit
    if USE_EXISTING_RESULTS and os.path.exists(out_file):
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

    ligand_coords = np.array(ligand_coords) if ligand_coords is not None and len(ligand_coords) > 0 else np.array([])
    result["ligand_coords"]  = ligand_coords
    result["protein_coords"] = protein_coords.tolist()
    result["protein_atoms"]  = protein_atoms

    if np.unique(protein_coords, axis=0).shape[0] < 4:
        logging.error(f"Too few unique atoms in {pdb_path}. Skipping.")
        return result

    # Delaunay triangulation
    logging.info("Running Delaunay triangulation...")
    try:
        simplices = (
            cgal_weighted_delaunay(protein_coords, radii)[0]
            if USE_CGAL
            else Delaunay(protein_coords).simplices
        )
    except Exception as e:
        logging.error(f"Delaunay failed for {pdb_path}: {e}")
        return result
    tetra_positions = protein_coords[simplices]

    # Alpha shape
    logging.info("Computing alpha shape...")
    alpha = get_alpha()
    logging.info(f"Alpha = {alpha}")
    alpha_mask, radii, tetra_positions, simplices = compute_alpha_complex_from_tetrahedra(
        simplices, tetra_positions, alpha, base,
        protein_coords, FILTER_ALPHA_BY_SASA, SASA_CONTACT_THRESHOLD, pdb_path,
    )

    # Pocket detection
    logging.info("Detecting pockets...")
    pocket_info = detect_pockets(
        protein_coords, simplices, tetra_positions, alpha_mask,
        min_volume_threshold=MIN_POCKET_VOLUME,
        flow_params=FLOW_KWARGS,
        molecule_name=base,
        radii=radii,
    )
    ranked_pockets      = pocket_info.get("ranked_pockets", [])
    representative_pts  = pocket_info.get("representative_points", [])
    ranking_scores      = pocket_info.get("ranking_scores", [])
    result.update({"ranked_pockets": ranked_pockets, "representative_points": representative_pts, "ranking_scores": ranking_scores})

    # Volumes and depths
    pocket_volumes, pocket_depths = [], []
    for idx, pocket in enumerate(ranked_pockets):
        valid_idx = [t for t in pocket if 0 <= t < len(protein_coords)]
        if valid_idx:
            volume = compute_analytic_pocket_volume(tetra_positions, valid_idx)
            depth  = float(np.max([np.linalg.norm(protein_coords[t] - protein_coords[valid_idx[0]]) for t in valid_idx]))
            logging.info(f"Pocket {idx+1}: volume={volume:.2f}, depth={depth:.2f}")
        else:
            volume, depth = 0.0, 0.0
            logging.warning(f"Pocket {idx+1}: no valid tetrahedra.")
        pocket_volumes.append(volume)
        pocket_depths.append(depth)
    result["pocket_volumes"] = pocket_volumes
    result["pocket_depths"]  = pocket_depths

    # Mouth parameters
    mouth_params = compute_mouth_parameters(ranked_pockets, simplices, protein_coords, ligand_coords)
    result["mouth_area"]      = [m.get("mouth_area", 0)         for m in mouth_params]
    result["mouth_perimeter"] = [m.get("mouth_circumference", 0) for m in mouth_params]

    # Validation
    validation_methods      = []
    ligand_containment_mesh = []
    ligand_mesh_distances   = []
    ligand_validation       = []
    method = VALIDATION_METHOD.lower()

    if method == "mesh_extrusion":
        for i, mouth in enumerate(mouth_params):
            rim_atoms   = mouth.get("rim_atoms", [])
            mouth_coords = np.array(protein_coords)[rim_atoms] if rim_atoms else np.array([])
            in_contact, mesh, min_dist = validate_ligand_in_extruded_mesh(
                mouth_coords, ligand_coords, extrusion_distance=MESH_EXTRUSION_DISTANCE
            )
            validation_methods.append("Mesh" if in_contact else "None")
            ligand_containment_mesh.append(in_contact)
            ligand_mesh_distances.append(min_dist)
            ligand_validation.append(in_contact)
            if SAVE_EXTRUDED_MESHES and mesh:
                mesh_dir = os.path.join(OUTPUT_DIR, VERSION_TAG, "extruded_meshes", base)
                os.makedirs(mesh_dir, exist_ok=True)
                mesh.export(os.path.join(mesh_dir, f"{base}_pocket_{i+1}.ply"))

    elif method == "sasa":
        sasa_before = compute_sasa(pdb_path, SASA_METHOD)
        sasa_after  = compute_sasa(pdb_path, SASA_METHOD, remove_ligand=True)
        tree = KDTree(protein_coords)
        for pocket in ranked_pockets:
            tetra_pts     = np.concatenate([tetra_positions[i] for i in pocket], axis=0)
            _, atom_idx   = tree.query(tetra_pts, k=1)
            valid = evaluate_sasa_contact(
                list(set(atom_idx)), sasa_before, sasa_after,
                ligand_coords, protein_coords, SASA_CONTACT_THRESHOLD,
            )
            validation_methods.append("SASA" if valid else "None")
            ligand_validation.append(valid)

    elif method == "fake_ball":
        for mouth in mouth_params:
            rim_atoms    = mouth.get("rim_atoms", [])
            mouth_coords = np.array(protein_coords)[rim_atoms] if rim_atoms else np.array([])
            valid = validate_ligand_in_fake_sphere(mouth_coords, ligand_coords, FAKE_SPHERE_RADIUS)
            validation_methods.append("FakeBall" if valid else "None")
            ligand_validation.append(valid)

    else:
        validation_methods = ["None"] * len(ranked_pockets)
        ligand_validation  = [False]  * len(ranked_pockets)

    result["validation_methods"]      = validation_methods
    result["ligand_validation"]        = ligand_validation
    result["ligand_containment_mesh"]  = ligand_containment_mesh
    result["ligand_mesh_distances"]    = ligand_mesh_distances
    result["step_to_ligand"]           = next((i + 1 for i, m in enumerate(validation_methods) if m != "None"), None)
    result["step_to_ligand_mesh"]      = result["step_to_ligand"] if method == "mesh_extrusion" else None

    # Save pocket PDBs
    protein_df     = ppdb.df["ATOM"]
    output_pdb_dir = os.path.join(OUTPUT_DIR, VERSION_TAG, "pocket_pdbs", base)
    os.makedirs(output_pdb_dir, exist_ok=True)
    tree = KDTree(protein_coords)
    for rank, pocket in enumerate(ranked_pockets[:NUM_POCKETS_TO_SAVE], 1):
        pocket_atoms = np.unique(np.concatenate([tetra_positions[i] for i in pocket]), axis=0)
        _, idx       = tree.query(pocket_atoms, k=1)
        from biopandas.pdb import PandasPdb
        pocket_df           = protein_df.iloc[idx].copy()
        pocket_df["b_factor"] = rank
        pdb_out             = PandasPdb()
        pdb_out.df["ATOM"]  = pocket_df
        pdb_out.to_pdb(os.path.join(output_pdb_dir, f"{base}_ranked_{rank}.pdb"))
        logging.info(f"Pocket {rank} saved.")

    save_pocket_properties(result, base, OUTPUT_DIR)
    logging.info(f"=== Done: {pdb_path} | pockets={len(ranked_pockets)}, top={result['step_to_ligand']} ===")
    return convert_numpy_types(result)


# ── Analysis modes ────────────────────────────────────────────────────────────

def run_analysis_mode(pdb_sources, analysis_type="single", correspondence_df=None):
    results = []
    out_dir = os.path.join(OUTPUT_DIR, VERSION_TAG)
    os.makedirs(out_dir, exist_ok=True)

    if analysis_type == "paired" and correspondence_df is not None:
        for _, row in correspondence_df.iterrows():
            b_name = row["bounded"]
            u_name = row["unbounded"]
            b_path = os.path.join(BOUNDED_DIR,   f"{b_name}.pdb")
            u_path = os.path.join(UNBOUNDED_DIR,  f"{u_name}.pdb")
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
                b_res["protein_atoms"], u_res["protein_atoms"], allowed_atom_names=ALIGN_ATOMS
            )
            _, rmsd, (R, t) = align_using_superimposer(fixed_c, moving_c, return_transform=True)
            u_res["protein_coords"] = (np.dot(np.array(u_res["protein_coords"]), R) + t).tolist()
            logging.info(f"Pair {b_name}/{u_name} RMSD={rmsd:.3f} Å")

            b_rank, _ = analyze_top_pockets(b_res)
            u_rank, _ = analyze_top_pockets(u_res)
            results.append({
                "bounded_file":          b_name,
                "unbounded_file":        u_name,
                "bounded_pocket_rank":   b_rank,
                "unbounded_pocket_rank": u_rank,
                "is_top1_bounded":       b_rank == 1 if b_rank else False,
                "is_top3_bounded":       b_rank is not None and b_rank <= 3,
                "is_top5_bounded":       b_rank is not None and b_rank <= 5,
                "is_top1_unbounded":     u_rank == 1 if u_rank else False,
                "is_top3_unbounded":     u_rank is not None and u_rank <= 3,
                "is_top5_unbounded":     u_rank is not None and u_rank <= 5,
                "alignment_rmsd":        rmsd,
                "step_to_ligand":        b_res.get("step_to_ligand"),
                "step_to_ligand_mesh":   b_res.get("step_to_ligand_mesh"),
                "ligand_to_pocket_distances": b_res.get("ligand_to_pocket_distances", []),
                "ligand_mesh_distances":      b_res.get("ligand_mesh_distances", []),
                "ligand_containment_mesh":    b_res.get("ligand_containment_mesh", []),
                "ligand_containment_strict":  b_res.get("ligand_containment_strict", []),
                "pocket_volumes":   b_res.get("pocket_volumes", []),
                "pocket_depths":    b_res.get("pocket_depths", []),
                "mouth_area":       b_res.get("mouth_area", []),
                "mouth_perimeter":  b_res.get("mouth_perimeter", []),
                "ligand_coords":    b_res.get("ligand_coords", []),
            })
            save_json(b_res, os.path.join(out_dir, f"bounded_{b_name}.json"))
            save_json(u_res, os.path.join(out_dir, f"unbounded_{u_name}.json"))

    else:
        for pdb_file in pdb_sources:
            logging.info(f"Processing {pdb_file}")
            result = process_pdb(pdb_file)
            if not result or not result.get("ranked_pockets"):
                logging.warning(f"Skipping {pdb_file}: no pockets.")
                continue
            pdb_name  = os.path.splitext(os.path.basename(pdb_file))[0]
            save_json(result, os.path.join(out_dir, f"{pdb_name}_pockets.json"))
            rank_sasa = result.get("step_to_ligand")
            rank_mesh = result.get("step_to_ligand_mesh")
            chosen    = rank_mesh if VALIDATION_METHOD == "mesh_extrusion" else rank_sasa
            results.append({
                "pdb_path":    pdb_file,
                "file":        os.path.basename(pdb_file),
                "pocket_rank": chosen,
                "step_to_ligand":      rank_sasa,
                "step_to_ligand_mesh": rank_mesh,
                "is_top1": chosen == 1 if isinstance(chosen, int) else False,
                "is_top3": chosen <= 3  if isinstance(chosen, int) else False,
                "is_top5": chosen <= 5  if isinstance(chosen, int) else False,
                "pocket_volumes":  result.get("pocket_volumes", []),
                "pocket_depths":   result.get("pocket_depths", []),
                "mouth_areas":     result.get("mouth_area", []),
                "mouth_perimeters":result.get("mouth_perimeter", []),
                "ligand_coords":   result.get("ligand_coords", []),
                "ligand_to_pocket_distances": result.get("ligand_to_pocket_distances", []),
                "num_pockets": len(result.get("ranked_pockets", [])),
            })
    return results


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    logging.info(f"=== PyCASTA analysis start | version={VERSION_TAG} ===")
    out_dir = os.path.join(OUTPUT_DIR, VERSION_TAG)
    os.makedirs(out_dir, exist_ok=True)

    bounded_files = glob.glob(os.path.join(BOUNDED_DIR, "*.pdb"))
    if not bounded_files:
        logging.error("No PDB files found in BOUNDED_DIR.")
        sys.exit(1)

    if UNBOUNDED_DIR and os.path.exists(UNBOUNDED_DIR) and CORRESPONDENCE_FILE and os.path.exists(CORRESPONDENCE_FILE):
        logging.info("Mode: paired")
        results = run_analysis_mode(
            pdb_sources=bounded_files, analysis_type="paired",
            correspondence_df=pd.read_excel(CORRESPONDENCE_FILE),
        )
        summarize_and_print_results(
            results, analysis_type="paired",
            output_file=os.path.join(out_dir, "paired_summary.json"),
            csv_output_file=os.path.join(out_dir, "paired_summary.csv"),
        )
    elif not UNBOUNDED_DIR or not os.path.exists(UNBOUNDED_DIR):
        logging.info("Mode: single (bounded)")
        results = run_analysis_mode(pdb_sources=bounded_files, analysis_type="single")
        summarize_and_print_results(
            results, analysis_type="single",
            output_file=os.path.join(out_dir, "single_summary.json"),
            csv_output_file=os.path.join(out_dir, "single_summary.csv"),
        )
    elif UNBOUNDED_DIR and os.path.exists(UNBOUNDED_DIR) and not CORRESPONDENCE_FILE:
        logging.info("Mode: unbounded")
        unbounded_files = glob.glob(os.path.join(UNBOUNDED_DIR, "*.pdb"))
        if not unbounded_files:
            logging.error("No PDB files found in UNBOUNDED_DIR.")
            sys.exit(1)
        results = run_analysis_mode(pdb_sources=unbounded_files, analysis_type="unbounded")
        summarize_and_print_results(
            results, analysis_type="unbounded",
            output_file=os.path.join(out_dir, "unbounded_summary.json"),
            csv_output_file=os.path.join(out_dir, "unbounded_summary.csv"),
        )
    else:
        logging.error("Invalid configuration: check UNBOUNDED_DIR and CORRESPONDENCE_FILE.")
        sys.exit(1)

    logging.info(f"=== Analysis complete. Results in: {out_dir} ===")


if __name__ == "__main__":
    main()
