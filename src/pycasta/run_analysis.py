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


print("Working directory:", os.getcwd())


import sys
from datetime import datetime

import numpy as np
import pandas as pd
from Bio.PDB import Superimposer
from biopandas.pdb import PandasPdb
from Bio.PDB.Atom import Atom
from scipy.spatial import Delaunay, KDTree
from scipy.spatial.distance import cdist
from tabulate import tabulate

from config import (
    BOUNDED_DIR,
    UNBOUNDED_DIR,
    DATASET_DIR,
    OUTPUT_DIR,
    INCLUDE_WATER,
    MIN_POCKET_VOLUME,
    FLOW_KWARGS,
    MERGE_CLUSTERS,
    MERGE_THRESHOLD,
    USE_EXISTING_RESULTS,
    VERSION_TAG,
    CORRESPONDENCE_FILE,
    VOLUME_METHOD,
    SIMPLE_SPLIT,
    USE_CGAL,
    USE_SASA,
    FILTER_ALPHA_BY_SASA,
    SASA_THRESHOLD,
    SASA_METHOD,
    LIGAND_CONTACT_METHOD,
    LIGAND_CONTACT_THRESHOLD,
    NUM_POCKETS_TO_SAVE,
    STRICT_DISTANCE_THRESHOLD,
    USE_SASA_CONTACT_VALIDATION,
    SASA_CONTACT_THRESHOLD,
    MESH_EXTRUSION_DISTANCE,
    NUM_POCKETS_TO_VALIDATE,
    VALIDATION_METHOD,
    FAKE_SPHERE_RADIUS,
    SAVE_EXTRUDED_MESHES,
)

print("DEBUG: Looking for PDBs in BOUNDED_DIR:", os.path.abspath(BOUNDED_DIR))
if os.path.isdir(BOUNDED_DIR):
    print("DEBUG: Files found in directory:", os.listdir(BOUNDED_DIR))
else:
    print("DEBUG: Directory does not exist:", BOUNDED_DIR)


# Core modules
from pocket_detection import detect_pockets, compute_analytic_pocket_volume
from alpha_shape import compute_alpha_complex_from_tetrahedra
from cgal_wdelaunay import cgal_weighted_delaunay

# Utility modules
from utils.preprocessing_utils import load_and_separate_pdb, calculate_atomic_radii
from utils.data_utils import (
    save_pocket_results_csv,
    save_json,
    load_json,
    convert_numpy_types,
    summarize_and_print_results,
)
from utils.geometry_utils import (
    validate_ligand_in_extruded_mesh,
    validate_ligand_in_fake_sphere,
)
from utils.sasa_utils import (
    compute_sasa,
    compute_ligand_contact_sasa,
    evaluate_sasa_contact,
)
from utils.pocket_utils import (
    apply_splitting_to_shared_atoms,
    compute_volume_for_tetra_group,
    merge_nearby_clusters,
    compute_mouth_parameters,
    save_pocket_properties,
)
from ranking import compute_ranking_scores


print("Working directory:", os.getcwd())


# Set up logging configuration
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

logging.info(f"🔍 Using Validation Method: {VALIDATION_METHOD}")


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


# DummyAtom class for compatibility with Superimposer
class DummyAtom:
    def __init__(self, coord):
        self.coord = np.array(coord, dtype=float)

    def get_coord(self):
        return self.coord


def match_atoms_by_residue_and_name(
    fixed_atoms, moving_atoms, allowed_atom_names=["CA"]
):
    # Filter atoms based on allowed names (corrected for dicts)
    fixed_filtered = [a for a in fixed_atoms if a["atom_name"] in allowed_atom_names]
    moving_filtered = [a for a in moving_atoms if a["atom_name"] in allowed_atom_names]

    # Build dictionaries for fast lookup by (resid, atom_name)
    fixed_dict = {(a["resid"], a["atom_name"]): a["coord"] for a in fixed_filtered}
    moving_dict = {(a["resid"], a["atom_name"]): a["coord"] for a in moving_filtered}

    # Find common keys and extract coordinates
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
        fullname=" X  ",
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
    else:
        return aligned, rmsd


def analyze_top_pockets(result):
    """
    Determines the first validated pocket based on the chosen validation method.
    Returns a tuple: (top_index, validation_tags)
    where top_index is the 1-based index of the first pocket that validated (or None),
    and validation_tags is a list (one per pocket) indicating the chosen method if validated, or "None".
    This function now checks the unified field 'ligand_validation'.
    """
    from config import VALIDATION_METHOD

    method = VALIDATION_METHOD.lower()
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


def run_analysis_mode(pdb_sources, analysis_type="single", correspondence_df=None):
    results = []
    versioned_output_dir = os.path.join(OUTPUT_DIR, VERSION_TAG)
    os.makedirs(versioned_output_dir, exist_ok=True)

    # Paired mode
    if analysis_type == "paired" and correspondence_df is not None:
        for _, row in correspondence_df.iterrows():
            bounded_name = row["bounded"]
            unbounded_name = row["unbounded"]
            bounded_path = os.path.join(BOUNDED_DIR, f"{bounded_name}.pdb")
            unbounded_path = os.path.join(UNBOUNDED_DIR, f"{unbounded_name}.pdb")

            if not os.path.exists(bounded_path) or not os.path.exists(unbounded_path):
                logging.warning(
                    f"Skipping pair {bounded_name} & {unbounded_name}: Missing file(s)."
                )
                continue

            logging.info(f"Processing pair: {bounded_name} & {unbounded_name}")
            bounded_result = process_pdb(bounded_path)
            if not bounded_result or not bounded_result.get("ranked_pockets"):
                logging.warning(f"Skipping {bounded_name}: No pockets found.")
                continue

            unbounded_result = process_pdb(unbounded_path)
            if not unbounded_result or not unbounded_result.get("ranked_pockets"):
                logging.warning(f"Skipping {unbounded_name}: No pockets found.")
                continue

            # Align Proteins
            # Optional: choose backbone or CA atoms for alignment
            ALIGN_ATOM_NAMES = ["CA"]  # or ["N", "CA", "C", "O"]

            # Match atoms by residue ID and atom name
            fixed_coords, moving_coords = match_atoms_by_residue_and_name(
                bounded_result["protein_atoms"],
                unbounded_result["protein_atoms"],
                allowed_atom_names=ALIGN_ATOM_NAMES,
            )

            # Perform alignment and retrieve transformation
            aligned_coords, rmsd, (R, t) = align_using_superimposer(
                fixed_coords, moving_coords, return_transform=True
            )

            # Apply rotation/translation to all unbounded protein coords
            all_coords = np.array(unbounded_result["protein_coords"])
            aligned_all = np.dot(all_coords, R) + t
            unbounded_result["protein_coords"] = aligned_all.tolist()

            logging.info(
                f"Alignment RMSD (filtered on {ALIGN_ATOM_NAMES}): {rmsd:.3f} Å"
            )

            # ✅ Unpack rank + validation info correctly
            bounded_rank, bounded_validation = analyze_top_pockets(bounded_result)
            unbounded_rank, unbounded_validation = analyze_top_pockets(unbounded_result)

            paired_data = {
                "bounded_file": bounded_name,
                "unbounded_file": unbounded_name,
                "bounded_pocket_rank": bounded_rank,
                "unbounded_pocket_rank": unbounded_rank,
                "is_top1_bounded": (bounded_rank == 1) if bounded_rank else False,
                "is_top3_bounded": (bounded_rank and bounded_rank <= 3),
                "is_top5_bounded": (bounded_rank and bounded_rank <= 5),
                "is_top1_unbounded": (unbounded_rank == 1) if unbounded_rank else False,
                "is_top3_unbounded": (unbounded_rank and unbounded_rank <= 3),
                "is_top5_unbounded": (unbounded_rank and unbounded_rank <= 5),
                "alignment_rmsd": rmsd,
                "step_to_ligand": bounded_result.get("step_to_ligand", None),
                "step_to_ligand_mesh": bounded_result.get("step_to_ligand_mesh", None),
                "ligand_to_pocket_distances": bounded_result.get(
                    "ligand_to_pocket_distances", []
                ),
                "ligand_mesh_distances": bounded_result.get(
                    "ligand_mesh_distances", []
                ),
                "ligand_containment_mesh": bounded_result.get(
                    "ligand_containment_mesh", []
                ),
                "ligand_containment_strict": bounded_result.get(
                    "ligand_containment_strict", []
                ),
                "pocket_volumes": bounded_result.get("pocket_volume", []),
                "pocket_depths": bounded_result.get("pocket_depths", []),
                "mouth_area": bounded_result.get("mouth_area", []),
                "mouth_perimeter": bounded_result.get("mouth_perimeter", []),
                "ligand_coords": bounded_result.get("ligand_coords", []),
            }

            results.append(paired_data)

            save_json(
                bounded_result,
                os.path.join(versioned_output_dir, f"bounded_{bounded_name}.json"),
            )
            save_json(
                unbounded_result,
                os.path.join(versioned_output_dir, f"unbounded_{unbounded_name}.json"),
            )

    # Single (bounded) mode
    else:
        for pdb_file in pdb_sources:
            logging.info(f"Processing {pdb_file}")
            result = process_pdb(pdb_file)
            if (
                not result
                or "ranked_pockets" not in result
                or len(result["ranked_pockets"]) == 0
            ):
                logging.warning(f"Skipping {pdb_file}: No pockets found.")
                continue
            pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
            json_path = os.path.join(versioned_output_dir, f"{pdb_name}_pockets.json")
            save_json(result, json_path)
            rank_sasa = result.get("step_to_ligand", None)
            rank_mesh = result.get("step_to_ligand_mesh", None)
            chosen_rank = (
                rank_mesh if VALIDATION_METHOD == "mesh_extrusion" else rank_sasa
            )
            result_entry = {
                "pdb_path": pdb_file,
                "file": os.path.basename(pdb_file),
                "pocket_rank": chosen_rank,
                "step_to_ligand": rank_sasa,
                "step_to_ligand_mesh": rank_mesh,
                "is_top1": (
                    (chosen_rank == 1) if isinstance(chosen_rank, int) else False
                ),
                "is_top3": (
                    (chosen_rank <= 3) if isinstance(chosen_rank, int) else False
                ),
                "is_top5": (
                    (chosen_rank <= 5) if isinstance(chosen_rank, int) else False
                ),
                "pocket_volumes": result.get("pocket_volumes", []),
                "pocket_depths": result.get("pocket_depths", []),
                "mouth_areas": result.get("mouth_area", []),
                "mouth_perimeters": result.get("mouth_perimeter", []),
                "ligand_coords": result.get("ligand_coords", []),
                "ligand_to_pocket_distances": result.get(
                    "ligand_to_pocket_distances", []
                ),
                "num_pockets": len(result.get("ranked_pockets", [])),
            }
            results.append(result_entry)
    return results


def get_output_filename(pdb_path, version_tag):
    base = os.path.splitext(os.path.basename(pdb_path))[0]
    out_dir = os.path.join("results", version_tag)
    os.makedirs(out_dir, exist_ok=True)  # Crea la directory se non esiste
    return os.path.normpath(os.path.join(out_dir, f"{base}_pockets.json"))


def process_pdb(pdb_path):
    logging.info(f"=== Starting processing of PDB: {pdb_path} ===")

    result = {
        "pdb_path": pdb_path,
        "ranked_pockets": [],
        "ranking_scores": [],
        "ligand_coords": [],
        "protein_coords": [],
        "step_to_ligand": None,
        "validation_methods": [],
        "pocket_volumes": [],
        "pocket_depths": [],
        "mouth_area": [],
        "mouth_perimeter": [],
        "ligand_to_pocket_distances": [],
    }

    base = os.path.splitext(os.path.basename(pdb_path))[0]
    out_dir = os.path.join(OUTPUT_DIR, VERSION_TAG)
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, f"{base}_pockets.json")

    print(f"CHECKING IF OUTPUT EXISTS: {out_file}")
    if USE_EXISTING_RESULTS and os.path.exists(out_file):
        print(f"EXISTING RESULT FOUND: {out_file}")
        cached_result = load_json(out_file)
        print(
            f"Loaded cached_result keys: {list(cached_result.keys()) if cached_result else cached_result}"
        )
        if cached_result:
            print("USING EXISTING RESULT: skipping calculation")
            return cached_result
        print("EXISTING RESULT INVALID OR EMPTY: recalculating")

    try:
        protein_coords, ligand_coords, protein_atoms, ppdb = load_and_separate_pdb(
            pdb_path
        )
        radii = calculate_atomic_radii(ppdb)
    except Exception as e:
        logging.error(f"Failed to process PDB {pdb_path}: {e}")
        return result

    ligand_coords = (
        np.array(ligand_coords)
        if ligand_coords is not None and len(ligand_coords) > 0
        else np.array([])
    )

    result["ligand_coords"] = ligand_coords
    result["protein_coords"] = protein_coords.tolist()
    result["protein_atoms"] = protein_atoms

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
        logging.error(f"Delaunay triangulation failed for {pdb_path}: {e}")
        return result
    tetra_positions = protein_coords[simplices]

    # Alpha shape
    logging.info("Computing alpha shape...")
    from config import get_alpha

    alpha = get_alpha()
    logging.info(f"Using alpha shape with alpha = {alpha}")
    alpha_mask, radii, tetra_positions, simplices = (
        compute_alpha_complex_from_tetrahedra(
            simplices,
            tetra_positions,
            alpha,
            base,
            protein_coords,
            FILTER_ALPHA_BY_SASA,
            SASA_CONTACT_THRESHOLD,
            pdb_path,
        )
    )

    # Pocket detection
    logging.info("Detecting pockets...")
    pocket_info = detect_pockets(
        protein_coords,
        simplices,
        tetra_positions,
        alpha_mask,
        min_volume_threshold=MIN_POCKET_VOLUME,
        flow_params=FLOW_KWARGS,
        molecule_name=base,
        radii=radii,
    )

    ranked_pockets = pocket_info.get("ranked_pockets", [])
    representative_points = pocket_info.get("representative_points", [])
    ranking_scores = pocket_info.get("ranking_scores", [])
    result.update(
        {
            "ranked_pockets": ranked_pockets,
            "representative_points": representative_points,
            "ranking_scores": ranking_scores,
        }
    )

    # Compute pocket distances, volumes, and depths
    ligand_to_pocket_distances = []

    pocket_volumes = []
    pocket_depths = []

    for idx, pocket in enumerate(ranked_pockets):
        valid_indices = [t for t in pocket if 0 <= t < len(protein_coords)]
        if valid_indices:
            volume = compute_analytic_pocket_volume(tetra_positions, valid_indices)
            depth = np.max(
                [
                    np.linalg.norm(protein_coords[t] - protein_coords[valid_indices[0]])
                    for t in valid_indices
                ]
            )
            pocket_volumes.append(volume)
            pocket_depths.append(depth)
            logging.info(
                f"📦 Pocket {idx+1}: Volume = {volume:.2f}, Depth = {depth:.2f}"
            )
        else:
            pocket_volumes.append(0)
            pocket_depths.append(0)
            logging.warning(
                f"⚠️ Pocket {idx+1}: No valid tetrahedra. Volume and depth set to 0."
            )

    result["ligand_to_pocket_distances"] = ligand_to_pocket_distances
    result["pocket_volumes"] = pocket_volumes
    result["pocket_depths"] = pocket_depths

    # Mouth calculations
    mouth_params = compute_mouth_parameters(
        ranked_pockets, simplices, protein_coords, ligand_coords
    )
    result["mouth_area"] = [m.get("mouth_area", 0) for m in mouth_params]
    result["mouth_perimeter"] = [m.get("mouth_circumference", 0) for m in mouth_params]

    # Pocket validation based on selected method
    validation_methods = []
    method = VALIDATION_METHOD.lower()

    ligand_containment_mesh = []
    ligand_mesh_distances = []

    if method == "mesh_extrusion":
        for mouth in mouth_params:
            rim_atoms = mouth.get("rim_atoms", [])
            mouth_coords = np.array(protein_coords)[rim_atoms] if rim_atoms else []
            in_contact, mesh, min_dist = validate_ligand_in_extruded_mesh(
                mouth_coords, ligand_coords, extrusion_distance=MESH_EXTRUSION_DISTANCE
            )
            validation_methods.append("Mesh" if in_contact else "None")
            ligand_containment_mesh.append(in_contact)
            ligand_mesh_distances.append(min_dist)

            if SAVE_EXTRUDED_MESHES and mesh:
                mesh_dir = os.path.join(
                    OUTPUT_DIR, VERSION_TAG, "extruded_meshes", base
                )
                os.makedirs(mesh_dir, exist_ok=True)
                mesh.export(
                    os.path.join(
                        mesh_dir, f"{base}_pocket_{len(validation_methods)}.ply"
                    )
                )

    elif method == "sasa":
        sasa_before = compute_sasa(pdb_path, SASA_METHOD)
        sasa_after = compute_sasa(pdb_path, SASA_METHOD, remove_ligand=True)
        tree = KDTree(protein_coords)

        for pocket in ranked_pockets:
            tetra_pts = np.concatenate([tetra_positions[i] for i in pocket], axis=0)
            _, atom_indices = tree.query(tetra_pts, k=1)
            unique_atoms = list(set(atom_indices))
            valid = evaluate_sasa_contact(
                unique_atoms,
                sasa_before,
                sasa_after,
                ligand_coords,
                protein_coords,
                SASA_CONTACT_THRESHOLD,
            )
            validation_methods.append("SASA" if valid else "None")

    elif method == "fake_ball":
        for mouth in mouth_params:
            rim_atoms = mouth.get("rim_atoms", [])
            mouth_coords = np.array(protein_coords)[rim_atoms] if rim_atoms else []
            valid = validate_ligand_in_fake_sphere(
                mouth_coords, ligand_coords, FAKE_SPHERE_RADIUS
            )
            validation_methods.append("FakeBall" if valid else "None")

    else:
        validation_methods = ["None"] * len(ranked_pockets)

    result["validation_methods"] = validation_methods
    first_valid = next(
        (i + 1 for i, m in enumerate(validation_methods) if m != "None"), None
    )
    result["step_to_ligand"] = first_valid
    result["ligand_containment_mesh"] = ligand_containment_mesh
    print(ligand_containment_mesh)
    result["ligand_mesh_distances"] = ligand_mesh_distances
    print(ligand_mesh_distances)
    # Save top pockets as PDBs
    protein_df = ppdb.df["ATOM"]
    output_pdb_dir = os.path.join(OUTPUT_DIR, VERSION_TAG, "pocket_pdbs", base)
    os.makedirs(output_pdb_dir, exist_ok=True)

    for rank, pocket in enumerate(ranked_pockets[:NUM_POCKETS_TO_SAVE], 1):
        pocket_atoms = np.unique(
            np.concatenate([tetra_positions[i] for i in pocket]), axis=0
        )
        tree = KDTree(protein_coords)
        _, idx = tree.query(pocket_atoms, k=1)
        pocket_df = protein_df.iloc[idx].copy()
        pocket_df["b_factor"] = rank
        pdb_out = PandasPdb()
        pdb_out.df["ATOM"] = pocket_df
        pdb_out.to_pdb(os.path.join(output_pdb_dir, f"{base}_ranked_{rank}.pdb"))

        logging.info(f"✅ Pocket {rank} saved.")

    # Save properties
    save_pocket_properties(result, base, OUTPUT_DIR)

    logging.info(f"Total pockets: {len(ranked_pockets)}")
    logging.info(f"=== Finished processing {pdb_path} ===")

    return convert_numpy_types(result)


def main():
    logging.info(f"=== Starting analysis. Version: {VERSION_TAG} ===")
    versioned_output_dir = os.path.join(OUTPUT_DIR, VERSION_TAG)
    os.makedirs(versioned_output_dir, exist_ok=True)
    bounded_files = glob.glob(os.path.join(BOUNDED_DIR, "*.pdb"))
    if not bounded_files:
        logging.error("No PDB files found in BOUNDED_DIR.")
        sys.exit(1)
    if (
        UNBOUNDED_DIR
        and os.path.exists(UNBOUNDED_DIR)
        and CORRESPONDENCE_FILE
        and os.path.exists(CORRESPONDENCE_FILE)
    ):
        logging.info("Running paired analysis...")
        correspondence_df = pd.read_excel(CORRESPONDENCE_FILE)
        paired_results = run_analysis_mode(
            pdb_sources=bounded_files,
            analysis_type="paired",
            correspondence_df=correspondence_df,
        )
        summarize_and_print_results(
            paired_results,
            analysis_type="paired",
            output_file=os.path.join(versioned_output_dir, "paired_summary.json"),
            csv_output_file=os.path.join(versioned_output_dir, "paired_summary.csv"),
        )
    elif not UNBOUNDED_DIR or not os.path.exists(UNBOUNDED_DIR):
        logging.info("Running single (bounded only) analysis...")
        results = run_analysis_mode(pdb_sources=bounded_files, analysis_type="single")
        summarize_and_print_results(
            results,
            analysis_type="single",
            output_file=os.path.join(versioned_output_dir, "single_summary.json"),
            csv_output_file=os.path.join(versioned_output_dir, "single_summary.csv"),
        )
    elif UNBOUNDED_DIR and os.path.exists(UNBOUNDED_DIR) and not CORRESPONDENCE_FILE:
        logging.info("Running unbounded analysis...")
        unbounded_files = glob.glob(os.path.join(UNBOUNDED_DIR, "*.pdb"))
        if not unbounded_files:
            logging.error("No PDB files found in UNBOUNDED_DIR.")
            sys.exit(1)
        results = run_analysis_mode(
            pdb_sources=unbounded_files, analysis_type="unbounded"
        )
        summarize_and_print_results(
            results,
            analysis_type="unbounded",
            output_file=os.path.join(versioned_output_dir, "unbounded_summary.json"),
            csv_output_file=os.path.join(versioned_output_dir, "unbounded_summary.csv"),
        )
    else:
        logging.error(
            "Invalid configuration: check UNBOUNDED_DIR and CORRESPONDENCE_FILE settings."
        )
        sys.exit(1)
    logging.info(
        f"=== Analysis completed. Results saved in: {versioned_output_dir} ==="
    )


if __name__ == "__main__":
    main()
