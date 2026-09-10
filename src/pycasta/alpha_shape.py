#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
alpha_shape.py

Module for computing the alpha shape and for the flow_detection function.
"""

import os
import numpy as np
import logging
from scipy.spatial import KDTree
from biopandas.pdb import PandasPdb

from .alpha_constraints import apply_alpha_constraints

from .config import (
   
    SIMPLE_SPLIT,
    ATOM_SPLITTING,
    ATOM_SPLITTING_MODE,
    SAVE_TRIANGULATION_DEBUG,
    SAVE_ALPHA,
    VERSION_TAG,
    OUTPUT_DIR,
    FILTER_AFTER_MERGE,
    USE_GRAPH_FILTERED_ALPHA,
    ALPHA_MAX_EDGE_FACTOR,
    ALPHA_LOCAL_SCALE_K,
)

from .utils.geometry_utils import compute_all_circumsphere_radii, compute_power_radii,log_radii_statistics
from .utils.pocket_utils import apply_splitting_to_shared_atoms
from .utils.sasa_utils import compute_sasa


def save_triangulation_debug_data(
    simplices,
    tetra_positions,
    protein_coords,
    molecule_name,
    radii=None,
):
    """
    Save raw triangulation data before alpha filtering, for debugging/visualization.

    Saved fields:
        protein_coords: (N, 3)
        simplices: (M, 4) indices into protein_coords
        tetrahedra: (M, 4, 3) tetra vertex coordinates
        circumsphere_radii: (M,) optional
    """
    if molecule_name is None or protein_coords is None:
        logging.warning(
            "[triangulation_debug] molecule_name or protein_coords missing; skip save."
        )
        return None

    version_dir = os.path.join(OUTPUT_DIR, VERSION_TAG)
    os.makedirs(version_dir, exist_ok=True)

    out_file = os.path.join(version_dir, f"{molecule_name}.triangulation_debug.npz")

    payload = {
        "protein_coords": np.asarray(protein_coords, dtype=float),
        "simplices": np.asarray(simplices, dtype=np.int32),
        "tetrahedra": np.asarray(tetra_positions, dtype=float),
    }

    if radii is not None:
        payload["circumsphere_radii"] = np.asarray(radii, dtype=float)

    np.savez(out_file, **payload)
    logging.info(f"[triangulation_debug] saved to {out_file}")
    return out_file


def filter_alpha_complex_by_sasa(
    tetrahedra, alpha_values, protein_coords, pdb_file, sasa_threshold
):
    logging.info("Computing SASA for filtering the alpha complex...")
    sasa_dict = compute_sasa(pdb_file, method="pymol", remove_ligand=True)
    if sasa_dict is None:
        logging.error("SASA computation failed; no filtering applied.")
        return tetrahedra, alpha_values, np.ones_like(alpha_values, dtype=bool)
    ppdb = PandasPdb().read_pdb(pdb_file)
    atoms_df = ppdb.df["ATOM"]
    atom_residues = atoms_df["residue_number"].values
    

    tree = KDTree(protein_coords)
    filtered_indices = []
    total_tetra = len(alpha_values)
    for i, tet in enumerate(tetrahedra):
        centroid = np.mean(tet, axis=0)
        dist, idx = tree.query(centroid)
        residue = atom_residues[idx]
        sasa = sasa_dict.get(str(residue), 0)
        if sasa < sasa_threshold:
            filtered_indices.append(i)
    logging.info(
        f"Tetrahedra: {total_tetra} total, {len(filtered_indices)} pass the SASA filter (threshold={sasa_threshold})."
    )
    filtered_tetrahedra = tetrahedra[filtered_indices]
    filtered_alpha_values = alpha_values[filtered_indices]
    filtered_mask = np.zeros_like(alpha_values, dtype=bool)
    filtered_mask[filtered_indices] = True
    return filtered_tetrahedra, filtered_alpha_values, filtered_mask

def compute_alpha_complex_from_tetrahedra(
    simplices,
    tetra_positions,
    alpha_value,
    molecule_name=None,
    protein_coords=None,
    radii=None,  # NEW: atomic radii for weighted Delaunay (power radii)
    filter_by_sasa=False,
    pdb_file=None,
    save_alpha_data=True,
):
    """
    Compute the alpha complex for a set of tetrahedra.
    
    This function supports both standard (unweighted) and weighted (regular) Delaunay 
    triangulations:
    
    
    
    All existing features are preserved:
    - Graph-based alpha constraints (USE_GRAPH_FILTERED_ALPHA)
    - SASA filtering (filter_by_sasa)
    - Debug data saving (SAVE_TRIANGULATION_DEBUG, SAVE_ALPHA)
    
    Parameters
    ----------
    simplices : np.ndarray
        Array of shape (N, 4) containing vertex indices for each tetrahedron 
        in the original triangulation.
    tetra_positions : np.ndarray
        Array of shape (N, 4, 3) containing 3D coordinates of the 4 vertices 
        for each tetrahedron.
    alpha_value : float
        Alpha threshold for filtering tetrahedra (tetrahedra with radius < alpha_value are kept).
    molecule_name : str, optional
        Name of the molecule (used for logging and file naming).
    protein_coords : np.ndarray, optional
        Array of shape (M, 3) containing coordinates of all atoms in the protein.
        Required for weighted Delaunay (power radii) and graph-based constraints.
    radii : np.ndarray, optional
        Array of shape (M,) containing atomic radii (e.g., van der Waals radii).
        Required for weighted Delaunay (power radii). If None, falls back to 
        standard circumradii.
    filter_by_sasa : bool, optional
        If True, apply SASA-based filtering to the alpha complex.
    sasa_threshold : float, optional
        SASA threshold for filtering (default: SASA_CONTACT_THRESHOLD).
    pdb_file : str, optional
        Path to the PDB file (required for SASA filtering).
    save_alpha_data : bool, optional
        If True, save alpha shape data to NPZ files for debugging.
    
    Returns
    -------
    alpha_mask : np.ndarray
        Boolean mask of shape (N,) indicating which tetrahedra belong to the 
        alpha complex (radius < alpha_value).
    radii : np.ndarray
        Array of shape (N,) containing the radii used for filtering 
        (either circumradii or power radii).
    tetra_positions : np.ndarray
        Possibly updated tetra_positions after SASA filtering.
    simplices : np.ndarray
        Filtered simplices array (aligned with tetra_positions after SASA filtering).
    """
    import os
    import logging
    import numpy as np
    from .config import (
        SAVE_TRIANGULATION_DEBUG,
        SAVE_ALPHA,
        OUTPUT_DIR,
        VERSION_TAG,
        USE_GRAPH_FILTERED_ALPHA,
        ALPHA_MAX_EDGE_FACTOR,
        ALPHA_LOCAL_SCALE_K,
    )
    from .utils.geometry_utils import compute_all_circumsphere_radii, compute_power_radii, log_radii_statistics
    from .alpha_constraints import apply_alpha_constraints
    #from .utils.sasa_utils import filter_alpha_complex_by_sasa
    
    # Validate input shape
    assert tetra_positions.ndim == 3 and tetra_positions.shape[1:] == (
        4,
        3,
    ), f"Invalid tetra_positions shape: {tetra_positions.shape}"
    
    # =========================================================================
    # STEP 1: Compute radii (circumradii or power radii)
    # =========================================================================
    

   

    protein_coords = np.asarray(protein_coords, dtype=float)
    radii_array = np.asarray(radii, dtype=float)
       
    radii = compute_power_radii(simplices, protein_coords, radii_array)
    
    log_radii_statistics(radii)
    
    # Compute initial alpha mask based on radius threshold
    alpha_mask = (radii < alpha_value).astype(bool)
    
    # =========================================================================
    # STEP 2: Save raw triangulation debug data (before alpha filtering)
    # =========================================================================
    
    if SAVE_TRIANGULATION_DEBUG:
        try:
            save_triangulation_debug_data(
                simplices=simplices,
                tetra_positions=tetra_positions,
                protein_coords=protein_coords,
                molecule_name=molecule_name,
                radii=radii,
            )
        except Exception as e:
            logging.warning(f"[triangulation_debug] failed to save: {e}")
    
    # =========================================================================
    # STEP 3: Apply graph-based alpha constraints (if enabled)
    # =========================================================================
    
    constraint_stats = {
        "graph_filter_enabled": False,
        "num_alpha_before_graph_filter": int(np.sum(alpha_mask)),
        "num_alpha_after_graph_filter": int(np.sum(alpha_mask)),
        "num_removed_from_alpha_by_graph_filter": 0,
    }
    
    if USE_GRAPH_FILTERED_ALPHA:
        if protein_coords is None:
            logging.warning(
                "USE_GRAPH_FILTERED_ALPHA=True but protein_coords is None; "
                "skipping graph-based alpha constraints."
            )
        else:
            alpha_mask, constraint_stats = apply_alpha_constraints(
                alpha_mask=alpha_mask,
                tetra_positions=tetra_positions,
                protein_coords=protein_coords,
                use_graph_filter=True,
                max_edge_factor=ALPHA_MAX_EDGE_FACTOR,
                local_scale_k=ALPHA_LOCAL_SCALE_K,
            )
            
            logging.info(
                "[alpha_constraints] before=%d after=%d removed=%d edge_factor=%.3f k=%d",
                constraint_stats["num_alpha_before_graph_filter"],
                constraint_stats["num_alpha_after_graph_filter"],
                constraint_stats["num_removed_from_alpha_by_graph_filter"],
                ALPHA_MAX_EDGE_FACTOR,
                ALPHA_LOCAL_SCALE_K,
            )
    
    # =========================================================================
    # STEP 4: Save alpha shape data (before SASA filtering)
    # =========================================================================
    
    if save_alpha_data and SAVE_ALPHA and molecule_name and protein_coords is not None:
        version_dir = os.path.join(OUTPUT_DIR, VERSION_TAG)
        os.makedirs(version_dir, exist_ok=True)
        
        out_file = os.path.join(version_dir, f"{molecule_name}.alpha.npz")
        
        np.savez(
            out_file,
            alpha_value=float(alpha_value),
            alpha_values=radii,
            alpha_mask=alpha_mask,
            protein_coords=protein_coords,
            tetrahedra=tetra_positions,
        )
        
        logging.info(f"Alpha shape data saved to {out_file}")
    
    # =========================================================================
    # STEP 5: Apply SASA filtering (if enabled)
    # =========================================================================
    
    if filter_by_sasa:
        if pdb_file is None:
            logging.error("pdb_file must be provided when filtering by SASA.")
        else:
            # Filter alpha complex by SASA
            filtered_tetra, filtered_alpha_values, filtered_mask = (
                filter_alpha_complex_by_sasa(
                    tetra_positions, radii, protein_coords, pdb_file, sasa_threshold
                )
            )
            
            # Update tetra_positions and radii based on the SASA filter
            tetra_positions = filtered_tetra
            radii = filtered_alpha_values
            
            # Recalculate alpha_mask from the filtered radii
            alpha_mask = (radii < alpha_value).astype(bool)
            
            # Reset constraint stats for post-SASA filtering
            constraint_stats = {
                "graph_filter_enabled": False,
                "num_alpha_before_graph_filter": int(np.sum(alpha_mask)),
                "num_alpha_after_graph_filter": int(np.sum(alpha_mask)),
                "num_removed_from_alpha_by_graph_filter": 0,
            }
            
            # Re-apply graph-based constraints after SASA filtering
            if USE_GRAPH_FILTERED_ALPHA:
                if protein_coords is None:
                    logging.warning(
                        "USE_GRAPH_FILTERED_ALPHA=True but protein_coords is None; "
                        "skipping graph-based alpha constraints after SASA filtering."
                    )
                else:
                    alpha_mask, constraint_stats = apply_alpha_constraints(
                        alpha_mask=alpha_mask,
                        tetra_positions=tetra_positions,
                        protein_coords=protein_coords,
                        use_graph_filter=True,
                        max_edge_factor=ALPHA_MAX_EDGE_FACTOR,
                        local_scale_k=ALPHA_LOCAL_SCALE_K,
                    )
                    
                    logging.info(
                        "[alpha_constraints][post_sasa] before=%d after=%d removed=%d edge_factor=%.3f k=%d",
                        constraint_stats["num_alpha_before_graph_filter"],
                        constraint_stats["num_alpha_after_graph_filter"],
                        constraint_stats["num_removed_from_alpha_by_graph_filter"],
                        ALPHA_MAX_EDGE_FACTOR,
                        ALPHA_LOCAL_SCALE_K,
                    )
            
            # Update simplices using the same boolean mask from SASA filtering
            simplices = simplices[filtered_mask]
            
            # Save filtered alpha shape data
            version_dir = os.path.join(OUTPUT_DIR, VERSION_TAG)
            os.makedirs(version_dir, exist_ok=True)
            out_file_filtered = os.path.join(
                version_dir, f"{molecule_name}.alpha_filtered.npz"
            )
            
            np.savez(
                out_file_filtered,
                alpha_values=radii,
                alpha_mask=alpha_mask,
                protein_coords=protein_coords,
                tetrahedra=tetra_positions,
            )
            logging.info(f"Filtered alpha shape data saved to {out_file_filtered}")
    
    # =========================================================================
    # STEP 6: Return results
    # =========================================================================
    
    return alpha_mask, radii, tetra_positions, simplices




def flow_detection(
    simplices,
    alpha_mask,
    tetra_positions,
    protein_coords,
    flow_params=None,
    proxy_values=None,  # NEW: precomputed radii (circum or power)
):
    """
    Apply discrete flow on the tetrahedra to grow pockets.
    This version uses a relative tolerance (tol_fraction) to decide if a neighbor is lower.
    
    Parameters
    ----------
    proxy_values : np.ndarray, optional
        Precomputed proxy values (e.g., circumradii or power radii) for all tetrahedra.
        If None, defaults to computing circumradii from tetra_positions.
    """
    if flow_params is None:
        flow_params = {}
    sigma_p = flow_params.get("sigma_p", 1.0)
    # Use a relative tolerance fraction (e.g., 0.01 means a 1% drop)
    tol_fraction = flow_params.get("tol_fraction", 0.01)
    max_steps = flow_params.get("max_steps", 100)
    adaptive = flow_params.get("adaptive", False)
    adaptive_factor = flow_params.get("adaptive_factor", 0.9)
    min_steps = flow_params.get("min_steps", 0)

    # Compute proxy values for all tetrahedra (e.g., based on circumsphere radii)
    if proxy_values is None:
        # Fallback: compute circumradii (classic mode)
        proxies = compute_all_circumsphere_radii(tetra_positions)
        proxy_values = sigma_p * proxies
    else:
        # Use precomputed proxy values (e.g., power radii from weighted Delaunay)
        # sigma_p is already applied if needed by the caller
        proxy_values = np.asarray(proxy_values, dtype=float)

    # Get indices of tetrahedra that are outside the alpha complex
    empty_indices = np.where(~alpha_mask)[0]

    # Build a face-to-tetrahedron mapping (neighbors sharing a face)
    face_to_tetra = {}
    for i, simplex in enumerate(simplices):
        faces = [
            tuple(sorted(simplex[[a, b, c]]))
            for (a, b, c) in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
        ]
        for face in faces:
            face_to_tetra.setdefault(face, []).append(i)

    # Build neighbor dictionary for each empty tetrahedron
    neighbor_dict = {}
    for i in empty_indices:
        faces = [
            tuple(sorted(simplices[i, [a, b, c]]))
            for (a, b, c) in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
        ]
        neighbors = set()
        for face in faces:
            neighbors.update(face_to_tetra.get(face, []))
        # Exclude the tetrahedron itself and those that are not in the empty region
        neighbor_dict[i] = [j for j in neighbors if j != i and not alpha_mask[j]]

    # Flow: for each tetrahedron in empty_indices, try to follow a gradient descent
    flow_target = {}
    flow_steps = {}
    total_steps = 0

    for i in empty_indices:
        current = i
        steps = 0
        current_proxy = proxy_values[current]
        while steps < max_steps:
            # Use relative tolerance: neighbor qualifies if its proxy is at least tol_fraction lower than current
            lower_neighbors = [
                j
                for j in neighbor_dict[current]
                if (current_proxy - proxy_values[j]) > tol_fraction * current_proxy
            ]
            if not lower_neighbors:
                if adaptive:
                    # Optionally, if no neighbor is found, increase tolerance (or adapt in some way)
                    tol_fraction *= adaptive_factor
                    lower_neighbors = [
                        j
                        for j in neighbor_dict[current]
                        if (current_proxy - proxy_values[j])
                        > tol_fraction * current_proxy
                    ]
                    if not lower_neighbors:
                        break
                else:
                    break
            next_current = min(lower_neighbors, key=lambda j: proxy_values[j])
            current = next_current
            current_proxy = proxy_values[current]
            steps += 1
        if steps >= min_steps:
            flow_target[i] = current
            flow_steps[i] = steps
            total_steps += steps

    average_steps = total_steps / max(len(flow_target), 1)
    logging.info(f"Average number of flow steps per tetrahedron: {average_steps:.2f}")

    # Group tetrahedra by their final sink (flow target)
    sink_groups = {}
    for i, sink in flow_target.items():
        sink_groups.setdefault(sink, []).append(i)

    final_pockets = []
    connectivity_graph = {}
    # For each sink group, build a connectivity graph among tetrahedra
    for sink, tetra_indices in sink_groups.items():
        graph = {i: [] for i in tetra_indices}
        for idx, i in enumerate(tetra_indices):
            for j in tetra_indices[idx + 1 :]:
                if len(set(simplices[i]) & set(simplices[j])) == 3:
                    graph[i].append(j)
                    graph[j].append(i)
        visited = set()
        for node in tetra_indices:
            if node not in visited:
                component = []
                stack = [node]
                while stack:
                    cur = stack.pop()
                    if cur not in visited:
                        visited.add(cur)
                        component.append(cur)
                        stack.extend(graph[cur])
                final_pockets.append(component)
                connectivity_graph[len(final_pockets) - 1] = graph

    
    final_pockets = apply_splitting_to_shared_atoms(
        final_pockets, protein_coords, flow_params.get("atom_radii")
    )
    return final_pockets, flow_steps, connectivity_graph
