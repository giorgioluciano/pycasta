#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
refactor_weighted_only.py

Script di refactoring per rimuovere tutto il supporto non pesato (standard Delaunay)
e fissare i bug principali nel flusso weighted-only.

Esegui da C:\github\pycasta:
    python refactor_weighted_only.py
"""

from pathlib import Path
import re

ROOT = Path("src/pycasta")

def fix_alpha_shape_py():
    """
    Fix alpha_shape.py:
    - Rimuovi supporto circumradii fallback
    - Fix proxy_values nel flow_detection call
    - Fix tol_fraction globale
    - Pulisci docstring
    """
    path = ROOT / "alpha_shape.py"
    if not path.exists():
        print(f"[SKIP] {path} non esiste")
        return

    content = path.read_text(encoding="utf-8")

    # 1. Fix docstring: rimuovi riferimento a standard Delaunay
    old_doc = """Compute the alpha complex for a set of tetrahedra.
    
    This function supports both standard (unweighted) and weighted (regular) Delaunay 
    triangulations:
    
    
    All existing features are preserved:
    - Graph-based alpha constraints (USE_GRAPH_FILTERED_ALPHA)
    - SASA filtering (filter_by_sasa)
    - Debug data saving (SAVE_TRIANGULATION_DEBUG, SAVE_ALPHA)"""

    new_doc = """Compute the alpha complex for a set of tetrahedra (weighted Delaunay only).
    
    This function assumes weighted (regular) Delaunay triangulation.
    The `radii` parameter (atomic van der Waals radii) must be provided.
    
    Features:
    - Graph-based alpha constraints (USE_GRAPH_FILTERED_ALPHA)
    - SASA filtering (filter_by_sasa)
    - Debug data saving (SAVE_TRIANGULATION_DEBUG, SAVE_ALPHA)"""

    content = content.replace(old_doc, new_doc)

    # 2. Fix docstring: rimuovi "If None, falls back to standard circumradii"
    content = re.sub(
        r"radii : np\.ndarray, optional\n.*?Required for weighted Delaunay.*?If None, falls back to standard circumradii\.",
        "radii : np.ndarray\n    Array of shape (M,) containing atomic radii (e.g., van der Waals radii).\n    Required for weighted Delaunay (power radii).",
        content,
        flags=re.DOTALL,
    )

    # 3. Fix docstring: rimuovi "either circumradii or power radii"
    content = re.sub(
        r"Array of shape \(N,\) containing the radii used for filtering \(either circumradii or power radii\)\.",
        "Array of shape (N,) containing the power radii used for filtering.",
        content,
    )

    # 4. Fix compute_alpha_complex_from_tetrahedra: rimuovi fallback a circumradii
    # Cerca il blocco:
    # protein_coords = np.asarray(protein_coords, dtype=float)
    # radii_array = np.asarray(radii, dtype=float)
    # radii = compute_power_radii(...)
    # e assicurati che non ci sia fallback

    # 5. Fix flow_detection call: usa local_radii invece di radii
    # Cerca:
    # pockets, flow_steps, connectivity_graph = flow_detection(
    #     simplices,
    #     alpha_mask,
    #     tetra_positions,
    #     protein_coords,
    #     flow_params,
    #     proxy_values=radii,
    # )
    # e sostituisci con local_radii

    old_flow_call = """pockets, flow_steps, connectivity_graph = flow_detection(
        simplices,
        alpha_mask,
        tetra_positions,
        protein_coords,
        flow_params,
        proxy_values=radii,
    )"""

    new_flow_call = """pockets, flow_steps, connectivity_graph = flow_detection(
        local_simplices,
        alpha_mask,
        local_tetra_positions,
        protein_coords,
        flow_params,
        proxy_values=local_radii,
    )"""

    content = content.replace(old_flow_call, new_flow_call)

    # 6. Fix flow_detection: tol_fraction locale per ogni traiettoria
    # Cerca il blocco while e sostituisci con versione corretta
    old_flow_while = """current_proxy = proxy_values[current]
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
                    tol_fraction *= adaptive_factor"""

    new_flow_while = """current_proxy = proxy_values[current]
        local_tol_fraction = tol_fraction  # Local copy for this trajectory
        while steps < max_steps:
            # Use relative tolerance: neighbor qualifies if its proxy is at least tol_fraction lower than current
            lower_neighbors = [
                j
                for j in neighbor_dict[current]
                if (current_proxy - proxy_values[j]) > local_tol_fraction * current_proxy
            ]
            if not lower_neighbors:
                if adaptive:
                    # Optionally, if no neighbor is found, increase tolerance (or adapt in some way)
                    local_tol_fraction *= adaptive_factor"""

    content = content.replace(old_flow_while, new_flow_while)

    path.write_text(content, encoding="utf-8")
    print(f"[FIX] {path}")


def fix_alpha_stability_py():
    """
    Fix alpha_stability.py:
    - Fix bug: 'radii' non definito, usa local_radii
    """
    path = ROOT / "alpha_stability.py"
    if not path.exists():
        print(f"[SKIP] {path} non esiste")
        return

    content = path.read_text(encoding="utf-8")

    # Fix: proxy_values=radii -> proxy_values=local_radii
    old_call = """pockets, flow_steps, connectivity_graph = flow_detection(
        simplices,
        alpha_mask,
        tetra_positions,
        protein_coords,
        flow_params,
        proxy_values=radii,
    )"""

    new_call = """pockets, flow_steps, connectivity_graph = flow_detection(
        local_simplices,
        alpha_mask,
        local_tetra_positions,
        protein_coords,
        flow_params,
        proxy_values=local_radii,
    )"""

    content = content.replace(old_call, new_call)

    path.write_text(content, encoding="utf-8")
    print(f"[FIX] {path}")


def fix_alpha_diagnostics_py():
    """
    Fix alpha_diagnostics.py:
    - Fix proxy_values nel flow_detection call
    - Rimuovi calcolo raw_radii come circumsphere (non serve)
    """
    path = ROOT / "alpha_diagnostics.py"
    if not path.exists():
        print(f"[SKIP] {path} non esiste")
        return

    content = path.read_text(encoding="utf-8")

    # 1. Fix flow_detection call: usa local_radii
    old_call = """pockets, flow_steps, connectivity_graph = flow_detection(
        simplices,
        alpha_mask,
        tetra_positions,
        protein_coords,
        flow_params,
        proxy_values=radii,
    )"""

    new_call = """pockets, flow_steps, connectivity_graph = flow_detection(
        local_simplices,
        alpha_mask,
        local_tetra_positions,
        protein_coords,
        flow_params,
        proxy_values=local_radii,
    )"""

    content = content.replace(old_call, new_call)

    # 2. Rimuovi calcolo raw_radii (non serve più per weighted-only)
    old_raw_radii = """try:
        from .utils.geometry_utils import compute_all_circumsphere_radii
        raw_radii = compute_all_circumsphere_radii(tetra_positions)
    except Exception as e:
        logging.warning(f"[alpha_diagnostics] Could not compute raw radii separately: {e}")
        raw_radii = np.array([], dtype=float)

    stats = _compute_histogram_stats(raw_radii)
    alpha_grid = _alpha_grid_from_base(
        base_alpha=base_alpha,
        radii=raw_radii,
        alpha_values=alpha_values,
        alpha_offsets=alpha_offsets,
    )"""

    new_raw_radii = """# For weighted-only, skip raw_radii computation (not needed)
    raw_radii = np.array([], dtype=float)
    stats = _compute_histogram_stats(raw_radii)
    alpha_grid = _alpha_grid_from_base(
        base_alpha=base_alpha,
        radii=raw_radii,
        alpha_values=alpha_values,
        alpha_offsets=alpha_offsets,
    )"""

    content = content.replace(old_raw_radii, new_raw_radii)

    path.write_text(content, encoding="utf-8")
    print(f"[FIX] {path}")


def fix_geometry_utils_py():
    """
    Fix geometry_utils.py:
    - Depreca compute_all_circumsphere_radii (non serve più)
    - Pulisci commenti
    """
    path = ROOT / "utils" / "geometry_utils.py"
    if not path.exists():
        print(f"[SKIP] {path} non esiste")
        return

    content = path.read_text(encoding="utf-8")

    # 1. Fix docstring del modulo
    old_mod_doc = """Module containing functions for:
- Computing circumsphere radii of tetrahedra.
- Estimating pocket volume using Monte Carlo integration.
- Validating pockets via mesh extrusion and fake sphere tests."""

    new_mod_doc = """Module containing functions for:
- Computing power radii for weighted Delaunay triangulation.
- Estimating pocket volume using analytic formulas.
- Validating pockets via mesh extrusion and fake sphere tests."""

    content = content.replace(old_mod_doc, new_mod_doc)

    # 2. Aggiungi deprecazione a compute_all_circumsphere_radii
    old_circum_func = """def compute_all_circumsphere_radii(tetra_positions):
    radii = np.array([compute_circumsphere_radius(tetra) for tetra in tetra_positions])
    log_radii_statistics(radii)
    return radii"""

    new_circum_func = """def compute_all_circumsphere_radii(tetra_positions):
    \"\"\"
    DEPRECATED: This function is kept for backward compatibility only.
    For weighted Delaunay, use compute_power_radii() instead.
    \"\"\"
    import warnings
    warnings.warn(
        "compute_all_circumsphere_radii() is deprecated for weighted-only mode. "
        "Use compute_power_radii() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    radii = np.array([compute_circumsphere_radius(tetra) for tetra in tetra_positions])
    log_radii_statistics(radii)
    return radii"""

    content = content.replace(old_circum_func, new_circum_func)

    # 3. Fix log_radii_statistics: cambia messaggio
    old_log_msg = """logging.info(f"Circumsphere radii stats: {stats}")"""
    new_log_msg = """logging.info(f"Power radii stats: {stats}")"""

    content = content.replace(old_log_msg, new_log_msg)

    path.write_text(content, encoding="utf-8")
    print(f"[FIX] {path}")


def fix_cgal_wdelaunay_py():
    """
    Fix cgal_wdelaunay.py:
    - Rimuovi fallback a SciPy standard Delaunay
    - Rinomina funzione in weighted_delaunay
    """
    path = ROOT / "cgal_wdelaunay.py"
    if not path.exists():
        print(f"[SKIP] {path} non esiste")
        return

    content = path.read_text(encoding="utf-8")

    # 1. Rinomina funzione
    content = content.replace(
        "def cgal_weighted_delaunay(protein_coords, radii):",
        "def weighted_delaunay(protein_coords, radii):",
    )

    # 2. Rimuovi fallback a SciPy
    old_fallback = """simplices = Delaunay(protein_coords).simplices
    logging.info(f"[SciPy] Standard Delaunay built: {len(simplices)} simplices.")
    return simplices, None"""

    new_fallback = """logging.error(
        "[weighted_delaunay] CGAL/pocketlab is required for weighted Delaunay. "
        "Please install pocketlab or enable CGAL."
    )
    raise RuntimeError(
        "Weighted Delaunay triangulation failed: CGAL/pocketlab not available."
    )"""

    content = content.replace(old_fallback, new_fallback)

    # 3. Fix docstring
    old_doc = """Weighted Delaunay triangulation wrapper.
Uses CGAL via pocketlab when enabled, otherwise falls back to SciPy Delaunay."""

    new_doc = """Weighted (regular) Delaunay triangulation.
Requires CGAL via pocketlab. No fallback to standard Delaunay."""

    content = content.replace(old_doc, new_doc)

    path.write_text(content, encoding="utf-8")
    print(f"[FIX] {path}")


def fix_config_py():
    """
    Fix config.py:
    - Rimuovi commenti che menzionano standard vs weighted
    - Fissa flag weighted-only
    """
    path = ROOT / "config.py"
    if not path.exists():
        print(f"[SKIP] {path} non esiste")
        return

    content = path.read_text(encoding="utf-8")

    # 1. Aggiungi commento weighted-only
    old_top = "# config.py - versione pulita"
    new_top = "# config.py - weighted Delaunay only\n# Standard (unweighted) Delaunay support has been removed."

    content = content.replace(old_top, new_top)

    path.write_text(content, encoding="utf-8")
    print(f"[FIX] {path}")


def main():
    print("=== Refactoring weighted-only ===\n")

    fix_alpha_shape_py()
    fix_alpha_stability_py()
    fix_alpha_diagnostics_py()
    fix_geometry_utils_py()
    fix_cgal_wdelaunay_py()
    fix_config_py()

    print("\n=== Refactoring completato ===")
    print("\nProssimi passi:")
    print("1. Verifica i file modificati con git diff")
    print("2. Esegui un test su una singola struttura")
    print("3. Se tutto OK, procedi con il benchmark")


if __name__ == "__main__":
    main()