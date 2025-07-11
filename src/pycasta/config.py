#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Configuration file for the pocket detection and evaluation pipeline.
Adjust the paths and parameters as needed for your environment.
"""
_current_alpha = None
# =============================================================================
# Directories and Files
# =============================================================================
# Directories (relative or absolute paths)

BOUNDED_DIR = "data/bounded"
UNBOUNDED_DIR = "data/unbounded"  # None

DATASET_DIR = "data"  # Directory containing additional data such as tables
OUTPUT_DIR = "results"

# Correspondence file (for paired analysis)
# CORRESPONDENCE_FILE = "data/tables/correspondence.xlsx"
CORRESPONDENCE_FILE = None

# =============================================================================
# General Processing Settings
# =============================================================================
INCLUDE_WATER = False  # Whether to include water molecules during processing
CONNECTIVITY_THRESHOLD = 1

# =============================================================================
# PyMOL and SASA Settings
# =============================================================================
# Use PyMOL2 if True; otherwise, use standard PyMOL API
USE_PYMOL2 = False
USE_SASA = True
SASA_METHOD = "pymol"  # Options: "pymol2" or "pymol"

# =============================================================================
# Delaunay and CGAL Settings
# =============================================================================
USE_CGAL = True  # Enable weighted Delaunay triangulation
WEIGHTED_DELAUNAY = True  # (Currently a dummy alternative)

# =============================================================================
# Alpha Shape Filtering
# =============================================================================
FILTER_ALPHA_BY_SASA = False  # Whether to filter the alpha complex using SASA
SASA_THRESHOLD = 0  # SASA threshold for filtering

# =============================================================================
# Splitting Settings
# =============================================================================
ATOM_SPLITTING = True  # Enable atom splitting (True/False)
ATOM_SPLITTING_MODE = (
    "CAST"  # Options: "SIMPLE" for equal splitting, "CAST" for a CAST-like split
)
SIMPLE_SPLIT = False

# =============================================================================
# Pocket Calculation Settings
# =============================================================================
CALCULATE_MOUTH_PARAMETERS = True  # Optionally compute pocket mouth parameters
VOLUME_METHOD = "analytic"  # Options: "analytic" or "montecarlo"

# =============================================================================
# Ligand Contact Settings
# =============================================================================
LIGAND_CONTACT_METHOD = "sasa"  # Options: "sasa" or "distance"
LIGAND_CONTACT_THRESHOLD = 4.0

LIGAND_THRESHOLD = 4.0
STRICT_DISTANCE_THRESHOLD = 4.0

# =============================================================================
# Alpha Shape and Flow Settings
# ============================================================================

SEMIAUTO_ALPHA = False
MANUAL_ALPHA = 2.0
ALPHA_PREFERENCE = None


def set_alpha(value):
    global _current_alpha
    _current_alpha = value


def get_alpha():
    if _current_alpha is not None:
        print(f"[CONFIG] get_alpha(): using CURRENT_ALPHA = {_current_alpha}")
        return _current_alpha
    elif ALPHA_PREFERENCE is not None:
        print(f"[CONFIG] get_alpha(): using ALPHA_PREFERENCE = {ALPHA_PREFERENCE}")
        return ALPHA_PREFERENCE
    elif SEMIAUTO_ALPHA:
        print(f"[CONFIG] get_alpha(): using SEMIAUTO_ALPHA = {MANUAL_ALPHA}")
        return MANUAL_ALPHA
    else:
        print(f"[CONFIG] get_alpha(): using MANUAL_ALPHA = {MANUAL_ALPHA}")
        return MANUAL_ALPHA


SAVE_ALPHA = True

MIN_POCKET_VOLUME = 50  # Minimum volume threshold for valid pockets
FILTER_AFTER_MERGE = True

FLOW_METHOD = "discrete_approx"
FLOW_KWARGS = {
    "max_steps": 100,
    "tol": 0.01,
    "sigma_p": 1.4,
    "adaptive": True,  # Enable adaptive tolerance
    "adaptive_factor": 0.5,  # Factor to reduce tolerance when no lower neighbor is found
    "base_step_size": 3,
    "min_steps": 3,
    "debug": True,
    "write_steps": True,
}

# =============================================================================
# Merging Settings
# =============================================================================
MERGE_CLUSTERS = True
MERGE_THRESHOLD = 18


def set_merge_threshold(value):
    global MERGE_THRESHOLD
    MERGE_THRESHOLD = value
    print(
        f"[CONFIG] set_merge_threshold(): updated MERGE_THRESHOLD = {MERGE_THRESHOLD}"
    )


def get_merge_threshold():
    return MERGE_THRESHOLD


# =============================================================================
# Pocket Validation Settings
# =============================================================================
USE_SASA_CONTACT_VALIDATION = False
SASA_CONTACT_THRESHOLD = 4.0

VALIDATION_METHOD = "mesh_extrusion"  # Options: "sasa", "fake_ball", "mesh_extrusion"
FAKE_SPHERE_RADIUS = 50
MESH_EXTRUSION_DISTANCE = 4
SAVE_EXTRUDED_MESHES = True  # Save extruded pocket meshes as .ply or .obj files

# =============================================================================
# Output and Version Settings
# =============================================================================
NUM_POCKETS_TO_VALIDATE = 20
ALPHA_RANKING = 1
BETA_RANKING = 1
NUM_POCKETS_TO_SAVE = 5
USE_EXISTING_RESULTS = True
DEBUG = True

VERSION_TAG = "bov1"
