#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path

"""
Configuration file for the pocket detection and evaluation pipeline.
Paths are resolved from the project root so the code is path-independent.
"""

_current_alpha = None

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"
TESTSHORT_DIR = DATA_DIR / "testShort"
OUTPUT_DIR = PROJECT_ROOT / "output"

BOUNDED_DIR = TESTSHORT_DIR / "bounded"
UNBOUNDED_DIR = TESTSHORT_DIR / "unbounded"
DATASET_DIR = DATA_DIR
CORRESPONDENCE_FILE = TESTSHORT_DIR / "tables" / "correspondence.xlsx"

INCLUDE_WATER = False
CONNECTIVITY_THRESHOLD = 1

USE_SASA = True
SASA_METHOD = "pymol"

USE_CGAL = False
WEIGHTED_DELAUNAY = False

FILTER_ALPHA_BY_SASA = False
SASA_THRESHOLD = 0

ATOM_SPLITTING = True
ATOM_SPLITTING_MODE = "CAST"
SIMPLE_SPLIT = False

CALCULATE_MOUTH_PARAMETERS = True
VOLUME_METHOD = "analytic"

LIGAND_CONTACT_METHOD = "sasa"
LIGAND_CONTACT_THRESHOLD = 4.0
LIGAND_THRESHOLD = 4.0
STRICT_DISTANCE_THRESHOLD = 4.0

SEMIAUTO_ALPHA = False
MANUAL_ALPHA = 2.0
ALPHA_PREFERENCE = None

SAVE_ALPHA = True
MIN_POCKET_VOLUME = 50
FILTER_AFTER_MERGE = True

FLOW_METHOD = "discrete_approx"
FLOW_KWARGS = {
    "max_steps": 100,
    "tol": 0.01,
    "sigma_p": 1.4,
    "adaptive": True,
    "adaptive_factor": 0.5,
    "base_step_size": 3,
    "min_steps": 3,
    "debug": True,
    "write_steps": True,
}

MERGE_CLUSTERS = True
MERGE_THRESHOLD = 18

USE_SASA_CONTACT_VALIDATION = False
SASA_CONTACT_THRESHOLD = 4.0

VALIDATION_METHOD = "mesh_extrusion"
FAKE_SPHERE_RADIUS = 50
MESH_EXTRUSION_DISTANCE = 4
SAVE_EXTRUDED_MESHES = True

NUM_POCKETS_TO_VALIDATE = 20
ALPHA_RANKING = 1
BETA_RANKING = 1
NUM_POCKETS_TO_SAVE = 5
USE_EXISTING_RESULTS = True
DEBUG = True

VERSION_TAG = "smoke"
_current_alpha = None
_current_merge_threshold = None


def get_alpha():
    return MANUAL_ALPHA if _current_alpha is None else _current_alpha


def set_alpha(value):
    global _current_alpha
    _current_alpha = float(value)


def reset_alpha():
    global _current_alpha
    _current_alpha = None


def get_merge_threshold():
    return MERGE_THRESHOLD if _current_merge_threshold is None else _current_merge_threshold


def set_merge_threshold(value):
    global _current_merge_threshold
    _current_merge_threshold = float(value)


def reset_merge_threshold():
    global _current_merge_threshold
    _current_merge_threshold = None