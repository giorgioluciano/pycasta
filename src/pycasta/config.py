# config.py - weighted Delaunay only
# Standard (unweighted) Delaunay support has been removed.
from pathlib import Path

# Alpha e merge (parametri runtime modificabili)
_current_alpha = None
_current_merge_threshold = None
_current_min_pocket_volume = None

# Paths
PROJECT_ROOT = Path(r"C:\github\pycasta")
DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR = r"C:\github\pycasta\results"
VERSION_TAG = "pocketpicker_classic_alpha28"


BOUNDED_DIR = r"C:\github\pycasta\data\pocketpicker\holo"
UNBOUNDED_DIR = None
DATASET_DIR = DATA_DIR
CORRESPONDENCE_FILE = None

# Geometria
INCLUDE_WATER = False
CONNECTIVITY_THRESHOLD = 1
USE_CGAL = False


# ============================================
# SASA Configuration (Richards surface, 1977)
# ============================================
SASA_PROBE_RADIUS = 1.4  # Å, standard water probe
SASA_METHOD = "freesasa"

# Atom splitting
ATOM_SPLITTING = True
ATOM_SPLITTING_MODE = "CAST"
SIMPLE_SPLIT = False

# Mouth parameters
CALCULATE_MOUTH_PARAMETERS = True
VOLUME_METHOD = "analytic"

# Ligand contact
LIGAND_CONTACT_METHOD = "sasa"
LIGAND_CONTACT_THRESHOLD = 4.0
LIGAND_THRESHOLD = 4.0
STRICT_DISTANCE_THRESHOLD = 4.0

# Alpha selection
SEMIAUTO_ALPHA = False
MANUAL_ALPHA = 2.8
ALPHA_PREFERENCE = None

# Pocket filtering
SAVE_ALPHA = True
MIN_POCKET_VOLUME = 50
FILTER_AFTER_MERGE = True

# Flow detection
FLOW_METHOD = "discrete_approx"
FLOW_KWARGS = {
    "max_steps": 100,
    "tol_fraction": 0.01,
    "sigma_p": 1.4,
    "adaptive": True,
    "adaptive_factor": 0.5,
    "base_step_size": 3,
    "min_steps": 3,
    "debug": False,  # disattiva per batch
    "write_steps": False,
}

# Merge clusters
MERGE_CLUSTERS = True
MERGE_THRESHOLD = 10

# Validation
VALIDATION_METHOD = "mesh_extrusion"
FAKE_SPHERE_RADIUS = 8
MESH_EXTRUSION_DISTANCE = 10
SAVE_EXTRUDED_MESHES = True
NUM_POCKETS_TO_VALIDATE = 8
USE_EXISTING_RESULTS = False
DEBUG = True

# Ranking
RANKING_METHOD = "percentile_log"
ALPHA_RANKING = 1.0
BETA_RANKING = 1.0
NUM_POCKETS_TO_SAVE = 5

# Graph filtering (opzionale)
USE_GRAPH_FILTERED_ALPHA = False
ALPHA_MAX_EDGE_FACTOR = 2.5
ALPHA_LOCAL_SCALE_K = 6

# Debug triangulation
SAVE_TRIANGULATION_DEBUG = False  # disattiva per batch

# Helper functions
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

def get_min_pocket_volume():
    return MIN_POCKET_VOLUME if _current_min_pocket_volume is None else _current_min_pocket_volume

def set_min_pocket_volume(value):
    global _current_min_pocket_volume
    _current_min_pocket_volume = float(value)

def reset_min_pocket_volume():
    global _current_min_pocket_volume
    _current_min_pocket_volume = None

# Filtri per AlphaShape
FILTER_ALPHA_BY_SASA = False


# =============================================================================
# Ranking weights
# =============================================================================
VOLUME_WEIGHT = 0.4
FLOW_WEIGHT = 0.3
CONNECTIVITY_WEIGHT = 0.2
SASA_WEIGHT = 0.1

# SASA parameters
SASA_PROBE_RADIUS = 1.4  # Angstrom
SASA_N_SPHERES = 100
SASA_THRESHOLD = 0.0  # Soglia minima SASA per pocket

