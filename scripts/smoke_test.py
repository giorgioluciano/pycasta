from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pycasta import config
from pycasta.run_analysis import main as run_main

DATA = ROOT / "data" / "testshort"
OUT = ROOT / "output" / "smoke_test"
OUT.mkdir(parents=True, exist_ok=True)

config.BOUNDED_DIR = DATA / "bounded"
config.UNBOUNDED_DIR = DATA / "unbounded"
config.CORRESPONDENCE_FILE = DATA / "correspondence.xlsx"
config.OUTPUT_DIR = OUT
config.USE_SASA = True
config.SASA_METHOD = "freesasa"
config.USE_CGAL = False
config.FILTER_ALPHA_BY_SASA = False
config.VALIDATION_METHOD = "mesh_extrusion"
config.USE_EXISTING_RESULTS = False
config.DEBUG = True
config.VERSION_TAG = "smoke"

if __name__ == "__main__":
    run_main()