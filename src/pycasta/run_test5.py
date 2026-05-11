import config, glob, os

config.BOUNDED_DIR   = r"C:\github\pycasta\data\coach420\pdbs"
config.UNBOUNDED_DIR = None
config.OUTPUT_DIR    = r"C:\github\pycasta\src\pycasta\results\coach420_test"
config.VERSION_TAG   = "coach420_test"
config.CORRESPONDENCE_FILE = None

# Solo prime 5 PDB
import glob
all_pdbs = sorted(glob.glob(r"C:\github\pycasta\data\coach420\pdbs\*.pdb"))[:5]
print("Test su:", [os.path.basename(p) for p in all_pdbs])

from run_analysis import run_analysis_mode, summarize_and_print_results
results = run_analysis_mode(pdb_sources=all_pdbs, analysis_type="single")
summary = summarize_and_print_results(results, analysis_type="single",
                                      output_file=None, csv_output_file=None)
print(summary)