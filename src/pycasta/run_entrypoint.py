import glob
import logging
import os
import pandas as pd

from .config import BOUNDED_DIR, UNBOUNDED_DIR, CORRESPONDENCE_FILE, OUTPUT_DIR, VERSION_TAG


def run_all_analysis():
    versioned_output_dir = os.path.join(str(OUTPUT_DIR), VERSION_TAG)
    os.makedirs(versioned_output_dir, exist_ok=True)

    bounded_files = glob.glob(os.path.join(str(BOUNDED_DIR), "*.pdb"))
    if not bounded_files:
        raise FileNotFoundError("No PDB files found in BOUNDED_DIR.")

    unbounded_exists = UNBOUNDED_DIR is not None and os.path.exists(str(UNBOUNDED_DIR))
    correspondence_exists = (
        CORRESPONDENCE_FILE is not None and os.path.exists(str(CORRESPONDENCE_FILE))
    )

    if unbounded_exists and correspondence_exists:
        correspondence_df = pd.read_excel(str(CORRESPONDENCE_FILE))
        paired_results = run_analysis_mode(
            pdb_sources=bounded_files,
            analysis_type="paired",
            correspondence_df=correspondence_df,
        )
        return summarize_and_print_results(
            paired_results,
            analysis_type="paired",
            output_file=None,
            csv_output_file=None,
        )

    elif unbounded_exists and not correspondence_exists:
        unbounded_files = glob.glob(os.path.join(str(UNBOUNDED_DIR), "*.pdb"))
        if not unbounded_files:
            raise FileNotFoundError("No PDB files found in UNBOUNDED_DIR.")
        results = run_analysis_mode(
            pdb_sources=unbounded_files,
            analysis_type="unbounded",
        )
        return summarize_and_print_results(
            results,
            analysis_type="unbounded",
            output_file=None,
            csv_output_file=None,
        )

    else:
        results = run_analysis_mode(
            pdb_sources=bounded_files,
            analysis_type="single",
        )
        return summarize_and_print_results(
            results,
            analysis_type="single",
            output_file=None,
            csv_output_file=None,
        )
        
        
        