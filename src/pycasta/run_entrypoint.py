# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 21:13:13 2025

@author: giorgio
"""

# run_entrypoint.py

# Patch to make run_analysis callable as a function with minimal changes

import os
import glob
import pandas as pd

from run_analysis import run_analysis_mode, summarize_and_print_results
from config import (
    OUTPUT_DIR,
    VERSION_TAG,
    CORRESPONDENCE_FILE,
    UNBOUNDED_DIR,
    BOUNDED_DIR,
)


def run_all_analysis():
    """
    Runs the full analysis and returns the result summary (dict) without writing to file unless needed.
    Can be reused in batch optimization or notebooks.
    """
    versioned_output_dir = os.path.join(OUTPUT_DIR, VERSION_TAG)
    os.makedirs(versioned_output_dir, exist_ok=True)

    bounded_files = glob.glob(os.path.join(BOUNDED_DIR, "*.pdb"))
    if not bounded_files:
        raise FileNotFoundError("No PDB files found in BOUNDED_DIR.")

    # Paired mode
    if (
        UNBOUNDED_DIR
        and os.path.exists(UNBOUNDED_DIR)
        and CORRESPONDENCE_FILE
        and os.path.exists(CORRESPONDENCE_FILE)
    ):
        correspondence_df = pd.read_excel(CORRESPONDENCE_FILE)
        paired_results = run_analysis_mode(
            pdb_sources=bounded_files,
            analysis_type="paired",
            correspondence_df=correspondence_df,
        )
        summary = summarize_and_print_results(
            paired_results,
            analysis_type="paired",
            output_file=None,
            csv_output_file=None,
        )
        return summary

    # Single (bounded) mode
    elif not UNBOUNDED_DIR or not os.path.exists(UNBOUNDED_DIR):
        results = run_analysis_mode(pdb_sources=bounded_files, analysis_type="single")
        summary = summarize_and_print_results(
            results, analysis_type="single", output_file=None, csv_output_file=None
        )
        return summary

    # Unbounded mode
    elif UNBOUNDED_DIR and os.path.exists(UNBOUNDED_DIR) and not CORRESPONDENCE_FILE:
        unbounded_files = glob.glob(os.path.join(UNBOUNDED_DIR, "*.pdb"))
        if not unbounded_files:
            raise FileNotFoundError("No PDB files found in UNBOUNDED_DIR.")
        results = run_analysis_mode(
            pdb_sources=unbounded_files, analysis_type="unbounded"
        )
        summary = summarize_and_print_results(
            results, analysis_type="unbounded", output_file=None, csv_output_file=None
        )
        return summary

    else:
        raise ValueError(
            "Invalid configuration: check UNBOUNDED_DIR and CORRESPONDENCE_FILE settings."
        )
