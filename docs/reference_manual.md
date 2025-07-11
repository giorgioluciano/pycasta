# pycasta Reference Manual

A technical reference for the main modules and functions in the pycasta project.

---

## `alpha_shape.py`

### Purpose
Functions for computing alpha shapes from protein structures and performing discrete flow pocket detection.

### Functions

#### `filter_alpha_complex_by_sasa(tetrahedra, alpha_values, protein_coords, pdb_file, sasa_threshold)`
- **Purpose:**  
  Filter tetrahedra in the alpha complex by Solvent Accessible Surface Area (SASA) of residues. Only tetrahedra whose centroid maps to residues with SASA below `sasa_threshold` are retained.
- **Parameters:**  
    - `tetrahedra`: ndarray (N,4,3) – Cartesian coordinates of tetrahedra.
    - `alpha_values`: ndarray (N,) – Circumsphere radii for each tetrahedron.
    - `protein_coords`: ndarray (M,3) – Coordinates of all protein atoms.
    - `pdb_file`: str – Path to the original PDB.
    - `sasa_threshold`: float – SASA cutoff.
- **Returns:**  
    - Filtered tetrahedra, alpha_values, and mask (boolean array of kept indices).

#### `compute_alpha_complex_from_tetrahedra(simplices, tetra_positions, alpha_value, molecule_name=None, protein_coords=None, filter_by_sasa=False, sasa_threshold=SASA_CONTACT_THRESHOLD, pdb_file=None)`
- **Purpose:**  
  Compute the boolean mask for the alpha complex from a set of tetrahedra, with optional SASA filtering and export of alpha data.
- **Parameters:**  
    - `simplices`: ndarray (K,4) – Tetrahedron indices.
    - `tetra_positions`: ndarray (K,4,3) – Vertex positions.
    - `alpha_value`: float – Alpha parameter (Å).
    - `molecule_name`: str (optional) – For output filenames.
    - `protein_coords`: ndarray (M,3) (optional).
    - `filter_by_sasa`: bool (optional).
    - `sasa_threshold`: float.
    - `pdb_file`: str.
- **Returns:**  
    - `alpha_mask`: Boolean mask of tetrahedra kept.
    - `radii`: Circumsphere radii.
    - `tetra_positions`, `simplices`: Possibly filtered.

#### `flow_detection(simplices, alpha_mask, tetra_positions, protein_coords, flow_params=None)`
- **Purpose:**  
  Run discrete flow (iterative neighbor descent) to group tetrahedra into pockets based on geometric/topological features.
- **Parameters:**  
    - `simplices`: ndarray (K,4).
    - `alpha_mask`: bool array.
    - `tetra_positions`: ndarray (K,4,3).
    - `protein_coords`: ndarray.
    - `flow_params`: dict (optional) – Hyperparameters (`sigma_p`, `tol_fraction`, etc.).
- **Returns:**  
    - `final_pockets`: list of lists of tetrahedron indices (each pocket).
    - `flow_steps`: dict (tetrahedron idx → number of flow steps).
    - `connectivity_graph`: dict for each pocket.

---

## `atomic_radii.py`

### Purpose
Dictionary and default for atomic radii (in Å) for elements, used to compute weighted Delaunay/alpha shapes.

### Constants

#### `ATOMIC_RADII_DICT`
- **Purpose:**  
  Mapping from element symbol (str) to radius (float, Å).

#### `DEFAULT_RADIUS`
- **Purpose:**  
  Default radius to use if element is not found.

---

## `cgal_wdelaunay.py`

### Purpose
(Placeholder) for weighted Delaunay triangulation via CGAL.  
Currently falls back to `scipy.spatial.Delaunay`.

### Functions

#### `cgal_weighted_delaunay(protein_coords, radii)`
- **Purpose:**  
  Perform (optionally weighted) Delaunay triangulation of atomic coordinates.
- **Parameters:**  
    - `protein_coords`: ndarray (N,3).
    - `radii`: ndarray (N,) – atomic radii.
- **Returns:**  
    - `simplices`: ndarray (M,4) – Indices of vertices of each tetrahedron.
    - None (for future use).

---

## `check_pocket_stats.py`

### Purpose
Script for batch evaluation of detection performance across JSON result files.

### Main Logic

- **Scans** all `.json` files in a results directory.
- **Tallies** the number of files where the first, first three, or any pockets are validated (using `ligand_containment_mesh`).
- **Prints** summary statistics (top1/top3/topALL found).

---

## `optimize_alpha.py`

### Purpose
Tools for hyperparameter optimization: scan alpha and flow parameters, record pocket detection/validation results.

### Functions

#### `set_flow_params(min_steps=None, min_volume=None, adaptive_factor=None, tol=None)`
- **Purpose:**  
  Dynamically update global flow parameters for batch runs.

#### `optimize_alpha_and_flow(alpha_values, min_steps_values, min_volume_values, adaptive_factor_values, tol_values, merge_thresholds, output_csv="alpha_flow_optimization_results.csv")`
- **Purpose:**  
  Batch test all combinations of parameters, recording performance metrics (top1/top3/step-to-ligand, etc.) to CSV.
- **Parameters:**  
    - `alpha_values`: list of float.
    - `min_steps_values`, etc.: lists of candidate hyperparameter values.
    - `output_csv`: str – CSV output.
- **Returns:**  
    - `df`: pandas DataFrame of results.

---

## `pocket_detection.py`

### Purpose
Core logic for pocket detection and ranking using alpha shapes and flow clustering.

### Functions

#### `compute_analytic_pocket_volume(tetra_positions, pocket_indices)`
- **Purpose:**  
  Compute analytic (3D) volume for a set of tetrahedra in a pocket.
- **Parameters:**  
    - `tetra_positions`: ndarray (K,4,3).
    - `pocket_indices`: list of int.
- **Returns:**  
    - `float`: total volume.

#### `detect_pockets(protein_coords, simplices=None, tetra_positions=None, alpha_mask=None, min_volume_threshold=None, flow_params=None, merge_clusters=None, merge_threshold=None, molecule_name=None, radii=None)`
- **Purpose:**  
  Complete pipeline for pocket detection, including merging and ranking.
- **Parameters:**  
    - `protein_coords`, `simplices`, `alpha_mask`, `flow_params`, etc.
- **Returns:**  
    - Dict with keys: `ranked_pockets`, `representative_points`, `ranking_scores`, `dual_sets_info`, `mouths_info`.

---

## `run_analysis.py`

### Purpose
Main script to run the entire pipeline on PDB files in single or paired mode, saving results and summary statistics.

### Functions

#### `compute_pocket_distances(pocket_points, ligand_coords)`
- **Purpose:**  
  Compute min/max/average distances between pocket and ligand points.
- **Parameters:**  
    - `pocket_points`: ndarray.
    - `ligand_coords`: ndarray.
- **Returns:**  
    - dict: `{min, max, avg}` distances.

#### `match_atoms_by_residue_and_name(fixed_atoms, moving_atoms, allowed_atom_names=["CA"])`
- **Purpose:**  
  Utility for matching atoms between two protein structures (for superposition).
- **Parameters:**  
    - `fixed_atoms`, `moving_atoms`: list of dicts.
    - `allowed_atom_names`: list of str.
- **Returns:**  
    - Lists of matched coordinates.

#### `align_using_superimposer(fixed_coords, moving_coords, return_transform=False)`
- **Purpose:**  
  Superimpose two sets of atomic coordinates (usually for alignment of bound/unbound proteins).
- **Parameters:**  
    - `fixed_coords`, `moving_coords`: list/array.
    - `return_transform`: bool.
- **Returns:**  
    - `aligned`, `rmsd`, (`R`, `t`) if requested.

#### `analyze_top_pockets(result)`
- **Purpose:**  
  Determine the top validated pocket index and validation methods used for each pocket.
- **Parameters:**  
    - `result`: dict from process_pdb.
- **Returns:**  
    - `(top_index, validation_tags)`.

#### `run_analysis_mode(pdb_sources, analysis_type="single", correspondence_df=None)`
- **Purpose:**  
  High-level function to process all files (single, paired, unbounded modes).
- **Parameters:**  
    - `pdb_sources`: list of str.
    - `analysis_type`: "single", "paired", "unbounded".
    - `correspondence_df`: DataFrame (for paired).
- **Returns:**  
    - List of results/summary dicts.

#### `get_output_filename(pdb_path, version_tag)`
- **Purpose:**  
  Compute output JSON filename for a PDB and version.

#### `process_pdb(pdb_path)`
- **Purpose:**  
  Complete pipeline for a single PDB: separation, alpha, pocket detection, validation, output.
- **Returns:**  
    - Dict of results (with all keys).

#### `main()`
- **Purpose:**  
  Entrypoint: run the whole analysis (single or paired).

---

# Notes

- For utility/internal functions, see the `utils/` directory.
- All output files are saved according to config parameters (`OUTPUT_DIR`, `VERSION_TAG`).
- See `config_parameters.md` for configuration options.
- Some modules expect specific PDB formatting (e.g., ATOM/HETATM naming).

---

_Last updated: July 2025_

