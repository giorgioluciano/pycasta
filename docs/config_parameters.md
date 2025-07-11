# pycasta Configuration Parameters Reference

This document explains the advanced parameters available in `config.py` for pycasta.
These settings allow you to fine-tune the behavior of cavity detection, merging, validation, and output management.

---

## Merging Settings

- **`MERGE_CLUSTERS`**  
  *(bool)*  
  If `True`, nearby or overlapping pockets (clusters) are merged into a single pocket if within a set threshold.
  *Reduces redundant pockets representing the same region.*

- **`MERGE_THRESHOLD`**  
  *(int)*  
  Cutoff for merging clusters. Controls how close two pockets need to be (distance or mesh units) to be merged together.

- **`set_merge_threshold(value)` / `get_merge_threshold()`**  
  Utility functions to set/get the threshold dynamically from scripts.

---

## Pocket Validation Settings

- **`USE_SASA_CONTACT_VALIDATION`**  
  *(bool)*  
  Enables validation step based on Solvent Accessible Surface Area (SASA) contact between ligand and pocket.

- **`SASA_CONTACT_THRESHOLD`**  
  *(float, Å²)*  
  Minimum SASA contact required between ligand and pocket to be considered valid.

- **`VALIDATION_METHOD`**  
  *(str: `"sasa"`, `"fake_ball"`, `"mesh_extrusion"`)*  
  Selects pocket validation method:
    - `"sasa"`: Validate by SASA contact area.
    - `"fake_ball"`: Use an artificial sphere for validation (benchmarking/testing).
    - `"mesh_extrusion"`: Default. Extrudes pocket mesh toward ligand for geometric validation.

- **`FAKE_SPHERE_RADIUS`**  
  *(float, Å)*  
  Radius of fake sphere for `"fake_ball"` validation.

- **`MESH_EXTRUSION_DISTANCE`**  
  *(float, Å)*  
  Distance to extrude mesh for validation in `"mesh_extrusion"` mode.

- **`SAVE_EXTRUDED_MESHES`**  
  *(bool)*  
  If `True`, saves extruded pocket meshes for further visualization/inspection.

---

## Output and Version Settings

- **`NUM_POCKETS_TO_VALIDATE`**  
  *(int)*  
  Max number of top-ranked pockets to consider for validation.

- **`ALPHA_RANKING`**  
  *(int)*  
  Which alpha rank to use for pocket selection (in multi-alpha workflows).

- **`BETA_RANKING`**  
  *(int)*  
  For workflows using a secondary ranking variable.

- **`NUM_POCKETS_TO_SAVE`**  
  *(int)*  
  Number of top pockets to save to output (PDB, PLY, etc).

- **`USE_EXISTING_RESULTS`**  
  *(bool)*  
  If `True`, skip recalculation if result files already exist.

- **`DEBUG`**  
  *(bool)*  
  Enables detailed log/debug output during analysis.

- **`VERSION_TAG`**  
  *(str)*  
  String to tag this run’s output folders/files (e.g. `"bov1"`).

---

## Tips

- Adjust merging and validation thresholds to tune sensitivity/specificity of pocket detection.
- Set `USE_EXISTING_RESULTS` to `False` to force reprocessing all inputs.
- Use unique `VERSION_TAG`s for each parameter sweep or batch run.
- Always check SASA/contact-based validation settings for biologically meaningful pocket selection.

---

*For any parameter, you can override values directly in your script or notebook for custom runs:*

```python
import config
config.VERSION_TAG = "experiment2025"
config.NUM_POCKETS_TO_SAVE = 10
