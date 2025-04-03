# PyCASTa

**PyCASTa** (Pocket and Cavity Analysis with Shape-based Techniques) is a scientific pipeline for detecting and validating protein binding pockets using Delaunay triangulation, alpha shapes, and ligand-based validation.  
It supports both single (bounded) and paired protein structure analysis.

---

## 🚀 Features

- Alpha shape-based pocket detection with Delaunay triangulation.
- Customizable pocket ranking and volume estimation.
- Ligand validation with:
  - **SASA contact**
  - **Mesh extrusion** (geometry-based)
  - **Fake sphere** contact
- Support for paired structure alignment and RMSD calculation.
- Batch optimization of detection parameters (alpha, flow, merging).

---

## 📦 Installation

Clone the repository and install locally using `pip`:

```bash
git clone https://github.com/your-username/pycast.git
cd pycast
pip install -e .

