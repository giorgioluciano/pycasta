
# pycasta

A Python package for the detection and analysis of protein cavities.

---

## Installation

To set up the dependencies for pycasta, you can use either pip (with requirements.txt) or conda (with environment.yml).

### Using pip

Make sure you are in the root directory of the repository (where `requirements.txt` is located), then run:

```bash
pip install -r requirements.txt
```

### Using conda (recommended for scientific workflows)

Make sure you are in the root directory of the repository (where `environment.yml` is located), then run:

```bash
conda env create -f environment.yml
conda activate pycasta-env
```

> **Note:**  
> Some dependencies such as `freesasa` and `pymol2` may require additional system libraries or specific installation steps, especially on Windows.  
> If you encounter installation issues, please consult the official documentation of those packages.

---

## Quick Start

After installing the dependencies, you can run the included example analysis with:

```bash
cd src
python run_analysis.py
```

This will execute a demo analysis with the default settings.  
Check the output and log messages for results or any errors.

---

## Example Data

A small example dataset is included in the `data/` directory for testing and demonstration purposes.

The `data/` folder contains three subdirectories:

- `bounded/` – for analysis of molecules in the **bound** state (paired analysis)
- `unbounded/` – for analysis of molecules in the **unbound** state (paired analysis)
- `tables/` – contains Excel files that map the correspondence between each bound and unbound molecule for paired analyses

### Paired Analysis

If you want to perform a paired analysis (comparing bound and unbound forms),  
you **must provide an Excel file** (for example: `mapping.xlsx`) inside the `tables/` directory.  
This file should specify the correspondence between each bound and unbound molecule.

The Excel file should contain at least two columns:

| bound_molecule     | unbound_molecule   |
|:------------------:|:-----------------:|
| 1abc_bound.pdb     | 1abc_unbound.pdb  |
| 2xyz_bound.pdb     | 2xyz_unbound.pdb  |
| ...                | ...               |

> **Place the mapping file inside the `tables/` directory.**

### Single Analysis

For a single-molecule analysis (not paired), simply use a structure where the heteroatom is already present.  
You can place such files directly in the appropriate folder.

---

Feel free to use or modify the example data to test different types of analyses.  
For larger or custom datasets, please follow the same folder structure and file naming conventions.

---





