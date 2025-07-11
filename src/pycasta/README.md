# pycasta

Python package for detection and analysis of protein cavities.

------------------------------------------------------------------------

## Installation

-   **Recommended:** Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/).

-   From the repository root, create and activate the environment:

    ``` bash
    conda env create -f environment.yml
    conda activate pycasta-env
    ```

    *or, alternatively:*

    ``` bash
    pip install -r requirements.txt
    ```

    > Note: Some libraries (`freesasa`, `pymol2`, etc.) may require extra steps on Windows.

------------------------------------------------------------------------

## Example Usage

### 1. Using Jupyter Notebook (Recommended)

-   From your environment, launch:

    ``` bash
    jupyter notebook
    ```

-   Open the provided example notebook for step-by-step instructions and visualization.

### 2. Command-line Usage

-   To run from the console:

    ``` bash
    cd source
    python run_analysis.py
    ```

-   Analysis options can be configured in `config.py`.

-   For paired analysis (bound/unbound), ensure that you have set up the directories and correspondence file as below.

------------------------------------------------------------------------

## Example Data

A sample dataset is included for demonstration:

-   `data/bounded/` — contains bound-state PDBs:\
    `1a6w.pdb`, `1qif.pdb`, `3app.pdb`
-   `data/unbounded/` — contains unbound-state PDBs:\
    `1a6u.pdb`, `1acj.pdb`, `apu.pdb`
-   `data/tables/correspondence.xlsx` — Excel file mapping each bound to its unbound structure, with columns: \| bound_molecule \| unbound_molecule \| \|:--------------:\|:----------------:\| \| 1a6w.pdb \| 1a6u.pdb \| \| 1qif.pdb \| 1acj.pdb \| \| 3app.pdb \| apu.pdb \|

For single-molecule analysis, simply use a PDB file (with heteroatom) in `data/bounded/`.

------------------------------------------------------------------------

## File Placement

-   Place `requirements.txt` and `environment.yml` in the root of the repository.

------------------------------------------------------------------------

## Need Help?

-   Check the example notebook for a guided workflow.
-   For support, open an issue on GitHub.
