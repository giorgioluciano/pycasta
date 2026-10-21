# -*- coding: utf-8 -*-
"""
sasa_utils.py

Solvent Accessible Surface Area (SASA) calculations using FreeSASA.
Implements Richards' rolling ball algorithm (1977) via FreeSASA library.

References
----------
Richards FM. Areas, volumes, packing and protein structure. 
Annu Rev Biophys Bioeng. 1977;6:151-176.

FreeSASA: http://freesasa.github.io/
"""

import logging
import numpy as np
import freesasa
from biopandas.pdb import PandasPdb


def compute_sasa_freesasa(pdb_path, probe_radius=1.4):
    """
    Compute SASA for a PDB file using FreeSASA via calcCoord.
    """
    print(f"DEBUG: compute_sasa_freesasa called with pdb_path={pdb_path}")
    import freesasa
    from Bio.PDB import PDBParser
    import numpy as np
    print("DEBUG: imports done")
    
    # Parse PDB
    parser = PDBParser(QUIET=True)
    print("DEBUG: parser created")
    structure = parser.get_structure('protein', pdb_path)
    print("DEBUG: structure loaded")
    
    # Parse PDB
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_path)
    
    # Estrai coordinate e radii
    coords = []
    radii = []
    print(f"DEBUG: starting atom extraction, total atoms: {sum(1 for _ in structure.get_atoms())}")

    # Radii standard (Å)
    atom_radii = {
        'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8,
        'H': 1.2, 'P': 1.8, 'Se': 1.9, 'F': 1.47,
        'Cl': 1.75, 'Br': 1.85, 'I': 1.98,
    }

    for i, atom in enumerate(structure.get_atoms()):
        if i % 500 == 0:
            print(f"DEBUG: processed {i} atoms")
        coord = atom.get_coord()
        coords.extend(coord.tolist())
        element = atom.element.strip().upper()
        radius = atom_radii.get(element, 1.7)
        radii.append(radius)

    print(f"DEBUG: extracted {len(coords)} coords, {len(radii)} radii")
    
    coords = np.array(coords, dtype=float)
    radii = np.array(radii, dtype=float)
    print(f"DEBUG: arrays created, coords shape={coords.shape}, radii shape={radii.shape}")

    # Modifica sasa_utils.py righe 69-79:

    # Calcola SASA
    print("DEBUG: calling freesasa.calcCoord...", flush=True)
    result, _ = freesasa.calcBioPDB(structure)
    print(f"DEBUG: calcBioPDB returned", flush=True)

    total_sasa = result.totalArea()
    residue_sasa = result.residueAreas()
    atom_sasa = None
    return total_sasa, residue_sasa, atom_sasa
    atom_sasa = None


def compute_pocket_sasa(
    pocket_tetra_indices,
    tetra_positions,
    protein_coords=None,
    pdb_file=None,
    probe_radius=1.4,
):
    """Compute SASA for a pocket using FreeSASA."""
    from scipy.spatial import KDTree
    from Bio.PDB import PDBParser
    
    # Get pocket vertices
    pocket_tetra = tetra_positions[pocket_tetra_indices]
    pocket_vertices = pocket_tetra.reshape(-1, 3)
    pocket_vertices = np.unique(pocket_vertices, axis=0)
    
    # Build KD-tree for pocket vertices
    pocket_tree = KDTree(pocket_vertices)
    
    # Find protein atoms within threshold of pocket
    threshold = 5.0
    distances, indices = pocket_tree.query(protein_coords, k=1)
    pocket_atom_indices = set(indices[distances < threshold])
    
    # Parse PDB to get residue mapping
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    # Map atom index to (chain, resid_str)
    atom_to_residue = {}
    for i, atom in enumerate(structure.get_atoms()):
        residue = atom.get_parent()
        chain = residue.get_parent().id
        resid = str(residue.id[1])
        atom_to_residue[i] = (chain, resid)
    
    # Find unique residues in pocket
    pocket_residues = set()
    for idx in pocket_atom_indices:
        if idx in atom_to_residue:
            pocket_residues.add(atom_to_residue[idx])
    
    # Compute SASA for entire protein
    total_sasa, residue_sasa, atom_sasa = compute_sasa_freesasa(pdb_file, probe_radius)
    
    # Sum SASA for pocket-lining residues
    # residue_sasa format: {chain: {resid_str: ResidueArea}}
    pocket_sasa = 0.0
    for chain, res_dict in residue_sasa.items():
        for resid_str, res_area in res_dict.items():
            if (chain, resid_str) in pocket_residues:
                pocket_sasa += res_area.total
    
    return float(pocket_sasa)






def compute_pocket_surface_area_analytic(
    pocket_tetra_indices,
    tetra_positions,
):
    """
    Compute analytic surface area of pocket from tetrahedra.
    
    This computes the surface area of the tetrahedral mesh itself,
    not the SASA.
    
    Parameters
    ----------
    pocket_tetra_indices : list of int
        Indices of tetrahedra forming the pocket.
    tetra_positions : np.ndarray, shape (N, 4, 3)
        Vertex coordinates of all tetrahedra.
    
    Returns
    -------
    surface_area : float
        Surface area in Å².
    """
    pocket_tetra = tetra_positions[pocket_tetra_indices]
    
    # Surface area of each tetrahedron (4 faces)
    total_area = 0.0
    
    for tetra in pocket_tetra:
        # 4 faces
        faces = [
            [tetra[0], tetra[1], tetra[2]],
            [tetra[0], tetra[1], tetra[3]],
            [tetra[0], tetra[2], tetra[3]],
            [tetra[1], tetra[2], tetra[3]],
        ]
        
        for face in faces:
            # Area of triangle
            v0, v1, v2 = face
            edge1 = v1 - v0
            edge2 = v2 - v0
            cross = np.cross(edge1, edge2)
            area = 0.5 * np.linalg.norm(cross)
            total_area += area
    
    # Divide by 2 because each face is shared by 2 tetrahedra (approximately)
    # This is a rough approximation
    return float(total_area) / 2.0

# =============================================================================
# Wrapper per compatibilitàª½ con ranking.py
# =============================================================================

def compute_pocket_surface_area(
    pocket_tetra_indices,
    tetra_positions,
    protein_coords=None,
    probe_radius=1.4,
    n_spheres=100,
):
    """
    Wrapper per compute_pocket_surface_area_analytic per compatibilitàª½.
    
    Parametri:
        pocket_tetra_indices: indici dei tetraedri nella pocket
        tetra_positions: posizioni dei vertici dei tetraedri
        protein_coords: (opzionale, non usato nella versione analitica)
        probe_radius: (opzionale, non usato nella versione analitica)
        n_spheres: (opzionale, non usato nella versione analitica)
    
    Returns:
        float: area superficiale approssimata in ÅÂ²
    """
    return compute_pocket_surface_area_analytic(
        pocket_tetra_indices,
        tetra_positions,
    )

