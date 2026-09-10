# tools/occlusion_viewer.py
"""
Minimal occlusion viewer for PyCASTA.
Shows protein atoms colored by occlusion fraction, with ligand overlaid.
"""

import json
import os
import sys
import numpy as np
import plotly.graph_objects as go
from pathlib import Path


def load_pdb_atoms(pdb_file):
    """Load atoms from PDB file (simple parser)."""
    atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atom_name = line[12:16].strip()
                residue = line[17:20].strip()
                chain = line[21]
                resid = int(line[22:26])
                element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                is_hetero = line.startswith('HETATM')
                atoms.append({
                    'coord': [x, y, z],
                    'atom_name': atom_name,
                    'residue': residue,
                    'chain': chain,
                    'resid': resid,
                    'element': element,
                    'is_hetero': is_hetero,
                })
    return atoms


def compute_occlusion_for_atoms(atoms, num_rays=64, max_distance=15.0):
    """
    Compute occlusion fraction for each atom.
    Occlusion = fraction of rays that hit other atoms.
    """
    coords = np.array([a['coord'] for a in atoms], dtype=float)
    n_atoms = len(atoms)
    occlusion_values = np.zeros(n_atoms)
    
    # Simple neighbor-based occlusion (faster than full raycasting for visualization)
    from scipy.spatial import KDTree
    tree = KDTree(coords)
    
    for i in range(n_atoms):
        center = coords[i]
        # Find neighbors within max_distance
        neighbors = tree.query_ball_point(center, r=max_distance)
        # Occlusion ~ fraction of solid angle covered by neighbors
        # Simplified: count neighbors and their relative sizes
        if len(neighbors) > 1:
            # More neighbors = more occluded
            occlusion_values[i] = min(1.0, (len(neighbors) - 1) / 20.0)  # Normalize
    
    return occlusion_values


def make_occlusion_figure(pdb_file, ligand_residues=None, output_html=None):
    """
    Create an interactive HTML visualization of occlusion.
    
    Parameters
    ----------
    pdb_file : str
        Path to PDB file.
    ligand_residues : list, optional
        List of residue names/IDs that are ligands (e.g., ['ATP', 'HEM']).
        If None, all HETATM are considered ligands.
    output_html : str, optional
        Output HTML file path. If None, returns figure object.
    """
    # Load atoms
    atoms = load_pdb_atoms(pdb_file)
    coords = np.array([a['coord'] for a in atoms], dtype=float)
    
    # Separate protein and ligand
    protein_atoms = [a for a in atoms if not a['is_hetero']]
    ligand_atoms = [a for a in atoms if a['is_hetero']]
    
    if ligand_residues is not None:
        ligand_atoms = [a for a in ligand_atoms if a['residue'] in ligand_residues]
    
    protein_coords = np.array([a['coord'] for a in protein_atoms], dtype=float)
    ligand_coords = np.array([a['coord'] for a in ligand_atoms], dtype=float)
    
    # Compute occlusion for protein atoms
    if len(protein_coords) > 0:
        occlusion_values = compute_occlusion_for_atoms(protein_atoms)
    else:
        occlusion_values = np.array([])
    
    # Create figure
    fig = go.Figure()
    
    # Protein atoms (colored by occlusion)
    if len(protein_coords) > 0:
        fig.add_trace(go.Scatter3d(
            x=protein_coords[:, 0],
            y=protein_coords[:, 1],
            z=protein_coords[:, 2],
            mode='markers',
            marker=dict(
                size=3,
                color=occlusion_values,
                colorscale='RdYlBu_r',  # Red = occluded, Blue = exposed
                colorbar=dict(title='Occlusion'),
                opacity=0.8,
            ),
            name='Protein (occlusion)',
            text=[f"{a['residue']} {a['resid']}" for a in protein_atoms],
        ))
    
    # Ligand atoms (as spheres)
    if len(ligand_coords) > 0:
        fig.add_trace(go.Scatter3d(
            x=ligand_coords[:, 0],
            y=ligand_coords[:, 1],
            z=ligand_coords[:, 2],
            mode='markers',
            marker=dict(
                size=6,
                color='red',
                opacity=0.9,
            ),
            name='Ligand',
        ))
    
    # Layout
    fig.update_layout(
        title=f"Occlusion Viewer - {Path(pdb_file).stem}",
        scene=dict(
            aspectmode='data',
            xaxis_title='X (Å)',
            yaxis_title='Y (Å)',
            zaxis_title='Z (Å)',
        ),
        margin=dict(l=0, r=0, t=40, b=0),
        height=700,
    )
    
    # Save or return
    if output_html:
        fig.write_html(output_html, include_plotlyjs='cdn')
        print(f"Occlusion viewer saved to: {output_html}")
        return output_html
    else:
        return fig


def main():
    if len(sys.argv) < 2:
        print("Usage: python tools/occlusion_viewer.py <pdb_file> [output_html]")
        print("Example: python tools/occlusion_viewer.py data/bounded/1a6w.pdb occlusion_1a6w.html")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    output_html = sys.argv[2] if len(sys.argv) > 2 else None
    
    if output_html is None:
        base_name = Path(pdb_file).stem
        output_html = f"occlusion_{base_name}.html"
    
    make_occlusion_figure(pdb_file, output_html=output_html)


if __name__ == "__main__":
    main()