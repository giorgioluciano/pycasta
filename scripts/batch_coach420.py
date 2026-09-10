#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Batch grid search su tutto COACH420.
"""

import sys
from pathlib import Path
import numpy as np
import json
import logging
from scipy.spatial import KDTree
from biopandas.pdb import PandasPdb

# Setup paths
PROJECT_ROOT = Path(r"C:\github\pycasta")
SRC_DIR = PROJECT_ROOT / "src"
PDB_DIR = PROJECT_ROOT / "data" / "coach420" / "pdbs"
RESULTS_DIR = PROJECT_ROOT / "results" / "coach420_full420"

sys.path.insert(0, str(SRC_DIR))

from pycasta.run_analysis import process_pdb
from pycasta import config as cfg

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def load_pdb_ligand(pdb_path):
    """Carica PDB e estrae ligando."""
    ppdb = PandasPdb().read_pdb(str(pdb_path))
    
    protein_df = ppdb.df['ATOM']
    protein_coords = protein_df[['x_coord', 'y_coord', 'z_coord']].values
    
    hetatm_df = ppdb.df['HETATM']
    if len(hetatm_df) > 0:
        ligand_df = hetatm_df[hetatm_df['residue_name'] != 'HOH']
        if len(ligand_df) > 0:
            ligand_coords = ligand_df[['x_coord', 'y_coord', 'z_coord']].values
        else:
            ligand_coords = None
    else:
        ligand_coords = None
    
    pdb_id = pdb_path.stem
    return protein_coords, ligand_coords, pdb_id


def compute_pocket_ligand_distances(pockets, tetra_positions, ligand_coords):
    """Calcola distanze tasca-ligando."""
    if ligand_coords is None or len(ligand_coords) == 0:
        return None
    
    ligand_tree = KDTree(ligand_coords)
    pocket_distances = []
    
    for pocket_idx, pocket_tetra_indices in enumerate(pockets):
        vertices = []
        for ti in pocket_tetra_indices:
            tetra = tetra_positions[ti]
            vertices.extend(tetra)
        vertices = np.array(vertices)
        
        distances, _ = ligand_tree.query(vertices)
        
        pocket_distances.append({
            'pocket_idx': pocket_idx,
            'min_dist': float(distances.min()),
            'mean_dist': float(distances.mean()),
            'n_within_5A': int(np.sum(distances < 5.0)),
            'n_within_3A': int(np.sum(distances < 3.0)),
        })
    
    return pocket_distances


def run_alpha_grid(pdb_path, alpha_values=None):
    """Esegue grid search per una proteina."""
    if alpha_values is None:
        alpha_values = [1.68, 2.0, 2.35, 3.27, 4.0, 5.0]
    
    protein_coords, ligand_coords, pdb_id = load_pdb_ligand(pdb_path)
    
    logging.info(f"Processing {pdb_id} (ligand: {len(ligand_coords) if ligand_coords is not None else 0} atoms)")
    
    results = []
    
    for alpha in alpha_values:
        logging.info(f"  Alpha = {alpha:.2f}")
        
        original_alpha = cfg.MANUAL_ALPHA
        cfg.MANUAL_ALPHA = alpha
        
        try:
            result = process_pdb(str(pdb_path))
            
            if result is None:
                logging.warning(f"  process_pdb returned None")
                results.append({'alpha': alpha, 'error': 'process_pdb returned None'})
                cfg.MANUAL_ALPHA = original_alpha
                continue
            
            pockets = result.get('ranked_pockets', [])
            
            # Carica tetra_positions
            import os
            triangulation_file = f"C:\\github\\pycasta\\src\\pycasta\\results\\coach420_test\\{pdb_id}.triangulation_debug.npz"
            
            if os.path.exists(triangulation_file):
                triangulation_data = np.load(triangulation_file, allow_pickle=True)
                tetra_positions = triangulation_data['tetrahedra']
            else:
                logging.warning(f"  Triangulation file not found: {triangulation_file}")
                tetra_positions = None
            
            # Calcola distanze
            if tetra_positions is not None and ligand_coords is not None and len(pockets) > 0:
                pocket_distances = compute_pocket_ligand_distances(pockets, tetra_positions, ligand_coords)
            else:
                pocket_distances = None
            
            # Trova best pocket
            if pocket_distances:
                best_pocket = min(pocket_distances, key=lambda x: x['min_dist'])
                best_dist = best_pocket['min_dist']
                best_pocket_idx = best_pocket['pocket_idx']
            else:
                best_dist = None
                best_pocket_idx = None
            
            results.append({
                'alpha': alpha,
                'n_pockets': len(pockets),
                'best_pocket_idx': best_pocket_idx,
                'best_min_dist': best_dist,
                'pocket_distances': pocket_distances,
            })
            
            logging.info(f"  N pockets: {len(pockets)}, Best dist: {best_dist:.2f} Å" if best_dist else f"  N pockets: {len(pockets)}")
            
        except Exception as e:
            logging.error(f"  ERROR: {e}")
            results.append({'alpha': alpha, 'error': str(e)})
        
        finally:
            cfg.MANUAL_ALPHA = original_alpha
    
    return results, pdb_id, ligand_coords


def save_results(results, pdb_id, output_dir=None):
    """Salva risultati JSON."""
    if output_dir is None:
        output_dir = RESULTS_DIR
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = output_dir / f"{pdb_id}_alpha_grid.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    return output_file


def run_batch(pdb_dir=None, results_dir=None):
    """Esegue batch su tutte le proteine (resumable)."""
    if pdb_dir is None:
        pdb_dir = PDB_DIR
    if results_dir is None:
        results_dir = RESULTS_DIR
    
    # Trova tutti i PDB
    pdb_files = list(Path(pdb_dir).glob("*.pdb"))
    logging.info(f"Found {len(pdb_files)} PDB files")
    
    # Controlla quali sono già stati processati
    already_done = []
    to_process = []
    
    for pdb_file in pdb_files:
        pdb_id = pdb_file.stem
        output_file = Path(results_dir) / f"{pdb_id}_alpha_grid.json"
        
        if output_file.exists():
            # Controlla se il file è valido (non vuoto)
            try:
                with open(output_file, 'r') as f:
                    data = json.load(f)
                if len(data) > 0:
                    already_done.append(pdb_id)
                    logging.info(f"Skipping {pdb_id} (already done)")
                    continue
            except:
                pass
        
        to_process.append(pdb_file)
    
    logging.info(f"\nAlready done: {len(already_done)}/{len(pdb_files)}")
    logging.info(f"To process: {len(to_process)}/{len(pdb_files)}")
    
    if len(to_process) == 0:
        logging.info("All proteins already processed!")
        return
    
    # Summary (carica se esiste già)
    summary_file = Path(results_dir) / "batch_summary.json"
    if summary_file.exists():
        with open(summary_file, 'r') as f:
            summary = json.load(f)
        logging.info(f"Loaded existing summary with {len(summary)} entries")
    else:
        summary = []
    
    # Processa
    for i, pdb_file in enumerate(to_process):
        logging.info(f"\n{'='*60}")
        logging.info(f"Processing {i+1}/{len(to_process)}: {pdb_file.name}")
        logging.info(f"{'='*60}")
        
        try:
            results, pdb_id, ligand_coords = run_alpha_grid(pdb_file)
            save_results(results, pdb_id, results_dir)
            
            # Estrai best result
            valid_results = [r for r in results if 'error' not in r and r['best_min_dist'] is not None]
            if valid_results:
                best = min(valid_results, key=lambda x: x['best_min_dist'])
                summary.append({
                    'pdb_id': pdb_id,
                    'has_ligand': ligand_coords is not None,
                    'best_alpha': best['alpha'],
                    'best_min_dist': best['best_min_dist'],
                    'n_pockets_at_best': best['n_pockets'],
                })
            else:
                summary.append({
                    'pdb_id': pdb_id,
                    'has_ligand': ligand_coords is not None,
                    'error': 'No valid results',
                })
            
        except Exception as e:
            logging.error(f"Failed {pdb_file.name}: {e}")
            summary.append({
                'pdb_id': pdb_file.stem,
                'error': str(e),
            })
        
        # Salva summary dopo ogni proteina (così se si pianta, non perdi tutto)
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
    
    logging.info(f"\nBatch completed! Summary saved to {summary_file}")
    
    # Stats
    valid = [s for s in summary if 'error' not in s]
    logging.info(f"Valid results: {len(valid)}/{len(summary)}")
    
    if valid:
        best_alphas = [s['best_alpha'] for s in valid]
        logging.info(f"Best alpha distribution:")
        for alpha in sorted(set(best_alphas)):
            count = best_alphas.count(alpha)
            logging.info(f"  Alpha {alpha:.2f}: {count}/{len(valid)} ({100*count/len(valid):.1f}%)")

if __name__ == "__main__":
    run_batch()