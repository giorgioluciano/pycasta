from pathlib import Path
from run_analysis import process_pdb

pdb_dir = Path(r"C:\github\pycasta\data\coach420\pdbs")
pdbs = sorted(pdb_dir.glob("*.pdb"))[:3]

for pdb in pdbs:
    print(f"\n=== {pdb.name} ===")
    result = process_pdb(str(pdb))
    print("step_to_ligand:", result.get("step_to_ligand"))
    print("backend:", result.get("triangulation_backend", "no backend field"))
    print("num_pockets:", len(result.get("ranked_pockets", [])))