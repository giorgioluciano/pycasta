
import sys
sys.path.insert(0, r"C:\github\pocketlab\src")
sys.path.insert(0, r"C:\github\pycasta\src\pycasta")

from run_analysis import process_pdb

result = process_pdb(r"C:\github\pycasta\data\coach420\pdbs\1a26A.pdb")
print(result["step_to_ligand"])
print(result["triangulation_backend"] if "triangulation_backend" in result else "no backend field")