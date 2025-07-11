import os
import json

# Path alla directory dove ci sono i tuoi file json
json_dir = "C:/my_test/pycas/results/"  # <-- CAMBIA con il tuo percorso

total_files = 0
top1_found = 0
top3_found = 0
topall_found = 0

print("Processing files...", end="", flush=True)

for i, (root, dirs, files) in enumerate(os.walk(json_dir)):
    for file in files:
        if file.endswith(".json"):
            total_files += 1
            file_path = os.path.join(root, file)
            with open(file_path, "r") as f:
                data = json.load(f)
            mesh = data.get("ligand_containment_mesh", [])
            if mesh:
                if mesh[0] is True:
                    top1_found += 1
                if any(mesh[:3]):
                    top3_found += 1
                if any(mesh):
                    topall_found += 1

            # Progresso ogni 100 file
            if total_files % 100 == 0:
                print(".", end="", flush=True)

print(" done!\n")

# Calcola percentuali
if total_files > 0:
    pct_top1 = 100 * top1_found / total_files
    pct_top3 = 100 * top3_found / total_files
    pct_topall = 100 * topall_found / total_files

    print(f"Total files: {total_files}")
    print(f"Top1 found: {top1_found} ({pct_top1:.1f}%)")
    print(f"Top3 found: {top3_found} ({pct_top3:.1f}%)")
    print(f"Top ALL found: {topall_found} ({pct_topall:.1f}%)")
else:
    print("No json files found.")
