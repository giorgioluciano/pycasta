import os
import urllib.request

base = "https://raw.githubusercontent.com/rdk/p2rank-datasets/master/coach420/"
list_url = "https://raw.githubusercontent.com/rdk/p2rank-datasets/master/coach420.ds"

os.makedirs("data/coach420/pdbs", exist_ok=True)

# Scarica la lista dei file
print("Scarico lista...")
urllib.request.urlretrieve(list_url, "coach420.ds")

with open("coach420.ds") as f:
    lines = [l.strip() for l in f if l.strip() and not l.startswith("#")]

print(f"Trovati {len(lines)} file PDB")
for i, line in enumerate(lines):
    fname = line.split()[0]
    url = f"https://raw.githubusercontent.com/rdk/p2rank-datasets/master/coach420/{fname}"
    dest = f"data/coach420/pdbs/{os.path.basename(fname)}"
    if not os.path.exists(dest):
        try:
            urllib.request.urlretrieve(url, dest)
            print(f"[{i+1}/{len(lines)}] {fname}")
        except Exception as e:
            print(f"ERRORE {fname}: {e}")

print("Download completato!")