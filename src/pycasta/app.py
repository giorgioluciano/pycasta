from flask import Flask, jsonify, send_file, abort
import os, json, glob
import numpy as np

app = Flask(__name__)

BASE_DIR    = r"C:\github\pycasta\src\pycasta"
RESULTS_DIR = os.path.join(BASE_DIR, "results")
DATA_DIR    = r"C:\github\pycasta\data\coach420\pdbs"

def load_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)

@app.route("/")
def index():
    return send_file(os.path.join(BASE_DIR, "viewer.html"))

@app.route("/api/datasets")
def datasets():
    ds = [d for d in os.listdir(RESULTS_DIR)
          if os.path.isdir(os.path.join(RESULTS_DIR, d))]
    return jsonify(sorted(ds))

@app.route("/api/molecules/<dataset>")
def molecules(dataset):
    folder = os.path.join(RESULTS_DIR, dataset)
    if not os.path.isdir(folder):
        return jsonify([])
    mols = [f.replace("_pockets.json","")
            for f in os.listdir(folder)
            if f.endswith("_pockets.json")]
    return jsonify(sorted(mols))

@app.route("/api/pockets/<dataset>/<molecule>")
def pockets(dataset, molecule):
    path = os.path.join(RESULTS_DIR, dataset, f"{molecule}_pockets.json")
    if not os.path.exists(path):
        abort(404)
    return jsonify(load_json(path))

@app.route("/api/alpha/<dataset>/<molecule>")
def alpha(dataset, molecule):
    # Cerca il file .alpha.npz
    path = os.path.join(RESULTS_DIR, dataset, f"{molecule}.alpha.npz")
    if not os.path.exists(path):
        abort(404)
    try:
        data = np.load(path, allow_pickle=True)
        # Converti tetrahedra e alpha_mask in liste JSON-serializzabili
        result = {
            "tetrahedra":  data["tetrahedra"].tolist(),
            "alpha_mask":  data["alpha_mask"].tolist(),
        }
        # Aggiungi alpha_value se presente
        if "alpha_value" in data:
            result["alpha_value"] = float(data["alpha_value"])
        return jsonify(result)
    except Exception as e:
        print(f"Error loading alpha: {e}")
        abort(500)
       

@app.route("/api/pdb/<molecule>")
def pdb(molecule):
    path = os.path.join(DATA_DIR, f"{molecule}.pdb")
    if not os.path.exists(path):
        abort(404)
    return send_file(path, mimetype="text/plain")

@app.route("/api/pocket_pdb/<dataset>/<molecule>/<int:rank>")
def pocket_pdb(dataset, molecule, rank):
    path = os.path.join(RESULTS_DIR, dataset, molecule, f"pocket_{rank}.pdb")
    if not os.path.exists(path):
        abort(404)
    return send_file(path, mimetype="text/plain")

@app.route("/api/mesh/<dataset>/<molecule>/<int:rank>")
def mesh(dataset, molecule, rank):
    path = os.path.join(RESULTS_DIR, dataset, molecule, f"pocket_{rank}.ply")
    if not os.path.exists(path):
        abort(404)
    return send_file(path, mimetype="application/octet-stream")

@app.route("/api/export/<dataset>/<molecule>")
def export(dataset, molecule):
    path = os.path.join(RESULTS_DIR, dataset, f"{molecule}_pockets.json")
    if not os.path.exists(path):
        abort(404)
    return send_file(path, as_attachment=True,
                     download_name=f"{molecule}_pockets.json")

if __name__ == "__main__":
    app.run(debug=True, port=5000)