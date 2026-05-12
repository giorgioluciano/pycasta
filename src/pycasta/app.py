from flask import Flask, jsonify, send_file, abort, request
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

# ── Lab / Step-1 routes ──────────────────────────────────────────
current_lab_dir = DATA_DIR

@app.route("/lab/step1")
def lab_step1():
    return send_file(os.path.join(BASE_DIR, "alpha_lab_step1.html"))

@app.route("/api/set_data_dir", methods=["POST"])
def set_data_dir():
    global current_lab_dir
    body = request.get_json(force=True)
    d = body.get("dir", "").strip()
    if not d or not os.path.isdir(d):
        return jsonify({"ok": False, "error": "Directory not found: " + d})
    current_lab_dir = d
    return jsonify({"ok": True})

@app.route("/api/lab/pdb_list")
def lab_pdb_list():
    try:
        files = [os.path.splitext(f)[0] for f in os.listdir(current_lab_dir)
                 if f.lower().endswith(".pdb")]
        return jsonify(sorted(files))
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route("/api/lab/compute_alpha", methods=["POST"])
def lab_compute_alpha():
    from scipy.spatial import Delaunay
    body    = request.get_json(force=True)
    mol     = body.get("mol", "")
    alpha   = float(body.get("alpha", 1.4))
    maxtet  = int(body.get("maxtet", 50000))
    mode    = body.get("mode", "standard")
    pdb_path = os.path.join(current_lab_dir, mol + ".pdb")
    if not os.path.exists(pdb_path):
        return jsonify({"error": "PDB not found"}), 404
    try:
        # Parse CA atoms
        coords, radii_w = [], []
        with open(pdb_path) as f:
            for line in f:
                if line.startswith(("ATOM","HETATM")):
                    try:
                        x=float(line[30:38]); y=float(line[38:46]); z=float(line[46:54])
                        coords.append([x,y,z])
                        radii_w.append(1.8)
                    except: pass
        if len(coords) < 4:
            return jsonify({"error": "Not enough atoms"}), 400
        pts = np.array(coords)
        n_atoms = len(pts)
        if mode == "weighted":
            rw = np.array(radii_w)
            pts_w = np.hstack([pts, (rw**2 - alpha**2).reshape(-1,1)])
            tri = Delaunay(pts_w[:, :3])
        else:
            tri = Delaunay(pts)
        tets = tri.simplices
        if len(tets) > maxtet:
            tets = tets[:maxtet]
        # Circumradii
        def circumradius(p):
            a,b,c,d = p
            A=b-a; B=c-a; C=d-a
            M=np.array([[np.dot(A,A),np.dot(A,B),np.dot(A,C)],
                        [np.dot(B,A),np.dot(B,B),np.dot(B,C)],
                        [np.dot(C,A),np.dot(C,B),np.dot(C,C)]])
            rhs=0.5*np.array([np.dot(A,A),np.dot(B,B),np.dot(C,C)])
            try:
                u=np.linalg.solve(M,rhs)
                return float(np.linalg.norm(u))
            except: return float('inf')
        radii_arr = np.array([circumradius(pts[t]) for t in tets])
        finite_mask = np.isfinite(radii_arr)
        tets = tets[finite_mask]
        radii_arr = radii_arr[finite_mask]
        alpha_mask = radii_arr <= alpha
        # Build flattened line segments
        EDGES = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
        lines_in, lines_out = [], []
        for i, (tet, inside) in enumerate(zip(tets, alpha_mask)):
            lst = lines_in if inside else lines_out
            for e0,e1 in EDGES:
                p0=pts[tet[e0]]; p1=pts[tet[e1]]
                lst.extend([float(p0[0]),float(p0[1]),float(p0[2]),
                             float(p1[0]),float(p1[1]),float(p1[2])])
        # Read PDB string
        with open(pdb_path) as f:
            pdb_str = f.read()
        return jsonify({
            "mol": mol, "mode": mode, "alpha": alpha,
            "n_total": int(len(tets)),
            "n_in":    int(alpha_mask.sum()),
            "n_out":   int((~alpha_mask).sum()),
            "n_atoms": n_atoms,
            "radii":   [round(float(r),4) for r in radii_arr[:5000]],
            "lines_in":  lines_in,
            "lines_out": lines_out,
            "pdb_str": pdb_str
        })
    except Exception as e:
        import traceback; traceback.print_exc()
        return jsonify({"error": str(e)}), 500

    app.run(debug=True, port=5000)
