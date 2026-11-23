#!/usr/bin/env python3
"""Viewer HTML per le regioni CAST di un run PyCASTa (solo libreria standard).

Uso:
  python3 cast_region_viewer.py --pdb 1fkf.pdb --out OUTDIR RUN_DIR [RUN_DIR ...]
RUN_DIR contiene manifest.json e alpha_results/cast_alpha_*/results.json.
Il browser carica 3Dmol.js da internet (CDN).
"""
import argparse, json, math
from pathlib import Path

COL = {"ligand": "gold", "second": "dodgerblue", "cavity": "magenta",
       "pocket": "limegreen", "other": "orange"}

TEMPLATE = r"""<!doctype html><html><head><meta charset="utf-8"><title>__TITLE__</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
body{margin:0;font-family:sans-serif;display:flex;height:100vh}
#side{width:430px;overflow:auto;padding:10px;font-size:13px;border-right:1px solid #ccc}
#v{flex:1;position:relative}
.row{display:flex;align-items:center;gap:6px;margin:3px 0}
.sw{width:14px;height:14px;border:1px solid #444;display:inline-block}
button{margin-right:4px}
</style></head><body>
<div id="side"><h3>__TITLE__</h3>
<div>Rosso = ligando FK5 (non usato nel calcolo). Giallo = regione con piu' contatti col ligando.
Blu = seconda tasca con bocche. Magenta = cavita'. Verde = altre tasche. Arancio = altre.</div>
<p><button onclick="setAll(true)">tutte</button><button onclick="setAll(false)">nessuna</button>
<label><input type="checkbox" id="mesh" checked onchange="render()"> tetraedri</label>
<label><input type="checkbox" id="atoms" checked onchange="render()"> atomi</label>
<label><input type="checkbox" id="labels" onchange="render()"> residui</label></p>
<div id="list"></div></div><div id="v"></div>
<script>
const DATA=__DATA__; const PDB=__PDB__;
const viewer=$3Dmol.createViewer('v',{backgroundColor:'white'});
viewer.addModel(PDB,'pdb');
viewer.setStyle({},{cartoon:{color:'lightgrey',opacity:0.6}});
viewer.setStyle({resn:'FK5'},{stick:{color:'red',radius:0.25}});
const list=document.getElementById('list');
DATA.regions.forEach(r=>{
  const d=document.createElement('div'); d.className='row';
  d.innerHTML='<input type="checkbox" id="c'+r.id+'" '+(r.show?'checked':'')+' onchange="render()">'
   +'<span class="sw" style="background:'+r.color+'"></span>'
   +'<span>#'+r.id+' '+r.cls+' | celle '+r.ncells+' | bocche '+r.nm+' | atomi '+r.atoms.length
   +' | contatti&lt;=4A '+r.near4+' | dmin '+r.dmin.toFixed(2)+'</span>';
  list.appendChild(d);});
function setAll(v){DATA.regions.forEach(r=>document.getElementById('c'+r.id).checked=v);render();}
function render(){
  viewer.removeAllShapes(); viewer.removeAllLabels();
  const showMesh=document.getElementById('mesh').checked, showAt=document.getElementById('atoms').checked,
        showLab=document.getElementById('labels').checked;
  DATA.regions.forEach(r=>{
    if(!document.getElementById('c'+r.id).checked) return;
    if(showMesh){
      const V=[],N=[],F=[];
      r.tets.forEach(t=>{
        const P=t.map(i=>DATA.xyz[i]);
        [[0,1,2],[0,1,3],[0,2,3],[1,2,3]].forEach(f=>{
          const a=P[f[0]],b=P[f[1]],c=P[f[2]];
          const u=[b[0]-a[0],b[1]-a[1],b[2]-a[2]], w=[c[0]-a[0],c[1]-a[1],c[2]-a[2]];
          let n=[u[1]*w[2]-u[2]*w[1],u[2]*w[0]-u[0]*w[2],u[0]*w[1]-u[1]*w[0]];
          const l=Math.hypot(n[0],n[1],n[2])||1; n=n.map(x=>x/l);
          const s=V.length;
          [a,b,c].forEach(p=>{V.push({x:p[0],y:p[1],z:p[2]});N.push({x:n[0],y:n[1],z:n[2]});});
          F.push(s,s+1,s+2,s,s+2,s+1);
        });
      });
      viewer.addCustom({vertexArr:V,normalArr:N,faceArr:F,color:r.color,alpha:0.45});
    }
    if(showAt) r.atoms.forEach(i=>{const p=DATA.xyz[i];
      viewer.addSphere({center:{x:p[0],y:p[1],z:p[2]},radius:0.45,color:r.color,alpha:0.9});});
    if(showLab){const seen={}; r.atoms.forEach(i=>{const k=DATA.res[i]; if(seen[k])return; seen[k]=1;
      const p=DATA.xyz[i]; viewer.addLabel(k,{position:{x:p[0],y:p[1],z:p[2]},fontSize:11,
        backgroundOpacity:0.5,fontColor:'black'});});}
  });
  viewer.render();
}
render(); viewer.zoomTo({resn:'FK5'}); viewer.zoom(0.7); viewer.render();
</script></body></html>
"""

def read_pdb(path):
    atoms, lig, text = [], [], Path(path).read_text()
    for line in text.splitlines():
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        xyz = [float(line[i:i + 8]) for i in (30, 38, 46)]
        if rec == "ATOM":
            atoms.append((xyz, line[17:20].strip() + line[22:26].strip()))
        elif line[17:20].strip() == "FK5":
            lig.append(xyz)
    return atoms, lig, text

def build(run, pdb, out_dir):
    run = Path(run)
    m = json.loads((run / "manifest.json").read_text())
    res = json.loads(next((run / "alpha_results").glob("cast_alpha_*/results.json")).read_text())
    atoms, lig, text = read_pdb(pdb)
    xyz = m["coordinates"]
    if len(atoms) != len(xyz):
        raise SystemExit(f"{run.name}: {len(atoms)} atomi PDB vs {len(xyz)} nel manifest")
    for i, (a, b) in enumerate(zip(atoms, xyz)):
        if max(abs(p - q) for p, q in zip(a[0], b)) > 1e-3:
            raise SystemExit(f"{run.name}: coordinate diverse all'atomo {i}")
    resl = [a[1] for a in atoms]
    dmin_atom = [min(math.dist(p, l) for l in lig) for p in xyz]

    regs = []
    for rec in res["regions"]:
        tets = [list(t) for t in rec["dual_set"]["tetrahedra"]]
        verts = sorted({a for t in tets for a in t})
        regs.append({"id": rec["region_id"], "cls": rec["classification"],
                     "nm": rec["n_mouths"], "ncells": len(tets), "tets": tets,
                     "atoms": verts, "near4": sum(dmin_atom[i] <= 4.0 for i in verts),
                     "dmin": min(dmin_atom[i] for i in verts)})
    site = max(regs, key=lambda r: (r["near4"], r["ncells"]))
    cand = [r for r in regs if r is not site and r["nm"] > 0]
    second = max(cand, key=lambda r: r["ncells"]) if cand else None
    for r in regs:
        if r is site: r["color"], r["show"] = COL["ligand"], True
        elif r is second: r["color"], r["show"] = COL["second"], True
        elif r["cls"] == "cavity": r["color"], r["show"] = COL["cavity"], False
        elif r["nm"] > 0: r["color"], r["show"] = COL["pocket"], False
        else: r["color"], r["show"] = COL["other"], False
    regs.sort(key=lambda r: -r["ncells"])
    html = (TEMPLATE.replace("__TITLE__", f"FKBP12 1FKF - {run.name}")
            .replace("__DATA__", json.dumps({"xyz": xyz, "res": resl, "regions": regs}))
            .replace("__PDB__", json.dumps(text)))
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    f = out / f"{run.name}.html"
    f.write_text(html, encoding="utf-8")
    print(f"{f}  (regioni: {len(regs)}, sito: #{site['id']}, seconda: "
          f"{'#' + str(second['id']) if second else 'nessuna'})")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("runs", nargs="+")
    a = ap.parse_args()
    for r in a.runs:
        build(r, a.pdb, a.out)
