import argparse, hashlib, json, subprocess, sys, math
from pathlib import Path
from datetime import datetime

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser(); p.add_argument('pdb',type=Path); p.add_argument('baseline',type=Path); p.add_argument('--ca-radius',type=float,required=True); a=p.parse_args()
    m=json.loads(a.baseline.read_text(encoding='utf-8-sig'))['provenance']['manifest']
    if sha(a.pdb)!=m['input_sha256'] or m['alpha_squared']!=0: p.error('Baseline PDB hash or alpha mismatch')
    if not math.isfinite(a.ca_radius) or a.ca_radius<=0: p.error('Invalid radius')
    lines=[l for l in a.pdb.read_text().splitlines() if l.startswith('HETATM') and l[17:20].strip()=='CA' and l[76:78].strip()=='CA']
    if len(lines)!=1: p.error('Expected exactly one calcium record')
    xyz=[float(lines[0][i:i+8]) for i in (30,38,46)]
    cli=Path(__file__).with_name('cast_alpha.py');
    if not cli.is_file(): p.error('Place next to cast_alpha.py')
    import pocketlab.core as core
    if sha(Path(core.__file__))!=m['backend_sha256']: p.error('Backend binary differs from baseline')
    m.update(coordinates=m['coordinates']+[xyz],atomic_radii=m['atomic_radii']+[a.ca_radius],n_atoms=len(m['coordinates'])+1,backend=core.__file__,baseline_report_sha256=sha(a.baseline),calcium_record=lines[0],calcium_radius_A=a.ca_radius,model_note='Calcium sensitivity trial; THP and waters excluded; no probe inflation')
    g=core.build_cast_complex_exact(m['coordinates'],m['atomic_radii'],0.0)
    if not g.get('ok') or not g.get('alpha_all_dimensions_exported'): p.error('Exact backend export failed; no fallback')
    out=Path('output_isolated_tests')/('cast_1snc_calcium_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f')); out.mkdir(parents=True)
    gp=out/'cast_complex_all_dims.json'; mp=out/'manifest.json'
    for f,obj in ((gp,g),(mp,m)): f.write_text(json.dumps(obj,allow_nan=False),encoding='utf-8')
    print('CGAL EXPORT:',out.resolve(),flush=True)
    return subprocess.call([sys.executable,str(cli),'--complex',str(gp),'--manifest',str(mp),'--samples','0','--out',str(out/'alpha_results')])
if __name__=='__main__': raise SystemExit(main())
