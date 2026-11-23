import argparse, hashlib, json, math, subprocess, sys
from pathlib import Path
from datetime import datetime

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def main():
    a=argparse.ArgumentParser()
    a.add_argument('manifest',type=Path);a.add_argument('--radii-file',type=Path,required=True)
    a.add_argument('--out',type=Path,default=Path('output_isolated_tests'))
    a.add_argument('--samples',type=int,default=0); args=a.parse_args()
    m=read(args.manifest); b=read(args.radii_file)
    xyz=m['coordinates']; old=m['atomic_radii']; new=b['atomic_radii']; alpha=m['alpha_squared']
    if alpha != 0 or len(xyz)!=len(old) or len(new)!=len(old): a.error('Expected alpha_squared=0 and equal atom counts')
    if not m.get('input_sha256') or b.get('input_sha256')!=m['input_sha256'] or b.get('coordinates')!=xyz:
        a.error('PDB hash or ordered coordinates differ')
    if not isinstance(b.get('source'),str) or not b['source'].strip(): a.error('Document the radius source')
    if not all(type(r) in (float,int) and math.isfinite(r) and r>0 for r in new): a.error('Invalid radii')
    changed=sum(x!=y for x,y in zip(old,new))
    if not changed: a.error('No radius changed')
    if args.samples<0 or args.samples==1: a.error('Use samples=0 or >=2')
    cli=Path(__file__).with_name('cast_alpha.py')
    if not cli.is_file(): a.error('Place this script next to cast_alpha.py')
    import pocketlab.core as core
    g=core.build_cast_complex_exact(xyz,new,alpha)
    if not g.get('ok') or not g.get('alpha_all_dimensions_exported'): a.error('Exact backend export failed; no fallback')
    d=args.out/('cast_radii_B_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'));d.mkdir(parents=True,exist_ok=False)
    m.update(atomic_radii=new, backend=core.__file__,backend_sha256=sha(Path(core.__file__)),
             model_note='Explicit alternative radii; no additional probe inflation',radius_source=b['source'],
             baseline_manifest_sha256=sha(args.manifest), radius_file_sha256=sha(args.radii_file), changed_radii=changed)
    gp=d/'cast_complex_all_dims.json';mp=d/'manifest.json'
    for p,data in [(gp,g),(mp,m)]: p.write_text(json.dumps(data,indent=2,allow_nan=False),encoding='utf-8')
    print(f'Changed radii: {changed}/{len(old)}; export: {d.resolve()}',flush=True)
    return subprocess.call([sys.executable,str(cli),'--complex',str(gp),'--manifest',str(mp),'--samples',str(args.samples),'--out',str(d/'alpha_results')])
if __name__=='__main__': raise SystemExit(main())
