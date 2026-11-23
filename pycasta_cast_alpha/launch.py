"""Discover a single frozen run and invoke the consolidated alpha."""
import argparse,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--run-dir',required=True,type=Path)
p.add_argument('--samples',type=int,default=128)
p.add_argument('--seed',type=int,default=1729)
a=p.parse_args()
root=a.run_dir.resolve()
def find(name,required=True):
    hits=list(root.rglob(name))
    if not hits and not required:return None
    if len(hits)!=1:raise SystemExit(f'Expected one {name} under {root}, found {len(hits)}. Use cast_alpha.py with explicit paths.')
    return hits[0]
geo=find('cast_complex_all_dims.json');manifest=find('manifest.json');ref=find('topology.json',False)
cmd=[sys.executable,str(Path(__file__).with_name('cast_alpha.py')),'--complex',str(geo),
     '--manifest',str(manifest),'--samples',str(a.samples),'--seed',str(a.seed),'--out',str(root/'alpha_results')]
if ref:cmd+=['--reference-topology',str(ref)]
raise SystemExit(subprocess.call(cmd))
