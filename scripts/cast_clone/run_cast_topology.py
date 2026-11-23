import argparse
from pathlib import Path
from datetime import datetime
import json
from collections import Counter
from pocketlab.core import build_cast_complex_exact
from cast_topology import build_topology

p=argparse.ArgumentParser()
p.add_argument('manifest',type=Path)
a=p.parse_args()
m=json.loads(a.manifest.read_text(encoding='utf-8'))
g=build_cast_complex_exact(m['coordinates'],m['atomic_radii'],m['alpha_squared'])
if not g['ok']:raise SystemExit('Geometry degenerate; stop: '+g['status'])
t=build_topology(g)
d=a.manifest.resolve().parent/('topology_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
d.mkdir(exist_ok=False)
(d/'cast_complex_all_dims.json').write_text(json.dumps(g),encoding='utf-8')
(d/'topology.json').write_text(json.dumps(t),encoding='utf-8')
(d/'manifest.json').write_text(json.dumps(m,indent=2),encoding='utf-8')
summary={'status':t['status'],'n_cells':len(g['cells_idx']),
         'n_regions':len(t['regions']),
         'classes':dict(Counter(r['classification'] for r in t['regions'])),
         'regions':[{'id':r['region_id'],'tetrahedra':len(r['tetra_indices']),
                     'mouths':r['n_mouths'],
                     'mouth_triangle_counts':[len(x['triangles']) for x in r['mouths']],
                     'boundary_edge_manifold_check':r['boundary_edge_manifold_check']}
                    for r in t['regions']]}
(d/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2));print('OUTPUT:',d)
