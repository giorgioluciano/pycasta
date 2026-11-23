"""CAST research alpha: consume exact CGAL exports, never silently retriangulate."""
import argparse,csv,hashlib,html,json,platform,sys,time
from collections import Counter,defaultdict,deque
from datetime import datetime
from itertools import combinations
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
VERSION='0.1.0'
NAMES=('vertices','edges','triangles','tetrahedra')
def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def dump(path,data):
    Path(path).write_text(json.dumps(data,indent=2,allow_nan=False),encoding='utf-8')
def components(nodes,adj):
    left=set(nodes);out=[]
    while left:
        seed=min(left);left.remove(seed);stack=[seed];part=[]
        while stack:
            i=stack.pop();part.append(i)
            for j in adj(i):
                if j in left:left.remove(j);stack.append(j)
        out.append(sorted(part))
    return sorted(out,key=lambda p:(-len(p),p[0]))
def serialize(sets):
    return {name:[list(x) for x in sorted(s)] for name,s in zip(NAMES,sets)}
def simplex_components(sets):
    nodes=set().union(*sets);adj=defaultdict(set)
    for s in nodes:
        for size in range(1,len(s)):
            for sub in combinations(s,size):
                if sub in nodes:adj[s].add(sub);adj[sub].add(s)
    return [serialize([{s for s in part if len(s)==d+1} for d in range(4)])
            for part in components(nodes,lambda i:adj[i])]
def alpha_sets(g):
    cells=[tuple(sorted(x)) for x in g['cells_idx']]
    K=[{(a,) for a,v in g['alpha_vertex_rows'] if v},
       {tuple(sorted((a,b))) for a,b,v in g['alpha_edge_rows'] if v},set(),
       {cells[i] for i,v in enumerate(g['tetra_in_alpha']) if v}]
    for i,cell in enumerate(g['cells_idx']):
        for k,v in enumerate(g['face_in_alpha'][i]):
            if v:K[2].add(tuple(sorted(cell[:k]+cell[k+1:])))
    for d in range(1,4):
        for s in K[d]:
            if any(f not in K[d-1] for f in combinations(s,d)):
                raise ValueError('Alpha complex not closed under faces')
    return cells,K

def reconstruct(g):
    if not g.get('ok') or not g.get('alpha_all_dimensions_exported'):
        raise ValueError('Need successful all-dimensional CGAL export')
    if g.get('index_space')!='input_atoms':raise ValueError('Input-atom indexing required')
    if g.get('zero_face_predicates') or g.get('equal_power_flow_faces'):
        raise ValueError('Degenerate flow predicates: symbolic handling not implemented')
    cells,K=alpha_sets(g);n=len(cells);ext=g['exterior']
    if ext!=n or len(g['successors'])!=n+1:raise ValueError('Invalid exterior/successor schema')
    owners=defaultdict(list)
    for i,cell in enumerate(cells):
        if len(set(cell))!=4:raise ValueError('Repeated vertex in tetrahedron')
        for f in combinations(cell,3):owners[f].append(i)
    if len(set(cells))!=n:raise ValueError('Duplicate tetrahedra')
    adj=[set() for _ in cells]
    for f,ids in owners.items():
        if len(ids)>2:raise ValueError('Face has more than two tetrahedra')
        if len(ids)==2:
            a,b=ids;adj[a].add(b);adj[b].add(a)
    for i,cell in enumerate(g['cells_idx']):
        if len(g['neighbors'][i])!=4:raise ValueError('Invalid neighbors')
        for k,j in enumerate(g['neighbors'][i]):
            f=tuple(sorted(cell[:k]+cell[k+1:]));expected=set(owners[f])-{i}
            if expected:
                if j not in expected:raise ValueError('Opposite-face neighbor mismatch')
            elif j not in (-1,ext):raise ValueError('Invalid hull neighbor')
    reverse=[[] for _ in range(n+1)]
    for i,targets in enumerate(g['successors']):
        for j in targets:
            if not 0<=j<=n:raise ValueError('Invalid flow target')
            if i<n and j<n and j not in adj[i]:raise ValueError('Flow between nonneighbors')
            reverse[j].append(i)
    reached={ext};queue=deque([ext])
    while queue:
        for i in reverse[queue.popleft()]:
            if i not in reached:reached.add(i);queue.append(i)
    outside=reached-{ext}
    if outside!=set(g['outside_ancestors']):raise ValueError('Outside reachability mismatch')
    eligible={i for i in range(n) if i not in outside and not g['tetra_in_alpha'][i]}
    parts=components(eligible,lambda i:(j for j in adj[i] if tuple(sorted(set(cells[i]) & set(cells[j]))) not in K[2]))
    if 'pocket_tetra_components' in g:
        if {frozenset(x) for x in parts}!={frozenset(x) for x in g['pocket_tetra_components']}:
            raise ValueError('Component rule differs from backend export')
    regions=[]
    for rid,part in enumerate(parts):
        C=[set() for _ in range(4)];counts=Counter()
        for i in part:
            for d in range(4):C[d].update(combinations(cells[i],d+1))
            counts.update(combinations(cells[i],3))
        B=[set() for _ in range(4)];B[2]={f for f,nf in counts.items() if nf==1}
        for f in B[2]:
            B[0].update(combinations(f,1));B[1].update(combinations(f,2))
        P=[C[d]-B[d]-K[d] for d in range(4)]
        mouth=[B[d]-K[d] for d in range(4)]
        mouths=simplex_components(mouth)
        degree=Counter(e for f in B[2] for e in combinations(f,2))
        bad=[{'edge':list(e),'degree':v,'in_alpha':e in K[1]}
             for e,v in sorted(degree.items()) if v!=2]
        regions.append(dict(region_id=rid,tetra_indices=part,dual_set=serialize(P),
                            closure_boundary=serialize(B),mouths=mouths,n_mouths=len(mouths),
                            classification='cavity' if not mouths else 'pocket' if len(mouths)==1 else 'multiple_connector',
                            raw_boundary_edge_contacts=bad,
                            boundary_edge_manifold_check=not bad))
    return dict(status='research_alpha',regions=regions,
                physical_analytic_measures_computed=False,
                original_CAST_equivalence_certified=False,
                flow_predicate_implementation_reverified=False,
                checks={'alpha_closure':True,'neighbors':True,'outside_reachability':True,'backend_partition':True},
                limitations=['Consumes existing exact backend output; does not regenerate or prove flow predicates.',
                             'No symbolic perturbation; detected ambiguous predicates rejected.',
                             'Raw manifold warnings retained; vertex-link checks not implemented here.',
                             'No analytic molecular surface areas, mouth areas/perimeters or atom attribution.',
                             'Numerical volume is optional and is not CAST analytic volume.',
                             'Inherited radii/preprocessing not certified equivalent to 1998 CAST.'])

def free_fraction(points,centers,radii,tree):
    hits=tree.query_ball_point(points,float(radii.max()))
    lengths=np.fromiter((len(x) for x in hits),dtype=int,count=len(points))
    if not lengths.sum():return np.ones(len(points),dtype=bool)
    atom=np.concatenate([np.asarray(x,dtype=int) for x in hits if len(x)])
    sample=np.repeat(np.arange(len(points)),lengths)
    inside=np.sum((points[sample]-centers[atom])**2,axis=1)<=radii[atom]**2
    covered=np.bincount(sample[inside],minlength=len(points))>0
    return ~covered

def measures(g,m,report,samples,seed):
    xyz=np.asarray(m['coordinates'],dtype=float);radii=np.asarray(m['atomic_radii'],dtype=float)
    if xyz.shape!=(len(radii),3) or not np.isfinite(xyz).all() or not np.isfinite(radii).all() or (radii<=0).any():
        raise ValueError('Invalid coordinates/radii')
    cells=np.asarray(g['cells_idx'],dtype=int)
    if cells.min()<0 or cells.max()>=len(xyz):raise ValueError('Atom index outside manifest')
    if m.get('alpha_squared')!=g.get('alpha_squared'):raise ValueError('Manifest/complex alpha mismatch')
    corners=xyz[cells];vol=np.abs(np.linalg.det(corners[:,1:]-corners[:,:1]))/6
    if (vol<=0).any():raise ValueError('Zero geometric tetra volume')
    # This is a diagnostic proxy, deliberately not the physical molecular volume.
    for r in report['regions']:
        ids=np.asarray(r['tetra_indices']);v=vol[ids]
        r['tetrahedral_support_volume_A3']=float(v.sum())
        centroid=corners[ids].mean(axis=1)
        r['support_centroid_A']=np.average(centroid,axis=0,weights=v).tolist()
        r['free_volume_mc_A3']=None;r['free_volume_mc_standard_error_A3']=None
    if not samples:return xyz
    if samples<2:raise ValueError('At least two MC samples per tetrahedron')
    rng=np.random.default_rng(seed);tree=cKDTree(xyz)
    for r in report['regions']:
        total=0.;variance=0.
        ids=r['tetra_indices']
        for start in range(0,len(ids),64):
            ix=np.asarray(ids[start:start+64]);nt=len(ix)
            w=rng.exponential(size=(nt,samples,4));w/=w.sum(axis=2,keepdims=True)
            pts=np.einsum('bsi,bij->bsj',w,corners[ix]).reshape(-1,3)
            free=free_fraction(pts,xyz,radii,tree).reshape(nt,samples)
            p=free.mean(axis=1);v=vol[ix]
            total+=float(np.sum(v*p))
            variance+=float(np.sum(v*v*p*(1-p)/(samples-1)))
        r['free_volume_mc_A3']=total;r['free_volume_mc_standard_error_A3']=variance**.5
        r['mc_samples_per_tetrahedron']=samples
    return xyz

def export_obj(path,faces,xyz):
    vertices=sorted({i for f in faces for i in f});mapping={v:k+1 for k,v in enumerate(vertices)}
    with Path(path).open('w',encoding='ascii') as f:
        f.write('# Diagnostic triangle support; not a molecular surface; winding not guaranteed.\n')
        for i in vertices:f.write('v '+' '.join(format(v,'.9g') for v in xyz[i])+'\n')
        for face in faces:f.write('f '+' '.join(str(mapping[i]) for i in face)+'\n')

def compare(report,old):
    bypart={frozenset(r['tetra_indices']):r for r in old['regions']};errors=[]
    if len(bypart)!=len(report['regions']):errors.append('Region count mismatch')
    norm=lambda x: {tuple(sorted(y)) for y in x}
    mouthnorm=lambda m: tuple(frozenset(norm(m.get(name,[]))) for name in NAMES)
    for r in report['regions']:
        prev=bypart.get(frozenset(r['tetra_indices']))
        if prev is None:errors.append('Missing corresponding region');continue
        for section in ('dual_set','closure_boundary'):
            for name in NAMES:
                if norm(r[section][name])!=norm(prev[section][name]):errors.append(f'{section}/{name} differs')
        if {mouthnorm(m) for m in r['mouths']}!={mouthnorm(m) for m in prev['mouths']}:errors.append('Mouths differ')
    return {'matches_reference_topology':not errors,'differences':errors}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--complex',required=True,type=Path);p.add_argument('--manifest',required=True,type=Path)
    p.add_argument('--reference-topology',type=Path);p.add_argument('--out',type=Path,default=Path('output_isolated_tests'))
    p.add_argument('--samples',type=int,default=0,help='MC samples per tetra; 0 disables numerical integration')
    p.add_argument('--seed',type=int,default=1729);a=p.parse_args()
    if a.samples<0 or a.samples==1:p.error('--samples must be 0 or >=2')
    started=time.perf_counter();g=load(a.complex);m=load(a.manifest);report=reconstruct(g)
    xyz=measures(g,m,report,a.samples,a.seed)
    if a.reference_topology:
        report['reference_comparison']=compare(report,load(a.reference_topology))
        if not report['reference_comparison']['matches_reference_topology']:
            raise ValueError(str(report['reference_comparison']))
    dest=a.out/('cast_alpha_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'));dest.mkdir(parents=True,exist_ok=False)
    report['version']=VERSION;report['sampling']={'seed':a.seed,'samples_per_tetrahedron':a.samples,
        'method':'independent uniform volume sampling via Dirichlet barycentric weights',
        'uncertainty_note':'MC standard error only; zero observed variance is not proof of exact volume.'}
    report['provenance']={'python':sys.version,'platform':platform.platform(),
        'input_sha256':{str(x.resolve()):hashlib.sha256(x.read_bytes()).hexdigest()
                        for x in [a.complex,a.manifest]+([a.reference_topology] if a.reference_topology else [])},
        'manifest':m,'kernel':g.get('kernel'),'alpha_squared':g.get('alpha_squared')}
    ranked=sorted(report['regions'],key=lambda r:-(r['free_volume_mc_A3'] if a.samples else r['tetrahedral_support_volume_A3']))
    rows=[]
    for rank,r in enumerate(ranked,1):
        r['provisional_rank']=rank
        rows.append({'rank':rank,'region_id':r['region_id'],'classification':r['classification'],
                     'tetrahedra':len(r['tetra_indices']),'mouths':r['n_mouths'],
                     'support_volume_A3':r['tetrahedral_support_volume_A3'],
                     'free_volume_mc_A3':r['free_volume_mc_A3'],
                     'mc_standard_error_A3':r['free_volume_mc_standard_error_A3'],
                     'raw_bad_edges':len(r['raw_boundary_edge_contacts'])})
        folder=dest/f"region_{r['region_id']:03d}";folder.mkdir()
        export_obj(folder/'boundary_support.obj',r['closure_boundary']['triangles'],xyz)
        for j,mouth in enumerate(r['mouths']):export_obj(folder/f'mouth_{j:03d}_support.obj',mouth['triangles'],xyz)
    with (dest/'summary.csv').open('w',newline='',encoding='utf-8') as f:
        fields=list(rows[0]) if rows else ['region_id'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    report['elapsed_seconds']=time.perf_counter()-started;dump(dest/'results.json',report)
    columns=list(rows[0]) if rows else []
    table='<table border="1"><tr>'+''.join('<th>'+html.escape(k)+'</th>' for k in columns)+'</tr>'
    for row in rows:table+='<tr>'+''.join('<td>'+html.escape(str(row[k]))+'</td>' for k in columns)+'</tr>'
    (dest/'report.html').write_text('<!doctype html><meta charset="utf-8"><title>CAST research alpha</title><h1>CAST research alpha 0.1</h1><p>NOT certified original CAST. OBJ files are flat diagnostic supports, not molecular surfaces.</p><ul>'+''.join('<li>'+html.escape(x)+'</li>' for x in report['limitations'])+'</ul>'+table+'</table>',encoding='utf-8')
    print(f"Regions: {len(rows)}; mouths: {sum(r['n_mouths'] for r in report['regions'])}")
    print('Reference:',report.get('reference_comparison','not supplied'))
    print('Ranking:', 'numerical free-volume estimate' if a.samples else 'tetrahedral support proxy (NOT physical volume)')
    print('OUTPUT:',dest.resolve())
if __name__=='__main__':main()
