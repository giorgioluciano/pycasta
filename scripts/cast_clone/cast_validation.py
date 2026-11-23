"""Validation increment 4. Reads saved data; never modifies geometry or inputs."""
import argparse
from collections import Counter, defaultdict
from datetime import datetime
from itertools import combinations
import csv
import hashlib
import json
from pathlib import Path

NAMES=('vertices','edges','triangles','tetrahedra')
def unpack(d):
    return {tuple(sorted(s)) for name in NAMES for s in d.get(name,[])}
def closure(S):
    return {f for s in S for k in range(1,len(s)+1) for f in combinations(s,k)}
def component_count(S):
    if not S:return 0
    adj={s:set() for s in S}
    for s in S:
        if len(s)>1:
            for f in combinations(s,len(s)-1):
                if f in S:adj[s].add(f);adj[f].add(s)
    remaining=set(S);count=0
    while remaining:
        count+=1;stack=[remaining.pop()]
        while stack:
            for f in adj[stack.pop()]:
                if f in remaining:remaining.remove(f);stack.append(f)
    return count

def surface_checks(S):
    triangles=[s for s in S if len(s)==3]
    edge_faces=defaultdict(list);vertex_faces=defaultdict(list)
    for f in triangles:
        for e in combinations(f,2):edge_faces[e].append(f)
        for v in f:vertex_faces[v].append(f)
    bad_edges=[{'edge':list(s),'triangle_degree':len(edge_faces[s])}
               for s in sorted(S) if len(s)==2 and len(edge_faces[s])!=2]
    bad_vertices=[]
    for s in sorted(S):
        if len(s)!=1:continue
        v=s[0];adj=defaultdict(set)
        for f in vertex_faces[v]:
            a,b=[x for x in f if x!=v];adj[a].add(b);adj[b].add(a)
        unseen=set(adj);count=0
        while unseen:
            count+=1;stack=[unseen.pop()]
            while stack:
                for w in adj[stack.pop()]:
                    if w in unseen:unseen.remove(w);stack.append(w)
        if count!=1 or any(len(x)!=2 for x in adj.values()):
            bad_vertices.append({'vertex':v,'link_components':count,
                                 'link_degree_histogram':dict(Counter(map(len,adj.values())))})
    return {'included_edge_failures':bad_edges,'included_vertex_link_failures':bad_vertices,
            'local_surface_checks_passed':not bad_edges and not bad_vertices,
            'components':component_count(S)}

def audit(g,t):
    if not g.get('ok') or not g.get('alpha_all_dimensions_exported'):
        raise ValueError('Need nondegenerate geometry with all alpha dimensions')
    cells=[tuple(sorted(x)) for x in g['cells_idx']]
    face_owners=defaultdict(list);face_state={};K=set()
    for i,raw in enumerate(g['cells_idx']):
        if g['tetra_in_alpha'][i]:K.add(cells[i])
        for k in range(4):
            f=tuple(sorted(raw[j] for j in range(4) if j!=k))
            face_owners[f].append(i)
            state=bool(g['face_in_alpha'][i][k])
            if f in face_state and face_state[f]!=state:raise ValueError('Inconsistent face alpha')
            face_state[f]=state
            if state:K.add(f)
    for a,b,inside in g['alpha_edge_rows']:
        if inside:K.add(tuple(sorted((a,b))))
    for a,inside in g['alpha_vertex_rows']:
        if inside:K.add((a,))
    D=closure(cells)
    errors=[]
    if closure(K)!=K or not K<=D:errors.append('invalid_alpha_complex')
    outside=set(g['outside_ancestors'])
    expected={i for i in range(len(cells)) if i not in outside and not g['tetra_in_alpha'][i]}
    flat=[i for r in t['regions'] for i in r['tetra_indices']]
    if len(flat)!=len(set(flat)) or set(flat)!=expected:errors.append('invalid_region_partition')
    rows=[]
    for r in t['regions']:
        local=[];ids=set(r['tetra_indices'])
        if not ids<=set(range(len(cells))):raise ValueError('Cell index out of range')
        C=closure(cells[i] for i in ids)
        faces={f for f,owners in face_owners.items() if any(i in ids for i in owners)
               and (len(owners)==1 or any(i not in ids for i in owners))}
        B=closure(faces);P=(C-B)-K;actual=unpack(r['dual_set'])
        if actual!=P:local.append('dual_set_mismatch')
        if unpack(r['closure_boundary'])!=B:local.append('boundary_mismatch')
        if component_count(P)!=1:local.append('dual_set_not_connected')
        if any(f in P and s not in P for s in D if len(s)>1 for f in combinations(s,len(s)-1)):
            local.append('dual_set_not_open')
        expected_mouth=B-K;mouth_sets=[unpack(m) for m in r['mouths']]
        merged=set().union(*mouth_sets) if mouth_sets else set()
        if merged!=expected_mouth or sum(map(len,mouth_sets))!=len(merged):
            local.append('mouth_partition_mismatch')
        if component_count(expected_mouth)!=len(mouth_sets):local.append('mouth_component_count_mismatch')
        mouth_checks=[]
        for j,S in enumerate(mouth_sets):
            check=surface_checks(S);check['mouth_id']=j
            check['n_triangles']=sum(len(s)==3 for s in S)
            if check['components']!=1:local.append('disconnected_mouth_'+str(j))
            if not check['n_triangles']:local.append('mouth_without_triangles_'+str(j))
            if not check['local_surface_checks_passed']:local.append('mouth_local_surface_failure_'+str(j))
            mouth_checks.append(check)
        raw=surface_checks(B)
        bad_in=sum(tuple(x['edge']) in K for x in raw['included_edge_failures'])
        if local:errors.append({'region_id':r['region_id'],'errors':local})
        rows.append({'region_id':r['region_id'],'n_tetrahedra':len(ids),'n_mouths':len(mouth_sets),
                     'combinatorial_errors':local,'raw_boundary':raw,
                     'bad_raw_edges_in_alpha':bad_in,
                     'bad_raw_edges_outside_alpha':len(raw['included_edge_failures'])-bad_in,
                     'mouth_checks':mouth_checks})
    warnings=[]
    if any(not x['raw_boundary']['local_surface_checks_passed'] for x in rows):
        warnings.append('Raw closure boundary is not locally a surface; relation to original construction requires review.')
    return {'execution_status':'completed','validation_status':'failed' if errors else 'passed_checks_with_warnings' if warnings else 'passed_implemented_checks',
            'original_CAST_equivalence_certified':False,'errors':errors,'warnings':warnings,
            'checks_scope':'Saved simplex topology only. No recomputation of exact geometry, symbolic perturbation or physical measures.',
            'regions':rows}

def main():
    p=argparse.ArgumentParser();p.add_argument('geometry',type=Path);p.add_argument('topology',type=Path)
    a=p.parse_args();g=json.loads(a.geometry.read_text(encoding='utf-8'));t=json.loads(a.topology.read_text(encoding='utf-8'))
    report=audit(g,t)
    report['input_sha256']={str(f.resolve()):hashlib.sha256(f.read_bytes()).hexdigest() for f in (a.geometry,a.topology)}
    dest=a.topology.resolve().parent/('validation_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'));dest.mkdir(exist_ok=False)
    (dest/'validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    columns=['region_id','n_tetrahedra','n_mouths','raw_bad_edges','bad_raw_edges_in_alpha','bad_raw_edges_outside_alpha','raw_bad_vertex_links','mouth_bad_edges','mouth_bad_vertex_links','combinatorial_errors']
    summaries=[]
    for r in report['regions']:
        row={k:r[k] for k in columns[:3]}
        row.update(raw_bad_edges=len(r['raw_boundary']['included_edge_failures']),
                   bad_raw_edges_in_alpha=r['bad_raw_edges_in_alpha'],bad_raw_edges_outside_alpha=r['bad_raw_edges_outside_alpha'],
                   raw_bad_vertex_links=len(r['raw_boundary']['included_vertex_link_failures']),
                   mouth_bad_edges=sum(len(m['included_edge_failures']) for m in r['mouth_checks']),
                   mouth_bad_vertex_links=sum(len(m['included_vertex_link_failures']) for m in r['mouth_checks']),
                   combinatorial_errors=len(r['combinatorial_errors']))
        summaries.append(row)
    with (dest/'validation_summary.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(summaries)
    print(json.dumps({k:report[k] for k in ('execution_status','validation_status','original_CAST_equivalence_certified','warnings','errors')},indent=2))
    print(json.dumps(summaries,indent=2));print('OUTPUT:',dest)
    if report['errors']:raise SystemExit(1)
if __name__=='__main__':main()
