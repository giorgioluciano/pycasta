"""Classify boundary contacts without changing geometry or topology."""
import argparse
import csv
import hashlib
import json
from collections import Counter,defaultdict
from datetime import datetime
from itertools import combinations
from pathlib import Path

def classify_contacts(g,t):
    if not g.get('ok') or not g.get('alpha_all_dimensions_exported'):
        raise ValueError('Need classified, nondegenerate, all-dimensional complex')
    n=len(g['cells_idx']);outside=set(g['outside_ancestors'])
    region_of={}
    for r in t['regions']:
        for i in r['tetra_indices']:
            if i in region_of:raise ValueError('Overlapping regions')
            region_of[i]=r['region_id']
    eligible={i for i in range(n) if i not in outside and not g['tetra_in_alpha'][i]}
    if set(region_of)!=eligible:raise ValueError('Region partition mismatch')
    alpha_edges={tuple(sorted((a,b))):bool(v) for a,b,v in g['alpha_edge_rows']}
    face_cells=defaultdict(list);edge_cells=defaultdict(list)
    for i,cell in enumerate(g['cells_idx']):
        for face in combinations(sorted(cell),3):face_cells[face].append(i)
        for edge in combinations(sorted(cell),2):edge_cells[edge].append(i)
    if any(len(owners)>2 for owners in face_cells.values()):
        raise ValueError('Nonmanifold triangulation face')
    def boundary(selected):
        return [f for f,owners in face_cells.items()
                if any(i in selected for i in owners)
                and (len(owners)==1 or any(i not in selected for i in owners))]
    def bad_edges(faces):
        degree=Counter(e for f in faces for e in combinations(f,2))
        return {e:d for e,d in degree.items() if d!=2}
    def sectors(edge,selected):
        ids=set(edge_cells[edge]);adj={i:set() for i in ids}
        for i in ids:
            for k,j in enumerate(g['neighbors'][i]):
                if j in ids:
                    face=set(g['cells_idx'][i]);face.remove(g['cells_idx'][i][k])
                    if set(edge)<=face:adj[i].add(j)
        remaining=ids&selected;count=0
        while remaining:
            count+=1;stack=[remaining.pop()]
            while stack:
                for j in adj[stack.pop()]:
                    if j in remaining:remaining.remove(j);stack.append(j)
        return count
    records=[]
    pre=set(range(n))-outside
    for e,d in sorted(bad_edges(boundary(pre)).items()):
        regions=sorted({region_of[i] for i in edge_cells[e] if i in region_of})
        in_alpha=alpha_edges[e]
        category=('removed_alpha_edge' if in_alpha else
                  'cross_region_contact' if len(regions)>1 else
                  'requires_review')
        records.append(dict(scope='global_pre_alpha',region_id=None,edge=list(e),
                            boundary_degree=d,in_alpha=in_alpha,regions=regions,
                            selected_sectors=sectors(e,pre),category=category))
    for r in t['regions']:
        selected=set(r['tetra_indices'])
        for e,d in sorted(bad_edges(boundary(selected)).items()):
            in_alpha=alpha_edges[e]
            records.append(dict(scope='region_raw_boundary',region_id=r['region_id'],
                                edge=list(e),boundary_degree=d,in_alpha=in_alpha,
                                regions=sorted({region_of[i] for i in edge_cells[e] if i in region_of}),
                                selected_sectors=sectors(e,selected),
                                category='removed_alpha_edge' if in_alpha else 'same_region_edge_requires_review'))
    mouth_failures=[]
    for r in t['regions']:
        for j,m in enumerate(r['mouths']):
            inc=Counter(e for f in m['triangles'] for e in combinations(sorted(f),2))
            for raw in m['edges']:
                e=tuple(sorted(raw))
                if alpha_edges[e]:raise ValueError('Mouth includes an alpha edge')
                if inc[e]!=2:
                    mouth_failures.append(dict(region_id=r['region_id'],mouth_id=j,
                                               edge=list(e),triangle_degree=inc[e]))
    counts={scope:dict(Counter(r['category'] for r in records if r['scope']==scope))
            for scope in ('global_pre_alpha','region_raw_boundary')}
    unresolved=[r for r in records if r['category'] in
                ('requires_review','same_region_edge_requires_review')]
    return dict(execution_status='completed',
                contact_status='review_required' if unresolved or mouth_failures else 'all_edge_contacts_categorized',
                original_CAST_equivalence_certified=False,
                scope_note='Edge contacts only. Does not replace vertex-link checks or resolve original boundary-manifold assumptions.',
                counts=counts,unresolved_contact_count=len(unresolved),
                included_mouth_edge_failure_count=len(mouth_failures),
                included_mouth_edge_failures=mouth_failures,records=records)

def main():
    p=argparse.ArgumentParser();p.add_argument('geometry',type=Path);p.add_argument('topology',type=Path)
    a=p.parse_args();g=json.loads(a.geometry.read_text(encoding='utf-8'));t=json.loads(a.topology.read_text(encoding='utf-8'))
    report=classify_contacts(g,t)
    report['input_sha256']={str(x.resolve()):hashlib.sha256(x.read_bytes()).hexdigest() for x in (a.geometry,a.topology)}
    dest=a.topology.resolve().parent/('contacts_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'));dest.mkdir(exist_ok=False)
    (dest/'contacts.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    with (dest/'contacts.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['scope','region_id','edge','boundary_degree','in_alpha','regions','selected_sectors','category']);w.writeheader();w.writerows(report['records'])
    print(json.dumps({k:v for k,v in report.items() if k not in ('records','input_sha256','included_mouth_edge_failures')},indent=2))
    print('OUTPUT:',dest)
    if report['unresolved_contact_count'] or report['included_mouth_edge_failure_count']:raise SystemExit(1)
if __name__=='__main__':main()
