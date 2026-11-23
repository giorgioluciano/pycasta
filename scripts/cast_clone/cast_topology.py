"""Combinatorial dual sets and mouths; no geometric floating-point tests."""
from itertools import combinations
from collections import defaultdict

def closure(simplices):
    result=set()
    for simplex in simplices:
        for size in range(1,len(simplex)+1):
            result.update(combinations(simplex,size))
    return result

def packed(simplices):
    return {name:[list(s) for s in sorted(simplices) if len(s)==dim+1]
            for dim,name in enumerate(('vertices','edges','triangles','tetrahedra'))}

def incidence_components(simplices):
    adj={s:set() for s in simplices}
    for s in simplices:
        if len(s)>1:
            for f in combinations(s,len(s)-1):
                if f in simplices:
                    adj[s].add(f);adj[f].add(s)
    remaining=set(simplices);groups=[]
    while remaining:
        seed=min(remaining);remaining.remove(seed);stack=[seed];group=set()
        while stack:
            s=stack.pop();group.add(s)
            for f in adj[s]:
                if f in remaining:
                    remaining.remove(f);stack.append(f)
        groups.append(group)
    return sorted(groups,key=lambda g:(-sum(len(s)==3 for s in g),min(g)))

def build_topology(data):
    if not data.get('ok'):
        raise ValueError('Flow not classified: degenerate or failed geometry')
    if not data.get('alpha_all_dimensions_exported'):
        raise ValueError('Rebuild increment 3: edge/vertex alpha classifications missing')
    cells=[tuple(sorted(map(int,t))) for t in data['cells_idx']]
    if len(set(cells))!=len(cells):
        raise ValueError('Duplicate tetrahedra')
    face_owners=defaultdict(list);K=set()
    for i,raw in enumerate(data['cells_idx']):
        if data['tetra_in_alpha'][i]:K.add(cells[i])
        for k in range(4):
            f=tuple(sorted(raw[j] for j in range(4) if j!=k))
            face_owners[f].append(i)
            if data['face_in_alpha'][i][k]:K.add(f)
    for a,b,inside in data['alpha_edge_rows']:
        if inside:K.add(tuple(sorted((int(a),int(b)))))
    for a,inside in data['alpha_vertex_rows']:
        if inside:K.add((int(a),))
    all_simplices=closure(cells)
    if not K<=all_simplices or closure(K)!=K:
        raise ValueError('Exported alpha complex is not closed under faces')
    if any(len(owners)>2 for owners in face_owners.values()):
        raise ValueError('Nonmanifold triangulation face')
    outside=set(data['outside_ancestors'])
    eligible={i for i in range(len(cells)) if i not in outside and not data['tetra_in_alpha'][i]}
    groups=data['pocket_tetra_components']
    flat=[i for g in groups for i in g]
    if len(flat)!=len(set(flat)) or set(flat)!=eligible:
        raise ValueError('Invalid tetrahedral partition')
    regions=[]
    for region_id,group in enumerate(groups):
        ids=set(group);C=closure(cells[i] for i in ids)
        boundary_faces={f for f,owners in face_owners.items()
                        if any(i in ids for i in owners)
                        and (len(owners)==1 or any(i not in ids for i in owners))}
        B=closure(boundary_faces);I=C-B;P=I-K
        if {s for s in P if len(s)==4}!={cells[i] for i in ids}:
            raise ValueError('Lost eligible tetrahedra')
        # Verify the open-set condition in the finite Delaunay complex.
        for t in all_simplices:
            if len(t)>1:
                for s in combinations(t,len(t)-1):
                    if s in P and t not in P:
                        raise ValueError('Dual set is not open in Delaunay complex')
        mouth_open=B-K
        mouths=incidence_components(mouth_open)
        for mouth in mouths:
            if not any(len(s)==3 for s in mouth):
                raise ValueError('Mouth component without triangles: inspect topology')
        edge_degrees=defaultdict(int)
        for f in boundary_faces:
            for e in combinations(f,2):edge_degrees[e]+=1
        nonmanifold=[list(e) for e,count in edge_degrees.items() if count!=2]
        regions.append({
            'region_id':region_id,'tetra_indices':sorted(ids),
            'dual_set':packed(P),'closure_boundary':packed(B),
            'mouths':[packed(g) for g in mouths],
            'n_mouths':len(mouths),
            'classification':('void' if not mouths else 'normal_pocket' if len(mouths)==1
                              else 'simple_connector' if len(mouths)==2 else 'multiple_connector'),
            'boundary_edges_degree_not_two':nonmanifold,
            'boundary_edge_manifold_check':not nonmanifold,
        })
    return {'status':'ok','regions':regions,
            'physical_volumes_computed':False,
            'symbolic_perturbation_implemented':False,
            'boundary_vertex_manifold_check_performed':False}
