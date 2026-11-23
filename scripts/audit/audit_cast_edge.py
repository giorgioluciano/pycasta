"""Local link audit around an input-atom edge. No geometry changes."""
import argparse,json,csv
from pathlib import Path
from datetime import datetime
from collections import defaultdict

def ordered_star(g,edge):
    edge=set(edge);link=defaultdict(list);cells={}
    for i,t in enumerate(g['cells_idx']):
        if edge<=set(t):
            pair=tuple(sorted(set(t)-edge))
            if len(pair)!=2:raise ValueError('Invalid tetrahedron')
            cells[i]=pair
            for v in pair:link[v].append(i)
    if not cells:raise ValueError('Edge absent')
    if any(len(x)!=2 for x in link.values()):
        raise ValueError('Link is not a closed cycle: hull or invalid configuration')
    start=min(link);v=start;previous=None;ordered=[];used=set()
    while True:
        choices=[i for i in link[v] if i!=previous]
        i=min(choices) if previous is None else choices[0]
        if i in used:raise ValueError('Repeated tetrahedron before cycle completion')
        used.add(i);a,b=cells[i];w=b if a==v else a
        ordered.append((i,v,w));previous=i;v=w
        if v==start:break
    if used!=set(cells):raise ValueError('Disconnected edge link')
    return ordered

def runs(values):
    if not any(values):return 0
    if all(values):return 1
    return sum(v and not values[i-1] for i,v in enumerate(values))

def inspect(g,t,edge,region):
    star=ordered_star(g,edge);selected=set(t['regions'][region]['tetra_indices'])
    outside=set(g['outside_ancestors']);rows=[]
    for i,a,b in star:
        inside=bool(g['tetra_in_alpha'][i]);out=i in outside
        rows.append(dict(cell_id=i,link_from=a,link_to=b,in_region=i in selected,
                         in_alpha=inside,reaches_infinity=out,
                         retained_before_alpha=not out,
                         eligible_after_alpha=not out and not inside))
    alpha_edge=next(bool(v) for a,b,v in g['alpha_edge_rows'] if tuple(sorted((a,b)))==tuple(sorted(edge)))
    return dict(edge=list(edge),region_id=region,edge_in_alpha=alpha_edge,
                cyclic_star=rows,
                pre_alpha_selected_runs=runs([x['retained_before_alpha'] for x in rows]),
                post_alpha_selected_runs=runs([x['eligible_after_alpha'] for x in rows]),
                region_runs=runs([x['in_region'] for x in rows]),
                note='Cycle order is combinatorial, up to reversal. This diagnoses sectors; it does not repair topology.')

def main():
    p=argparse.ArgumentParser();p.add_argument('geometry',type=Path);p.add_argument('topology',type=Path)
    p.add_argument('--edge',nargs=2,type=int,default=[5,6]);p.add_argument('--region',type=int,default=0)
    a=p.parse_args();g=json.loads(a.geometry.read_text());t=json.loads(a.topology.read_text())
    report=inspect(g,t,a.edge,a.region)
    dest=a.topology.resolve().parent/('edge_audit_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'));dest.mkdir(exist_ok=False)
    (dest/'edge_audit.json').write_text(json.dumps(report,indent=2))
    with (dest/'edge_star.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(report['cyclic_star'][0]));w.writeheader();w.writerows(report['cyclic_star'])
    print(json.dumps(report,indent=2));print('OUTPUT:',dest)
if __name__=='__main__':main()
