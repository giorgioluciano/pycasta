import argparse,csv,json,math
from pathlib import Path
import numpy as np
from scipy.integrate import quad

def area(t,x,r):
    a,b,c=t; u=b-a; u=u/np.linalg.norm(u); n=np.cross(b-a,c-a); total=np.linalg.norm(n)/2
    if total==0: raise ValueError('Degenerate triangle')
    n=n/(2*total); v=np.cross(n,u); xy=(t-a)@np.array([u,v]).T
    lo,hi=min(xy[:,0]),max(xy[:,0]); xx=x-a; h=xx@n; rr=r*r-h*h
    q=xx@np.array([u,v]).T; rad=np.sqrt(np.maximum(rr,0)); good=(rr>0)&np.all(q+rad[:,None]>xy.min(0),axis=1)&np.all(q-rad[:,None]<xy.max(0),axis=1)
    disks=[(s,y,R) for (s,y),R in zip(q[good],rad[good])]; edges=[]
    for P,Q in zip(xy,np.roll(xy,1,axis=0)):
        if P[0]!=Q[0]: edges.append((min(P[0],Q[0]),max(P[0],Q[0]),P,(Q[1]-P[1])/(Q[0]-P[0])))
    def width(s):
        ys=[P[1]+k*(s-P[0]) for l,h,P,k in edges if l<=s<=h]; l,h=min(ys),max(ys); seg=[]
        for cx,cy,R in disks:
            z=R*R-(s-cx)**2
            if z>0:
                z=math.sqrt(z); L,H=max(l,cy-z),min(h,cy+z)
                if H>L:seg.append((L,H))
        covered=0.; end=l
        for L,H in sorted(seg):
            covered+=max(0.,H-max(end,L));end=max(end,H)
        return max(0.,h-l-covered)
    pts=set(xy[:,0]); pts.update(s+sign*R for s,y,R in disks for sign in [-1,1]); pts.update(np.linspace(lo,hi,17)); pts=sorted(p for p in pts if lo<p<hi)
    value,err=quad(width,lo,hi,points=pts,epsabs=1e-7,limit=max(200,len(pts)*4));return float(total),float(np.clip(value,0,total)),err

def main():
    p=argparse.ArgumentParser();p.add_argument('results',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args();d=json.loads(a.results.read_text(encoding='utf-8-sig'));m=d['provenance']['manifest'];x=np.asarray(m['coordinates']);r=np.asarray(m['atomic_radii']);rows=[]; cache={}
    for reg in d['regions']:
        for i,mo in enumerate(reg['mouths']):
            faces={tuple(sorted(f)) for f in mo['triangles']}; vals=[]
            for f in sorted(faces):
                if f not in cache:cache[f]=area(x[list(f)],x,r)
                vals.append(cache[f])
            t,v,e=np.sum(vals,axis=0) if vals else (0.,0.,0.)
            rows.append(dict(region_id=reg['region_id'],mouth_index=i,triangles=len(faces),support_area_A2=t,free_planar_area_A2=v,quadrature_error_estimate_A2=e))
    with a.out.open('x',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['region_id','mouth_index','triangles','support_area_A2','free_planar_area_A2','quadrature_error_estimate_A2']);w.writeheader();w.writerows(rows)
    print(len(rows),'mouths;',sum(q['free_planar_area_A2'] for q in rows),'A^2 free planar area;',a.out)
if __name__=='__main__':main()
