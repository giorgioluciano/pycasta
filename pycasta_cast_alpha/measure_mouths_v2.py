import argparse,csv,json,math
from itertools import combinations
from pathlib import Path
import numpy as np
from scipy.integrate import quad

def measure(t,x,r):
    origin=t[0]; u=t[1]-origin; normal=np.cross(u,t[2]-origin); area=np.linalg.norm(normal)/2
    if not area>0: raise ValueError('Degenerate triangle')
    u/=np.linalg.norm(u); n=normal/(2*area); basis=np.array([u,np.cross(n,u)]).T
    tri=(t-origin)@basis; q=(x-origin)@basis; rr=r*r-((x-origin)@n)**2; rad=np.sqrt(np.maximum(rr,0))
    mask=(rr>0)&np.all(q+rad[:,None]>=tri.min(0),axis=1)&np.all(q-rad[:,None]<=tri.max(0),axis=1)
    disks=list(zip(q[mask],rad[mask])); edges=list(zip(tri,np.roll(tri,-1,axis=0))); cuts=list(tri[:,0]); lo=min(cuts); hi=max(cuts)
    for c,R in disks:
        cuts.extend([c[0]-R,c[0]+R])
        for a,b in edges:
            v=b-a; z=a-c; A=v@v; B=z@v; disc=B*B-A*(z@z-R*R)
            if disc>=0:
                for s in ((-B-math.sqrt(disc))/A,(-B+math.sqrt(disc))/A):
                    if 0<s<1:cuts.append(a[0]+s*v[0])
    for (a,R),(b,S) in combinations(disks,2):
        v=b-a; d=np.linalg.norm(v)
        if abs(R-S)<d<R+S:
            h=(R*R-S*S+d*d)/(2*d); off=math.sqrt(max(0,R*R-h*h)); mid=a[0]+h*v[0]/d; dx=off*v[1]/d;cuts.extend([mid-dx,mid+dx])
    cuts=sorted(set(c for c in cuts if lo<=c<=hi));
    def width(s):
        ys=[a[1]+(s-a[0])*(b[1]-a[1])/(b[0]-a[0]) for a,b in edges if min(a[0],b[0])<s<max(a[0],b[0])];
        if len(ys)<2:return 0.
        l,h=min(ys),max(ys); seg=[]
        for (cx,cy),R in disks:
            z=R*R-(s-cx)**2
            if z>0:
                z=math.sqrt(z); L=max(l,cy-z); H=min(h,cy+z)
                if H>L:seg.append((L,H))
        end=l; free=0.
        for L,H in sorted(seg):free+=max(0.,L-end);end=max(end,H)
        return free+h-end
    val=err=0.; notes=[]
    for a,b in zip(cuts,cuts[1:]):
        if b-a<=1e-13:continue
        ans=quad(width,a,b,epsabs=1e-9,epsrel=1e-8,full_output=True);val+=ans[0];err+=ans[1]
        if len(ans)>3:notes.append(ans[3])
    if not -1e-8<=val<=area+1e-8:notes.append('Area outside physical bounds')
    return area,val,err,' | '.join(dict.fromkeys(notes))

def main():
    p=argparse.ArgumentParser();p.add_argument('results',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args(); data=json.loads(a.results.read_text(encoding='utf-8-sig'));m=data['provenance']['manifest'];x=np.array(m['coordinates']);r=np.array(m['atomic_radii']); a.out.mkdir(parents=True,exist_ok=False); cache={}; rows=[]; faces=[]
    for reg in data['regions']:
        for i,mo in enumerate(reg['mouths']):
            group=[]
            for f in sorted({tuple(sorted(f)) for f in mo['triangles']}):
                if f not in cache:cache[f]=measure(x[list(f)],x,r)
                v=cache[f];group.append(v);faces.append((reg['region_id'],i,*f,*v))
            rows.append((reg['region_id'],i,len(group),*(sum(v[k] for v in group) for k in range(3)),sum(bool(v[3]) for v in group)))
    for name,header,body in [('faces.csv','region,mouth,a,b,c,support_A2,free_A2,error_estimate_A2,note',faces),('mouths.csv','region,mouth,triangles,support_A2,free_A2,error_estimate_A2,flagged_faces',rows)]:
        with (a.out/name).open('x',newline='',encoding='utf-8') as f:w=csv.writer(f);w.writerow(header.split(','));w.writerows(body)
    print(len(rows),sum(row[4] for row in rows),len(cache),sum(bool(v[3]) for v in cache.values()),sorted([row[4] for row in rows],reverse=True)[:3])
if __name__=='__main__':main()
