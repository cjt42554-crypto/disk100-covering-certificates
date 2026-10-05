# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Integer reconstruction of the incumbent's restricted Delaunay complex.

Consumes the saved binary64 center file as exact rationals and the numerical DT
as a proposal.  Every face orientation and all empty-circle signs are checked
with integers.  Dual edges/rays certify circle crossings without sampling.
This is a graph identity certificate, not a covering-radius certificate.
"""
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[n]='1'
import sys,json,hashlib,time
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np
from scipy.spatial import Delaunay
from graph_tools import canonical_disk
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parent/'best147_boundary_transfer/frozen_proposals/grid0025_R0.11544695488558258.npy'
if not SOURCE.exists():SOURCE=ROOT/'inputs/current.npy'


def dot(a,b):return a[0]*b[0]+a[1]*b[1]
def orient(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def main():
    tic=time.monotonic();C=np.load(SOURCE)
    ratios=[[float(t).as_integer_ratio() for t in row] for row in C]
    Q=max(d for row in ratios for n,d in row)
    assert all(Q%d==0 for row in ratios for n,d in row)
    P=[tuple(n*(Q//d) for n,d in row) for row in ratios]
    assert all(dot(p,p)<Q*Q for p in P),'center_not_strictly_inside_unit_disk'
    faces=[]
    for f in Delaunay(C).simplices:
        f=list(map(int,f))
        if orient(*(P[x] for x in f))<0:f[1],f[2]=f[2],f[1]
        assert orient(*(P[x] for x in f))>0
        faces.append(f)
    directed={};edges=defaultdict(list);circles=[];inside=[];checks=0;ties=0
    for fi,f in enumerate(faces):
        a,b,c=(P[x] for x in f);u=(b[0]-a[0],b[1]-a[1]);v=(c[0]-a[0],c[1]-a[1])
        D=2*orient(a,b,c);uu=dot(u,u);vv=dot(v,v)
        pr=(uu*v[1]-vv*u[1],u[0]*vv-v[0]*uu)
        N=(a[0]*D+pr[0],a[1]*D+pr[1]);rnum2=dot(pr,pr)
        for j,p in enumerate(P):
            if j in f:continue
            delta=(p[0]*D-N[0],p[1]*D-N[1]);sign=dot(delta,delta)-rnum2
            assert sign>=0,(fi,j,'nonempty_circle')
            ties+=sign==0;checks+=1
        side=(Q*D)**2-dot(N,N)
        assert side!=0,'circumcenter_exactly_on_unit_circle'
        inside.append(side>0);circles.append((N,D))
        for j in range(3):
            e=(f[j],f[(j+1)%3]);assert e not in directed
            directed[e]=f[(j+2)%3];edges[tuple(sorted(e))].append((fi,e))
    # Validate full face complex and exact hull orientation/convexity.
    whole=canonical_disk(faces);cycle=whole['boundary_cycle']
    assert len(whole['vertex_to_canonical'])==100
    for j,a in enumerate(cycle):
        b=cycle[(j+1)%len(cycle)]
        assert all(orient(P[a],P[b],p)>=0 for p in P),'boundary_not_convex_hull'
    cross=[];outout=[];hulloutside=[]
    for e,adj in edges.items():
        assert len(adj) in (1,2)
        if len(adj)==2:
            (f,ef),(g,eg)=adj;assert ef==eg[::-1]
            if inside[f]!=inside[g]:cross.append(e);continue
            if inside[f]:continue
            # Rational common-denominator dual segment.  Both endpoints outside.
            N1,D1=circles[f];N2,D2=circles[g]
            A=(N1[0]*D2,N1[1]*D2);B=(N2[0]*D1,N2[1]*D1);S=Q*D1*D2
            w=(B[0]-A[0],B[1]-A[1]);ww=dot(w,w);t=-dot(A,w)
            if 0<t<ww:
                assert dot(A,A)*ww-t*t>S*S*ww,('extra_outside_dual_crossing',e)
            outout.append(e)
        else:
            f,ef=adj[0]
            if inside[f]:cross.append(e);continue
            # Voronoi ray points to the right of the CCW hull edge.
            a,b=(P[x] for x in ef);w=(b[1]-a[1],a[0]-b[0]);N,D=circles[f]
            ww=dot(w,w);t=-dot(N,w)
            if t>0:assert dot(N,N)*ww-t*t>(Q*D)**2*ww,('extra_hull_ray_crossing',e)
            hulloutside.append(e)
    restricted=[f for f,keep in zip(faces,inside) if keep]
    ca=canonical_disk(restricted)
    rc=ca['boundary_cycle'];expected={tuple(sorted((a,rc[(j+1)%len(rc)]))) for j,a in enumerate(rc)}
    assert set(cross)==expected and len(cross)==len(expected),'dual_crossing_boundary_mismatch'
    numerical=json.loads((ROOT/'current_graph_numerical.json').read_text())
    assert ca['fingerprint']==numerical['fingerprint']
    row=dict(status='PASS',source=str(SOURCE),source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
             rational_common_denominator=str(Q),all_centers=100,full_DT_faces=len(faces),restricted_DT_faces=len(restricted),
             all_centers_exact_strictly_inside_unitdisk=True,
             hull_boundary=len(cycle),restricted_boundary=len(rc),emptycircle_integer_checks=checks,exact_emptycircle_other_site_ties=ties,
             all_exact_orientations_positive=True,all_circumcenter_unitdisk_signs_strict=True,
             exact_boundary_dual_crossings=len(cross),outside_outside_segments_certified_miss_disk=len(outout),
             outside_hull_rays_certified_miss_disk=len(hulloutside),fingerprint=ca['fingerprint'],
             boundary_cycle=rc,faces=restricted,vertex_to_canonical=ca['vertex_to_canonical'],
             seconds=time.monotonic()-tic,scope='Exact restricted-DT combinatorial identity of the saved binary64 coordinates; no radius or catalog-completeness certificate')
    target=Path(sys.argv[1]) if len(sys.argv)>1 else ROOT/'current_graph_integer_verified.json'
    target.write_text(json.dumps(row,indent=2)+'\n')
    print(json.dumps({k:v for k,v in row.items() if k not in ('faces','vertex_to_canonical','boundary_cycle')},indent=2))

if __name__=='__main__':main()
