# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Integer full/restricted-DT identity audit for an arbitrary seed proposal."""
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[n]='1'
from collections import defaultdict
import numpy as np
from scipy.spatial import Delaunay
from graph_tools import extract,canonical_disk


def dot(a,b):return a[0]*b[0]+a[1]*b[1]
def orient(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def verify(C):
    C=np.asarray(C,float);ratios=[[float(x).as_integer_ratio() for x in row] for row in C]
    Q=max(d for row in ratios for _,d in row);P=[tuple(n*(Q//d) for n,d in row) for row in ratios]
    assert len(set(P))==len(P) and all(dot(p,p)<Q*Q for p in P)
    faces=[];circles=[];inside=[];edges=defaultdict(list);checks=0
    for tri in Delaunay(C).simplices:
        f=list(map(int,tri))
        if orient(*(P[j] for j in f))<0:f[1],f[2]=f[2],f[1]
        a,b,c=(P[j] for j in f);D=2*orient(a,b,c);assert D>0
        u=(b[0]-a[0],b[1]-a[1]);v=(c[0]-a[0],c[1]-a[1]);uu=dot(u,u);vv=dot(v,v)
        pr=(uu*v[1]-vv*u[1],u[0]*vv-v[0]*uu);N=(a[0]*D+pr[0],a[1]*D+pr[1]);r2=dot(pr,pr)
        for j,p in enumerate(P):
            if j in f:continue
            dp=(p[0]*D-N[0],p[1]*D-N[1]);sign=dot(dp,dp)-r2
            assert sign>0,('emptycircle_not_strict',len(faces),j,'TIE' if sign==0 else 'NEGATIVE',f)
            checks+=1
        side=(Q*D)**2-dot(N,N);assert side!=0
        fi=len(faces);faces.append(f);circles.append((N,D));inside.append(side>0)
        for j in range(3):edges[tuple(sorted((f[j],f[(j+1)%3])))].append((fi,(f[j],f[(j+1)%3])))
    full=canonical_disk(faces);bc=full['boundary_cycle'];assert len(full['vertex_to_canonical'])==len(C)
    for j,a in enumerate(bc):assert all(orient(P[a],P[bc[(j+1)%len(bc)]],p)>=0 for p in P)
    crossed=[];outside_segments=0;outside_rays=0
    for e,ad in edges.items():
        assert len(ad) in (1,2)
        if len(ad)==2:
            (f,ef),(g,eg)=ad;assert ef==eg[::-1]
            if inside[f]!=inside[g]:crossed.append(e);continue
            if inside[f]:continue
            N1,D1=circles[f];N2,D2=circles[g];A=(N1[0]*D2,N1[1]*D2);B=(N2[0]*D1,N2[1]*D1);S=Q*D1*D2
            w=(B[0]-A[0],B[1]-A[1]);ww=dot(w,w);t=-dot(A,w)
            if 0<t<ww:assert dot(A,A)*ww-t*t>S*S*ww,('outside_dual_segment_hits_disk',e)
            outside_segments+=1
        else:
            f,ef=ad[0]
            if inside[f]:crossed.append(e);continue
            a,b=(P[j] for j in ef);w=(b[1]-a[1],a[0]-b[0]);N,D=circles[f];ww=dot(w,w);t=-dot(N,w)
            if t>0:assert dot(N,N)*ww-t*t>(Q*D)**2*ww,('outside_dual_ray_hits_disk',e)
            outside_rays+=1
    restricted=[f for f,k in zip(faces,inside) if k];ca=canonical_disk(restricted);rc=ca['boundary_cycle']
    assert set(crossed)=={tuple(sorted((a,rc[(j+1)%len(rc)]))) for j,a in enumerate(rc)} and len(crossed)==len(rc)
    numerical=extract(C,True);assert ca['fingerprint']==numerical['fingerprint']
    return dict(status='PASS',common_denominator=str(Q),all_centers=len(C),all_centers_exact_strictly_inside_unitdisk=True,
                full_DT_faces=len(faces),restricted_DT_faces=len(restricted),full_hull_boundary=len(bc),restricted_boundary=len(rc),
                integer_emptycircle_checks=checks,emptycircle_other_site_ties=0,all_face_orientations_strict=True,
                all_circumcenter_unitdisk_signs_strict=True,exact_boundary_dual_disk_crossings=len(crossed),
                outside_dual_segments_certified_miss_disk=outside_segments,outside_hull_rays_certified_miss_disk=outside_rays,
                exact_restricted_fingerprint=ca['fingerprint'],faces=restricted,boundary_cycle=rc,vertex_to_canonical=ca['vertex_to_canonical'],
                scope='Integer full/restricted-DT and disk-dual crossing identity of the exact binary64 seed; no radius/completeness/optimality assertion')
