# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Numerical restricted-Delaunay graph and boundary-rooted disk canonicalization.

No optimization.  The extraction follows the archived geometry.topology rules
but records the numerical margins, circle arc adjacency and all failure reasons.
Canonicalization uses faces and the boundary cycle, never center matching.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
import hashlib,json
from collections import Counter,defaultdict,deque
import numpy as np
from scipy.spatial import Delaunay,cKDTree


def canonical_disk(faces, boundary=None):
    faces=[tuple(map(int,f)) for f in faces]
    directed={}
    vertices=set()
    for f in faces:
        if len(set(f))!=3: raise ValueError('degenerate_combinatorial_face')
        vertices.update(f)
        for j in range(3):
            edge=(f[j],f[(j+1)%3])
            if edge in directed: raise ValueError('duplicate_oriented_edge')
            directed[edge]=f[(j+2)%3]
    boundary_edges=[e for e in directed if e[::-1] not in directed]
    following={a:b for a,b in boundary_edges}
    if len(following)!=len(boundary_edges): raise ValueError('boundary_branches')
    start=min(following)
    cycle=[]
    v=start
    while v not in cycle:
        cycle.append(v)
        if v not in following: raise ValueError('open_boundary')
        v=following[v]
    if v!=start or len(cycle)!=len(following): raise ValueError('multiple_boundary_cycles')
    if boundary is not None and set(boundary)!=set(cycle): raise ValueError('boundary_vertex_set_mismatch')
    # E=(3F+B)/2, directed includes each internal edge twice.
    E=(len(directed)+len(boundary_edges))//2
    if len(vertices)-E+len(faces)!=1: raise ValueError('not_disk_Euler')
    links=defaultdict(lambda:defaultdict(set))
    for f in faces:
        for j in range(3):
            x,a,b=f[j],f[(j+1)%3],f[(j+2)%3]
            links[x][a].add(b);links[x][b].add(a)
    boundary_vertices=set(cycle)
    for v,link in links.items():
        seen={next(iter(link))};todo=list(seen)
        while todo:
            x=todo.pop()
            for y in link[x]-seen:seen.add(y);todo.append(y)
        if seen!=set(link):raise ValueError('disconnected_vertex_link')
        ds=Counter(len(a) for a in link.values())
        if v in boundary_vertices:
            if ds[1]!=2 or set(ds)-{1,2}:raise ValueError('boundary_link_not_path')
        elif set(ds)!={2}:raise ValueError('interior_link_not_cycle')
    candidates=[]
    for reflected in (False,True):
        edge_third={(b,a):c for (a,b),c in directed.items()} if reflected else directed
        seeds=[(b,a) for a,b in boundary_edges] if reflected else boundary_edges
        for seed in seeds:
            a,b=seed;c=edge_third[seed]
            labels={a:0,b:1,c:2};todo=deque([(a,b,c)])
            seen=set();mapped=[]
            while todo:
                f=todo.popleft();key=frozenset(f)
                if key in seen: continue
                seen.add(key)
                mapped.append(tuple(sorted(labels[x] for x in f)))
                for j in range(3):
                    x,y=f[j],f[(j+1)%3]
                    if (y,x) in edge_third:
                        z=edge_third[y,x]
                        if z not in labels: labels[z]=len(labels)
                        todo.append((y,x,z))
            if len(seen)!=len(faces) or set(labels)!=vertices: raise ValueError('disconnected_faces')
            normalized=tuple(sorted(mapped))
            candidates.append((normalized,labels,reflected,seed))
    normalized,labels,reflected,seed=min(candidates,key=lambda x:x[0])
    payload={'vertices':len(vertices),'boundary':len(cycle),'faces':[list(f) for f in normalized]}
    fingerprint=hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return dict(fingerprint=fingerprint,canonical_faces=payload['faces'],vertex_to_canonical=labels,
                boundary_cycle=cycle,reflection=reflected,root_boundary_edge=list(seed))


def extract(C, full_details=False):
    C=np.asarray(C,dtype=float)
    if C.ndim!=2 or C.shape[1]!=2: raise ValueError('not_Nx2_array')
    if not np.isfinite(C).all(): raise ValueError('nonfinite_coordinates')
    DT=Delaunay(C).simplices.copy()
    P=C[DT];a=P[:,0];u=P[:,1]-a;v=P[:,2]-a
    det=u[:,0]*v[:,1]-u[:,1]*v[:,0]
    negative=det<0
    DT[negative]=DT[negative][:,[0,2,1]]
    P=C[DT];a=P[:,0];u=P[:,1]-a;v=P[:,2]-a
    det=u[:,0]*v[:,1]-u[:,1]*v[:,0]
    uu=np.sum(u*u,axis=1);vv=np.sum(v*v,axis=1)
    p=a+np.stack((uu*v[:,1]-vv*u[:,1],u[:,0]*vv-v[:,0]*uu),1)/(2*det[:,None])
    pnorm2=np.sum(p*p,axis=1)
    keep=pnorm2<=1+1e-12
    T=DT[keep]
    E=np.column_stack(np.triu_indices(len(C),1))
    uu=C[E[:,1]]-C[E[:,0]];ss=np.sum(uu*uu,axis=1)
    h=(np.sum(C[E[:,1]]**2,axis=1)-np.sum(C[E[:,0]]**2,axis=1))/2
    ok=(ss>1e-24)&(h*h<=ss+1e-14)
    E=E[ok];uu=uu[ok];ss=ss[ok];h=h[ok]
    E=np.repeat(E,2,axis=0);uu=np.repeat(uu,2,axis=0);ss=np.repeat(ss,2);h=np.repeat(h,2)
    sg=np.tile(np.array([-1,1]),len(E)//2)
    z=np.sqrt(np.maximum(1e-20,ss-h*h));up=np.stack((-uu[:,1],uu[:,0]),1)
    q=(h[:,None]*uu+(sg*z)[:,None]*up)/ss[:,None]
    g=np.sum((q-C[E[:,0]])**2,axis=1);nearest=cKDTree(C).query(q)[0]**2
    valid=g<=nearest+2e-11
    boundary_pairs=E[valid];boundary_signs=sg[valid];anchors=q[valid]
    norms=np.linalg.norm(C,axis=1);F=np.flatnonzero(norms>1e-12)
    qq=-C[F]/norms[F,None];gg=np.sum((qq-C[F])**2,axis=1);nn=cKDTree(C).query(qq)[0]**2
    singleton=F[gg<=nn+2e-11]
    canon=canonical_disk(T)
    cycle=canon['boundary_cycle'];bs=set(cycle)
    cycle_edges={tuple(sorted((cycle[j],cycle[(j+1)%len(cycle)]))) for j in range(len(cycle))}
    arc_pairs=[tuple(sorted(map(int,e))) for e in boundary_pairs]
    if Counter(arc_pairs)!=Counter(cycle_edges): raise ValueError('circle_arcs_and_face_boundary_mismatch')
    if len(singleton): raise ValueError('boundary_singleton_present')
    adjacency=defaultdict(set);qcounts=Counter();pc=Counter()
    for f in T:
        f=list(map(int,f));qf=sum(tuple(sorted((f[j],f[(j+1)%3]))) in cycle_edges for j in range(3))
        qcounts[qf]+=1
        for j in range(3):adjacency[f[j]].update((f[(j+1)%3],f[(j+2)%3]))
        if qf==1:
            ap=[x for x in f if x not in bs]
            if len(ap)==1:pc[ap[0]]+=1
    interior=sorted(set(range(len(C)))-bs)
    counts=Counter((pc[i],len(adjacency[i])) for i in interior)
    result=dict(status='PASS',numerical_only=True,n=len(C),total_DT_faces=len(DT),restricted_faces=len(T),
                boundary=len(cycle),interior=len(interior),arc_pairs=len(arc_pairs),singletons=len(singleton),
                q_face_histogram=dict(sorted(qcounts.items())),interior_p_degree_histogram={f'{p},{d}':n for (p,d),n in sorted(counts.items())},
                max_center_norm=float(norms.max()),inside_outside_circumcenter_margin=float(np.min(np.abs(pnorm2-1))),
                min_face_positive_determinant=float(det.min()),fingerprint=canon['fingerprint'])
    if full_details:
        result.update(canon)
        result.update(faces=T.tolist(),boundary_pairs=boundary_pairs.tolist(),boundary_signs=boundary_signs.tolist(),boundary_anchors=anchors.tolist())
    return result
