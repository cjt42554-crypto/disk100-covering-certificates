# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Exhaustive all-triple numerical replay, independent of Delaunay. Not exact certification."""
import itertools
import numpy as np
def all_triple_cover(C):
    C=np.asarray(C,dtype=float);n=len(C);T=np.array(list(itertools.combinations(range(n),3)),dtype=int);P=C[T];u=P[:,1]-P[:,0];v=P[:,2]-P[:,0];cross=u[:,0]*v[:,1]-u[:,1]*v[:,0];ok=np.abs(cross)>1e-14;a=P[ok,0];u=u[ok];v=v[ok];cross=cross[ok];uu=np.sum(u*u,axis=1);vv=np.sum(v*v,axis=1);q=a+np.column_stack((uu*v[:,1]-vv*u[:,1],u[:,0]*vv-v[:,0]*uu))/(2*cross[:,None]);q=q[np.sum(q*q,axis=1)<=1]
    ii,jj=np.triu_indices(n,1);u=C[jj]-C[ii];ss=np.sum(u*u,axis=1);h=(np.sum(C[jj]**2,axis=1)-np.sum(C[ii]**2,axis=1))/2;ok=(ss>1e-20)&(ss-h*h>=0);u=u[ok];ss=ss[ok];h=h[ok];z=np.sqrt(ss-h*h);perp=np.column_stack((-u[:,1],u[:,0]));mid=h[:,None]*u/ss[:,None];offset=z[:,None]*perp/ss[:,None];bp=np.r_[mid-offset,mid+offset];norms=np.linalg.norm(C,axis=1);nz=norms>1e-14;ap=-C[nz]/norms[nz,None];points=np.r_[q,bp,ap,np.eye(2),-np.eye(2)];best=-1.;worst=None
    for j in range(0,len(points),512):
        chunk=points[j:j+512];nearest=np.min(np.sum((chunk[:,None,:]-C[None,:,:])**2,axis=2),axis=1);at=int(np.argmax(nearest))
        if nearest[at]>best:best=float(nearest[at]);worst=chunk[at].tolist()
    return float(np.sqrt(best)),dict(sites=n,triples=len(T),nondegenerate_triples=len(a),inside_triple_candidates=len(q),pairs=len(ii),boundary_pair_candidates=len(bp),antipodes=len(ap),worst_point=worst,numerical_only=True,method='all-triples+all-pair-circle+antipodes+axes actual-nearest-distance, independent of Delaunay')
