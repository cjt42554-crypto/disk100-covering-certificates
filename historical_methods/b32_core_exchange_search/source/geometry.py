# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Complete numerical covering and dynamic SLSQP engine for the mixed-mode branch.

Keep the materialized input immutable. All 100 centers are optimization variables.
The original center order/orientation is retained throughout.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
import concurrent.futures as futures
import hashlib
import json
from pathlib import Path
import time
import traceback
import numpy as np
from scipy.spatial import Delaunay, cKDTree
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent
INPUT = ROOT.parent/'library_import/latest_modes/materialized/disk100_branch_1155018_search_modes_20261004'
CONFIG_SHA = '024276d929e0325c6ad7447bb55779bcb17da8aa33e8463cd8fa9e4869415461'
CURRENT_REFERENCE = .11550182419454742434
PAIR_INDEX = np.array([(i,j) for i in range(100) for j in range(i+1,100)], dtype=int)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def atomic_json(path,data):
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    tmp.replace(path)

def triangles(C,T):
    P=C[T];a=P[:,0];u=P[:,1]-a;v=P[:,2]-a
    d=u[:,0]*v[:,1]-u[:,1]*v[:,0]
    uu=np.einsum('ij,ij->i',u,u);vv=np.einsum('ij,ij->i',v,v)
    pr=np.stack((uu*v[:,1]-vv*u[:,1],u[:,0]*vv-v[:,0]*uu),1)/(2*d[:,None])
    p=a+pr;g=np.einsum('ij,ij->i',pr,pr)
    w1=(pr[:,0]*v[:,1]-pr[:,1]*v[:,0])/d
    w2=(u[:,0]*pr[:,1]-u[:,1]*pr[:,0])/d
    w=np.stack((1-w1-w2,w1,w2),1)
    return g,p,2*w[:,:,None]*(P-p[:,None,:])

def boundary(C,E,sg):
    a,b=C[E[:,0]],C[E[:,1]];u=b-a
    ss=np.einsum('ij,ij->i',u,u)
    h=(np.einsum('ij,ij->i',b,b)-np.einsum('ij,ij->i',a,a))/2
    z=np.sqrt(np.maximum(1e-20,ss-h*h))
    up=np.stack((-u[:,1],u[:,0]),1)
    q=(h[:,None]*u+(sg*z)[:,None]*up)/ss[:,None]
    tangent=np.stack((-q[:,1],q[:,0]),1)
    wi=np.einsum('ij,ij->i',b,tangent)/np.einsum('ij,ij->i',u,tangent)
    w=np.stack((wi,1-wi),1)
    return np.einsum('ij,ij->i',q-a,q-a),q,2*w[:,:,None]*(C[E]-q[:,None,:])

def single_boundary(C,F):
    A=C[F];norm=np.linalg.norm(A,axis=1)
    return (1+norm)**2, 2*(1+1/norm)[:,None]*A

def topology(C):
    # Complete candidate set for the maximum of nearest-site distance on unit disk:
    # interior Voronoi vertices, active pairwise bisector/unit-circle crossings,
    # and an active single site's antipodal boundary stationary point.
    T=Delaunay(C).simplices.copy()
    _,p,_=triangles(C,T)
    T=T[np.sum(p*p,axis=1)<=1+1e-12]
    E=np.column_stack(np.triu_indices(len(C),1))
    u=C[E[:,1]]-C[E[:,0]]
    ss=np.sum(u*u,axis=1)
    h=(np.sum(C[E[:,1]]**2,axis=1)-np.sum(C[E[:,0]]**2,axis=1))/2
    ok=(ss>1e-24)&(h*h<=ss+1e-14)
    E=E[ok]
    E=np.repeat(E,2,axis=0);sg=np.tile(np.array([-1,1]),len(E)//2)
    g,q,_=boundary(C,E,sg)
    nearest=cKDTree(C).query(q)[0]**2
    valid=g<=nearest+2e-11
    E,sg=E[valid],sg[valid]
    norms=np.linalg.norm(C,axis=1)
    F=np.flatnonzero(norms>1e-12)
    q=-C[F]/norms[F,None]
    g=np.sum((q-C[F])**2,axis=1)
    nearest=cKDTree(C).query(q)[0]**2
    F=F[g<=nearest+2e-11]
    return T,E,sg,F

def covering(C):
    top=topology(C);T,E,sg,F=top
    tg=triangles(C,T)[0] if len(T) else np.empty(0)
    bg=boundary(C,E,sg)[0] if len(E) else np.empty(0)
    fg=single_boundary(C,F)[0] if len(F) else np.empty(0)
    return float(np.sqrt(max(tg.max(initial=0),bg.max(initial=0),fg.max(initial=0)))),top

def top_key(top):
    T,E,sg,F=top
    return (tuple(sorted(tuple(sorted(t)) for t in T)),tuple(sorted((int(e[0]),int(e[1]),int(s)) for e,s in zip(E,sg))),tuple(F))

class Model:
    def __init__(self,C,top):
        self.T,self.E,self.sg,self.F=top
        self.n=len(C);self.mt=len(self.T);self.mb=len(self.E);self.m=self.mt+self.mb+len(self.F)
        self.last=None
        self.tcols=(2*self.T[:,:,None]+np.array([0,1])[None,None,:]).reshape(self.mt,6)
        self.ecols=(2*self.E[:,:,None]+np.array([0,1])[None,None,:]).reshape(self.mb,4)
    def ev(self,x):
        if self.last is not None and np.array_equal(x,self.last):return self.gv,self.gj
        C=x[:-1].reshape(-1,2);r=x[-1]
        tg,_,tj=triangles(C,self.T)
        bg,_,bj=boundary(C,self.E,self.sg)
        fg,fj=single_boundary(C,self.F) if len(self.F) else (np.empty(0),np.empty((0,2)))
        g=np.r_[tg,bg,fg];J=np.zeros((self.m+len(C),2*self.n+1))
        J[np.arange(self.mt)[:,None],self.tcols]=-100*tj.reshape(self.mt,6)
        J[self.mt+np.arange(self.mb)[:,None],self.ecols]=-100*bj.reshape(self.mb,4)
        if len(self.F):
            fc=2*self.F[:,None]+np.array([0,1])[None,:]
            J[self.mt+self.mb+np.arange(len(self.F))[:,None],fc]=-100*fj
        J[:self.m,-1]=200*r
        J[self.m+np.arange(len(C))[:,None],2*np.arange(len(C))[:,None]+np.array([0,1])[None,:]]=-2*C
        self.gv=np.r_[100*(r*r-g),1-np.sum(C*C,axis=1)]
        self.gj=J;self.last=x.copy()
        return self.gv,self.gj
    def con(self,x):return self.ev(x)[0]
    def jac(self,x):return self.ev(x)[1]

def refine(C,maxiter=400,rounds=5,trust=.25):
    C=C.copy();best_r,top=covering(C);best=C.copy();details=[]
    objgrad=np.r_[np.zeros(C.size),1.]
    for outer in range(rounds):
        radius,top=covering(C);before=top_key(top)
        model=Model(C,top);x0=np.r_[C.ravel(),radius+1e-10]
        lower=np.r_[np.maximum(C.ravel()-trust,-1),.09]
        upper=np.r_[np.minimum(C.ravel()+trust,1),max(.4,radius+.05)]
        tic=time.monotonic()
        res=minimize(lambda x:x[-1],x0,jac=lambda x:objgrad,method='SLSQP',
          constraints=[dict(type='ineq',fun=model.con,jac=model.jac)],bounds=list(zip(lower,upper)),
          options=dict(maxiter=maxiter,ftol=1e-13,disp=False))
        next_C=res.x[:-1].reshape(-1,2)
        if not np.isfinite(next_C).all():break
        rr,next_top=covering(next_C)
        changed=before!=top_key(next_top)
        details.append(dict(round=outer,initial_radius=radius,actual_radius=rr,objective=float(res.x[-1]),
          iterations=int(res.nit),success=bool(res.success),status=int(res.status),message=str(res.message),
          topology_changed=changed,interior=len(top[0]),boundary=len(top[1]),boundary_singletons=len(top[3]),
          seconds=time.monotonic()-tic,max_norm=float(np.linalg.norm(next_C,axis=1).max()),
          min_frozen_constraint=float(model.con(res.x).min())))
        if rr<best_r and np.linalg.norm(next_C,axis=1).max()<=1+1e-10:
            best_r=rr;best=next_C.copy()
        C=next_C
        if res.success and abs(rr-res.x[-1])<1e-10 and not changed:break
    return best_r,best,details
