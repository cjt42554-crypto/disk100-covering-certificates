# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Actual smooth Voronoi-hole branches, t=R^2; no artificial anchors."""
import numpy as np
import geometry as geom

def inventory(C, top=None):
    T,E,sg,F=geom.topology(C) if top is None else top
    assert not len(F), 'singleton branch needs a separate Hessian implementation'
    tg,p,tj=geom.triangles(C,T); bg,q,bj=geom.boundary(C,E,sg)
    J=np.zeros((len(T)+len(E),C.size))
    for k,(ids,grad) in enumerate(zip(T,tj)):
        J[k,(2*ids[:,None]+[0,1]).ravel()]=grad.ravel()
    for k,(ids,grad) in enumerate(zip(E,bj),len(T)):
        J[k,(2*ids[:,None]+[0,1]).ravel()]=grad.ravel()
    keys=[['t',*sorted(map(int,t))] for t in T]+[['b',int(e[0]),int(e[1]),int(s)] for e,s in zip(E,sg)]
    return np.r_[tg,bg],np.r_[p,q],J,keys,(T,E,sg,F)

class Branches:
    def __init__(self,C,active=None):
        g,p,J,keys,top=inventory(C)
        active=np.flatnonzero(g.max()-g<=5e-14) if active is None else np.asarray(active)
        T,E,sg,F=top;ti=active[active<len(T)];bi=active[active>=len(T)]-len(T)
        self.T,self.E,self.sg=T[ti],E[bi],sg[bi]
        self.keys=[keys[i] for i in active];self.indices=active;self.m=len(active)

    def evaluate(self,C):
        tg,p,tj=geom.triangles(C,self.T);bg,q,bj=geom.boundary(C,self.E,self.sg)
        J=np.zeros((self.m,C.size))
        for k,(ids,grad) in enumerate(zip(self.T,tj)):
            J[k,(2*ids[:,None]+[0,1]).ravel()]=grad.ravel()
        for k,(ids,grad) in enumerate(zip(self.E,bj),len(self.T)):
            J[k,(2*ids[:,None]+[0,1]).ravel()]=grad.ravel()
        return np.r_[tg,bg],np.r_[p,q],J

    def hessian(self,C,weights):
        _,points,_=self.evaluate(C);H=np.zeros((C.size,C.size));local=[]
        for k,(p,lam) in enumerate(zip(points,weights)):
            inside=k<len(self.T);ids=self.T[k] if inside else self.E[k-len(self.T)]
            P=C[ids];n=len(ids);rhs=np.zeros((2,2*n))
            if inside:
                B=P[1:]-P[0]
                for j in (1,2):
                    rhs[j-1,:2]=-(P[0]-p);rhs[j-1,2*j:2*j+2]=P[j]-p
                w=np.linalg.solve(np.vstack((P.T,np.ones(3))),np.r_[p,1]);beta=0.
            else:
                B=np.vstack((P[1]-P[0],p));rhs[0,:2]=-(P[0]-p);rhs[0,2:4]=P[1]-p
                tangent=np.array([-p[1],p[0]]);w0=P[1]@tangent/((P[1]-P[0])@tangent)
                w=np.array([w0,1-w0]);beta=-lam*(1-(w@P)@p)
            Q=np.linalg.solve(B,rhs);hk=2*beta*Q.T@Q
            for j in range(n):
                D=Q.copy();D[:,2*j:2*j+2]-=np.eye(2);hk+=2*lam*w[j]*D.T@D
            cols=(2*ids[:,None]+[0,1]).ravel();H[np.ix_(cols,cols)]+=hk
            local.append(dict(key=self.keys[k],multiplier=float(lam),barycentric=w.tolist(),unit_beta=float(beta),IFT_det=float(np.linalg.det(B))))
        return H,local
