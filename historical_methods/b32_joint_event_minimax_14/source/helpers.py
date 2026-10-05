# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
from pathlib import Path
import json,hashlib,numpy as np
import geometry

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def csha(C):return hashlib.sha256(np.asarray(C,dtype='<f8').tobytes()).hexdigest()

def clean(v):
    if isinstance(v,np.ndarray):return clean(v.tolist())
    if isinstance(v,(np.integer,np.bool_)):return v.item()
    if isinstance(v,(np.floating,float)):return float(v) if np.isfinite(v) else None
    if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    return v

def write(p,v):
    p=Path(p);q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(clean(v),indent=2,allow_nan=False)+'\n');q.replace(p)

def inventory(C,top=None):
    T,E,sg,F=geometry.topology(C) if top is None else top
    tg,p,tj=geometry.triangles(C,T);bg,q,bj=geometry.boundary(C,E,sg);fg,fj=geometry.single_boundary(C,F) if len(F) else (np.empty(0),np.empty((0,2)))
    J=np.zeros((len(T)+len(E)+len(F),C.size))
    for k,(ids,gradient) in enumerate(zip(T,tj)):J[k,(2*ids[:,None]+[0,1]).ravel()]=gradient.ravel()
    for k,(ids,gradient) in enumerate(zip(E,bj),len(T)):J[k,(2*ids[:,None]+[0,1]).ravel()]=gradient.ravel()
    for k,(i,gradient) in enumerate(zip(F,fj),len(T)+len(E)):J[k,2*i:2*i+2]=gradient
    keys=[['t',*sorted(map(int,t))] for t in T]+[['b',int(e[0]),int(e[1]),int(s)] for e,s in zip(E,sg)]+[['f',int(i)] for i in F]
    return np.r_[tg,bg,fg],J,keys,(T,E,sg,F)
