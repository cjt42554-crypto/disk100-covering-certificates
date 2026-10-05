# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Paired transport engine: full100 nooffset epigraph; dynamic whole-disk scoring."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import numpy as np,json,hashlib,time,traceback,argparse,multiprocessing as mp,itertools
from scipy.linalg import null_space
from scipy.optimize import nnls,minimize
import geometry as geom
from whole_disk_oracle import true_cover
ROOT=Path(__file__).resolve().parent.parent;BASE=ROOT.parent
REFERENCE=.11544695488558258
THRESHOLD=.1154469548855
MATERIAL=1e-10
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def csha(C):return hashlib.sha256(np.asarray(C,dtype='<f8').tobytes()).hexdigest()
def clean(v):
    if isinstance(v,np.ndarray):return clean(v.tolist())
    if isinstance(v,(np.integer,np.bool_)):return v.item()
    if isinstance(v,(np.floating,float)):return float(v) if np.isfinite(v) else None
    if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    return v
def write(p,value):
    p=Path(p);q=Path(str(p)+'.tmp');q.write_text(json.dumps(clean(value),ensure_ascii=False,indent=2,allow_nan=False)+'\n');q.replace(p)
def inventory(C,top=None):
    if top is None:top=geom.topology(C)
    T,E,sg,F=top;tg,P,tj=geom.triangles(C,T);bg,Q,bj=geom.boundary(C,E,sg);fg,fj=geom.single_boundary(C,F) if len(F) else (np.empty(0),np.empty((0,2)))
    J=np.zeros((len(T)+len(E)+len(F),C.size),dtype=np.result_type(C));tc=(2*T[:,:,None]+np.array([0,1])[None,None,:]).reshape(len(T),6);ec=(2*E[:,:,None]+np.array([0,1])[None,None,:]).reshape(len(E),4)
    J[np.arange(len(T))[:,None],tc]=tj.reshape(len(T),6);J[len(T)+np.arange(len(E))[:,None],ec]=bj.reshape(len(E),4)
    if len(F):J[len(T)+len(E)+np.arange(len(F))[:,None],2*F[:,None]+np.array([0,1])[None,:]]=fj
    keys=[['t',*sorted(map(int,t))] for t in T]+[['b',int(e[0]),int(e[1]),int(s)] for e,s in zip(E,sg)]+[['f',int(i)] for i in F]
    sites=[list(map(int,t)) for t in T]+[list(map(int,e)) for e in E]+[[int(i)] for i in F]
    points=np.r_[P,Q,-C[F]/np.linalg.norm(C[F],axis=1)[:,None]]
    return np.r_[tg,bg,fg],points,J,keys,sites,top
def signature(keys):return hashlib.sha256(json.dumps(sorted(keys),separators=(',',':')).encode()).hexdigest()
def diagnose(C):
    g,p,J,kk,sites,top=inventory(C);active=np.flatnonzero(g.max()-g<1e-9);used=sorted({i for j in active for i in sites[j]});inactive=[i for i in range(len(g)) if i not in set(active)]
    return dict(active_count=len(active),binding_count=len(used),binding_indices=used,nonbinding_indices=sorted(set(range(len(C)))-set(used)),active_keys=[kk[i] for i in active],inactive_keys=[kk[i] for i in inactive],active_signature=signature([kk[i] for i in active]),inactive_signature=signature([kk[i] for i in inactive]),interior_count=len(top[0]),boundary_count=len(top[1]),singleton_count=len(top[3]),active_tolerance_squared=1e-9)
def known99(C):
    K=np.load(ROOT/'inputs/known99_reference.npy');ref=diagnose(K);dd=diagnose(C);used=ref['binding_indices'];ang=np.arctan2(np.sum(C[used,0]*K[used,1]-C[used,1]*K[used,0]),np.sum(C[used]*K[used]));rot=np.array([[np.cos(ang),-np.sin(ang)],[np.sin(ang),np.cos(ang)]]);dist=np.linalg.norm(C-K@rot,axis=1);R,detail=true_cover(C)
    match=dist[used].max()<2e-4 and dd['binding_indices']==used and dd['active_signature']==ref['active_signature'] and abs(R-.11551610349101131)<1e-8
    return dict(match=bool(match),binding_skeleton_max_distance=float(dist[used].max()),all_center_max_distance=float(dist.max()),active_keys_match=dd['active_signature']==ref['active_signature'],inactive_keys_match=dd['inactive_signature']==ref['inactive_signature'],whole_disk_radius=R,free_indices=dd['nonbinding_indices'],scope='Alignment of actual99 binding centers plus full active-key identity and whole-disk oracle; never radius-only detection.')
def optimize(C,folder,row):
    best=C.copy();bestR,initial_detail=true_cover(C);rounds=[];objgrad=np.r_[np.zeros(C.size),1.];index=row['grid_id']
    for outer in range(2):
        radius,top=geom.covering(C);before=geom.top_key(top);model=geom.Model(C,top);x0=np.r_[C.ravel(),radius+1e-10];trust=.18;lower=np.r_[np.maximum(C.ravel()-trust,-1),.09];upper=np.r_[np.minimum(C.ravel()+trust,1),max(.4,radius+.05)];tic=time.monotonic()
        res=minimize(lambda x:x[-1],x0,jac=lambda x:objgrad,method='SLSQP',constraints=[dict(type='ineq',fun=model.con,jac=model.jac)],bounds=list(zip(lower,upper)),options=dict(maxiter=450,ftol=1e-13,disp=False));nextC=res.x[:-1].reshape(-1,2);path=folder/'stages'/f'{index:04d}_{outer:02d}.npy';np.save(path,nextC)
        if not np.isfinite(nextC).all():raise ValueError('nonfinite_coordinates')
        R,detail=true_cover(nextC);rr,nexttop=geom.covering(nextC);changed=before!=geom.top_key(nexttop);mincon=float(model.con(res.x).min());feasible=abs(R-res.x[-1])<1e-10 and mincon>-2e-8 and np.linalg.norm(nextC,axis=1).max()<=1+1e-10
        oldkeys=inventory(C,top)[3];newg,newp,newJ,newkeys,newsites,newtop=inventory(nextC,nexttop);newactive=np.flatnonzero(newg.max()-newg<1e-9);missing=sorted(set(map(tuple,newkeys))-set(map(tuple,oldkeys)));removed=sorted(set(map(tuple,oldkeys))-set(map(tuple,newkeys)));newbinding=sorted(set(tuple(newkeys[j]) for j in newactive)-set(map(tuple,oldkeys)))
        stage=dict(round=outer,coordinate_file=str(path.relative_to(ROOT)),file_sha256=sha(path),coordinate_sha256=csha(nextC),radius=R,whole_disk_score=detail,model_radius=rr,objective=float(res.x[-1]),success=bool(res.success),status=int(res.status),message=str(res.message),iterations=int(res.nit),topology_changed=changed,new_actual_keys=missing,removed_actual_keys=removed,new_actual_binding_keys=newbinding,true_radius_minus_frozen_objective=float(R-res.x[-1]),minimum_frozen_constraint=mincon,feasible=bool(feasible),max_center_norm=float(np.linalg.norm(nextC,axis=1).max()),seconds=time.monotonic()-tic,offsets=0)
        if R<bestR and stage['max_center_norm']<=1+1e-10:bestR=R;best=nextC.copy()
        rounds.append(stage);C=nextC
        if R<THRESHOLD and REFERENCE-R>MATERIAL:
            q=folder/'proposals'/f'grid{index:04d}_round{outer:02d}_R{R:.17f}.npy';np.save(q,C);write(q.with_suffix('.json'),dict(row,stage=stage,proposal_file_sha256=sha(q),proposal_coordinate_sha256=csha(C),scope='Numerical new proposal, exact integer/surd certification required.'));print('POTENTIAL_IMPROVEMENT',row['worker'],index,R,str(q),flush=True)
        if feasible and R>.11551 and diagnose(C)['binding_count']==99:
            kb=known99(C);stage['known99']=kb
            if kb['match']:stage['second_round_skipped_known99_family']=True;break
        if feasible and res.success and not changed:break
    return bestR,best,initial_detail,rounds
