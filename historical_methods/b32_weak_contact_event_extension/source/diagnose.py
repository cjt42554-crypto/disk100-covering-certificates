# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Fresh weak-14 cone/Hessian diagnosis, numerical only, no whole optimizer."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import null_space,eigh
from scipy.optimize import nnls,minimize
import branches
from whole_disk_oracle import true_cover
ROOT=Path(__file__).resolve().parent.parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def clean(v):
    if isinstance(v,np.ndarray):return clean(v.tolist())
    if isinstance(v,(np.integer,np.bool_)):return v.item()
    if isinstance(v,(np.floating,float)):return float(v) if np.isfinite(v) else None
    if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    return v
def write(p,v):Path(p).write_text(json.dumps(clean(v),indent=2,allow_nan=False)+'\n')
def main():
    path=ROOT/'inputs/index14_grid0061.npy';C=np.load(path)
    assert sha(path)=='78d5a0ae08b07010da7c52e9634d1a14ecdcf12cf0e04ff5d1895895231452c0'
    b=branches.Branches(C);g,p,J=b.evaluate(C)
    lam,res=nnls(np.vstack((J.T,np.ones(b.m))),np.r_[np.zeros(200),1.],maxiter=10000);lam/=lam.sum()
    H,local=b.hessian(C,lam);gauge=np.column_stack((-C[:,1],C[:,0])).ravel();gauge/=np.linalg.norm(gauge)
    weak=np.flatnonzero(lam<1e-7);strong=np.flatnonzero(lam>=1e-7)
    Nall=null_space(np.vstack((J,gauge)),rcond=1e-9);Nstrong=null_space(np.vstack((J[strong],gauge)),rcond=1e-9)
    Hall=Nall.T@H@Nall;Hs=Nstrong.T@H@Nstrong;evall=eigh(Hall,eigvals_only=True);es,Vs=eigh(Hs)
    W=J[weak]@Nstrong
    # Ker(W) equals all-active tangent; isolate the actual weak release quotient.
    U,sv,Vh=np.linalg.svd(W,full_matrices=True);wrank=int(np.sum(sv>sv[0]*1e-9))
    Z=null_space(W,rcond=1e-9);Y=Vh[:wrank].T
    Hzz=Z.T@Hs@Z;Hyy=Y.T@Hs@Y;Hzy=Z.T@Hs@Y
    hzz_e=eigh(Hzz,eigvals_only=True)
    if hzz_e[0]>1e-10:
        relax=Y-Z@np.linalg.solve(Hzz,Hzy)
        S=relax.T@Hs@relax;A=W@relax
        metric=relax.T@relax;se,sevec=eigh(S,metric)
    else:
        relax=Y;S=Hyy;A=W@Y;metric=np.eye(wrank);se,sevec=eigh(S)
    fd=[];fdmat=[]
    for h in (1e-5,1e-6,1e-7):
        M=np.empty_like(H)
        for k in range(200):
            Cp=C.copy();Cm=C.copy();Cp.ravel()[k]+=h;Cm.ravel()[k]-=h
            M[:,k]=(b.evaluate(Cp)[2].T@lam-b.evaluate(Cm)[2].T@lam)/(2*h)
        asym=np.max(np.abs(M-M.T));M=(M+M.T)/2;fdmat.append(M)
        fd.append(dict(h=h,asymmetry=asym,full_entry_error=np.max(np.abs(M-H)),strong_reduced_operator_error=np.linalg.norm(Nstrong.T@(M-H)@Nstrong,2),all_active_operator_error=np.linalg.norm(Nall.T@(M-H)@Nall,2)))
    allg,allp,allJ,allkeys,top=branches.inventory(C);R,detail=true_cover(C)
    cones=[]
    # Deterministic eigen-directions of a new full-cover Hessian Schur quotient.
    # This is a finite heuristic cone optimization, not global cone certification.
    for k in range(min(6,wrank)):
        for sign in (-1,1):
            y0=sign*sevec[:,k];y0/=np.sqrt(y0@metric@y0)
            result=minimize(lambda y:float(y@S@y),y0,jac=lambda y:2*S@y,method='SLSQP',constraints=[dict(type='eq',fun=lambda y:y@metric@y-1,jac=lambda y:2*metric@y),dict(type='ineq',fun=lambda y:-A@y,jac=lambda y:-A)],options=dict(maxiter=400,ftol=1e-13,disp=False))
            y=result.x;d=Nstrong@relax@y;norm=np.linalg.norm(d);d/=norm
            feasible=bool(np.max(J[weak]@d)<1e-8 and np.max(np.abs(J[strong]@d))<1e-8)
            key=f'cone{k:02d}_{"plus" if sign>0 else "minus"}'
            np.save(ROOT/'matrices'/f'{key}.npy',d.reshape(100,2))
            cones.append(dict(name=key,status=int(result.status),success=bool(result.success),message=str(result.message),iterations=int(result.nit),curvature=float(d@H@d),maximum_weak_derivative=float(np.max(J[weak]@d)),weak_derivatives=J[weak]@d,strong_derivative_max=np.max(np.abs(J[strong]@d)),rotation_gauge_residual=float(gauge@d),feasible_numeric=feasible,direction_file=f'matrices/{key}.npy',direction_sha256=sha(ROOT/'matrices'/f'{key}.npy')))
    np.savez(ROOT/'matrices/analysis.npz',C=C,J=J,weights=lam,H=H,gauge=gauge,strong=strong,weak=weak,N_all=Nall,N_strong=Nstrong,weak_map=W,release_relax=relax,release_S=S,release_metric=metric,release_A=A,H_fd_1e5=fdmat[0],H_fd_1e6=fdmat[1],H_fd_1e7=fdmat[2])
    report=dict(input_file=str(path.relative_to(ROOT)),input_sha256=sha(path),numerical_only=True,objective='t=R^2',whole_radius=R,whole_disk=detail,active_tolerance_squared=5e-14,active_count=b.m,active_keys=b.keys,multipliers=lam,stationarity_inf=np.max(np.abs(J.T@lam)),stationarity_norm=np.linalg.norm(J.T@lam),active_equality_range=np.ptp(g),tiny_lambda_cutoff_heuristic=1e-7,weak_count=len(weak),strong_count=len(strong),weak_keys=[b.keys[i] for i in weak],weak_positive_multipliers=lam[weak],all_positive_exact_critical_cone_requires_all_equalities=True,all_active_dimension=Nall.shape[1],strong_tangent_dimension=Nstrong.shape[1],weak_map_rank=wrank,all_active_eigenvalues=evall,strong_eigenvalues=es,all_active_minimum=evall[0],strong_minimum=es[0],release_generalized_eigenvalues=se,release_generalized_minimum=se[0],fd_checks=fd,contacts=local,cone_trials=cones,inventory_keys=allkeys,squared_slacks=allg.max()-allg,matrix_sha256=sha(ROOT/'matrices/analysis.npz'),scope='Tiny positive multipliers are kept; weak-inequality relaxation is a numerical exploration model. Finite cone trials and curvature are not a strict SOSC/global proof. No full100 covering optimizer was run.')
    write(ROOT/'DIAGNOSTIC.json',report)
    print(json.dumps(clean({k:report[k] for k in ['whole_radius','active_count','weak_count','strong_count','all_active_dimension','strong_tangent_dimension','weak_map_rank','all_active_minimum','strong_minimum','release_generalized_minimum','stationarity_inf','fd_checks']}),indent=2),flush=True)
    print('CONE_TRIALS',[(r['name'],r['status'],r['curvature'],r['feasible_numeric']) for r in cones],flush=True)
if __name__=='__main__':main()
