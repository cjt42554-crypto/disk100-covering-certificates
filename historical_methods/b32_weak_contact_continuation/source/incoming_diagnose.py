# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Audit raw rank conventions and finite-scale incoming-contact cone."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import numpy as np
from scipy.linalg import null_space,eigh
from scipy.optimize import minimize,linprog
import branches
from diagnose import ROOT,write,sha

def main():
    saved=np.load(ROOT/'matrices/analysis.npz');C=saved['C'];J=saved['J'];lam=saved['weights'];H=saved['H'];gauge=saved['gauge'];strong=saved['strong'];weak=saved['weak']
    rank=[]
    for name,A,weights in [('all_positive',J,lam),('strong',J[strong],lam[strong])]:
        U,s,Vh=np.linalg.svd(A,full_matrices=True);row=[]
        for rc in (1e-10,1e-9,1e-8,1e-7,1e-6):
            K=np.vstack((A,gauge));N=null_space(K,rcond=rc);eig=eigh(N.T@H@N,eigvals_only=True)
            row.append(dict(rcond=rc,augmented_largest_singular_value=np.linalg.svd(K,compute_uv=False)[0],tangent_dimension=N.shape[1],minimum_curvature=eig[0],null_residual_max=np.max(np.abs(A@N)),FD_h1e6_operator_error=np.linalg.norm(N.T@(saved['H_fd_1e6']-H)@N,2)))
        Nold=null_space(np.vstack((A,gauge)),rcond=1e-9);Nnew=null_space(np.vstack((A,gauge)),rcond=1e-8)
        P=Nnew-Nold@(Nold.T@Nnew);u,ss,v=np.linalg.svd(P,full_matrices=False);extra=u[:,0]
        np.save(ROOT/'matrices'/f'{name}_additional_near_null.npy',extra.reshape(100,2))
        epi=np.column_stack((A,-np.ones(len(A))));aug=np.vstack((epi,np.r_[gauge,0],np.r_[np.zeros(200),1]))
        epiN=null_space(aug,rcond=1e-9)
        rank.append(dict(name=name,spatial_singular_values=s,stationarity_sigma_min_upper_bound=np.linalg.norm(A.T@weights)/np.linalg.norm(weights),stationarity_residual_norm=np.linalg.norm(A.T@weights),lambda_norm=np.linalg.norm(weights),rank_tolerance_sweep=row,additional_near_null_curvature=extra@H@extra,additional_near_null_Aresidual=np.max(np.abs(A@extra)),additional_near_null_rotation_residual=gauge@extra,additional_subspace_singular_values=ss,epigraph_plus_objective_plus_gauge_dimension=epiN.shape[1],epigraph_objective_null_residual=np.max(np.abs(aug@epiN)),epigraph_min_curvature=eigh(epiN[:200].T@H@epiN[:200],eigvals_only=True)[0]))
    g,p,Q,keys,top=branches.inventory(C);slacks=g.max()-g;near=np.flatnonzero((slacks>5e-14)&(slacks<1e-8))
    inequalities=np.r_[J[weak],Q[near]];labels=[branches.Branches(C).keys[i] for i in weak]+[keys[i] for i in near]
    N=null_space(np.vstack((J[strong],gauge)),rcond=1e-8);Hs=N.T@H@N;W=inequalities@N
    u,s,vh=np.linalg.svd(W,full_matrices=True);wrank=int(np.sum(s>s[0]*1e-7));Z=vh[wrank:].T;Y=vh[:wrank].T
    Hzz=Z.T@Hs@Z;relax=Y-Z@np.linalg.solve(Hzz,Z.T@Hs@Y);S=relax.T@Hs@relax;M=relax.T@relax;A=W@relax
    ev,V=eigh(S,M);trials=[];old_traces=[]
    for fn in sorted((ROOT/'matrices').glob('cone*.npy')):
        d=np.load(fn).ravel();old_traces.append(dict(name=fn.stem,source_sha256=sha(fn),incoming_derivatives=Q[near]@d,maximum_incoming_derivative=np.max(Q[near]@d),near_slack_linear_first_collision=min([slacks[i]/v for i,v in zip(near,Q[near]@d) if v>1e-14],default=None),qualified_as_actual_cover_tangent=False))
    # Coordinate LP witnesses determine whether the conservative quotient cone
    # contains nonzero directions; finite eigen-ray QPs rank costs, not global minima.
    lps=[]
    for k in range(wrank):
        for sign in (-1,1):
            obj=np.zeros(wrank);obj[k]=sign
            r=linprog(obj,A_ub=A,b_ub=np.zeros(len(A)),bounds=[(-1,1)]*wrank,method='highs')
            lps.append(dict(coordinate=k,sign=sign,status=int(r.status),maximum_coordinate=float(-r.fun) if r.success else None,solution=r.x if r.success else None))
    for k in range(min(6,wrank)):
        for sign in (-1,1):
            y0=sign*V[:,k];y0/=np.sqrt(y0@M@y0)
            r=minimize(lambda y:float(y@S@y),y0,jac=lambda y:2*S@y,method='SLSQP',constraints=[dict(type='eq',fun=lambda y:y@M@y-1,jac=lambda y:2*M@y),dict(type='ineq',fun=lambda y:-A@y,jac=lambda y:-A)],options=dict(maxiter=500,ftol=1e-13,disp=False))
            d=N@relax@r.x;d/=np.linalg.norm(d);f=ROOT/'matrices'/f'incoming_cone{k:02d}_{"plus" if sign>0 else "minus"}.npy';np.save(f,d.reshape(100,2))
            trueineq=inequalities@d;slmax=np.max(np.abs(J[strong]@d));feasible=bool(trueineq.max()<1e-8 and slmax<1e-8)
            mode_matches=[]
            for old in sorted((ROOT/'inputs/owner_modes/modes').glob('*.npz')):
                z=np.load(old);D=z['V'].ravel();D/=np.linalg.norm(D);mode_matches.append(dict(mode=old.name,cosine=float(d@D),mode_file_sha256=sha(old)))
            trials.append(dict(name=f.stem,status=int(r.status),success=bool(r.success),message=str(r.message),iterations=int(r.nit),curvature=float(d@H@d),constraint_derivatives=trueineq,maximum_inequality_derivative=trueineq.max(),strong_derivative_max=slmax,rotation_residual=gauge@d,finite_conservative_model_feasible=feasible,direction_file=str(f.relative_to(ROOT)),direction_sha256=sha(f),maximum_absolute_cosine_to_owner16=max(abs(v['cosine']) for v in mode_matches),owner16_comparisons=mode_matches))
    np.savez(ROOT/'matrices/incoming_analysis.npz',C=C,H=H,J=J,lambda_values=lam,strong=strong,weak=weak,near=near,near_J=Q[near],near_slacks=slacks[near],inequality_J=inequalities,Nstrong=N,release_relax=relax,release_S=S,release_metric=M,release_A=A)
    report=dict(input_sha256=sha(ROOT/'inputs/index14_grid0061.npy'),numerical_only=True,raw_pose_not_exact_KKT=True,rank_review=rank,prior_146_model_directions=old_traces,near_inactive_cutoff_squared=1e-8,near_inactive_count=len(near),near_inactive_keys=[keys[i] for i in near],near_inactive_slacks=slacks[near],combined_inequality_keys=labels,combined_map_singular_values=s,combined_quotient_rank=wrank,strong_complete_dimension=N.shape[1],strong_complete_Hmin=eigh(Hs,eigvals_only=True)[0],kernel_dimension=Z.shape[1],kernel_Hmin=eigh(Hzz,eigvals_only=True)[0],quotient_eigenvalues=ev,LP_nonzero_cone_witnesses=lps,finite_cone_trials=trials,matrix_sha256=sha(ROOT/'matrices/incoming_analysis.npz'),scope='The 8 positive-slack contacts are incoming finite-scale constraints, not exact active constraints. Positive tiny lambdas are not exact zero. Augmented raw-pose rank is explicitly tolerance dependent; near-null extra mode is retained. Original12 directions are unqualified when incoming contacts grow. No covering optimizer ran.')
    write(ROOT/'INCOMING_CONTACT_REVIEW.json',report)
    print('COMPLETE_DIM',N.shape[1],'near',len(near),'rank',wrank,'Schurmin',ev[0],flush=True)
    print('CONES',[(x['name'],x['status'],x['curvature'],x['maximum_inequality_derivative'],x['finite_conservative_model_feasible'],x['maximum_absolute_cosine_to_owner16']) for x in trials],flush=True)
if __name__=='__main__':main()
