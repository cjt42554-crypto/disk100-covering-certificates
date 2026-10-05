# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Contact linear algebra and curvature directions; no nonlinear optimization."""
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[n]='1'
import sys,json,hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import nnls
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent
sys.dont_write_bytecode=True;sys.path.insert(0,str(BASE.parent/'b32_catalog_search/frozen_handoff'))
from graph_tools import extract

def inventory(C,T,E,sg):
    P=C[T];a=P[:,0];u=P[:,1]-a;v=P[:,2]-a;d=u[:,0]*v[:,1]-u[:,1]*v[:,0]
    uu=np.sum(u*u,1);vv=np.sum(v*v,1);pr=np.stack((uu*v[:,1]-vv*u[:,1],u[:,0]*vv-v[:,0]*uu),1)/(2*d[:,None]);p=a+pr
    w1=(pr[:,0]*v[:,1]-pr[:,1]*v[:,0])/d;w2=(u[:,0]*pr[:,1]-u[:,1]*pr[:,0])/d
    w=np.stack((1-w1-w2,w1,w2),1);tj=2*w[:,:,None]*(P-p[:,None,:]);tg=np.sum(pr*pr,1)
    a,b=C[E[:,0]],C[E[:,1]];u=b-a;ss=np.sum(u*u,1);h=(np.sum(b*b,1)-np.sum(a*a,1))/2;z=np.sqrt(ss-h*h)
    up=np.stack((-u[:,1],u[:,0]),1);q=(h[:,None]*u+(sg*z)[:,None]*up)/ss[:,None];tan=np.stack((-q[:,1],q[:,0]),1)
    wi=np.sum(b*tan,1)/np.sum(u*tan,1);w=np.stack((wi,1-wi),1);bj=2*w[:,:,None]*(C[E]-q[:,None,:]);bg=np.sum((q-a)**2,1)
    J=np.zeros((len(T)+len(E),C.size));tc=(2*T[:,:,None]+np.array([0,1])[None,None,:]).reshape(len(T),6);ec=(2*E[:,:,None]+np.array([0,1])[None,None,:]).reshape(len(E),4)
    J[np.arange(len(T))[:,None],tc]=tj.reshape(len(T),6);J[len(T)+np.arange(len(E))[:,None],ec]=bj.reshape(len(E),4)
    keys=[['t',*sorted(map(int,t))] for t in T]+[['b',int(e[0]),int(e[1]),int(s)] for e,s in zip(E,sg)]
    return np.r_[tg,bg],J,keys

def weighted_hessian(C,T,E,sg,active,lam):
    H=np.empty((C.size,C.size));step=2e-6
    for j in range(C.size):
        a=C.copy();b=C.copy();a.flat[j]+=step;b.flat[j]-=step
        H[:,j]=(inventory(a,T,E,sg)[1][active].T@lam-inventory(b,T,E,sg)[1][active].T@lam)/(2*step)
    asym=float(np.max(abs(H-H.T)));return (H+H.T)/2,asym,step

def main():
    source=json.loads((ROOT/'source_full_integer_graph_audit.json').read_text())['records'];rows=[];modes=[];(ROOT/'modes').mkdir(exist_ok=True)
    for s in source:
        i=s['independent_catalog_index'];assert s['integer_actual_graph']['status']=='PASS';p=Path(s['coordinate_file']);assert hashlib.sha256(p.read_bytes()).hexdigest()==s['file_sha256'];C=np.load(p)
        a=extract(C,True);T=np.array(a['faces']);E=np.array(a['boundary_pairs']);sg=np.array(a['boundary_signs']);g,J,keys=inventory(C,T,E,sg);R=float(np.sqrt(g.max()))
        active=np.flatnonzero(g.max()-g<5e-14);A=J[active];lam,res=nnls(np.vstack((A.T,np.ones(len(A)))),np.r_[np.zeros(C.size),1.],maxiter=10000)
        stresserr=float(np.max(abs(A.T@lam)));H,asym,h=weighted_hessian(C,T,E,sg,active,lam)
        rot=np.column_stack((-C[:,1],C[:,0])).ravel();rot/=np.linalg.norm(rot)
        N=null_space(np.vstack((A,rot)),rcond=1e-10);evals,evecs=np.linalg.eigh(N.T@H@N)
        weak=np.flatnonzero(lam<1e-7);strong=np.flatnonzero(lam>=1e-7)
        r=dict(index=i,source_file=str(p),source_file_sha256=s['file_sha256'],source_radius_numeric=R,active_squared_tolerance=5e-14,active_count=len(active),active_keys=[keys[j] for j in active],NNLS_lambda=lam.tolist(),NNLS_residual_l2=float(res),NNLS_stress_residual_inf=stresserr,
            all_positive_contact_equality_null_dimension=N.shape[1],lowest_reduced_curvature_eigenvalues=evals[:10].tolist(),hessian_central_difference_step=h,hessian_asymmetry_before_symmetrizing=asym,
            weak_split_threshold_for_heuristic_model_only=1e-7,nearzero_lambda_count=len(weak),nearzero_lambda_values=lam[weak].tolist(),strong_lambda_min=float(lam[strong].min()),
            nearzero_lambda_is_exact_zero=False,scope='Numerical contact differential/eigen analysis. NNLS and tolerance classify a heuristic split only; all observed positive multipliers retain equality in the primary model. No exact critical cone, stability or descent claim.')
        for mode in range(min(2,len(evals))):
            V=(N@evecs[:,mode]).reshape(C.shape);V/=np.max(np.linalg.norm(V,axis=1))
            for sign in (1,-1):
                v=sign*V;path=ROOT/'modes'/f'index{i:02d}_equal_mode{mode:02d}_sign{sign:+d}.npz';np.savez(path,C=C,V=v)
                modes.append(dict(index=i,mode=mode,sign=sign,file=str(path),file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),model='All numerically active contacts as equalities plus rotation gauge',max_active_squaredradius_directional_derivative=float(np.max(A@v.ravel())),max_absolute_active_squaredradius_directional_derivative=float(np.max(abs(A@v.ravel()))),curvature_eigenvalue=float(evals[mode]),largest_velocity_labels=np.argsort(np.linalg.norm(v,axis=1))[-20:][::-1].tolist()))
        if i==14:
            weaksites=sorted({j for k in weak for j in keys[active[k]][1:4] if j>=0})
            adj={j:set() for j in range(100)}
            for f in T:
                for j in f:adj[int(j)].update(map(int,f))
            patches=[weaksites,sorted(set(weaksites)|{q for j in weaksites for q in adj[j]})]
            r['weak_local_patches']=[]
            for patch in patches:
                cols=np.array([2*j+k for j in patch for k in (0,1)]);Ns=null_space(A[strong][:,cols],rcond=1e-10);Na=null_space(A[:,cols],rcond=1e-10)
                r['weak_local_patches'].append(dict(labels=patch,strong_contact_heuristic_null_dimension=Ns.shape[1],all_positive_contact_equality_null_dimension=Na.shape[1]))
                if Ns.shape[1]:
                    vals,vecs=np.linalg.eigh(Ns.T@H[np.ix_(cols,cols)]@Ns)
                    local=Ns@vecs[:,0];v=np.zeros(C.size);v[cols]=local;v=v.reshape(C.shape);v/=np.max(np.linalg.norm(v,axis=1))
                    for sign in (1,-1):
                        V=sign*v;path=ROOT/'modes'/f'index14_weakpatch{len(patch):02d}_sign{sign:+d}.npz';np.savez(path,C=C,V=V)
                        wd=A[weak]@V.ravel()
                        modes.append(dict(index=14,file=str(path),file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),model='Heuristic strong equality with every nearzero positive multiplier contact retained as tested inequality; not an exact zero multiplier assertion',patch_labels=patch,max_strong_squaredradius_directional_derivative=float(np.max(abs(A[strong]@V.ravel()))),max_weak_squaredradius_directional_derivative=float(wd.max()),weak_inequalities_numerically_nonpositive=bool(wd.max()<=1e-11),curvature_eigenvalue=float(vals[0])))
        rows.append(r);print('index',i,'active',len(active),'null',N.shape[1],'curvature',evals[:4],'weak',len(weak),flush=True)
    (ROOT/'continuation_linear_algebra.json').write_text(json.dumps(dict(status='COMPLETE_NUMERICAL_LINEAR_ALGEBRA',records=rows,modes=modes,no_covering_or_proxy_optimizer=True),indent=2)+'\n')

if __name__=='__main__':main()
