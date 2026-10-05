# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Independent all-triple acceptance, active inventory and numerical deletion debt."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.optimize import nnls,linear_sum_assignment
from concurrent.futures import ProcessPoolExecutor,as_completed
from independent_oracle import all_triple_cover
from whole_disk_oracle import true_cover
from proxy_driver import inventory,write,sha,csha
import graph_tools
ROOT=Path(__file__).resolve().parent.parent
def replay(item):
    path,expected=item;assert sha(ROOT/path)==expected;C=np.load(ROOT/path);R,details=all_triple_cover(C);rr,dd=true_cover(C)
    return dict(file=path,file_sha256=expected,coordinate_sha256=csha(C),all_triple_radius=R,dynamic_radius=rr,absolute_difference=abs(R-rr),passed=abs(R-rr)<2e-11,independent_oracle=details)
def inventory_diagnostic(C):
    R,detail=true_cover(C);g,J,keys,top=inventory(C);slack=R*R-g;ids=np.flatnonzero(slack<=5e-14);near=np.flatnonzero((slack>5e-14)&(slack<=1e-8));A=J[ids].T
    lam,res=nnls(np.vstack((A,np.ones((1,len(ids))))),np.r_[np.zeros(200),1.],maxiter=10000);used=set()
    for j in ids:used.update(keys[j][1:4] if keys[j][0]=='t' else keys[j][1:3] if keys[j][0]=='b' else keys[j][1:2])
    free=sorted(set(range(100))-used);deleted=[]
    for i in free:
        D=np.delete(C,i,axis=0);rd,details=all_triple_cover(D);dynamic,dd=true_cover(D)
        deleted.append(dict(removed_center_index=i,deleted99_all_triple_radius=rd,deleted99_dynamic_radius=dynamic,oracle_difference=abs(rd-dynamic),delete_cavity_radius_debt=rd-R,original_full100_radius=R,scope='Numerical deletion diagnostic. Nonbinding center still covers inactive holes; not permission to delete.',independent_oracle=details))
    return dict(radius=R,holes=len(keys),active_squared_slack_threshold=5e-14,active_count=len(ids),near_inactive_squared_slack_cutoff=1e-8,near_inactive_count=len(near),binding_center_count=len(used),nonbinding_center_indices=free,all_keys=keys,all_squared_slacks=slack.tolist(),active_indices=ids.tolist(),near_inactive_indices=near.tolist(),active_nnls_lambda=lam.tolist(),lambda_minimum=float(lam.min()),positive_lambda_count=int(np.sum(lam>1e-8)),nnls_augmented_stationarity_residual=float(res),spatial_stationarity_inf=float(np.max(np.abs(A@lam))),lambda_sum=float(lam.sum()),no_exact_stationarity_or_sosc_claim=True,deletion_debts=deleted,actual_graph=graph_tools.extract(C,True))
def alignment(P,Q):
    best=None
    # No translation: the unit disk center is fixed. Hungarian/Kabsch runs are
    # only numerical matching; the radial bound is a necessary lower bound.
    for reflect in (False,True):
        A=P.copy();A[:,1]*=-1 if reflect else 1
        for theta in np.linspace(0,2*np.pi,16,endpoint=False):
            O=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
            previous=None
            for k in range(12):
                dist=np.sum((A[:,None,:]@O-Q[None,:,:])**2,axis=2);ri,ci=linear_sum_assignment(dist)
                U,s,V=np.linalg.svd(A.T@Q[ci]);S=np.eye(2);S[1,1]=np.linalg.det(U@V);O=U@S@V
                if previous is not None and np.array_equal(previous,ci):break
                previous=ci.copy()
            error=np.linalg.norm(A@O-Q[ci],axis=1);row=dict(maximum_all100_center_distance=float(error.max()),rms_all100_center_distance=float(np.sqrt(np.mean(error**2))),reflection=reflect,rotation=O.tolist(),target_assignment=ci.tolist())
            if best is None or row['rms_all100_center_distance']<best['rms_all100_center_distance']:best=row
    best['all100_sorted_radial_assignment_lower_bound']=float(np.max(np.abs(np.sort(np.linalg.norm(P,axis=1))-np.sort(np.linalg.norm(Q,axis=1)))))
    return best
def main():
    audit=ROOT/'audit_original';audit.mkdir(exist_ok=True);aggregate=json.loads((ROOT/'ORIGINAL_RESTORE_RESULT.json').read_text());refs={};records=[]
    for r in aggregate['results']:
        p=ROOT/r['result_file'];assert sha(p)==r['result_sha256'];d=json.loads(p.read_text());refs[d['start_file']]=d['start_file_sha256'];refs[d['final_file']]=d['final_file_sha256']
        for s in d['stages']:refs[s['file']]=s['file_sha256']
        C=np.load(ROOT/d['final_file']);diagnostic=inventory_diagnostic(C);records.append(dict(restore_id=r['restore_id'],final_file=d['final_file'],final_sha256=d['final_file_sha256'],final_diagnostic=diagnostic,stages_status=[s['status'] for s in d['stages']],last_solver_success=d['stages'][-1]['success'],optimizer_exception=d['failure']))
    unique={}
    for p,h in refs.items():assert sha(ROOT/p)==h;unique.setdefault(csha(np.load(ROOT/p)),(p,h))
    replays=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for f in as_completed([pool.submit(replay,x) for x in unique.values()]):replays.append(f.result());write(audit/'CHECKPOINT.json',dict(completed=len(replays),total=len(unique),passed=all(r['passed'] for r in replays)))
    source=ROOT/'inputs/pinned14.npy';P=np.load(source)
    for r in records:
        C=np.load(ROOT/r['final_file']);r['pinned14_full100_alignment']=alignment(C,P);r['pinned14_source_sha256']=sha(source)
    report=dict(status='PASS' if all(r['passed'] for r in replays) else 'FAIL',unique_arrays=len(unique),fields=len(refs),all_triple_replays=replays,final_records=records,maximum_oracle_difference=max(r['absolute_difference'] for r in replays),proxy_solves=28,new_original_restores=4,cumulative_original_restores=23,remaining_authorized_original_budget=73,no_material_upper_improvement=True,no_global_claim=True,numerical_only=True)
    write(audit/'INDEPENDENT_RESTORE_AUDIT.json',report);print(json.dumps(dict(status=report['status'],arrays=len(unique),fields=len(refs),maxdifference=report['maximum_oracle_difference'],finals=[dict(id=r['restore_id'],R=r['final_diagnostic']['radius'],free=r['final_diagnostic']['nonbinding_center_indices'],delete_debts=r['final_diagnostic']['deletion_debts'],align=r['pinned14_full100_alignment']['maximum_all100_center_distance']) for r in records]),indent=2),flush=True)
if __name__=='__main__':main()
