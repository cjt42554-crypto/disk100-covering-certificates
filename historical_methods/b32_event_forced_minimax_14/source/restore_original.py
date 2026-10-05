# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""At most two zero-event original minimax restores per independently gated path."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,time,traceback,hashlib
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.optimize import minimize
from proxy_driver import inventory,write,sha,csha
from whole_disk_oracle import true_cover
import geometry,graph_tools
ROOT=Path(__file__).resolve().parent.parent
def run(row):
    folder=ROOT/'original_restores'/row['restore_id'];folder.mkdir(parents=True,exist_ok=True)
    source=ROOT/'first_pilot_handoff'/row['source_handoff_file'];assert sha(source)==row['source_file_sha256'];C=np.load(source);np.save(folder/'start.npy',C)
    R,detail=true_cover(C);bestR=R;bestC=C.copy();w=np.column_stack((-C[:,1],C[:,0])).ravel();w/=np.linalg.norm(w);C0=C.copy();stages=[];failure=None
    for k in range(6):
        try:
            top=geometry.topology(C);keys=inventory(C,top)[2];old=set(map(tuple,keys));tic=time.monotonic()
            def fun(x):return 100*(x[-1]-inventory(x[:200].reshape(100,2),top)[0])
            def jac(x):return np.column_stack((-100*inventory(x[:200].reshape(100,2),top)[1],np.full(len(keys),100.)))
            rr,_=true_cover(C);x0=np.r_[C.ravel(),rr*rr+1e-12]
            res=minimize(lambda x:x[-1],x0,jac=lambda x:np.r_[np.zeros(200),1.],method='SLSQP',constraints=[dict(type='ineq',fun=fun,jac=jac),dict(type='eq',fun=lambda x:w@(x[:200]-C0.ravel()),jac=lambda x:np.r_[w,0.])],bounds=list(zip(np.r_[C.ravel()-.15,.09**2],np.r_[C.ravel()+.15,.2**2])),options=dict(maxiter=700,ftol=1e-14,disp=False))
            C=res.x[:200].reshape(100,2);assert np.isfinite(C).all();sf=folder/f'round{k:02d}.npy';np.save(sf,C)
            rr,dd=true_cover(C);g,J,actualkeys,newtop=inventory(C);new=sorted(set(map(tuple,actualkeys))-old);removed=sorted(old-set(map(tuple,actualkeys)));gap=rr-np.sqrt(res.x[-1])
            try:graph=graph_tools.extract(C,True)
            except Exception as e:graph=dict(status='UNCLASSIFIED',exception=str(e))
            stage=dict(round=k,file=str(sf.relative_to(ROOT)),file_sha256=sha(sf),coordinate_sha256=csha(C),status=int(res.status),success=bool(res.success),message=str(res.message),iterations=int(res.nit),actual_radius=rr,actual_oracle=dd,epigraph_t=float(res.x[-1]),actual_minus_model_radius=float(gap),minimum_frozen_scaled_constraint=float(fun(res.x).min()),gauge_residual=float(w@(C.ravel()-C0.ravel())),new_actual_keys=new,removed_actual_keys=removed,actual_graph=graph,seconds=time.monotonic()-tic,no_event_constraint=True,no_strong_contact_equalities=True,no_unit_center_guard=True,component_step_guard=.15)
            stages.append(stage)
            if rr<bestR:bestR=rr;bestC=C.copy()
            if rr<.1154469548855829-1e-10:
                proposal=folder/f'proposal_round{k:02d}.npy';np.save(proposal,C);write(proposal.with_suffix('.json'),dict(row=row,stage=stage,file_sha256=sha(proposal),numerical_only=True,exact_radius_certificate_needed=True));print('POTENTIAL_MATERIAL_UPPER',rr,str(proposal),sha(proposal),flush=True)
            write(folder/'CHECKPOINT.json',dict(row=row,stages=stages,best_radius=bestR))
            if res.success and abs(gap)<2e-10 and not new and not removed:break
        except Exception as e:failure=dict(type=type(e).__name__,message=str(e),traceback=traceback.format_exc());break
    final=folder/'final.npy';np.save(final,bestC)
    report=dict(row=row,start_file=str((folder/'start.npy').relative_to(ROOT)),start_file_sha256=sha(folder/'start.npy'),stages=stages,final_file=str(final.relative_to(ROOT)),final_file_sha256=sha(final),final_coordinate_sha256=csha(bestC),best_actual_radius=bestR,failure=failure,whole_original_restore_budget_consumed=1,no_global_claim=True,independent_full_triple_acceptance_pending=True)
    write(folder/'RESULT.json',report);return dict(restore_id=row['restore_id'],result_file=str((folder/'RESULT.json').relative_to(ROOT)),result_sha256=sha(folder/'RESULT.json'),best_actual_radius=bestR,rounds=len(stages),failure=failure)
def main():
    hpath=ROOT/'first_pilot_handoff/HANDOFF.json';h=json.loads(hpath.read_text());assert h['status']=='READY_INDEPENDENT_NUMERICAL_PHYSICAL_GATE' and h['physical_geometry_dedup_status']=='PASS'
    ready=json.loads((hpath.parent/'READY.json').read_text());assert sha(hpath.parent/'MANIFEST_SHA256.json')==ready['manifest_sha256'];manifest=json.loads((hpath.parent/'MANIFEST_SHA256.json').read_text());assert all(sha(hpath.parent/p)==s for p,s in manifest.items())
    prior=set(json.loads((ROOT/'prior_ordered_hashes.json').read_text()));unord=set(json.loads((ROOT/'prior_unordered_hashes.json').read_text()));rows=[]
    for event in h['events']:
        assert event['numerical_related_physical_event_gate']
        for frame in event['frames']:
            if frame['fraction']==1:continue
            C=np.load(hpath.parent/frame['handoff_file']);assert csha(C) not in prior and csha(C[np.lexsort((C[:,1],C[:,0]))]) not in unord
            rows.append(dict(restore_id=event['target_id']+('_before' if frame['fraction']<1 else '_after'),target_id=event['target_id'],event_labels=event['labels'],fraction=frame['fraction'],source_handoff_file=frame['handoff_file'],source_file_sha256=frame['file_sha256'],source_baseline_sha256=event['source_sha256']))
    assert len(rows)==4 and len({r['source_file_sha256'] for r in rows})==4
    write(ROOT/'ORIGINAL_RESTORE_LAUNCH.json',dict(status='SOURCE_AND_EVENT_GATES_PASS',source_handoff_manifest_sha256=ready['manifest_sha256'],whole_original_restores_requested=4,earlier_original_restores=19,cumulative_budget_limit=96,remaining_after_requested=73,rows=rows,proxy_solves_separately_counted=28,independent_source_and_event_gate_passed=True))
    result=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for f in as_completed([pool.submit(run,row) for row in rows]):result.append(f.result());write(ROOT/'ORIGINAL_RESTORE_CHECKPOINT.json',dict(completed=len(result),total=4,results=result))
    write(ROOT/'ORIGINAL_RESTORE_RESULT.json',dict(status='OPTIMIZERS_FINISHED_INDEPENDENT_AUDIT_PENDING',results=result,whole_original_restores=4,earlier_original_restores=19,total_cumulative=23,remaining_authorized_budget=73))
    print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
