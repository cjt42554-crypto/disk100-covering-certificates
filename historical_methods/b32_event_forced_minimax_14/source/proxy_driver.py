# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Finite event-forced minimax paths, all200 center coordinates + t, one gauge.

Hard source/review hash gates precede execution. Proxy solves are counted
separately from zero-event original problem restores. This executable never
performs those restores. First physical event gates require independent audit.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,time,traceback,argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.optimize import minimize
import geometry,graph_tools,event_jets
from whole_disk_oracle import true_cover
ROOT=Path(__file__).resolve().parent.parent
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

def run_path(target):
    targetid=target['target_id'];folder=ROOT/'workers'/targetid
    for name in ['starts','attempts','accepted','diagnostics','proposals']:(folder/name).mkdir(parents=True,exist_ok=True)
    source=ROOT/target['source_file'];assert sha(source)==target['source_sha256'];C0=np.load(source);kind=target['event_kind'];labels=target['event_labels']
    reference=ROOT/target['pinned_reference_file'];assert sha(reference)==target['pinned_reference_sha256'];Cref=np.load(reference)
    jet0=event_jets.event_jet(C0,kind,labels);e0=jet0['raw_value'];normalizer=jet0['gradient_norm']
    assert abs(e0)>1e-16 and normalizer>1e-12
    assert abs(e0-target['source_raw_event_value'])<1e-11*max(abs(e0),1e-8)
    assert abs(normalizer-target['source_raw_event_gradient_norm'])<1e-10*max(normalizer,1e-8)
    w=np.column_stack((-C0[:,1],C0[:,0])).ravel();w/=np.linalg.norm(w)
    R0,detail0=true_cover(C0);currentC=C0.copy();currentR=R0;currentf=0.;step=.10;accepted=[];attempts=[];proxy_solves=0;failed_nodes=0;termination=None
    prior_hashes=set(json.loads((ROOT/'prior_ordered_hashes.json').read_text()));prior_unordered=set(json.loads((ROOT/'prior_unordered_hashes.json').read_text()));new_starts=set()
    sourcegraph=graph_tools.extract(C0,True)
    config=dict(source_sha256=target['source_sha256'],pinned_reference_sha256=target['pinned_reference_sha256'],objective='t=R^2',all_center_variables=200,gauge='Fixed source infinitesimal rotation orthogonality, one linear equality',other_anchors=0,event_level_semantics='e(C)/e(C0)=1-f',solver_event_row='[e(C)-(1-f)e0]/||gradient_e(source)||, fixed source scaling',event_gradient_norm=normalizer,source_event_value=e0,source_radius=R0,actual_radius_cap=.1165,pinned_reference_maximum_center_motion=.20,per_step_coordinate_trust=.05,unit_norm_guard=None,maximum_dynamic_inventory_rounds=8,accepted_node_budget=40,proxy_solve_budget=400,failed_node_budget=40,whole_original_restore_budget_consumed=0,include_all_actual_holes=True,expected_event_epigraph_row_dependencies_not_an_infeasibility_test=True)
    write(folder/'CONFIG.json',config)
    def event_fun(x,f):return (event_jets.event_jet(x[:200].reshape(100,2),kind,labels)['raw_value']-(1-f)*e0)/normalizer
    def event_jac(x):return np.r_[event_jets.event_jet(x[:200].reshape(100,2),kind,labels)['raw_jacobian']/normalizer,0.]
    while len(accepted)<40 and proxy_solves<400 and failed_nodes<40:
        if currentf>=1.0005-1e-13:termination='reached intended post-event level';break
        # Deliberate approach nodes keep resolution near the root, independently
        # of success/status. f=1 remains an approximate floating event only.
        landmarks=[.8,.9,.98,.995,.9995,1.,1.0005]
        goal=min(currentf+step,1.0005)
        for value in landmarks:
            if currentf+1e-13<value<goal:goal=value;break
        df=goal-currentf;C=currentC.copy();lastjet=event_jets.event_jet(C,kind,labels)
        if currentf==0. and 'ranking_direction_file' in target:
            jetfile=ROOT/target['ranking_direction_file'];assert sha(jetfile)==target['ranking_direction_file_sha256']
            with np.load(jetfile) as jets:direction=jets['ranking_jet_direction'][target['ranking_direction_index']]
            assert abs(lastjet['raw_jacobian']@direction/e0+1)<1e-8
            correction=df*direction;predictor_kind='Reviewed source full-space event jet, initial predictor only; strong contacts not fixed'
        else:
            correction=-df*e0*lastjet['raw_jacobian']/(lastjet['raw_jacobian']@lastjet['raw_jacobian']);predictor_kind='Current event-normal predictor; all centers optimize freely'
        C+=correction.reshape(100,2)
        startfile=folder/'starts'/f'node{len(accepted):03d}_try{failed_nodes:03d}.npy';np.save(startfile,C);startsha=csha(C);unordered=csha(C[np.lexsort((C[:,1],C[:,0]))]);duplicate=startsha in prior_hashes or unordered in prior_unordered or startsha in new_starts;new_starts.add(startsha)
        record=dict(node_index=len(accepted),trial_index=failed_nodes,current_fraction=currentf,target_fraction=goal,predictor_kind=predictor_kind,starts_file=str(startfile.relative_to(ROOT)),starts_file_sha256=sha(startfile),starts_coordinate_sha256=startsha,starts_unordered_sha256=unordered,old_exact_geometry_duplicate=duplicate,rounds=[])
        if duplicate:record['failure']='duplicate predictor array; no optimization run';attempts.append(record);failed_nodes+=1;step*=.5;continue
        node_ok=False
        try:
            for outer in range(8):
                if proxy_solves>=400:break
                top=geometry.topology(C);oldkeys=inventory(C,top)[2];oldkeyset=set(map(tuple,oldkeys));model_top=top
                def cfun(x):return x[-1]-inventory(x[:200].reshape(100,2),model_top)[0]
                def cjac(x):return np.column_stack((-inventory(x[:200].reshape(100,2),model_top)[1],np.ones(len(oldkeys))))
                R,_=true_cover(C);x0=np.r_[C.ravel(),R*R+1e-12]
                boundlo=np.maximum(C.ravel()-.05,Cref.ravel()-.20);boundhi=np.minimum(C.ravel()+.05,Cref.ravel()+.20)
                bounds=list(zip(np.r_[boundlo,.09**2],np.r_[boundhi,.14**2]))
                round_base=C.copy()
                def motion_guard(x):return np.r_[.20**2-np.sum((x[:200].reshape(100,2)-Cref)**2,axis=1),.05**2-np.sum((x[:200].reshape(100,2)-round_base)**2,axis=1)]
                def motion_guard_jac(x):
                    M=np.zeros((200,201));cols=(2*np.arange(100)[:,None]+[0,1]);M[np.arange(100)[:,None],cols]=-2*(x[:200].reshape(100,2)-Cref);M[100+np.arange(100)[:,None],cols]=-2*(x[:200].reshape(100,2)-round_base);return M
                tic=time.monotonic();proxy_solves+=1
                res=minimize(lambda x:x[-1],x0,jac=lambda x:np.r_[np.zeros(200),1.],method='SLSQP',constraints=[dict(type='ineq',fun=cfun,jac=cjac),dict(type='ineq',fun=motion_guard,jac=motion_guard_jac),dict(type='eq',fun=lambda x:event_fun(x,goal),jac=event_jac),dict(type='eq',fun=lambda x:w@(x[:200]-C0.ravel()),jac=lambda x:np.r_[w,0.])],bounds=bounds,options=dict(maxiter=500,ftol=1e-13,disp=False))
                C=res.x[:200].reshape(100,2);stagefile=folder/'attempts'/f'node{len(accepted):03d}_try{failed_nodes:03d}_round{outer:02d}.npy';np.save(stagefile,C)
                if not np.isfinite(C).all():raise ValueError('nonfinite proxy output')
                actualR,detail=true_cover(C);actualg,actualJ,actualkeys,newtop=inventory(C);eventraw=event_jets.event_jet(C,kind,labels)['raw_value'];normalized=eventraw/e0;gap=actualR-np.sqrt(res.x[-1]);eventerr=abs(normalized-(1-goal));motion=np.linalg.norm(C-Cref,axis=1).max();newkeys=sorted(set(map(tuple,actualkeys))-oldkeyset);removed=sorted(oldkeyset-set(map(tuple,actualkeys)));missingmax=max([float(actualg[j]-res.x[-1]) for j,k in enumerate(actualkeys) if tuple(k) not in oldkeyset],default=0.)
                if abs(goal-1)<1e-12:dependencies=event_jets.differential_dependency(C,kind,labels)
                else:dependencies=None
                stage=dict(round=outer,proxy_solve_number=proxy_solves,file=str(stagefile.relative_to(ROOT)),file_sha256=sha(stagefile),coordinate_sha256=csha(C),actual_radius=actualR,epigraph_t=float(res.x[-1]),model_radius=float(np.sqrt(res.x[-1])),actual_minus_model_radius=gap,actual_oracle=detail,event_raw_value=eventraw,event_normalized_value=normalized,event_normalized_residual=eventerr,event_equation_scaled_residual=abs(event_fun(res.x,goal)),rotation_gauge_residual=float(w@(C.ravel()-C0.ravel())),minimum_frozen_constraint=float(cfun(res.x).min()),source_maximum_center_motion=motion,max_center_norm=np.linalg.norm(C,axis=1).max(),status=int(res.status),success=bool(res.success),message=str(res.message),iterations=int(res.nit),seconds=time.monotonic()-tic,new_actual_keys=newkeys,removed_actual_keys=removed,new_hole_maximum_t_debt=missingmax,event_expected_differential_dependency=dependencies)
                record['rounds'].append(stage)
                if actualR<.1154469548855829-1e-10:
                    proposal=folder/'proposals'/f'node{len(accepted):03d}_try{failed_nodes:03d}_round{outer:02d}.npy';np.save(proposal,C);write(proposal.with_suffix('.json'),dict(target=target,stage=stage,proposal_file_sha256=sha(proposal),numerical_only=True,needs_independent_alltriple_and_exact_radius_certificate=True));print('POTENTIAL_MATERIAL_UPPER',targetid,actualR,str(proposal),sha(proposal),flush=True)
                feasible=eventerr<2e-7 and abs(gap)<2e-8 and cfun(res.x).min()>-2e-10 and abs(w@(C.ravel()-C0.ravel()))<1e-9 and motion<=.20+1e-9
                if actualR>.1165:record['failure']='actualR research guard .1165 exceeded; not a global infeasibility claim';break
                if feasible and missingmax<2e-10:node_ok=True;break
                if motion>.20+1e-9:record['failure']='source displacement research guard exceeded';break
            if node_ok:
                acceptedfile=folder/'accepted'/f'node{len(accepted):03d}_f{goal:.8f}.npy';np.save(acceptedfile,C)
                try:graph=graph_tools.extract(C,True)
                except Exception as error:graph=dict(status='NUMERICAL_GRAPH_UNCLASSIFIED',type=type(error).__name__,message=str(error),no_infeasibility_or_new_graph_claim=True)
                node=dict(accepted_index=len(accepted),fraction=goal,file=str(acceptedfile.relative_to(ROOT)),file_sha256=sha(acceptedfile),coordinate_sha256=csha(C),actual_radius=actualR,actual_oracle=detail,normalized_event_residual=eventerr,status=int(res.status),success=bool(res.success),finite_feasible_status8_accept=bool(res.status==8),actual_graph=graph,source_graph_changed=graph.get('fingerprint')!=sourcegraph['fingerprint'] if 'fingerprint' in graph else None,proxy_attempt_index=len(attempts),physical_event_endpoint_not_yet_qualified=True)
                accepted.append(node);currentC=C.copy();currentR=actualR;currentf=goal;step=min(.1,step*1.5)
            else:
                failed_nodes+=1;step*=.5
                if 'failure' not in record:record['failure']='dynamic actual-hole/event feasibility not reached within8 rounds'
        except Exception as e:
            record['failure']=dict(type=type(e).__name__,message=str(e),traceback=traceback.format_exc());failed_nodes+=1;step*=.5
        attempts.append(record);write(folder/'CHECKPOINT.json',dict(target=target,config=config,accepted=accepted,attempts=attempts,proxy_solves=proxy_solves,whole_original_restores=0,failed_nodes=failed_nodes))
        if step<.0001:termination='finite minimum fractional step';break
    if termination is None:termination='finite proxy/accepted/failure budget'
    result=dict(target=target,config=config,accepted=accepted,attempts=attempts,proxy_solves=proxy_solves,whole_original_restores=0,failed_nodes=failed_nodes,termination=termination,source_prior_array_is_reused_with_new_event_equations=True,no_exact_cocircular_claim=True,no_global_claim=True)
    write(folder/'RESULT.json',result);return dict(target_id=targetid,accepted_nodes=len(accepted),proxy_solves=proxy_solves,failed_nodes=failed_nodes,termination=termination,last_fraction=currentf,last_actual_radius=currentR,result_file=str((folder/'RESULT.json').relative_to(ROOT)),result_sha256=sha(folder/'RESULT.json'))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--launch',action='store_true');args=parser.parse_args()
    if not args.launch:raise SystemExit('Prepared executable only; use --launch after immutable review/source gates.')
    ready=ROOT/'PREPARED_TARGETS.json';assert ready.is_file(),'No targets prepared; not launching.'
    prepared=json.loads(ready.read_text());assert prepared['status']=='READY_AFTER_INDEPENDENT_TARGET_REVIEW'
    review=ROOT.parent/'b32_event_targeting_14_review';assert sha(review/'MANIFEST_SHA256.json')==prepared['review_manifest_sha256']
    selected=prepared['targets'];assert 1<=len(selected)<=2
    for row in selected:assert row['review_passed'] and sha(ROOT/row['source_file'])==row['source_sha256']
    summaries=[]
    with ProcessPoolExecutor(max_workers=min(2,len(selected))) as pool:
        for f in as_completed([pool.submit(run_path,row) for row in selected]):summaries.append(f.result());write(ROOT/'CHECKPOINT.json',dict(completed_paths=len(summaries),total_paths=len(selected),summaries=summaries))
    write(ROOT/'PROXY_BATCH_RESULT.json',dict(status='PROXY_PATHS_FINISHED_AWAITING_INDEPENDENT_PHYSICAL_GATES',paths=summaries,total_proxy_solves=sum(r['proxy_solves'] for r in summaries),whole_original_restores=0,worker_limit=2,no_automatic_original_restores=True))
    print(json.dumps(clean(summaries),indent=2),flush=True)
if __name__=='__main__':main()
