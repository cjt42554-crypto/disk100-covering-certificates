# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Two new all-variable joint-event proxies. No original-problem restores."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,time,traceback,argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import ConvexHull,Delaunay
import geometry,event_jets,graph_tools
from helpers import write,sha,csha,inventory
from whole_disk_oracle import true_cover
ROOT=Path(__file__).resolve().parent.parent
EPS=2.**-18
def joint_jets(C,target):
    labels=target['primary_triangle'];values=[];rows=[]
    for i in target['added_sites_in_control_order']:
        j=event_jets.event_jet(C,'incircle',labels+[i]);values.append(j['raw_value']);rows.append(j['raw_jacobian'])
    return np.array(values),np.array(rows)
def fresh_gate(C,target,joint):
    S=target['primary_triangle'];added=target['added_sites_in_control_order'];ids=S+added;P=C[ids];hull=ConvexHull(P);convex=len(hull.vertices)==5;e=event_jets.event_jet(C,'incircle',S+[added[0]]);q=e['center'];rho=e['radius_squared'];mask=np.ones(100,bool);mask[ids]=False;nearest95=np.sum((C[mask]-q)**2,axis=1)-rho;h,J=joint_jets(C,target);faces={tuple(sorted(map(int,f))) for f in Delaunay(C).simplices};W=[];facechecks=[]
    height={v:np.eye(2)[added.index(v)] if v in added else np.zeros(2) for v in ids}
    for f in target['desired_post_faces']:
        mat=np.vstack((C[f].T,np.ones(3)))
        for v in set(ids)-set(f):W.append(height[v]-np.linalg.solve(mat,np.r_[C[v],1])@np.array([height[i] for i in f]))
        j=event_jets.triangle_jet(C,f);qq=j['center'];rr=j['radius_squared'];others=[v for v in range(100) if v not in f];m=float(np.min(np.sum((C[others]-qq)**2,axis=1)-rr));facechecks.append(dict(face=f,orientation=float(j['orientation']),q_norm_squared=float(qq@qq),minimum_other97_circle_power=m,face_present_in_actual_fullDT=tuple(sorted(f)) in faces,strict_local_gate=bool(j['orientation']>1e-10 and qq@qq<1-1e-8 and m>1e-10 and tuple(sorted(f)) in faces)))
    eta=np.array(target['postjoint_raw_eta']);fanmargin=float(np.min(np.array(W)@eta));inside5=bool(convex and np.all(np.cross(np.roll(P[hull.vertices],-1,axis=0)-P[hull.vertices],q-P[hull.vertices])>=-1e-10));qualified=bool(convex and q@q<1-1e-8 and nearest95.min()>1e-8 and np.max(np.abs(h))<1e-8 and fanmargin>1e-8) if joint else bool(convex and all(f['strict_local_gate'] for f in facechecks))
    return dict(kind='floating_joint_gate' if joint else 'post_five_site_threeface_gate',source_fan_recomputed_on_actual_geometry=True,raw_powers=h.tolist(),q=q.tolist(),rho_squared=float(rho),q_norm_squared=float(q@q),convex_five=convex,q_inside_five_hull_diagnostic=inside5,minimum_other95_circle_power=float(nearest95.min()),preferred_raw_eta_actual_fan_minimum_margin=fanmargin,preferred_face_checks=facechecks,numerically_qualified=qualified,exact_joint_or_fullcover_certificate=False)
def run_path(target):
    name=target['target_id'];folder=ROOT/'workers'/name
    for sub in ['starts','attempts','accepted','proposals']:(folder/sub).mkdir(exist_ok=True)
    source=ROOT/target['source_file'];assert sha(source)==target['source_file_sha256'];C0=np.load(source);Cref=C0.copy();e0=np.array(target['raw_source_e0']);normalizer=np.array(target['fixed_source_gradient_norm']);eta=np.array(target['postjoint_raw_eta']);h,J=joint_jets(C0,target);assert np.max(np.abs(h-e0))<1e-14 and np.max(np.abs(np.linalg.norm(J,axis=1)-normalizer))<1e-12
    w=np.column_stack((-C0[:,1],C0[:,0])).ravel();w/=np.linalg.norm(w);Ccurrent=C0.copy();Rcurrent,source_detail=true_cover(Ccurrent);phase='prejoint';fraction=0.;epsilon=0.;step=.1;accepted=[];attempts=[];solves=0;failed=0;termination=None;gate=None
    prior=set(json.load(open(ROOT/'prior_ordered_hashes.json')));prioru=set(json.load(open(ROOT/'prior_unordered_hashes.json')));seen=set()
    config=dict(raw_e0=e0.tolist(),fixed_gradient_normalizer=normalizer.tolist(),prejoint_raw_levels='e0*(1-f)',postjoint_raw_levels='epsilon*raw_eta',raw_eta=eta.tolist(),postjoint_epsilon=EPS,variables=201,spatial_variables=200,rotation_gauge_count=1,unit_norm_guard=None,actual_radius_cap=.1165,source_circle_motion_cap=.20,round_motion_cap=.05,dynamic_round_budget=8,accepted_budget=40,proxy_solve_budget=400,failure_budget=40,original_restores=0)
    write(folder/'CONFIG.json',dict(target=target,config=config,source_radius=Rcurrent,source_oracle=source_detail))
    while len(accepted)<40 and solves<400 and failed<40:
        if phase=='prejoint' and fraction>=1-1e-13:
            gate=fresh_gate(Ccurrent,target,True);write(folder/'FLOATING_JOINT_GATE.json',gate)
            if not gate['numerically_qualified']:termination='Actual floating joint or fresh preferred fan gate failed; no restart';break
            phase='postjoint';epsilon=0.;step=EPS/4
        if phase=='postjoint' and epsilon>=EPS*(1-1e-10):termination='Reached preferred finite raw post control';break
        if phase=='prejoint':
            goal=min(1.,fraction+step)
            for landmark in [.8,.9,.98,.995,.9995,.9999,1.]:
                if fraction+1e-13<landmark<goal:goal=landmark;break
            levels=(1-goal)*e0
        else:goal=min(EPS,epsilon+step);levels=goal*eta
        currenth,currentJ=joint_jets(Ccurrent,target);A=np.vstack((currentJ,w));rhs=np.r_[levels-currenth,-w@(Ccurrent.ravel()-C0.ravel())];delta=np.linalg.lstsq(A,rhs,rcond=1e-12)[0];C=Ccurrent+delta.reshape(100,2);sf=folder/'starts'/f'node{len(accepted):03d}_try{failed:03d}_{phase}.npy';np.save(sf,C);cs=csha(C);cu=csha(C[np.lexsort((C[:,1],C[:,0]))]);duplicate=cs in prior or cu in prioru or cs in seen;seen.add(cs)
        record=dict(node_index=len(accepted),trial_index=failed,phase=phase,goal=goal,target_raw_levels=levels.tolist(),start_file=str(sf.relative_to(ROOT)),start_file_sha256=sha(sf),coordinate_sha256=cs,unordered_sha256=cu,old_exact_duplicate=duplicate,predictor='Least-norm two actual full-space event gradients plus gauge; no fixed strong-contact mechanism',rounds=[])
        if duplicate:record['failure']='Duplicate predictor, no solve';attempts.append(record);failed+=1;step*=.5;continue
        ok=False
        try:
            for outer in range(8):
                if solves>=400:break
                top=geometry.topology(C);oldkeys=inventory(C,top)[2];oldset=set(map(tuple,oldkeys));base=C.copy();rr,_=true_cover(C);x0=np.r_[C.ravel(),rr*rr+1e-12]
                def fun(x):return x[-1]-inventory(x[:200].reshape(100,2),top)[0]
                def jac(x):return np.column_stack((-inventory(x[:200].reshape(100,2),top)[1],np.ones(len(oldkeys))))
                def eq(x):return (joint_jets(x[:200].reshape(100,2),target)[0]-levels)/normalizer
                def eqjac(x):return np.column_stack((joint_jets(x[:200].reshape(100,2),target)[1]/normalizer[:,None],np.zeros(2)))
                def motion(x):D=x[:200].reshape(100,2);return np.r_[.2**2-np.sum((D-Cref)**2,axis=1),.05**2-np.sum((D-base)**2,axis=1)]
                def motionjac(x):
                    D=x[:200].reshape(100,2);M=np.zeros((200,201));cols=2*np.arange(100)[:,None]+[0,1];M[np.arange(100)[:,None],cols]=-2*(D-Cref);M[100+np.arange(100)[:,None],cols]=-2*(D-base);return M
                lo=np.maximum(C.ravel()-.05,Cref.ravel()-.2);hi=np.minimum(C.ravel()+.05,Cref.ravel()+.2);tic=time.monotonic();solves+=1
                res=minimize(lambda x:x[-1],x0,jac=lambda x:np.r_[np.zeros(200),1.],method='SLSQP',constraints=[dict(type='ineq',fun=fun,jac=jac),dict(type='ineq',fun=motion,jac=motionjac),dict(type='eq',fun=eq,jac=eqjac),dict(type='eq',fun=lambda x:w@(x[:200]-C0.ravel()),jac=lambda x:np.r_[w,0])],bounds=list(zip(np.r_[lo,.09**2],np.r_[hi,.14**2])),options=dict(maxiter=500,ftol=1e-13,disp=False))
                C=res.x[:200].reshape(100,2);assert np.isfinite(C).all();p=folder/'attempts'/f'node{len(accepted):03d}_try{failed:03d}_round{outer:02d}.npy';np.save(p,C);actualR,detail=true_cover(C);g,JJ,keys,newtop=inventory(C);h,J=joint_jets(C,target);gap=actualR-np.sqrt(res.x[-1]);err=float(np.max(np.abs((h-levels)/e0)));newkeys=sorted(set(map(tuple,keys))-oldset);removed=sorted(oldset-set(map(tuple,keys)));missing=max([float(g[j]-res.x[-1]) for j,k in enumerate(keys) if tuple(k) not in oldset],default=0.);movement=float(np.linalg.norm(C-Cref,axis=1).max());fmin=float(fun(res.x).min());gauge_res=float(w@(C.ravel()-C0.ravel()))
                stage=dict(round=outer,proxy_solve_number=solves,file=str(p.relative_to(ROOT)),file_sha256=sha(p),coordinate_sha256=csha(C),status=int(res.status),success=bool(res.success),message=str(res.message),iterations=int(res.nit),actual_radius=actualR,actual_oracle=detail,epigraph_t=float(res.x[-1]),actual_minus_model_radius=float(gap),raw_event_values=h.tolist(),raw_event_goals=levels.tolist(),maximum_source_normalized_event_residual=err,gradient_scaled_residuals=eq(res.x).tolist(),rotation_gauge_residual=gauge_res,minimum_frozen_constraint=fmin,source_max_center_motion=movement,max_center_norm=float(np.linalg.norm(C,axis=1).max()),new_actual_keys=newkeys,removed_actual_keys=removed,new_hole_maximum_t_debt=missing,seconds=time.monotonic()-tic)
                record['rounds'].append(stage)
                if actualR<.1154469548855829-1e-10:
                    proposal=folder/'proposals'/p.name;np.save(proposal,C);write(proposal.with_suffix('.json'),dict(target=target,stage=stage,file_sha256=sha(proposal),numerical_only=True,needs_independent_and_exact_certificate=True));print('POTENTIAL_MATERIAL_UPPER',name,actualR,str(proposal),sha(proposal),flush=True)
                if actualR>.1165:record['failure']='Actual R research guard .1165 exceeded; no global infeasibility claim';break
                if err<2e-7 and abs(gap)<2e-8 and fmin>-2e-10 and abs(gauge_res)<1e-9 and movement<=.2+1e-9 and missing<2e-10:ok=True;break
                if movement>.2+1e-9:record['failure']='Source movement guard exceeded';break
            if ok:
                ap=folder/'accepted'/f'node{len(accepted):03d}_{phase}_{goal:.12f}.npy';np.save(ap,C)
                try:graph=graph_tools.extract(C,True)
                except Exception as e:graph=dict(status='NUMERICAL_GRAPH_UNCLASSIFIED',type=type(e).__name__,message=str(e))
                node=dict(index=len(accepted),phase=phase,goal=goal,raw_goal=levels.tolist(),file=str(ap.relative_to(ROOT)),file_sha256=sha(ap),coordinate_sha256=csha(C),actual_radius=actualR,status=int(res.status),success=bool(res.success),finite_feasible_status8_accept=bool(res.status==8),actual_graph=graph,maximum_event_residual=err,actual_minus_model_radius=float(gap),original_restart_qualified=False)
                if phase=='prejoint' and abs(goal-1)<1e-13:node['fresh_joint_gate']=fresh_gate(C,target,True)
                if phase=='postjoint':node['fresh_post_gate']=fresh_gate(C,target,False)
                accepted.append(node);Ccurrent=C.copy();Rcurrent=actualR
                if phase=='prejoint':fraction=goal;step=min(.1,step*1.5)
                else:epsilon=goal;step=min(EPS/2,step*1.5)
            else:failed+=1;step*=.5;record.setdefault('failure','Finite dynamic inventory/event feasibility not reached')
        except Exception as e:record['failure']=dict(type=type(e).__name__,message=str(e),traceback=traceback.format_exc());failed+=1;step*=.5
        attempts.append(record);write(folder/'CHECKPOINT.json',dict(target=target,config=config,phase=phase,fraction=fraction,epsilon=epsilon,accepted=accepted,attempts=attempts,proxy_solves=solves,failed_nodes=failed,original_restores=0))
        if step<(.0001 if phase=='prejoint' else EPS/10000):termination='Finite minimum phase-control step';break
    if termination is None:termination='Finite accepted/proxy/failure budget'
    result=dict(target=target,config=config,accepted=accepted,attempts=attempts,proxy_solves=solves,failed_nodes=failed,termination=termination,original_restores=0,no_exact_joint_or_radius_certificate=True,no_global_claim=True);write(folder/'RESULT.json',result);return dict(target_id=name,accepted_nodes=len(accepted),proxy_solves=solves,failed_nodes=failed,last_phase=phase,last_fraction=fraction,last_epsilon=epsilon,last_radius=Rcurrent,termination=termination,result_file=str((folder/'RESULT.json').relative_to(ROOT)),result_sha256=sha(folder/'RESULT.json'))
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--launch',action='store_true');args=parser.parse_args();assert args.launch
    d=json.load(open(ROOT/'PREPARED_PATHS.json'));assert d['status']=='READY_IMMUTABLE_REVIEW_GATED';review=ROOT.parent/'b32_five_site_joint_review';assert sha(review/'MANIFEST_SHA256.json')==d['review_manifest_sha256'];assert sha(review/'QUALIFIED_CONSTRAINTS.json')==d['source_qualified_constraints_sha256'];assert len(d['paths'])==2;out=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for f in as_completed([pool.submit(run_path,t) for t in d['paths']]):out.append(f.result());write(ROOT/'CHECKPOINT.json',dict(completed_paths=len(out),total_paths=2,results=out))
    write(ROOT/'PROXY_BATCH_RESULT.json',dict(status='FINITE_JOINT_PROXIES_FINISHED_PENDING_INDEPENDENT_AUDIT',results=out,proxy_solves=sum(r['proxy_solves'] for r in out),original_restores=0,no_automatic_original_restores=True));print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':main()
