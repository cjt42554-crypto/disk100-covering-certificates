# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Three deduplicated incoming-cone curves; pure geometry, no cover objective.

Strong contacts are corrected to common variable t; weak and eight incoming
contacts are inequalities. Phase and rotation gauge select each finite curve.
Newton corrections use minimum Euclidean norm. No disk is removed.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,time,traceback,hashlib
import numpy as np
from scipy.spatial import Delaunay
import branches,graph_tools
from whole_disk_oracle import true_cover
from diagnose import ROOT,write,sha,clean

def csha(C):return hashlib.sha256(np.asarray(C,dtype='<f8').tobytes()).hexdigest()
def unordered(C):return csha(C[np.lexsort((C[:,1],C[:,0]))])
def facekeys(T):return sorted(tuple(sorted(map(int,t))) for t in T)

def main():
    report=json.loads((ROOT/'INCOMING_CONTACT_REVIEW.json').read_text());z=np.load(ROOT/'matrices/incoming_analysis.npz')
    C0=z['C'];strong=z['strong'];weak=z['weak'];near=z['near'];b=branches.Branches(C0);g0,p0,J0,keys0,top0=branches.inventory(C0)
    allb=branches.Branches(C0,np.arange(len(g0)));ineqids=np.r_[b.indices[weak],near];gauge=np.column_stack((-C0[:,1],C0[:,0])).ravel();gauge/=np.linalg.norm(gauge)
    selected=[]
    for r in sorted(report['finite_cone_trials'],key=lambda v:v['curvature']):
        if not r['success'] or not r['finite_conservative_model_feasible']:continue
        d=np.load(ROOT/r['direction_file']).ravel()
        if any(d@v[1]>.999999 for v in selected):continue
        assert r['maximum_absolute_cosine_to_owner16']<.999999
        selected.append((r,d))
    assert len(selected)==3
    previous=ROOT.parent/'b32_core_exchange_search';ordered=json.loads((previous/'batches/batch002/prior_ordered_hashes.json').read_text());unord=json.loads((previous/'batches/batch002/prior_unordered_hashes.json').read_text())
    old_manifest=json.loads((previous/'MANIFEST_SHA256.json').read_text());scanned=0
    for name,h in old_manifest.items():
        f=previous/name
        if f.suffix not in ('.npy','.npz'):continue
        assert sha(f)==h
        if f.suffix=='.npy':arrays=[np.load(f,allow_pickle=False)]
        else:
            with np.load(f,allow_pickle=False) as q:arrays=[q[k] for k in q.files]
        for C in arrays:
            if C.shape==(100,2):ordered.setdefault(csha(C),[]).append(str(f));unord.setdefault(unordered(C),[]).append(str(f));scanned+=1
    baseline_graph=graph_tools.extract(C0,True);baselineT=facekeys(top0[0]);baselineE=sorted((int(e[0]),int(e[1]),int(s)) for e,s in zip(top0[1],top0[2]));full0=facekeys(Delaunay(C0).simplices)
    branches_records=[];allnewhashes=set();start_time=time.monotonic()
    config=dict(input_sha256=sha(ROOT/'inputs/index14_grid0061.npy'),pure_geometry_no_cover_optimizer=True,curves=3,maximum_phase_parameter=.08,phase_step=.001,maximum_actual_radius=.11565,strong_equalities=len(strong),weak_and_incoming_inequalities=len(ineqids),newton_max_iterations=30,minimum_step=1e-6,old_ordered_count=len(ordered),old_unordered_count=len(unord),additional_frozen19_arrays_scanned=scanned,old_unordered_scope='Complete prior B32 profile inventory supplemented with all frozen87 and19 arrays; exact coordinate duplicate check, not metric basin proof.',selected_names=[r['name'] for r,d in selected])
    write(ROOT/'CONTINUATION_CONFIG.json',config)
    for branchid,(selected_row,d) in enumerate(selected):
        folder=ROOT/'stages'/f'curve{branchid:02d}';folder.mkdir(exist_ok=True)
        V=d.reshape(100,2).copy();V/=np.linalg.norm(V,axis=1).max();Vflat=V.ravel();speed=np.linalg.norm(Vflat);phase=d.copy()
        X=np.r_[C0.ravel(),g0.max()];current_s=0.;step=.001;rows=[];failures=[];first_full=None;first_actual_contact=None;termination=None
        strongids=b.indices[strong]
        def solve(goal,Xinit):
            x=Xinit.copy();activated=[];its=[]
            for iteration in range(30):
                C=x[:200].reshape(100,2);g,p,J=allb.evaluate(C);t=x[-1]
                for j in ineqids:
                    if g[j]-t>2e-14 and int(j) not in activated:activated.append(int(j))
                ids=np.r_[strongids,activated].astype(int)
                residual=np.r_[g[ids]-t,gauge@(x[:200]-C0.ravel()),phase@(x[:200]-C0.ravel())-goal*speed]
                violation=float(max(0.,np.max(g[ineqids]-t)));maximum=float(np.max(np.abs(residual)))
                its.append(dict(iteration=iteration,residual_max=maximum,maximum_inequality_violation=violation,activated_weak_keys=[keys0[j] for j in activated]))
                if maximum<5e-14 and violation<5e-14:return x,True,its,activated
                K=np.vstack((np.column_stack((J[ids],-np.ones(len(ids)))),np.r_[gauge,0],np.r_[phase,0]))
                correction=np.linalg.lstsq(K,-residual,rcond=1e-11)[0]
                scale=1.
                for backtrack in range(12):
                    trial=x+scale*correction
                    if trial[-1]<=0:scale*=.5;continue
                    gg=allb.evaluate(trial[:200].reshape(100,2))[0]
                    rr=np.r_[gg[ids]-trial[-1],gauge@(trial[:200]-C0.ravel()),phase@(trial[:200]-C0.ravel())-goal*speed]
                    if np.linalg.norm(rr)<=np.linalg.norm(residual)*(1-1e-4*scale)+1e-16:x=trial;its[-1].update(damping=scale,maximum_center_correction=np.linalg.norm((scale*correction[:200]).reshape(100,2),axis=1).max(),t_correction=float(scale*correction[-1]));break
                    scale*=.5
                else:its[-1]['failure']='Newton backtracking failed';return x,False,its,activated
            return x,False,its,activated
        try:
            while current_s<.08-1e-12:
                goal=current_s+step;predict=X.copy();predict[:200]+=step*Vflat
                candidate,success,iterations,activated=solve(goal,predict)
                if not success:
                    failures.append(dict(goal_phase=goal,step=step,newton_iterations=iterations,reason='pure_geometry_corrector_failed',candidate_file=f'stages/curve{branchid:02d}/failed{len(failures):03d}.npy'))
                    np.save(folder/f'failed{len(failures)-1:03d}.npy',candidate[:200].reshape(100,2));step*=.5
                    if step<1e-6:termination='Newton minimum step';break
                    continue
                C=candidate[:200].reshape(100,2);f=folder/f'stage{len(rows):03d}.npy';np.save(f,C)
                R,detail=true_cover(C);actualg,actualp,actualJ,actualkeys,actualtop=branches.inventory(C);actualgaps=actualg.max()-actualg
                try:graph=graph_tools.extract(C,True)
                except Exception as e:graph=dict(status='EXCEPTION',type=type(e).__name__,message=str(e))
                restricted_changed=facekeys(actualtop[0])!=baselineT;circle_changed=sorted((int(e[0]),int(e[1]),int(s)) for e,s in zip(actualtop[1],actualtop[2]))!=baselineE
                fullfaces=facekeys(Delaunay(C).simplices);fullchanged=fullfaces!=full0
                active=np.flatnonzero(actualgaps<=5e-14);oldactive=set(tuple(k) for k in b.keys);incoming_active=[actualkeys[j] for j in active if tuple(actualkeys[j]) not in oldactive]
                if fullchanged and first_full is None:first_full=len(rows)
                if incoming_active and first_actual_contact is None:first_actual_contact=len(rows)
                actual_minus_model=R-float(np.sqrt(candidate[-1]));safemodel=abs(actual_minus_model)<1e-9
                ch=csha(C);uh=unordered(C);duplicates=dict(ordered=ordered.get(ch,[]),unordered=unord.get(uh,[]),within_newbatch=ch in allnewhashes)
                allnewhashes.add(ch)
                row=dict(stage=len(rows),phase=goal,coordinate_file=str(f.relative_to(ROOT)),file_sha256=sha(f),coordinate_sha256=ch,unordered_sha256=uh,old_duplicates=duplicates,actual_radius=R,model_radius=np.sqrt(candidate[-1]),actual_minus_model_radius=actual_minus_model,all100_cover_oracle=detail,max_center_norm=np.linalg.norm(C,axis=1).max(),maximum_center_motion=np.linalg.norm(C-C0,axis=1).max(),strong_equality_residual=np.max(np.abs(allb.evaluate(C)[0][strongids]-candidate[-1])),weak_and_incoming_max_violation=np.max(allb.evaluate(C)[0][ineqids]-candidate[-1]),rotation_residual=gauge@(C.ravel()-C0.ravel()),incoming_actual_active_keys=incoming_active,strict_active_count=len(active),actual_graph=graph,full_DT_changed=fullchanged,restricted_DT_changed=restricted_changed,circle_inventory_changed=circle_changed,newton_iterations=iterations,activated_inequality_keys=[keys0[j] for j in activated],true_cover_matches_repaired_branch_model=safemodel,pure_geometry_no_objective_optimization=True)
                rows.append(row);X=candidate;current_s=goal
                write(folder/'CHECKPOINT.json',dict(selected_direction=selected_row,completed_stages=len(rows),last=row,failed_correctors=failures))
                if restricted_changed or circle_changed:termination='first actual restricted-DT or circle inventory event';break
                if not safemodel:termination='actual-cover takeover beyond repaired branch inventory';break
                if np.linalg.norm(C,axis=1).max()>1+1e-12:termination='center unit disk boundary';break
                if R>.11565:termination='finite actual radius cap .11565';break
                if step<.001:step=min(.001,step*1.5)
            if termination is None:termination='finite phase cap .08'
        except Exception as e:
            termination='exception';failures.append(dict(type=type(e).__name__,message=str(e),traceback=traceback.format_exc()))
        result=dict(curve=branchid,direction=selected_row,direction_norm_convention='max individual center speed1; phase is global unit direction projection',stages=rows,failures=failures,termination=termination,first_full_DT_change_stage=first_full,first_incoming_actual_contact_stage=first_actual_contact,new_cover_optimizer_runs=0)
        write(folder/'RESULT.json',result);branches_records.append(result)
        write(ROOT/'CONTINUATION_CHECKPOINT.json',dict(completed_curves=len(branches_records),elapsed_seconds=time.monotonic()-start_time,curves=branches_records))
        print('CURVE',branchid,selected_row['name'],'stages',len(rows),'lastR',rows[-1]['actual_radius'] if rows else None,'termination',termination,'fails',len(failures),flush=True)
    write(ROOT/'CONTINUATION_RESULT.json',dict(config=config,curves=branches_records,numerical_only=True,no_covering_optimizer=True,no_new_U=True,elapsed_seconds=time.monotonic()-start_time))
if __name__=='__main__':main()
