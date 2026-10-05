# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Independent whole-disk replay and binding-skeleton completion classification."""
import engine as run
import rattler_search as search
from alignment import compare
from independent_oracle import all_triple_cover
import os,json,time,collections,multiprocessing as mp
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import linear_sum_assignment,nnls
ROOT=run.ROOT;BASE=ROOT.parent
def initializer(cpus):os.sched_setaffinity(0,{cpus[(mp.current_process()._identity[0]-1)%len(cpus)]})
def replay(task):
    kind,index,path,stored=task;C=np.load(path);R,d=all_triple_cover(C);return dict(kind=kind,index=index,path=str(Path(path).relative_to(ROOT)),radius=R,stored_radius=stored,absolute_error=abs(R-stored),detail=d)

def contact_compare(K,C,wide=False):
    """Explicit binding geometry+real active graph, with inactive completion kept."""
    kd=run.diagnose(K);cd=run.diagnose(C);ku=np.array(kd['binding_indices']);cu=np.array(cd['binding_indices']);kf=np.array(kd['nonbinding_indices'],dtype=int);cf=np.array(cd['nonbinding_indices'],dtype=int)
    if len(ku)!=len(cu) or kd['active_count']!=cd['active_count']:return dict(matches=False,reason='binding_or_active_contact_cardinality_differs')
    if wide:
        sk=search.skeleton_compare(K[ku],C[cu],144);candidates=[sk]
    else:
        candidates=[]
        for reflection in (False,True):
            V=C[cu].copy()
            if reflection:V[:,1]*=-1
            angle=0.
            for it in range(20):
                Q=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]);rows,cols=linear_sum_assignment(np.sum((K[ku,None,:]-(V@Q)[None,:,:])**2,axis=2));A=K[ku[rows]];Z=V[cols];angle=float(np.arctan2(np.sum(A[:,0]*Z[:,1]-A[:,1]*Z[:,0]),np.sum(A*Z)))
            candidates.append(dict(reflection=reflection,rotation_angle=angle,target_selected_rows=rows,source_selected_columns=cols))
    best=None
    for sk in candidates:
        V=C.copy()
        if sk['reflection']:V[:,1]*=-1
        angle=sk['rotation_angle'];Q=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]);V=V@Q;W=np.empty_like(C);rows=np.array(sk['target_selected_rows']);cols=np.array(sk['source_selected_columns']);W[ku[rows]]=V[cu[cols]]
        if len(kf):
            r,j=linear_sum_assignment(np.sum((K[kf,None,:]-V[cf][None,:,:])**2,axis=2));W[kf[r]]=V[cf[j]]
        wd=run.diagnose(W);d=np.linalg.norm(K-W,axis=1);active=kd['active_signature']==wd['active_signature'];inactive=kd['inactive_signature']==wd['inactive_signature'];match=bool(d[ku].max()<2e-6 and active);result=dict(matches=match,binding_center_count=len(ku),maximum_binding_center_distance=float(d[ku].max()),all100_maximum_center_distance=float(d.max()),active_inventory_matches=active,inactive_inventory_matches=inactive,full100_geometry_match=bool(d.max()<2e-6 and active and inactive),nonbinding_completion_displacements={int(f):float(d[f]) for f in kf},reflection=sk['reflection'],rotation_angle=angle,binding_target_rows=rows,binding_source_columns=cols,search='144-angle rotation/reflection+binding assignment' if wide else 'explicit20-iteration assignment/procrustes near orientation+reflection',scope='A binding contact family requires successful binding geometry and actual active-key match. Nonbinding completion positions and inactive topology remain separate; no disk is treated as deletable.')
        if best is None or result['maximum_binding_center_distance']<best['maximum_binding_center_distance']:best=result
        if match:return result
    return best

def main():
    tic=time.monotonic();cp=json.loads((ROOT/'checkpoint.json').read_text());done=json.loads((ROOT/'completion.json').read_text());rows=cp['records'];grid=json.loads((ROOT/'grid.json').read_text());source=json.loads((ROOT/'source_inventory.json').read_text());runtime=json.loads((ROOT/'runtime_provenance.json').read_text());ordered,unordered,nold,olderrors=search.prior_hashes();failures=[];startchecks=[];stagechecks=[];tasks=[];seen=set();seenu=set();frozen=search.frozen_checks();runtimechecks={p:run.sha(ROOT/p)==h for p,h in runtime['source_sha256'].items()}
    if not all(v['unchanged'] for v in frozen.values()):failures.append(dict(reason='original_frozen_tree_changed',trees=frozen))
    if not all(runtimechecks.values()):failures.append(dict(reason='runtime_source_changed',checks=runtimechecks))
    if olderrors:failures.append(dict(reason='prior_array_read_errors',errors=olderrors))
    if done['completed']!=cp['requested'] or len(rows)!=cp['requested'] or len({r['grid_id']for r in rows})!=cp['requested']:failures.append(dict(reason='count_or_repeated_grid_ids'))
    templates=[];inputs={}
    for s in source:
        i=s['source_index'];C=np.load(ROOT/s['copied_file']);inputs[i]=C
        if run.sha(ROOT/s['copied_file'])!=s['input_file_sha256'] or run.sha(s['input_original_file'])!=s['input_file_sha256']:failures.append(dict(reason='source_hash',source=i))
        tasks.append(('source',i,str(ROOT/s['copied_file']),s['radius']));dp=ROOT/'inputs'/f'source{i:02d}_delete_free{s["free_index"]:03d}.npy'
        if not np.array_equal(np.load(dp),np.delete(C,s['free_index'],axis=0)):failures.append(dict(reason='source99_deleted_coordinates',source=i))
        tasks.append(('nonbinding_deletion',i,str(dp),s['without_nonbinding_radius']));templates.append(dict(label=f'input_source{i}_99',file=s['copied_file'],coordinate_file_sha256=s['input_file_sha256']))
    refs=[('old65_known99','inputs/known99_reference.npy'),('old001','inputs/old001_100.npy'),('old003','inputs/old003_100.npy'),('current146','inputs/current146_100.npy'),('known115460','inputs/audit_known115460_100.npy')]+[(f'known_source{i}','inputs/audit_known_source%02d_100.npy'%i)for i in (3,4,5)]
    for label,path in refs:
        C=np.load(ROOT/path);templates.append(dict(label=label,file=path,coordinate_file_sha256=run.sha(ROOT/path)));tasks.append(('old_reference',label,str(ROOT/path),run.true_cover(C)[0]))
    for row in rows:
        if 'error' in row:failures.append(dict(reason='optimizer_exception',grid_id=row['grid_id'],error=row['error']));continue
        i=row['grid_id'];C=inputs[row['source_index']];free=row['free_center_index'];folder=ROOT/'workers'/f'worker{row["worker"]:02d}';p=folder/'starts'/f'{i:04d}.npy';S=np.load(p);digest=run.csha(S);uh=search.unordered_hash(S);used=np.array([j for j in range(100)if j!=free]);rr,dd=run.true_cover(S);sc=dict(grid_id=i,other99_bitwise_preserved=np.array_equal(S[used],C[used]),target_exactly_matches=np.array_equal(S[free],np.array(row['target']['point'])),start_coordinate_hash_matches=digest==row['start_coordinate_sha256'],start_unordered_geometry_hash_matches=uh==row['start_unordered_geometry_sha256'],start_file_hash_matches=run.sha(p)==row['start_file_sha256'],prior_ordered_duplicate=digest in ordered,prior_unordered_duplicate=uh in unordered,new_ordered_duplicate=digest in seen,new_unordered_duplicate=uh in seenu,initial_radius_error=abs(rr-row['initial_radius']),target_ledger_score_error=abs(rr-row['target']['start_actual_whole_disk_radius']));startchecks.append(sc);seen.add(digest);seenu.add(uh)
        if not all(sc[k]for k in ('other99_bitwise_preserved','target_exactly_matches','start_coordinate_hash_matches','start_unordered_geometry_hash_matches','start_file_hash_matches')) or any(sc[k]for k in ('prior_ordered_duplicate','prior_unordered_duplicate','new_ordered_duplicate','new_unordered_duplicate')) or max(sc['initial_radius_error'],sc['target_ledger_score_error'])>2e-11:failures.append(sc)
        tasks.append(('start',i,str(p),row['initial_radius']));final=ROOT/row['final_coordinate_file'];F=np.load(final)
        if run.sha(final)!=row['final_file_sha256'] or run.csha(F)!=row['final_coordinate_sha256']:failures.append(dict(reason='final_hash',grid_id=i))
        tasks.append(('final',i,str(final),row['final_radius']))
        for stage in row['rounds']:
            p=ROOT/stage['coordinate_file'];Z=np.load(p);rr,dd=run.true_cover(Z);sc=dict(path=stage['coordinate_file'],file_hash_matches=run.sha(p)==stage['file_sha256'],coordinate_hash_matches=run.csha(Z)==stage['coordinate_sha256'],score_error=abs(rr-stage['radius']));stagechecks.append(sc)
            if not sc['file_hash_matches'] or not sc['coordinate_hash_matches'] or sc['score_error']>2e-11:failures.append(sc)
            tasks.append(('stage',f'{i}:{stage["round"]}',str(p),stage['radius']))
    results=[];cpus=sorted(os.sched_getaffinity(0))
    with ProcessPoolExecutor(max_workers=4,initializer=initializer,initargs=(cpus,))as pool:
        for result in pool.map(replay,tasks,chunksize=1):
            results.append(result)
            if result['absolute_error']>2e-11:failures.append(result)
            if len(results)%32==0:run.write(ROOT/'audit_checkpoint.json',dict(completed=len(results),requested=len(tasks),maximum_error=max(x['absolute_error']for x in results),seconds=time.monotonic()-tic,**search.pidinfo()));print('AUDIT_PROGRESS',len(results),len(tasks),flush=True)
    successful=[r for r in rows if'final_radius'in r];best=min(successful,key=lambda r:r['final_radius']);B=np.load(ROOT/best['final_coordinate_file']);bestcomparisons={label:compare(np.load(ROOT/path),B,144)for label,path in refs};run.write(ROOT/'best_full100_reference_comparisons.json',dict(final_file=best['final_coordinate_file'],final_file_sha256=best['final_file_sha256'],grid_id=best['grid_id'],comparisons=bestcomparisons));np.save(ROOT/'best_final.npy',B);run.write(ROOT/'best_final.json',best)
    groups=[];classification=[];obstacles=[]
    for row in sorted(successful,key=lambda r:r['grid_id']):
        C=np.load(ROOT/row['final_coordinate_file']);dd=run.diagnose(C);matched=None
        for rep in templates+groups:
            cmp=contact_compare(np.load(ROOT/rep['file']),C)
            if cmp.get('matches'):matched=dict(label=rep['label'],comparison=cmp);break
        if matched is None:
            for rep in templates+groups:
                cmp=contact_compare(np.load(ROOT/rep['file']),C,True)
                if cmp.get('matches'):matched=dict(label=rep['label'],comparison=cmp);break
        if matched is None:
            label=f'new_numerical_contact_family{len(groups):02d}';groups.append(dict(label=label,file=row['final_coordinate_file'],representative_grid_id=row['grid_id'],representative_radius=row['final_radius'],file_sha256=row['final_file_sha256'],binding_count=dd['binding_count'],active_count=dd['active_count'],nonbinding_indices=dd['nonbinding_indices'],active_signature=dd['active_signature'],inactive_signature=dd['inactive_signature'],scope='Unmatched after finite coordinate/contact searches, not proven geometrically inequivalent.'));matched=dict(label=label,comparison=None)
        classification.append(dict(grid_id=row['grid_id'],source_index=row['source_index'],final_file=row['final_coordinate_file'],final_file_sha256=row['final_file_sha256'],radius=row['final_radius'],binding_count=dd['binding_count'],nonbinding_indices=dd['nonbinding_indices'],active_count=dd['active_count'],inactive_signature=dd['inactive_signature'],**matched));st=row['rounds'];obstacles.append(dict(grid_id=row['grid_id'],source_index=row['source_index'],target_kind=row['target']['kind'],target_key=row['target']['key'],radius=row['final_radius'],binding_count=dd['binding_count'],nonbinding_indices=dd['nonbinding_indices'],active_count=dd['active_count'],topology_events=sum(z['topology_changed']for z in st),new_actual_binding_keys=[z['new_actual_binding_keys']for z in st if z['new_actual_binding_keys']],maximum_actual_minus_frozen_objective=max(z['true_radius_minus_frozen_objective']for z in st),minimum_actual_minus_frozen_objective=min(z['true_radius_minus_frozen_objective']for z in st),final_feasible=st[-1]['feasible'],final_status=st[-1]['status'],rounds=len(st)))
        if len(classification)%16==0:print('CLASSIFY_PROGRESS',len(classification),len(successful),'newfamilies',len(groups),flush=True)
    grouping=dict(contact_family_counts=dict(collections.Counter(r['label']for r in classification)),classifications=classification,new_numerical_contact_families=groups,scope='Geometry and actual active contact inventory define only binding-skeleton families. Inactive topology/full100 completion displacements remain explicit. Finite nonmatch is not a proof of inequivalence; radius is never used as the grouping criterion.');run.write(ROOT/'contact_family_classification.json',grouping)
    # New multi-rattler representatives receive an explicit deletion-responsibility
    # inventory. Even several nonbinding circles are never assumed deletable.
    new_multi=[];extra_results=[]
    for rep in groups:
        if rep['binding_count']>=99:continue
        C=np.load(ROOT/rep['file']);individual=[]
        for f in rep['nonbinding_indices']:
            D=np.delete(C,f,axis=0);r,detail=run.true_cover(D);rr,ind=all_triple_cover(D);p=ROOT/'inputs'/f'{rep["label"]}_delete{f:03d}.npy';np.save(p,D);sc=dict(kind='new_multi_individual_nonbinding_deletion',index=f'{rep["label"]}:{f}',path=str(p.relative_to(ROOT)),radius=rr,stored_radius=r,absolute_error=abs(rr-r),detail=ind);extra_results.append(sc);individual.append(dict(index=f,radius_without_center=r,deletion_debt=r-rep['representative_radius'],actual_score=detail,independent_result=sc))
        D=np.delete(C,rep['nonbinding_indices'],axis=0);r,detail=run.true_cover(D);rr,ind=all_triple_cover(D);p=ROOT/'inputs'/f'{rep["label"]}_delete_all_nonbinding.npy';np.save(p,D);sc=dict(kind='new_multi_all_nonbinding_deletion',index=rep['label'],path=str(p.relative_to(ROOT)),radius=rr,stored_radius=r,absolute_error=abs(rr-r),detail=ind);extra_results.append(sc);new_multi.append(dict(representative=rep,individual_nonbinding_deletions=individual,all_nonbinding_deleted_radius=r,total_deletion_debt=r-rep['representative_radius'],all_nonbinding_deleted_actual_score=detail,all_deleted_independent_result=sc))
    for sc in extra_results:
        if sc['absolute_error']>2e-11:failures.append(sc)
    run.write(ROOT/'new_multi_nonbinding_responsibility.json',dict(representatives=new_multi,independent_extra_results=extra_results,scope='True whole-disk individual and joint deletion debts; nonbinding centers are essential inactive-cover resources, not freely removable disks.'))
    g,p,J,keys,sites,top=run.inventory(B);ac=np.flatnonzero(g.max()-g<1e-9);M=np.r_[J[ac].T,np.ones((1,len(ac)))];rhs=np.r_[np.zeros(B.size),1.];lam,res=nnls(M,rhs,maxiter=20000);run.write(ROOT/'fresh_best_diagnostic.json',dict(branch=f'fresh_rattler grid{best["grid_id"]} final only',input_file=best['final_coordinate_file'],input_file_sha256=best['final_file_sha256'],radius=best['final_radius'],diagnostic=run.diagnose(B),NNLS_residual_inf=float(np.max(np.abs(M@lam-rhs))),NNLS_lambda_min=float(lam.min()),NNLS_lambda_max=float(lam.max()),active_keys=[keys[a]for a in ac],NNLS_lambda=lam,scope='Numerical final diagnostic, not exact root, stationarity/local-minimum or global-optimality proof.'))
    obs=dict(paths_with_topology_event=sum(r['topology_events']>0 for r in obstacles),paths_with_new_binding_holes=sum(bool(r['new_actual_binding_keys'])for r in obstacles),paths_with_actual_gap_above1e8=sum(r['maximum_actual_minus_frozen_objective']>1e-8 for r in obstacles),final_not_feasible=sum(not r['final_feasible']for r in obstacles),rows=obstacles,scope='Finite saved-path observations only; no global topology barrier.');run.write(ROOT/'topology_obstacles.json',obs);material=[r for r in successful if r['final_radius']<run.THRESHOLD and run.REFERENCE-r['final_radius']>run.MATERIAL];allresults=results+extra_results
    report=dict(status='PASS'if not failures else'FAIL',grid_budget=96,planned=len(grid),optimized=len(rows),skipped=sum(r['skip']for r in grid),optimizer_exceptions=sum('error'in r for r in rows),all_triple_oracle_calls=len(allresults),independent_oracle_counts=dict(collections.Counter(r['kind']for r in allresults)),maximum_all_triple_error=max(r['absolute_error']for r in allresults),stage_array_count=len(stagechecks),stage_checks=stagechecks,start_checks=startchecks,prior_array_count=nold,prior_unique_ordered=len(ordered),prior_unique_unordered=len(unordered),prior_read_errors=olderrors,original_frozen_checks=frozen,runtime_source_checks=runtimechecks,worker_exitcodes=done['workers_exitcode'],final_SLSQP_status_counts=dict(collections.Counter(r['final_status']for r in obstacles)),optimizer_round_counts=dict(collections.Counter(r['rounds']for r in obstacles)),final_binding_counts=dict(collections.Counter(r['binding_count']for r in classification)),contact_family_counts=grouping['contact_family_counts'],new_contact_family_count=len(groups),new_multi_nonbinding_family_count=len(new_multi),best_grid_id=best['grid_id'],best_source_index=best['source_index'],best_free_index=best['free_center_index'],best_target=best['target'],best_final_radius=best['final_radius'],best_file=best['final_coordinate_file'],best_file_sha256=best['final_file_sha256'],best_known_full100_matches=[k for k,v in bestcomparisons.items()if v['matches']],material_candidate_count=len(material),independent_results=allresults,failures=failures,seconds=time.monotonic()-tic,scope='Finite new99-skeleton transport with full100 completion and independent continuous-boundary/triple numerical replay; no exact covering or global claim.');run.write(ROOT/'replay_report.json',report);summary={k:report[k]for k in ('status','grid_budget','planned','optimized','skipped','optimizer_exceptions','all_triple_oracle_calls','independent_oracle_counts','maximum_all_triple_error','stage_array_count','worker_exitcodes','final_SLSQP_status_counts','optimizer_round_counts','final_binding_counts','contact_family_counts','new_contact_family_count','new_multi_nonbinding_family_count','best_grid_id','best_source_index','best_free_index','best_final_radius','best_file','best_file_sha256','best_known_full100_matches','material_candidate_count','seconds','scope')};summary['topology_obstacles']={k:v for k,v in obs.items()if k not in ('rows','scope')};summary['original_all_six_frozen_trees_unchanged']=all(v['unchanged']for v in frozen.values());run.write(ROOT/'RESULT.json',summary)
    (ROOT/'README.txt').write_text('Three physically new99 binding skeleton rattler transports\n\nInput sources0/1/2 and slots24/73/80 are exact copies from the frozen hybrid handoff. Full100 rotation/reflection/assignment and real Voronoi/contact inventory comparisons plus separate99 binding geometry comparisons to old65 and old001/003/current146 are preserved before dependent optimization. This is not radius-only basin identification. The100th circle still has positive whole-disk deletion debt and important inactive cavity responsibility.\n\nThe finite budget gives priority sources1/2 all32 actual boundary maxima and8 actual internal maxima each; source0 gets12 boundary and4 internal controls. Internal/control selection covers cavity prior gain and fresh NNLS high/middle/weak stress. Every start initially moves only its nonbinding circle; the other99 coordinates are bitwise identical to that source. Ordered and sorted100-coordinate geometry hashes exclude every inspected prior initial/trace and every duplicate in this batch.\n\nAll100 centers become nooffset epigraph variables for at most two dynamic inventory/SLSQP rounds. Each status, exception, source SHA, skeleton SHA, target, prior gain, new actual keys, new binding holes and actual-minus-frozen objective is stored. Four BLAS1 workers use independent output files.\n\nIndependent replay covers every source, source deletion, known reference, start, stage and final with every triple circumcenter, every pair-circle intersection, antipodes and axes, scoring nearest distance to all remaining centers. New97/98 binding families additionally receive individual/joint nonbinding deletion-responsibility replays. Binding-skeleton families require explicit geometry and actual active contacts; inactive completion topology and all100 displacement remain separate to prevent rattler drift appearing as a new topology.\n\nNo numerical finite-search result establishes exact covering or global optimality/infeasibility. Reproduce in a writable copy: python3 -B source/rattler_search.py; python3 -B source/audit_finalize.py; python3 -B source/seal_result.py.\n');run.write(ROOT/'audit_completion.json',dict(status=report['status'],audit_source_sha256=run.sha(Path(__file__)),independent_oracle_source_sha256=run.sha(ROOT/'source/independent_oracle.py'),seconds=time.monotonic()-tic,**search.pidinfo()));print('AUDIT_DONE',json.dumps(run.clean(summary)),flush=True)
if __name__=='__main__':main()
