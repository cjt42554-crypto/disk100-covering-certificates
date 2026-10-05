# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Independent finite B32 metric audit: all triples, all boundary extrema."""
import engine as run
import metric_search as search
from metric_validation import validate
from independent_oracle import all_triple_cover
from contact_geometry import contact_compare
from alignment import compare
import os,json,time,collections,multiprocessing as mp
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import nnls
ROOT=run.ROOT;BASE=ROOT.parent

def initializer(cpus):
    os.sched_setaffinity(0,{cpus[(mp.current_process()._identity[0]-1)%len(cpus)]})

def replay(task):
    kind,index,file,stored=task;R,detail=all_triple_cover(np.load(ROOT/file))
    return dict(kind=kind,index=index,file=file,radius=R,stored_radius=stored,absolute_error=abs(R-stored),detail=detail)

def short_metric(v):
    return {k:x for k,x in v.items()if k!='actual_graph'}

def main():
    tic=time.monotonic();cp=json.loads((ROOT/'checkpoint.json').read_text());rows=cp['records'];grid=json.loads((ROOT/'grid.json').read_text());complete=json.loads((ROOT/'completion.json').read_text());runtime=json.loads((ROOT/'runtime_provenance.json').read_text());checks=json.loads((ROOT/'input_checks.json').read_text());catalog=json.loads((search.HANDOFF/'catalog.json').read_text());graphs={r['independent_index']:r for r in catalog['graphs']};ledger=json.loads((ROOT/'screening/handoff_seed_validation.json').read_text());old=json.loads((ROOT/'prior_ordered_hashes.json').read_text());oldu=json.loads((ROOT/'prior_unordered_hashes.json').read_text());oldrefs=json.loads((ROOT/'inputs/old_metric_references.json').read_text());failures=[];startchecks=[];stagechecks=[];tasks=[];actualchecks=[];seen=set();seenu=set();frozen=search.frozen_checks();sourcechecks=[]
    if len(rows)!=87 or complete['completed']!=87 or len({r['grid_id']for r in rows})!=87:failures.append(dict(reason='count87_or_duplicate_grid_ids'))
    if any(x!=0 for x in complete['worker_exitcodes']):failures.append(dict(reason='worker_exitcodes',codes=complete['worker_exitcodes']))
    runtimechecks={p:run.sha(ROOT/p)==h for p,h in runtime['source_sha256'].items()}
    if not all(runtimechecks.values()):failures.append(dict(reason='runtime_source_changed',checks=runtimechecks))
    if not all(v['unchanged']for v in frozen.values())or not frozen['catalog_handoff']['owner_unchanged']:failures.append(dict(reason='old_or_handoff_changed',checks=frozen))
    hm=json.loads((ROOT/'inputs/history_snapshot/MANIFEST_SHA256.json').read_text());historychecks={p:run.sha(ROOT/'inputs/history_snapshot'/p)==h for p,h in hm.items()}
    if not all(historychecks.values()):failures.append(dict(reason='sealed_history_changed'))
    for r in ledger['records']:
        p=search.HANDOFF/r['coordinate_relative_file'];C=np.load(p);sc=dict(independent_index=r['independent_index'],file=str(p.relative_to(ROOT)),file_sha_match=run.sha(p)==r['file_sha256'],coordinate_sha_match=run.csha(C)==r['coordinate_sha256']);sourcechecks.append(sc)
        if not sc['file_sha_match']or not sc['coordinate_sha_match']:failures.append(sc)
        tasks.append(('raw_harmonic_seed',len(tasks),str(p.relative_to(ROOT)),r['actual_whole_disk_radius']))
    for r in oldrefs:
        p=ROOT/r['file'];C=np.load(p);sc=dict(old_independent_index=r['independent_index'],file=r['file'],file_sha_match=run.sha(p)==r['source_file_sha256'],coordinate_sha_match=run.csha(C)==r['source_coordinate_sha256']);sourcechecks.append(sc)
        if not sc['file_sha_match']or not sc['coordinate_sha_match']:failures.append(sc)
        tasks.append(('old_metric_reference',r['independent_index'],r['file'],run.true_cover(C)[0]))
    for r in rows:
        i=r['grid_id']
        if'error'in r:failures.append(dict(reason='optimizer_exception',grid_id=i,error=r['error']));continue
        p=ROOT/'workers'/f'worker{r["worker"]:02d}'/'starts'/f'{i:04d}.npy';C=np.load(p);cs=run.csha(C);us=search.unordered_hash(C);R,d=run.true_cover(C);v=validate(C,graphs[r['independent_index']]);sc=dict(grid_id=i,file=str(p.relative_to(ROOT)),owner_input_bitwise_equal=np.array_equal(C,np.load(search.HANDOFF/r['owner_seed_file'])),
            file_sha_match=run.sha(p)==r['start_file_sha256'],coordinate_sha_match=cs==r['initial_coordinate_sha256'],unordered_sha_match=us==r['initial_unordered_geometry_sha256'],
            old_ordered_duplicate=cs in old,old_unordered_duplicate=us in oldu,new_ordered_duplicate=cs in seen,new_unordered_duplicate=us in seenu,
            complete_historical_unordered_invariant_check=r['historical_unordered_geometry_check'],target_metric_realized=v['target_metric_realized'],start_score_error=abs(R-r['initial_radius']));seen.add(cs);seenu.add(us);startchecks.append(sc)
        if not all(sc[k]for k in('owner_input_bitwise_equal','file_sha_match','coordinate_sha_match','unordered_sha_match','target_metric_realized'))or any(sc[k]for k in('old_ordered_duplicate','old_unordered_duplicate','new_ordered_duplicate','new_unordered_duplicate'))or sc['start_score_error']>2e-11:failures.append(sc)
        p=ROOT/r['final_coordinate_file'];C=np.load(p);v=validate(C,graphs[r['independent_index']]);actualchecks.append(dict(grid_id=i,stage='final',validation=short_metric(v)))
        if run.sha(p)!=r['final_file_sha256']or run.csha(C)!=r['final_coordinate_sha256']:failures.append(dict(reason='final_hash',grid_id=i))
        if v['actual_fingerprint']!=r['final_target_and_actual_metric'].get('actual_fingerprint'):failures.append(dict(reason='final_graph_replay',grid_id=i))
        tasks.append(('final',i,r['final_coordinate_file'],r['final_radius']))
        for z in r['rounds']:
            p=ROOT/z['coordinate_file'];C=np.load(p);R,d=run.true_cover(C);v=validate(C,graphs[r['independent_index']]);actualchecks.append(dict(grid_id=i,stage=z['round'],validation=short_metric(v)));sc=dict(file=z['coordinate_file'],file_sha_match=run.sha(p)==z['file_sha256'],coordinate_sha_match=run.csha(C)==z['coordinate_sha256'],score_error=abs(R-z['radius']),actual_graph_replay_equal=v.get('actual_fingerprint')==z['target_and_actual_metric'].get('actual_fingerprint'));stagechecks.append(sc)
            if not sc['file_sha_match']or not sc['coordinate_sha_match']or not sc['actual_graph_replay_equal']or sc['score_error']>2e-11:failures.append(sc)
            tasks.append(('stage',f'{i}:{z["round"]}',z['coordinate_file'],z['radius']))
    results=[];cpus=sorted(os.sched_getaffinity(0))
    with ProcessPoolExecutor(max_workers=4,initializer=initializer,initargs=(cpus,))as pool:
        for result in pool.map(replay,tasks,chunksize=1):
            results.append(result)
            if result['absolute_error']>2e-11:failures.append(result)
            if len(results)%32==0:run.write(ROOT/'audit_checkpoint.json',dict(completed=len(results),requested=len(tasks),maximum_error=max(x['absolute_error']for x in results),seconds=time.monotonic()-tic,**search.pidinfo()));print('METRIC_AUDIT_PROGRESS',len(results),len(tasks),flush=True)
    successful=[r for r in rows if'final_radius'in r];best=min(successful,key=lambda r:r['final_radius']);C=np.load(ROOT/best['final_coordinate_file']);np.save(ROOT/'best_final.npy',C);run.write(ROOT/'best_final.json',best)
    bestbaseline=compare(np.load(ROOT/'inputs/current146_100.npy'),C,144);run.write(ROOT/'best_current146_alignment.json',dict(final_file=best['final_coordinate_file'],final_sha256=best['final_file_sha256'],comparison=bestbaseline))
    references=[dict(label=f'old_index{r["independent_index"]:02d}',file=r['file'])for r in oldrefs];groups=[];classified=[];obstacles=[]
    for r in sorted(successful,key=lambda r:r['grid_id']):
        C=np.load(ROOT/r['final_coordinate_file']);diag=run.diagnose(C);match=None
        for rep in references+groups:
            cmp=contact_compare(np.load(ROOT/rep['file']),C)
            if cmp.get('matches'):match=dict(label=rep['label'],comparison=cmp);break
        if match is None:
            for rep in references+groups:
                cmp=contact_compare(np.load(ROOT/rep['file']),C,True)
                if cmp.get('matches'):match=dict(label=rep['label'],comparison=cmp);break
        if match is None:
            label=f'finite_unmatched_metric_family{len(groups):02d}';groups.append(dict(label=label,file=r['final_coordinate_file'],representative_grid_id=r['grid_id'],independent_index=r['independent_index'],actual_fingerprint=r['final_target_and_actual_metric'].get('actual_fingerprint'),radius=r['final_radius'],binding_count=diag['binding_count'],active_count=diag['active_count'],scope='Finite metric/contact classification; no uniqueness or local-minimum proof.'));match=dict(label=label,comparison=None)
        classified.append(dict(grid_id=r['grid_id'],independent_index=r['independent_index'],final_file_sha256=r['final_file_sha256'],radius=r['final_radius'],binding_count=diag['binding_count'],active_count=diag['active_count'],nonbinding_indices=diag['nonbinding_indices'],actual_fingerprint=r['final_target_and_actual_metric'].get('actual_fingerprint'),target_graph_realized=r['final_target_and_actual_metric']['target_metric_realized'],**match))
        z=r['rounds'];obstacles.append(dict(grid_id=r['grid_id'],independent_index=r['independent_index'],rounds=len(z),initial_radius=r['initial_radius'],final_radius=r['final_radius'],final_status=z[-1]['status'],final_success=z[-1]['success'],final_feasible=z[-1]['feasible'],
            topology_events=sum(s['topology_changed']for s in z),new_binding_keys=[s['new_actual_binding_keys']for s in z if s['new_actual_binding_keys']],maximum_actual_minus_frozen=max(s['true_radius_minus_frozen_objective']for s in z),target_graph_retained_at_final=r['final_target_and_actual_metric']['target_metric_realized']))
    grouping=dict(counts=dict(collections.Counter(r['label']for r in classified)),classifications=classified,finite_unmatched_metric_families=groups,scope='Radius does not define a family: binding-center geometry and actual active/inactive inventories are compared; complete100 geometry is required for old-baseline identities. Finite failed matching proves no inequivalence.');run.write(ROOT/'contact_family_classification.json',grouping)
    graphstats=[]
    for i in sorted({r['independent_index']for r in successful}):
        rr=[r for r in successful if r['independent_index']==i];b=min(rr,key=lambda r:r['final_radius']);B=np.load(ROOT/b['final_coordinate_file']);g,p,J,kk,sites,top=run.inventory(B);active=np.flatnonzero(g.max()-g<1e-9);M=np.r_[J[active].T,np.ones((1,len(active)))];rhs=np.r_[np.zeros(B.size),1.];lam,res=nnls(M,rhs,maxiter=20000)
        info=dict(independent_index=i,target_fingerprint=b['target_fingerprint'],starts=len(rr),best_grid_id=b['grid_id'],best_radius=b['final_radius'],best_file=b['final_coordinate_file'],best_file_sha256=b['final_file_sha256'],max_radius=max(r['final_radius']for r in rr),target_retained_finals=sum(r['final_target_and_actual_metric']['target_metric_realized']for r in rr),binding_count=run.diagnose(B)['binding_count'],active_count=len(active),NNLS_residual_inf=float(np.max(np.abs(M@lam-rhs))),NNLS_lambda_min=float(lam.min()),NNLS_lambda_max=float(lam.max()),scope='Numeric fresh best-of-six metric diagnosis for this inputSHA, not all-metrics optimum or stationary/root certificate.');graphstats.append(info);run.write(ROOT/'screening'/f'graph{i:02d}_best_diagnostic.json',dict(info,active_keys=[kk[a]for a in active],lambda_values=lam))
    run.write(ROOT/'graph_metric_results.json',dict(graphs=graphstats,scope='Finite69+18 harmonic-original-problem attempts; no global or per-family optimum proof.'))
    obs=dict(paths_with_topology_event=sum(r['topology_events']>0 for r in obstacles),paths_with_new_binding=sum(bool(r['new_binding_keys'])for r in obstacles),paths_with_actual_gap_above1e8=sum(r['maximum_actual_minus_frozen']>1e-8 for r in obstacles),final_not_feasible=sum(not r['final_feasible']for r in obstacles),rows=obstacles);run.write(ROOT/'topology_obstacles.json',obs)
    material=[r for r in successful if r['final_radius']<run.THRESHOLD and run.REFERENCE-r['final_radius']>run.MATERIAL]
    report=dict(status='PASS'if not failures else'FAIL',budget=96,owner_seeds=96,optimized=len(rows),skipped=len(grid)-len(rows),skip_reasons=dict(collections.Counter(r['skip_reason']for r in grid if r['skip'])),optimizer_exceptions=sum('error'in r for r in rows),
        independent_oracle_calls=len(results),independent_counts=dict(collections.Counter(r['kind']for r in results)),maximum_all_triple_error=max(r['absolute_error']for r in results),source_checks=sourcechecks,start_checks=startchecks,stage_checks=stagechecks,actual_metric_replay_checks=actualchecks,
        runtime_source_checks=runtimechecks,frozen_checks=frozen,history_payload_SHA_checks=historychecks,worker_exitcodes=complete['worker_exitcodes'],final_SLSQP_status_counts=dict(collections.Counter(r['final_status']for r in obstacles)),round_counts=dict(collections.Counter(r['rounds']for r in obstacles)),
        final_binding_counts=dict(collections.Counter(r['binding_count']for r in classified)),final_target_graph_retained=sum(r['target_graph_realized']for r in classified),contact_family_counts=grouping['counts'],finite_unmatched_metric_family_count=len(groups),best_grid_id=best['grid_id'],best_independent_index=best['independent_index'],best_radius=best['final_radius'],best_file=best['final_coordinate_file'],best_file_sha256=best['final_file_sha256'],best_matches_current146=bestbaseline['matches'],material_candidate_count=len(material),independent_results=results,failures=failures,seconds=time.monotonic()-tic,scope='Finite numerical all-triple+continuous-boundary acceptance and actual restricted-DT/contact matching; no new exact upper certificate, graph completeness, all-metric minimum or global conclusion.')
    run.write(ROOT/'replay_report.json',report);summary={k:report[k]for k in('status','budget','owner_seeds','optimized','skipped','skip_reasons','optimizer_exceptions','independent_oracle_calls','independent_counts','maximum_all_triple_error','worker_exitcodes','final_SLSQP_status_counts','round_counts','final_binding_counts','final_target_graph_retained','contact_family_counts','finite_unmatched_metric_family_count','best_grid_id','best_independent_index','best_radius','best_file','best_file_sha256','best_matches_current146','material_candidate_count','seconds','scope')};summary['topology_obstacles']={k:v for k,v in obs.items()if k!='rows'};summary['old_ten_frozen_trees_unchanged']=all(v['unchanged']for v in frozen.values());run.write(ROOT/'RESULT.json',summary)
    run.write(ROOT/'audit_completion.json',dict(status=report['status'],audit_source_sha256=run.sha(Path(__file__)),seconds=time.monotonic()-tic,**search.pidinfo()));print('METRIC_AUDIT_DONE',json.dumps(run.clean(summary)),flush=True)

if __name__=='__main__':main()
