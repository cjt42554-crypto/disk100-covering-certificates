# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Independent paired-transport replay, complete geometry comparison and freeze."""
import engine as run
import pair_search as search
from independent_oracle import all_triple_cover
import os,json,time,collections,multiprocessing as mp
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import linear_sum_assignment,nnls
ROOT=run.ROOT;BASE=ROOT.parent

def initializer(cpus):
    os.sched_setaffinity(0,{cpus[(mp.current_process()._identity[0]-1)%len(cpus)]})

def replay(task):
    kind,index,path,stored=task;C=np.load(path);R,d=all_triple_cover(C)
    return dict(kind=kind,index=index,path=str(Path(path).relative_to(ROOT)),radius=R,stored_radius=stored,absolute_error=abs(R-stored),detail=d)

def compare(A,B,seeds=144):
    best=None
    for reflection in (False,True):
        V=B.copy()
        if reflection:V[:,1]*=-1
        for initial in np.linspace(0,2*np.pi,seeds,endpoint=False):
            angle=float(initial)
            for it in range(10):
                rot=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]);W=V@rot;_,j=linear_sum_assignment(np.sum((A[:,None,:]-W[None,:,:])**2,axis=2));Z=V[j];angle=float(np.arctan2(np.sum(A[:,0]*Z[:,1]-A[:,1]*Z[:,0]),np.sum(A*Z)))
            rot=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]);W=V[j]@rot;dist=np.linalg.norm(A-W,axis=1);error=float(np.sum(dist*dist))
            if best is None or error<best['sum_squared_error']:
                da=run.diagnose(A);db=run.diagnose(W);best=dict(sum_squared_error=error,maximum_center_distance=float(dist.max()),rotation_angle=angle,reflection=reflection,permutation=j,active_inventory_matches=da['active_signature']==db['active_signature'],inactive_inventory_matches=da['inactive_signature']==db['inactive_signature'])
    best.update(matches=bool(best['maximum_center_distance']<2e-6 and best['active_inventory_matches'] and best['inactive_inventory_matches']),rotation_seeds=seeds,scope='All100 coordinate alignment plus fresh active and inactive inventories. A failed finite search is not proof of inequivalence. Radius does not define a basin.')
    return best

def main():
    tic=time.monotonic();cp=json.loads((ROOT/'checkpoint.json').read_text());complete=json.loads((ROOT/'completion.json').read_text());rows=cp['records'];grid=json.loads((ROOT/'grid.json').read_text());baseline=np.load(ROOT/'inputs/baseline.npy');donors=json.loads((ROOT/'donor_inventory.json').read_text());screens=json.loads((ROOT/'pair_screening.json').read_text());runtime=json.loads((ROOT/'runtime_provenance.json').read_text());prior,prior_arrays,prior_errors=search.prior_start_hashes();failures=[];startchecks=[];stagechecks=[];tasks=[];seen=set();swapgroups={};frozenchecks={}
    for name in ('rattler_contact_transfer','best147_boundary_transfer','next146_boundary_transfer'):
        mf=json.loads((BASE/name/'MANIFEST_SHA256.json').read_text());checks={p:run.sha(BASE/name/p)==h for p,h in mf.items()};frozenchecks[name]=dict(count=len(checks),checks=checks,all_unchanged=all(checks.values()))
        if not all(checks.values()):failures.append(dict(reason='frozen_tree_changed',tree=name))
    sourcechecks={p:run.sha(ROOT/p)==h for p,h in runtime['source_sha256'].items()}
    if not all(sourcechecks.values()):failures.append(dict(reason='runtime_source_changed',checks=sourcechecks))
    if complete['completed']!=cp['requested'] or len(rows)!=cp['requested'] or len({r['grid_id'] for r in rows})!=cp['requested']:failures.append(dict(reason='grid_count_or_duplicate_grid_ids'))
    if prior_errors:failures.append(dict(reason='prior_start_read_errors',errors=prior_errors))
    tasks.append(('baseline',0,str(ROOT/'inputs/baseline.npy'),run.REFERENCE))
    for d in donors:
        p=ROOT/'inputs'/f'delete_{d["index"]:03d}.npy';C=np.load(p)
        if not np.array_equal(C,np.delete(baseline,d['index'],axis=0)):failures.append(dict(reason='single_deletion_coordinates',index=d['index']))
        tasks.append(('single_deletion',d['index'],str(p),d['deletion_radius']))
    for row in screens:
        p=ROOT/row['deletion_file'];C=np.load(p);labels=[i for i in range(100) if i not in row['pair']]
        if not np.array_equal(C,baseline[labels]) or run.sha(p)!=row['deletion_file_sha256']:failures.append(dict(reason='pair_deletion_coordinates_or_hash',index=row['screening_index']))
        tasks.append(('pair_deletion',row['screening_index'],str(p),row['pair_deletion_radius']))
    for row in rows:
        if 'error' in row:failures.append(dict(reason='optimizer_exception',grid_id=row['grid_id'],error=row['error']));continue
        i=row['grid_id'];folder=ROOT/'workers'/f'worker{row["worker"]:02d}';p=folder/'starts'/f'{i:04d}.npy';S=np.load(p);digest=run.csha(S);labels=[k for k in range(100) if k not in row['pair']];other=np.array_equal(S[labels],baseline[labels]);target=np.array_equal(S[row['pair']],np.array([h['point'] for h in row['targets']]));rr,dd=run.true_cover(S)
        sc=dict(grid_id=i,pair=row['pair'],assignment=row['assignment'],other98_exactly_preserved=other,targets_exactly_match=target,start_coordinate_hash_matches=digest==row['start_coordinate_sha256'],start_file_hash_matches=run.sha(p)==row['start_file_sha256'],prior_trace_duplicate=digest in prior,new_grid_duplicate=digest in seen,initial_radius_error=abs(rr-row['initial_radius']))
        startchecks.append(sc);seen.add(digest);swapgroups.setdefault((row['pair_id'],row['target_id']),{})[row['assignment']]=S
        if not all(sc[k] for k in ('other98_exactly_preserved','targets_exactly_match','start_coordinate_hash_matches','start_file_hash_matches')) or sc['prior_trace_duplicate'] or sc['new_grid_duplicate'] or sc['initial_radius_error']>2e-11:failures.append(sc)
        tasks.append(('start',i,str(p),row['initial_radius']))
        final=ROOT/row['final_coordinate_file'];C=np.load(final)
        if run.sha(final)!=row['final_file_sha256'] or run.csha(C)!=row['final_coordinate_sha256']:failures.append(dict(reason='final_hash',grid_id=i))
        tasks.append(('final',i,str(final),row['final_radius']))
        for stage in row['rounds']:
            p=ROOT/stage['coordinate_file'];C=np.load(p);rr,dd=run.true_cover(C);sc=dict(path=stage['coordinate_file'],file_hash_matches=run.sha(p)==stage['file_sha256'],coordinate_hash_matches=run.csha(C)==stage['coordinate_sha256'],radius_error=abs(rr-stage['radius']))
            stagechecks.append(sc)
            if not sc['file_hash_matches'] or not sc['coordinate_hash_matches'] or sc['radius_error']>2e-11:failures.append(sc)
    swapchecks=[]
    for (pair_id,target_id),group in swapgroups.items():
        selected=json.loads((ROOT/'selected_pairs.json').read_text())[pair_id];a,b=selected['pair'];ok=len(group)==2 and np.array_equal(group[0][[a,b]],group[1][[b,a]]) and run.csha(group[0])!=run.csha(group[1]);swapchecks.append(dict(pair_id=pair_id,target_id=target_id,both_assignments_present=len(group)==2,swapped_pair_centers_exact=ok))
        if not ok:failures.append(swapchecks[-1])
    cpus=sorted(os.sched_getaffinity(0));results=[]
    with ProcessPoolExecutor(max_workers=4,initializer=initializer,initargs=(cpus,)) as pool:
        for result in pool.map(replay,tasks,chunksize=1):
            results.append(result)
            if result['absolute_error']>2e-11:failures.append(result)
            if len(results)%32==0:
                run.write(ROOT/'audit_checkpoint.json',dict(completed=len(results),requested=len(tasks),maximum_error=max(x['absolute_error'] for x in results),seconds=time.monotonic()-tic,**search.pidinfo()));print('AUDIT_PROGRESS',len(results),len(tasks),flush=True)
    successful=[r for r in rows if 'final_radius' in r];best=min(successful,key=lambda r:r['final_radius']);C=np.load(ROOT/best['final_coordinate_file']);alignment=compare(baseline,C);run.write(ROOT/'best_baseline_alignment.json',dict(input_baseline_sha256=run.sha(ROOT/'inputs/baseline.npy'),final_file=best['final_coordinate_file'],final_file_sha256=best['final_file_sha256'],grid_id=best['grid_id'],comparison=alignment));np.save(ROOT/'best_final.npy',C);run.write(ROOT/'best_final.json',best)
    g,p,J,keys,sites,top=run.inventory(C);ac=np.flatnonzero(g.max()-g<1e-9);M=np.r_[J[ac].T,np.ones((1,len(ac)))];rhs=np.r_[np.zeros(C.size),1.];lam,res=nnls(M,rhs,maxiter=10000);inactive=[i for i in range(len(g)) if i not in set(ac)];run.write(ROOT/'fresh_best_diagnostic.json',dict(branch=f'paired_transport grid{best["grid_id"]} final only',input_file=best['final_coordinate_file'],input_file_sha256=best['final_file_sha256'],input_coordinate_sha256=best['final_coordinate_sha256'],radius=best['final_radius'],diagnostic=run.diagnose(C),NNLS_residual_inf=float(np.max(np.abs(M@lam-rhs))),NNLS_lambda_min=float(lam.min()),NNLS_lambda_max=float(lam.max()),NNLS_lambda_sum=float(lam.sum()),active_squared_gap_spread=float(np.ptp(g[ac])),minimum_inactive_squared_gap=float(g.max()-max(g[inactive])) if inactive else None,active_keys=[keys[i] for i in ac],lambda_values=lam,scope='Fresh numerical final diagnostic only, not stationarity, local-minimum, exact-cover, algebraic-root or global-optimality proof. Baseline stress statistics are not transferred to this final.'))
    obstacle_rows=[];baseinv=run.diagnose(baseline);bk=set(map(tuple,baseinv['active_keys']+baseinv['inactive_keys']))
    for r in successful:
        start=np.load(ROOT/'workers'/f'worker{r["worker"]:02d}'/'starts'/f'{r["grid_id"]:04d}.npy');sd=run.diagnose(start);sk=set(map(tuple,sd['active_keys']+sd['inactive_keys']));fd=r['final_diagnostic'];fk=set(map(tuple,fd['active_keys']+fd['inactive_keys']));st=r['rounds']
        obstacle_rows.append(dict(grid_id=r['grid_id'],pair_family=r['pair_family'],target_family=r['target_family'],assignment=r['assignment'],initial_radius=r['initial_radius'],final_radius=r['final_radius'],initial_inventory_added=sorted(sk-bk),initial_inventory_removed=sorted(bk-sk),final_inventory_added_to_seed=sorted(fk-sk),final_inventory_removed_from_seed=sorted(sk-fk),optimization_topology_events=sum(z['topology_changed'] for z in st),new_actual_binding_holes=sum(len(z['new_actual_binding_keys']) for z in st),maximum_actual_minus_frozen_objective=max(z['true_radius_minus_frozen_objective'] for z in st),minimum_actual_minus_frozen_objective=min(z['true_radius_minus_frozen_objective'] for z in st),final_feasible=st[-1]['feasible'],final_status=st[-1]['status'],optimizer_rounds=len(st)))
    obstacles=dict(starts_changing_actual_inventory=sum(bool(r['initial_inventory_added'] or r['initial_inventory_removed']) for r in obstacle_rows),paths_with_topology_event=sum(r['optimization_topology_events']>0 for r in obstacle_rows),paths_with_new_binding_holes=sum(r['new_actual_binding_holes']>0 for r in obstacle_rows),paths_with_actual_radius_gap_above_1e8=sum(r['maximum_actual_minus_frozen_objective']>1e-8 for r in obstacle_rows),paths_with_actual_below_frozen_by_1e8=sum(r['minimum_actual_minus_frozen_objective']<-1e-8 for r in obstacle_rows),final_not_feasible=sum(not r['final_feasible'] for r in obstacle_rows),rows=obstacle_rows,scope='Observed finite saved paths only; no global obstruction/infeasibility or optimality conclusion.')
    run.write(ROOT/'topology_obstacles.json',obstacles);material=[r for r in successful if r['final_radius']<run.THRESHOLD and run.REFERENCE-r['final_radius']>run.MATERIAL]
    report=dict(status='PASS' if not failures else 'FAIL',grid_budget=128,planned=len(grid),optimized=len(rows),skipped=sum(r['skip'] for r in grid),optimizer_exceptions=sum('error' in r for r in rows),all_triple_oracle_calls=len(results),independent_oracle_counts=dict(collections.Counter(r['kind'] for r in results)),maximum_all_triple_error=max(r['absolute_error'] for r in results),stage_array_count=len(stagechecks),stage_checks=stagechecks,start_checks=startchecks,swap_checks=swapchecks,prior_start_arrays=prior_arrays,prior_unique_hashes=len(prior),prior_read_errors=prior_errors,frozen_tree_checks=frozenchecks,runtime_source_checks=sourcechecks,worker_exitcodes=complete['workers_exitcode'],final_SLSQP_status_counts=dict(collections.Counter(r['final_status'] for r in obstacle_rows)),optimizer_round_counts=dict(collections.Counter(r['optimizer_rounds'] for r in obstacle_rows)),best_grid_id=best['grid_id'],best_pair=best['pair'],best_target_family=best['target_family'],best_final_radius=best['final_radius'],best_file=best['final_coordinate_file'],best_file_sha256=best['final_file_sha256'],material_candidate_count=len(material),best_full100_baseline_alignment=alignment,independent_results=results,failures=failures,seconds=time.monotonic()-tic,scope='Finite numerical pair relocation with independent full critical-point replay. No exact covering, global-infeasibility or optimality certificate.')
    run.write(ROOT/'replay_report.json',report);summary={k:report[k] for k in ('status','grid_budget','planned','optimized','skipped','optimizer_exceptions','all_triple_oracle_calls','independent_oracle_counts','maximum_all_triple_error','stage_array_count','worker_exitcodes','final_SLSQP_status_counts','optimizer_round_counts','best_grid_id','best_pair','best_target_family','best_final_radius','best_file','best_file_sha256','material_candidate_count','seconds','scope')};summary['baseline_file_sha256']=run.sha(ROOT/'inputs/baseline.npy');summary['best_matches_baseline_full100_geometry_and_inventories']=alignment['matches'];summary['original_three_frozen_trees_unchanged']=all(v['all_unchanged'] for v in frozenchecks.values());summary['topology_obstacles']={k:v for k,v in obstacles.items() if k not in ('rows','scope')};run.write(ROOT/'RESULT.json',summary)
    (ROOT/'README.txt').write_text('Paired remove-two/refill-two transport from de82 numerical baseline\n\nFresh all100 deletion debts and NNLS incident stress identify a finite34-donor pool. Sixty-four actual98 deletion configurations were whole-disk screened. Four shared-active-star adjacent pairs and four separated disjoint-star pairs were selected with actual pair debt, weak incident stress and mixed ranks. This is finite heuristic selection, not global donor optimality.\n\nEach pair uses eight distinct unordered target sets from actual98 deletion cavities, dynamically refreshed99 maxima after a first insertion, and actual boundary gaps chosen for input stress and spatial spread. Both label assignments are retained explicitly; their initial unordered geometric sets coincide, while bytewise100-center starts differ. All other98 coordinates are bitwise input-identical. No disk is assumed freely removable. Full start hashes are checked against prior start/seed/predictor traces, with explicit candidate-target rejects.\n\nAll100 centers are nooffset epigraph variables. At most two SLSQP rounds per start use dynamic full interior/boundary inventories and true whole-disk scoring. Four BLAS1 workers write separate files; the parent aggregator owns checkpoints and proposal freezes. Every actual optimizer status and error is recorded. Fresh numerical diagnostics and baseline grouping refer to explicit final SHA; grouping uses all100 assignment/rotation/reflection plus complete fresh active/inactive inventories.\n\nIndependent audit scores baseline, all100 single deletions,64 pair deletions,128 starts and128 finals using every triple circumcenter, every pair-circle crossing, antipodes and axes with true nearest-site distances. Saved stage arrays and all source/input/output hashes are checked separately. The three older frozen trees are checked unchanged. This is floating numerical verification, not a rigorous covering or global optimality proof.\n\nReproduce in a new writable copy: python3 -B source/pair_search.py; python3 -B source/audit_finalize.py. Existing checkpoints resume completed entries without repeating optimization.\n')
    run.write(ROOT/'audit_completion.json',dict(status=report['status'],audit_source_sha256=run.sha(Path(__file__)),independent_oracle_source_sha256=run.sha(ROOT/'source/independent_oracle.py'),seconds=time.monotonic()-tic,**search.pidinfo()))
    # Write manifest last; after this call the entire directory is a frozen result.
    manifest={str(p.relative_to(ROOT)):run.sha(p) for p in sorted(ROOT.rglob('*')) if p.is_file() and p.name not in ('MANIFEST_SHA256.json','audit.log')};run.write(ROOT/'MANIFEST_SHA256.json',manifest)
    print('AUDIT_DONE',json.dumps(run.clean(summary)),flush=True)

if __name__=='__main__':main()
