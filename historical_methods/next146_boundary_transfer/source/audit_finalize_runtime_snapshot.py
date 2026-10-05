# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Independent exhaustive all-triple replay of the completed donor grid."""
import grid_search as run
import os,itertools,json,time,hashlib,collections,multiprocessing as mp
from pathlib import Path
import numpy as np
from concurrent.futures import ProcessPoolExecutor
ROOT=run.ROOT;BASE=ROOT.parent
def all_triple_cover(C):
    C=np.asarray(C,dtype=float);n=len(C);T=np.array(list(itertools.combinations(range(n),3)),dtype=int);P=C[T];u=P[:,1]-P[:,0];v=P[:,2]-P[:,0];cross=u[:,0]*v[:,1]-u[:,1]*v[:,0];ok=np.abs(cross)>1e-14;a=P[ok,0];u=u[ok];v=v[ok];cross=cross[ok];uu=np.sum(u*u,axis=1);vv=np.sum(v*v,axis=1);q=a+np.column_stack((uu*v[:,1]-vv*u[:,1],u[:,0]*vv-v[:,0]*uu))/(2*cross[:,None]);q=q[np.sum(q*q,axis=1)<=1]
    ii,jj=np.triu_indices(n,1);u=C[jj]-C[ii];ss=np.sum(u*u,axis=1);h=(np.sum(C[jj]**2,axis=1)-np.sum(C[ii]**2,axis=1))/2;ok=(ss>1e-20)&(ss-h*h>=0);u=u[ok];ss=ss[ok];h=h[ok];z=np.sqrt(ss-h*h);perp=np.column_stack((-u[:,1],u[:,0]));mid=h[:,None]*u/ss[:,None];offset=z[:,None]*perp/ss[:,None];bp=np.r_[mid-offset,mid+offset];norms=np.linalg.norm(C,axis=1);nz=norms>1e-14;ap=-C[nz]/norms[nz,None];points=np.r_[q,bp,ap,np.eye(2),-np.eye(2)];best=-1.;worst=None
    for j in range(0,len(points),512):
        chunk=points[j:j+512];nearest=np.min(np.sum((chunk[:,None,:]-C[None,:,:])**2,axis=2),axis=1);at=int(np.argmax(nearest))
        if nearest[at]>best:best=float(nearest[at]);worst=chunk[at].tolist()
    return float(np.sqrt(best)),dict(sites=n,triples=len(T),nondegenerate_triples=len(a),inside_triple_candidates=len(q),pairs=len(ii),boundary_pair_candidates=len(bp),antipodes=len(ap),worst_point=worst,numerical_only=True,method='all-triples+all-pair-circle+antipodes+axes actual-nearest-distance, independent of Delaunay')
def initializer(cpus):
    if hasattr(os,'sched_setaffinity'):os.sched_setaffinity(0,{cpus[(mp.current_process()._identity[0]-1)%len(cpus)]})
def replay_job(task):
    kind,index,path,stored=task;C=np.load(path);R,detail=all_triple_cover(C);return dict(kind=kind,index=index,path=str(Path(path).relative_to(ROOT)),radius=R,stored_radius=stored,absolute_error=abs(R-stored),detail=detail)
def main():
    tic=time.monotonic();cp=json.loads((ROOT/'checkpoint.json').read_text());complete=json.loads((ROOT/'completion.json').read_text());rows=cp['records'];grid=json.loads((ROOT/'grid.json').read_text());donors=json.loads((ROOT/'donor_inventory.json').read_text());baseline=np.load(ROOT/'inputs/baseline.npy');manifest=json.loads((BASE/'rattler_contact_transfer/MANIFEST_SHA256.json').read_text());oldchecks={p:run.sha(BASE/'rattler_contact_transfer'/p)==h for p,h in manifest.items()}; old147=json.loads((BASE/'best147_boundary_transfer/MANIFEST_SHA256.json').read_text()); old147checks={p:run.sha(BASE/'best147_boundary_transfer'/p)==h for p,h in old147.items()};runtime=json.loads((ROOT/'runtime_provenance.json').read_text());sourcechecks={p:run.sha(ROOT/p)==h for p,h in runtime['source_sha256'].items()};failures=[];startchecks=[];stagechecks=[];tasks=[];seen=set();prior83=set()
    for sub in ('','max_stress_supplement','boundary_followup'):
        for p in (BASE/'rattler_contact_transfer'/sub/'starts').glob('*.npy'):prior83.add(run.csha(np.load(p)))
    if complete['completed']!=cp['requested'] or len(rows)!=cp['requested'] or len({r['grid_id'] for r in rows})!=cp['requested']:failures.append(dict(reason='grid_count_or_duplicate_grid_ids'))
    for path in (BASE/'best147_boundary_transfer/workers').glob('worker*/starts/*.npy'):prior83.add(run.csha(np.load(path)))
    for d in donors:
        path=ROOT/'inputs'/f'delete_{d["index"]:03d}.npy';np.save(path,np.delete(baseline,d['index'],axis=0));tasks.append(('donor_deletion',d['index'],str(path),d['without_donor_radius']))
    for row in rows:
        if 'error' in row:failures.append(dict(reason='optimizer_exception',grid=row['grid_id'],error=row['error']));continue
        i=row['grid_id'];folder=ROOT/'workers'/f'worker{row["worker"]:02d}';path=folder/'starts'/f'{i:04d}.npy';S=np.load(path);digest=run.csha(S);other99=np.array_equal(S[np.arange(100)!=row['donor']],baseline[np.arange(100)!=row['donor']]);target=np.array_equal(S[row['donor']],row['target']['point']);R,detail=run.true_cover(S);sc=dict(grid_id=i,start_sha_matches=digest==row['start_coordinate_sha256'],other99_exactly_preserved=other99,target_exactly_matches=target,old83_duplicate=digest in prior83,newgrid_duplicate=digest in seen,initial_radius=R,initial_score_error=abs(R-row['initial_radius']))
        startchecks.append(sc);seen.add(digest)
        if not other99 or not target or not sc['start_sha_matches'] or sc['old83_duplicate'] or sc['newgrid_duplicate'] or sc['initial_score_error']>2e-11:failures.append(sc)
        final=ROOT/row['final_coordinate_file'];C=np.load(final)
        if run.sha(final)!=row['final_file_sha256'] or run.csha(C)!=row['final_coordinate_sha256']:failures.append(dict(reason='final_hash',grid_id=i))
        tasks.append(('final',i,str(final),row['final_radius']))
        for stage in row['rounds']:
            p=ROOT/stage['coordinate_file'];C=np.load(p);R,detail=run.true_cover(C);sc=dict(path=stage['coordinate_file'],file_sha_matches=run.sha(p)==stage['file_sha256'],coordinate_sha_matches=run.csha(C)==stage['coordinate_sha256'],radius_error=abs(R-stage['radius']))
            stagechecks.append(sc)
            if not sc['file_sha_matches'] or not sc['coordinate_sha_matches'] or sc['radius_error']>2e-11:failures.append(sc)
    if not all(oldchecks.values()):failures.append(dict(reason='original429_frozen_files_changed'))
    if not all(old147checks.values()):failures.append(dict(reason='original1467_frozen_files_changed'))
    if not all(sourcechecks.values()):failures.append(dict(reason='runtime_source_changed'))
    cpus=sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else [0,1,2,3];results=[]
    with ProcessPoolExecutor(max_workers=4,initializer=initializer,initargs=(cpus,)) as pool:
        for result in pool.map(replay_job,tasks,chunksize=1):
            results.append(result)
            if result['absolute_error']>2e-11:failures.append(result)
            if len(results)%32==0:
                run.write(ROOT/'audit_checkpoint.json',dict(completed=len(results),requested=len(tasks),maximum_error=max(x['absolute_error'] for x in results),seconds=time.monotonic()-tic));print('AUDIT_PROGRESS',len(results),len(tasks),flush=True)
    best=min((r for r in rows if 'final_radius' in r),key=lambda r:r['final_radius']);new=[r for r in rows if r.get('final_radius',1)<run.THRESHOLD and run.REFERENCE-r['final_radius']>run.MATERIAL]
    obstacles=dict(paths_with_topology_event=sum(any(z['topology_changed'] for z in r.get('rounds',[])) for r in rows),paths_with_new_binding_holes=sum(any(z.get('new_actual_binding_keys') for z in r.get('rounds',[])) for r in rows),final_unresolved_actual_radius_gap=sum(r['rounds'][-1].get('true_radius_minus_frozen_objective',0)>1e-8 for r in rows if r.get('rounds')),final_not_feasible_count=sum(not r['rounds'][-1]['feasible'] for r in rows if r.get('rounds')),maximum_actual_radius_minus_frozen_objective=max(z.get('true_radius_minus_frozen_objective',0) for r in rows for z in r.get('rounds',[])),scope='Actual inventory/epigraph mismatch after at most two rounds is a numerical finite-search obstacle, not a global topology barrier or infeasibility theorem.');run.write(ROOT/'topology_obstacles.json',obstacles)
    report=dict(status='PASS' if not failures else 'FAIL',grid_budget=len(grid),optimized=len(rows),skipped=sum(r['skip'] for r in grid),exceptions=sum('error' in r for r in rows),all100_deletions_independently_replayed=sum(r['kind']=='donor_deletion' for r in results),all_final_arrays_independently_replayed=sum(r['kind']=='final' for r in results),all_triple_oracle_calls=len(results),maximum_all_triple_score_error=max(x['absolute_error'] for x in results),stage_array_count=len(stagechecks),stage_checks=stagechecks,start_checks=startchecks,original429_frozen_checks=oldchecks,original1467_frozen_checks=old147checks,runtime_source_checks=sourcechecks,independent_results=results,worker_exitcodes=complete['workers_exitcode'],final_SLSQP_status_counts=dict(collections.Counter(r['rounds'][-1]['status'] for r in rows if 'rounds' in r)),best_grid_id=best['grid_id'],best_donor=best['donor'],best_boundary=best['boundary'],best_final_radius=best['final_radius'],best_file=best['final_coordinate_file'],best_file_sha256=best['final_file_sha256'],material_candidate_count=len(new),topology_obstacles=obstacles,known99_stops=sum(any(z.get('known99',{}).get('match',False) for z in r.get('rounds',[])) for r in rows),failures=failures,seconds=time.monotonic()-tic,scope='Finite floating numerical construction search, whole-disk critical-point replay and file audit. No exact covering or global-infeasibility/optimality certificate.')
    run.write(ROOT/'replay_report.json',report)
    outcome={k:report[k] for k in ['status','grid_budget','optimized','skipped','exceptions','all100_deletions_independently_replayed','all_final_arrays_independently_replayed','maximum_all_triple_score_error','stage_array_count','best_grid_id','best_donor','best_target','best_final_radius','best_file','best_file_sha256','material_candidate_count','final_SLSQP_status_counts','known99_stops','seconds','scope']};run.write(ROOT/'RESULT.json',outcome)
    (ROOT/'README.txt').write_text('New grid25/de82 146-active donor transfer experiment\n\nThe prior429 and1467 file trees remain unchanged. Baseline is grid25 SHAde82a12038c6f9a705b7a985c7e25c53c8e94f5c506946aba9a6c6f03a02852d, radius0.11544695488558258. All preparation NNLS/deletion/response statistics belong to this baseline only.\n\nThe selected12 donors comprise5 low-debt,4 low-incident-stress and3 mixed-rank centers. All100 donors were deleted and whole-disk scored independently before selection; no nonbinding disk is assumed removable. The core grid covers all32 actual active boundary maxholes, plus3 nearest nonincident active internal maxholes for each of the5 lowest-debt donors. There are399 total planned entries; skips and actual count are explicit in grid.json/RESULT.json. Every start preserves the other99 centers bitwise and checks its hash against all prior465 starts and this grid.\n\nAll100 centers are zero-offset epigraph variables. Four BLAS1 workers use separate output files. Acceptance uses the complete numerical whole-disk critical-point oracle, never sampled coverage. Stage inventories, newly exposed binding holes and frozen-objective versus actual-radius gaps are retained in topology_obstacles.json and each stage record. Limits of one or two optimizer rounds are stated and status8 cases remain explicit.\n\nIndependent audit covers all100 deletion configurations and every final by every triple circumcenter, every pair-circle crossing, antipodes and axes. This floating replay is independent of Delaunay and is not an exact covering or optimality proof. Exact certification is performed by the parent.\n\nReproduction: python3 -B source/grid_search.py --prepare-only; python3 -B source/grid_search.py; python3 -B source/audit_finalize.py; python3 -B source/freeze_result.py.\n')
    print('AUDIT_DONE',json.dumps(outcome),flush=True)
if __name__=='__main__':main()
