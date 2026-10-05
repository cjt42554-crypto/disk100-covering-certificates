# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Finite full100 restoration of verified B32 harmonic metric realizations."""
import engine as run
import os,json,time,traceback,multiprocessing as mp,argparse,queue as qm,collections
from pathlib import Path
import numpy as np
from metric_validation import validate
ROOT=run.ROOT;BASE=ROOT.parent;HANDOFF=ROOT/'inputs/catalog_handoff'
FROZEN_NAMES=('rattler_contact_transfer','best147_boundary_transfer','next146_boundary_transfer',
    'paired_transport_146','topology_patch_hybrid_146','event_geometry_146',
    'fresh_rattler_transfer_99','fresh_rattler_geometry_bridge','cavity_repair_transfer_99',
    'fresh99_second_order_diagnostic')

def unordered_hash(C):
    C=np.asarray(C,dtype='<f8');return run.csha(C[np.lexsort((C[:,1],C[:,0]))])

def pidinfo():
    return dict(pid=os.getpid(),namespace_pid_status=[s for s in Path('/proc/self/status').read_text().splitlines()if s.startswith('NSpid:')],affinity=sorted(os.sched_getaffinity(0)))

def frozen_checks():
    out={}
    for n in FROZEN_NAMES:
        p=BASE/n/'MANIFEST_SHA256.json';m=json.loads(p.read_text())
        out[n]=dict(count=len(m),manifest_sha256=run.sha(p),unchanged=all(run.sha(BASE/n/q)==h for q,h in m.items()))
    m=json.loads((HANDOFF/'MANIFEST_SHA256.json').read_text())
    out['catalog_handoff']=dict(count=len(m),manifest_sha256=run.sha(HANDOFF/'MANIFEST_SHA256.json'),unchanged=all(run.sha(HANDOFF/q)==h for q,h in m.items()),owner_unchanged=all(run.sha(BASE/'b32_catalog_search/frozen_handoff'/q)==h for q,h in m.items()))
    return out

def prior_hashes():
    """Old trace hash cache plus every subsequent sealed mathematical array.

    Designated catalog harmonic preparations are excluded because they have
    never undergone original-problem optimization; this exception is explicit.
    """
    old=json.loads((BASE/'cavity_repair_transfer_99/prior_ordered_hashes.json').read_text())
    oldu=json.loads((BASE/'cavity_repair_transfer_99/prior_unordered_hashes.json').read_text())
    n=0;errors=[]
    for name in ('cavity_repair_transfer_99','fresh_rattler_geometry_bridge','fresh99_second_order_diagnostic'):
        m=json.loads((BASE/name/'MANIFEST_SHA256.json').read_text())
        for relative,h in m.items():
            p=BASE/name/relative
            if p.suffix not in ('.npy','.npz'):continue
            try:
                if p.suffix=='.npy':values=[('',np.load(p,allow_pickle=False))]
                else:
                    with np.load(p,allow_pickle=False)as z:values=[(k,z[k])for k in z.files]
                for k,C in values:
                    if C.shape==(201,):C=C[:-1].reshape(100,2)
                    elif C.shape==(200,):C=C.reshape(100,2)
                    if C.shape!=(100,2):continue
                    label=str(p.relative_to(BASE))+(':'+k if k else '')
                    old.setdefault(run.csha(C),[]).append(label);oldu.setdefault(unordered_hash(C),[]).append(label);n+=1
            except Exception as e:errors.append(dict(path=str(p),error=str(e)))
    return old,oldu,n,errors

def prepare():
    tic=time.monotonic();frozen=frozen_checks();assert all(v['unchanged']for v in frozen.values());assert frozen['catalog_handoff']['owner_unchanged']
    ready=json.loads((HANDOFF/'READY.json').read_text());assert ready['status']=='SEALED';assert ready['current_independent_index']==10
    cat=json.loads((HANDOFF/'catalog.json').read_text());graphs={g['independent_index']:g for g in cat['graphs']};ledger=json.loads((ROOT/'screening/handoff_seed_validation.json').read_text())
    history=json.loads((ROOT/'inputs/history_snapshot/READY.json').read_text());assert history['status']=='READY_FROZEN'
    old,oldu,nnew,errors=prior_hashes();assert not errors
    # The larger historical coordinate inventory supplies ordered binary64 SHA
    # records in the same convention, without reading arrays being written.
    hist=json.loads((ROOT/'inputs/history_snapshot/retained_coordinate_records.json').read_text())
    histrows=hist if isinstance(hist,list)else hist['records']
    history_count=0
    for r in histrows:
        cs=r.get('coordinate_sha256',r.get('coordinate_hash'))
        if cs:
            old.setdefault(cs,[]).append('historical_snapshot:'+str(r.get('path',r.get('file',r.get('index','unknown')))));history_count+=1
    run.write(ROOT/'prior_ordered_hashes.json',old);run.write(ROOT/'prior_unordered_hashes.json',oldu)
    historysummary=json.loads((ROOT/'inputs/history_snapshot/historical_matching_summary.json').read_text())
    previously_seen=set(map(int,historysummary['independent_catalog_matching_unique_coordinate_counts']))
    historical_matches=json.loads((ROOT/'inputs/history_snapshot/historical_b32_profile_matches.json').read_text())
    records=sorted(ledger['records'],key=lambda r:(r['independent_index']==10,r['independent_index']in previously_seen,r['independent_index'],r['boundary_mode'],r['core_edge_weight']))
    grid=[];seen=set();seenu=set();executable=0
    for i,item in enumerate(records):
        idx=item['independent_index'];path=HANDOFF/item['coordinate_relative_file'];assert run.sha(path)==item['file_sha256'];C=np.load(path);cs=run.csha(C);us=unordered_hash(C);assert cs==item['coordinate_sha256'];v=validate(C,graphs[idx]);v.pop('actual_graph',None);R,d=run.true_cover(C);skip=None
        norm=float(np.linalg.norm(C,axis=1).max())
        samegraph=[r for r in historical_matches if r['fingerprint']==graphs[idx]['fingerprint']]
        possibleunordered=[r for r in samegraph if abs(r['max_center_norm']-norm)<1e-10]
        # An unordered equality must have the same strict, label-invariant
        # restricted-DT canonical graph and the same maximum center norm.
        # The complete historical B32-profile extraction has484 arrays;
        # these harmonic radii(.94) exclude all old same-graph metrics(>.98).
        # This supplies the larger inventory check without reading other projects.
        assert not possibleunordered
        historical_unordered=dict(method='necessary canonical strict-DT graph and maximum center norm equality',
            samegraph_historical_unique_arrays=len(samegraph),possible_unordered_duplicates=len(possibleunordered),
            minimum_maxnorm_separation=min((abs(r['max_center_norm']-norm)for r in samegraph),default=None),
            scope='Finite34868-array retained history only; geometric invariants supplement prior sorted-coordinate hashes.')
        if idx==10:skip='current_best_graph_reference_only_not_repeated'
        elif not v['target_metric_realized']:skip='harmonic_target_DT_mismatch_not_catalog_family_realization'
        elif cs in old:skip='duplicate_prior_ordered_coordinates'
        elif us in oldu:skip='duplicate_prior_unordered_geometry'
        elif cs in seen or us in seenu:skip='duplicate_new_batch_geometry'
        worker=executable%4
        row=dict(grid_id=i,independent_index=idx,target_fingerprint=graphs[idx]['fingerprint'],target_graph=graphs[idx],family_name_status='independent index; parent exact catalog bijection pending',
            seed_method=item['method'],boundary_mode=item['boundary_mode'],core_edge_weight=item['core_edge_weight'],owner_seed_file=item['coordinate_relative_file'],owner_file_sha256=item['file_sha256'],
            initial_coordinate_sha256=cs,initial_unordered_geometry_sha256=us,preoptimization_target_metric=v,preoptimization_radius=R,preoptimization_whole_disk_score=d,
            target_graph_seen_in_retained_history=idx in previously_seen,historical_unordered_geometry_check=historical_unordered,history_scope=historysummary.get('scope'),worker=worker,skip=skip is not None,skip_reason=skip,
            offsets=0,all_center_variables=100,original_problem_constraints_only=True)
        if skip is None:
            folder=ROOT/'workers'/f'worker{worker:02d}'
            for n in ('starts','finals','stages','proposals'):(folder/n).mkdir(parents=True,exist_ok=True)
            p=folder/'starts'/f'{i:04d}.npy';np.save(p,C);row['start_file_sha256']=run.sha(p);seen.add(cs);seenu.add(us);executable+=1
        grid.append(row)
    assert executable<=96;assert all(r['preoptimization_target_metric']['target_metric_realized']for r in grid if not r['skip']);assert not any(r['independent_index']==10 and not r['skip']for r in grid)
    checks=dict(stage='prepared',budget=96,owner_planned=len(grid),executable=executable,skips=dict(collections.Counter(r['skip_reason']for r in grid if r['skip'])),
        graph_start_counts=dict(collections.Counter(r['independent_index']for r in grid if not r['skip'])),prior_ordered_unique=len(old),prior_unordered_unique=len(oldu),additional_sealed_array_count=nnew,
        historical_ordered_records_added=history_count,prior_read_errors=errors,previously_seen_target_indices=sorted(previously_seen),historical_snapshot_ready_sha256=run.sha(ROOT/'inputs/history_snapshot/READY.json'),
        source_handoff_manifest_sha256=run.sha(HANDOFF/'MANIFEST_SHA256.json'),frozen_trees=frozen,seconds=time.monotonic()-tic,**pidinfo(),
        scope='Finite heuristic16-graph metric seeds,87 planned original-problem starts; no completeness or global-optimality assertion. Designated harmonic preparations were never previously optimized. Same historical graph with a novel metric is explicitly labeled.')
    run.write(ROOT/'input_checks.json',checks);run.write(ROOT/'grid.json',grid);run.write(ROOT/'status.json',checks);print('METRIC_PREPARED',json.dumps(run.clean(checks)),flush=True);return grid

def worker_main(worker,rows,queue,cpus):
    os.sched_setaffinity(0,{cpus[worker%len(cpus)]});folder=ROOT/'workers'/f'worker{worker:02d}';records=[];tic=time.monotonic();run.write(folder/'runtime.json',dict(worker=worker,BLAS_threads=1,**pidinfo()))
    if(folder/'checkpoint.json').exists():records=json.loads((folder/'checkpoint.json').read_text())['records']
    done={r['grid_id']:r for r in records}
    with open(folder/'records.jsonl','a',buffering=1)as fp:
        for row in rows:
            if row['grid_id']in done:queue.put(dict(type='record',record=done[row['grid_id']]));continue
            tt=time.monotonic();i=row['grid_id']
            try:
                C=np.load(folder/'starts'/f'{i:04d}.npy');R,C,initial,rounds=run.optimize(C,folder,row);p=folder/'finals'/f'{i:04d}.npy';np.save(p,C);R,d=run.true_cover(C);v=validate(C,row['target_graph']);v.pop('actual_graph',None)
                result=dict(row,initial_radius=float(np.sqrt(initial['maximum_squared_distance'])),initial_whole_disk_score=initial,final_radius=R,final_whole_disk_score=d,final_diagnostic=run.diagnose(C),
                    final_target_and_actual_metric=v,final_coordinate_file=str(p.relative_to(ROOT)),final_file_sha256=run.sha(p),final_coordinate_sha256=run.csha(C),rounds=rounds,seconds=time.monotonic()-tt,numerical_only=True)
            except Exception as e:result=dict(row,error=str(e),traceback=traceback.format_exc(),seconds=time.monotonic()-tt)
            records.append(result);fp.write(json.dumps(run.clean(result),allow_nan=False)+'\n');run.write(folder/'checkpoint.json',dict(worker=worker,completed=len(records),requested=len(rows),records=records,best_radius=min((r.get('final_radius',1)for r in records),default=1),seconds=time.monotonic()-tic,BLAS_threads=1,**pidinfo()));queue.put(dict(type='record',record=result))
    queue.put(dict(type='done',worker=worker))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--prepare-only',action='store_true');args=ap.parse_args();grid=json.loads((ROOT/'grid.json').read_text())if(ROOT/'grid.json').exists()else prepare()
    if args.prepare_only:return
    rows=[r for r in grid if not r['skip']];cpus=sorted(os.sched_getaffinity(0));runtime=dict(source_sha256={str(p.relative_to(ROOT)):run.sha(p)for p in sorted((ROOT/'source').glob('*.py'))},worker_cpu_assignments=[cpus[w%len(cpus)]for w in range(4)],workers=4,BLAS_threads_per_worker=1,**pidinfo());run.write(ROOT/'runtime_provenance.json',runtime)
    ctx=mp.get_context('spawn');queue=ctx.Queue();procs=[ctx.Process(target=worker_main,args=(w,[r for r in rows if r['worker']==w],queue,cpus))for w in range(4)]
    for p in procs:p.start()
    runtime['worker_pids']=[p.pid for p in procs];run.write(ROOT/'runtime_provenance.json',runtime);print('METRIC_WORKERS',json.dumps(runtime),flush=True)
    records=[];done=set();best=float('inf');tic=time.monotonic()
    while len(done)<4:
        try:item=queue.get(timeout=20)
        except qm.Empty:
            dead=[p for p in procs if p.exitcode not in(None,0)]
            if dead:raise RuntimeError('worker failed: '+str([(p.pid,p.exitcode)for p in dead]))
            continue
        if item['type']=='done':done.add(item['worker']);continue
        r=item['record'];records.append(r);R=r.get('final_radius',1)
        if R<run.THRESHOLD and run.REFERENCE-R>run.MATERIAL:
            p=ROOT/'frozen_proposals'/f'grid{r["grid_id"]:04d}_R{R:.17f}.npy';np.save(p,np.load(ROOT/r['final_coordinate_file']));run.write(p.with_suffix('.json'),dict(r,frozen_file_sha256=run.sha(p),frozen_coordinate_sha256=run.csha(np.load(p))));print('FROZEN_PROPOSAL',r['grid_id'],R,str(p),run.sha(p),flush=True)
        if R<best:
            best=R;np.save(ROOT/'best.npy',np.load(ROOT/r['final_coordinate_file']));run.write(ROOT/'best.json',r);print('METRIC_BEST',r['grid_id'],r['independent_index'],best,flush=True)
        run.write(ROOT/'checkpoint.json',dict(completed=len(records),requested=len(rows),budget=96,skipped=len(grid)-len(rows),workers_completed=sorted(done),best_radius=best,records=records,seconds=time.monotonic()-tic))
        if len(records)%8==0:print('METRIC_PROGRESS',len(records),len(rows),'best',best,'seconds',time.monotonic()-tic,flush=True)
    for p in procs:p.join()
    run.write(ROOT/'completion.json',dict(completed=len(records),requested=len(rows),best_radius=best,errors=sum('error'in r for r in records),worker_exitcodes=[p.exitcode for p in procs],seconds=time.monotonic()-tic))
    print('METRIC_DONE',len(records),best,'seconds',time.monotonic()-tic,flush=True)

if __name__=='__main__':main()
