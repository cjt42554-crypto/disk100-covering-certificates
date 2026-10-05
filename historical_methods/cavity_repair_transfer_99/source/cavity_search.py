# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Original-problem restoration of immutable OR-cavity/QCQP repair inputs."""
import engine as run
import os,sys,json,itertools,time,traceback,multiprocessing as mp,argparse,queue as queue_mod,collections
from pathlib import Path
import numpy as np
ROOT=run.ROOT;BASE=ROOT.parent;OWNER=BASE/'fresh_rattler_geometry_bridge'

def unordered_hash(C):
    C=np.asarray(C,dtype='<f8');return run.csha(C[np.lexsort((C[:,1],C[:,0]))])
def frozen_checks():
    out={}
    for name in('rattler_contact_transfer','best147_boundary_transfer','next146_boundary_transfer','paired_transport_146','topology_patch_hybrid_146','event_geometry_146','fresh_rattler_transfer_99','fresh_rattler_geometry_bridge'):
        mf=json.loads((BASE/name/'MANIFEST_SHA256.json').read_text());out[name]=dict(count=len(mf),manifest_sha256=run.sha(BASE/name/'MANIFEST_SHA256.json'),unchanged=all(run.sha(BASE/name/p)==h for p,h in mf.items()))
    return out

def prepare():
    tic=time.monotonic();run.write(ROOT/'status.json',dict(stage='48 immutable OR-repair whole-disk starts',budget=64,**pidinfo()));handoff=json.loads((ROOT/'inputs/owner_HANDOFF.json').read_text());repairs=json.loads((ROOT/'inputs/owner_REPAIR_TARGETS.json').read_text());manifest=json.loads((ROOT/'inputs/owner_MANIFEST_SHA256.json').read_text());assert handoff['candidate_count']==48;assert run.sha(ROOT/'inputs/owner_HANDOFF.json')==manifest['HANDOFF.json'];assert run.sha(ROOT/'inputs/owner_REPAIR_TARGETS.json')==manifest['REPAIR_TARGETS.json'];records={r['name']:r for r in repairs['records']};sources={};sourcechecks=[]
    for frame in handoff['source_frames']:
        i=frame['source'];C=np.load(ROOT/'inputs'/f'source{i:02d}_raw.npy');assert run.sha(ROOT/'inputs'/f'source{i:02d}_raw.npy')==frame['raw_source_sha256'];V=C.copy()
        if frame['reflection_before_rotation']:V[:,1]*=-1
        aligned=V[frame['permutation']]@np.array(frame['row_coordinate_rotation_matrix']);target=BASE/frame['aligned_source'];sc=dict(source=i,raw_sha256=frame['raw_source_sha256'],aligned_sha256=frame['aligned_source_sha256'],aligned_transform_reproduces_bitwise=np.array_equal(aligned,np.load(target)),aligned_original_sha_matches=run.sha(target)==frame['aligned_source_sha256']);assert sc['aligned_transform_reproduces_bitwise']and sc['aligned_original_sha_matches'];sourcechecks.append(sc);sources[i]=C
    old,oldu,oldarrays,olderrors=prior_hashes();assert not olderrors;run.write(ROOT/'prior_ordered_hashes.json',old);run.write(ROOT/'prior_unordered_hashes.json',oldu);grid=[];seen=set();seenu=set()
    for index,item in enumerate(handoff['priority']):
        name=item['name'];record=records[name];source=item['source'];path=OWNER/item['file'];assert run.sha(path)==item['file_sha256']==manifest[item['file']];S=np.load(path);assert run.csha(S)==item['coordinate_sha256'];C=sources[source];free=record['free_label'];patch=record['patch_labels'];outside=sorted(set(range(100))-set(patch)-{free});outsideok=np.array_equal(S[outside],C[outside]);freeok=np.array_equal(S[free],np.array(record['new_free_center']));assert outsideok and freeok
        R,detail=run.true_cover(S);cs=run.csha(S);us=unordered_hash(S);reason=None
        if cs in old:reason='duplicate_prior_original_problem_input_sha'
        elif us in oldu:reason='duplicate_prior_original_problem_unordered_geometry'
        elif cs in seen or us in seenu:reason='duplicate_new_batch_ordered_or_unordered_geometry'
        elif not np.isfinite(S).all()or np.linalg.norm(S,axis=1).max()>1+1e-10:reason='invalid_circle_geometry'
        worker=index%4;row=dict(grid_id=index,source_index=source,parent_raw_source_file=f'inputs/source{source:02d}_raw.npy',parent_raw_source_file_sha256=item['source_sha256'],handoff_name=name,handoff_original_file=item['file'],handoff_input_file_sha256=item['file_sha256'],handoff_coordinate_sha256=item['coordinate_sha256'],free_center_index=free,patch_indices=patch,all_modified_indices=sorted(set(patch)|{free}),outside_patch_bitwise_preserved=outsideok,new_free_center=record['new_free_center'],owner_QCQP_success=record['QCQP_success'],owner_QCQP_status=record['QCQP_status'],owner_LP_status=record['LP_status'],owner_QCQP_record=record,OR_semantics=handoff['OR_semantics'],preoptimization_radius=R,preoptimization_whole_disk_score=detail,preoptimization_owner_radius_error=abs(R-item['actual_start_radius']),preoptimization_diagnostic=run.diagnose(S),start_coordinate_sha256=cs,start_unordered_geometry_sha256=us,worker=worker,skip=reason is not None,skip_reason=reason,offsets=0,all_center_variables=100,original_problem_constraints_only=True)
        assert row['preoptimization_owner_radius_error']<2e-11
        if reason is None:
            folder=ROOT/'workers'/f'worker{worker:02d}'
            for d in('starts','finals','stages','proposals'):(folder/d).mkdir(parents=True,exist_ok=True)
            p=folder/'starts'/f'{index:04d}.npy';np.save(p,S);row['start_file_sha256']=run.sha(p);seen.add(cs);seenu.add(us)
        grid.append(row)
    frozen=frozen_checks();assert all(v['unchanged']for v in frozen.values());checks=dict(stage='prepared',budget=64,planned=len(grid),executable=sum(not r['skip']for r in grid),skipped=sum(r['skip']for r in grid),source_counts=dict(collections.Counter(r['source_index']for r in grid)),QCQP_proxy_success=sum(r['owner_QCQP_success']for r in grid),QCQP_proxy_failure=sum(not r['owner_QCQP_success']for r in grid),source_frame_checks=sourcechecks,old_initial_arrays=oldarrays,old_unique_ordered=len(old),old_unique_unordered=len(oldu),old_read_errors=olderrors,excluded_preparation_namespace='fresh_rattler_geometry_bridge/warmstarts: these designated inputs had only local QCQP preparation and no original full100 optimizer; their old13798 novelty proof and this larger prior-original-problem hash check are retained.',frozen_trees=frozen,OR_semantics=handoff['OR_semantics'],seconds=time.monotonic()-tic,**pidinfo(),scope='Whole100 restoration of OR-owner branch cavity repairs; proxy/QCQP statuses are preserved and do not imply full covering or global improvement.');run.write(ROOT/'input_checks.json',checks);run.write(ROOT/'grid.json',grid);run.write(ROOT/'status.json',checks);print('PREPARED',json.dumps(run.clean(checks)),flush=True);return grid


def pidinfo():return dict(pid=os.getpid(),namespace_pid_status=[s for s in Path('/proc/self/status').read_text().splitlines() if s.startswith('NSpid:')],affinity=sorted(os.sched_getaffinity(0)))

def prior_hashes():
    ordered={};unordered={};arrays=0;errors=[];excluded={'.git','.aws','.agents','.codex','github_archive_staging','cavity_repair_transfer_99','fresh_rattler_geometry_bridge','__pycache__'}
    frozenmath={'rattler_contact_transfer','best147_boundary_transfer','next146_boundary_transfer','paired_transport_146','topology_patch_hybrid_146','event_geometry_146','fresh_rattler_transfer_99'}
    for folder,dirs,files in os.walk(BASE):
        dirs[:]=[d for d in dirs if d not in excluded]
        for name in files:
            p=Path(folder)/name;parts=p.relative_to(BASE).parts
            if p.suffix not in ('.npy','.npz'):continue
            trace=any(s in (name+' '+str(p.parent.relative_to(BASE))).lower() for s in ('start','seed','predictor','initial','paired_event_search'))
            if not trace and parts[0] not in frozenmath:continue
            try:
                if p.suffix=='.npy':values=[('',np.load(p,allow_pickle=False))]
                else:
                    with np.load(p,allow_pickle=False) as z:values=[(k,z[k]) for k in z.files]
                for key,a in values:
                    if a.shape==(201,):a=a[:-1].reshape(100,2)
                    elif a.shape==(200,):a=a.reshape(100,2)
                    if a.shape==(100,2):
                        label=str(p.relative_to(BASE))+(':'+key if key else '');ordered.setdefault(run.csha(a),[]).append(label);unordered.setdefault(unordered_hash(a),[]).append(label);arrays+=1
            except Exception as e:errors.append(dict(path=str(p.relative_to(BASE)),error=str(e)))
    return ordered,unordered,arrays,errors

def worker_main(worker,rows,queue,cpus):
    os.sched_setaffinity(0,{cpus[worker%len(cpus)]});folder=ROOT/'workers'/f'worker{worker:02d}';records=[];tic=time.monotonic();run.write(folder/'runtime.json',dict(worker=worker,BLAS_threads=1,**pidinfo()))
    if (folder/'checkpoint.json').exists():records=json.loads((folder/'checkpoint.json').read_text())['records']
    done={r['grid_id']:r for r in records}
    with open(folder/'records.jsonl','a',buffering=1) as fp:
        for row in rows:
            if row['grid_id'] in done:queue.put(dict(type='record',record=done[row['grid_id']]));continue
            tt=time.monotonic();i=row['grid_id']
            try:
                start=np.load(folder/'starts'/f'{i:04d}.npy');R,C,initial,rounds=run.optimize(start,folder,row);path=folder/'finals'/f'{i:04d}.npy';np.save(path,C);R,detail=run.true_cover(C)
                result=dict(row,initial_radius=float(np.sqrt(initial['maximum_squared_distance'])),initial_whole_disk_score=initial,final_radius=R,final_whole_disk_score=detail,final_diagnostic=run.diagnose(C),final_coordinate_file=str(path.relative_to(ROOT)),final_file_sha256=run.sha(path),final_coordinate_sha256=run.csha(C),rounds=rounds,seconds=time.monotonic()-tt,numerical_only=True)
            except Exception as e:result=dict(row,error=str(e),traceback=traceback.format_exc(),seconds=time.monotonic()-tt)
            records.append(result);fp.write(json.dumps(run.clean(result),ensure_ascii=False,allow_nan=False)+'\n');run.write(folder/'checkpoint.json',dict(worker=worker,completed=len(records),requested=len(rows),records=records,best_radius=min((r.get('final_radius',1) for r in records),default=1),seconds=time.monotonic()-tic,BLAS_threads=1,**pidinfo()));queue.put(dict(type='record',record=result))
    queue.put(dict(type='done',worker=worker))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--prepare-only',action='store_true');args=ap.parse_args();grid=json.loads((ROOT/'grid.json').read_text()) if (ROOT/'grid.json').exists() else prepare()
    if args.prepare_only:return
    executed=[r for r in grid if not r['skip']];cpus=sorted(os.sched_getaffinity(0));runtime=dict(source_sha256={str(p.relative_to(ROOT)):run.sha(p) for p in sorted((ROOT/'source').glob('*.py'))},worker_cpu_assignments=[cpus[w%len(cpus)] for w in range(4)],workers=4,BLAS_threads_per_worker=1,**pidinfo());run.write(ROOT/'runtime_provenance.json',runtime)
    ctx=mp.get_context('spawn');queue=ctx.Queue();processes=[ctx.Process(target=worker_main,args=(w,[r for r in executed if r['worker']==w],queue,cpus)) for w in range(4)]
    for p in processes:p.start()
    runtime['worker_pids']=[p.pid for p in processes];run.write(ROOT/'runtime_provenance.json',runtime);print('CAVITY_WORKERS',json.dumps(runtime),flush=True)
    records=[];done=set();tic=time.monotonic();best=run.REFERENCE
    while len(done)<4:
        try:item=queue.get(timeout=20)
        except queue_mod.Empty:
            dead=[p for p in processes if p.exitcode not in (None,0)]
            if dead:raise RuntimeError('worker exited without completion: '+str([(p.pid,p.exitcode) for p in dead]))
            continue
        if item['type']=='done':done.add(item['worker']);continue
        row=item['record'];records.append(row);R=row.get('final_radius',1)
        if R<run.THRESHOLD and run.REFERENCE-R>run.MATERIAL:
            q=ROOT/'frozen_proposals'/f'grid{row["grid_id"]:04d}_R{R:.17f}.npy';np.save(q,np.load(ROOT/row['final_coordinate_file']));run.write(q.with_suffix('.json'),dict(row,frozen_file_sha256=run.sha(q),frozen_coordinate_sha256=run.csha(np.load(q))));print('FROZEN_PROPOSAL',row['grid_id'],R,str(q),run.sha(q),flush=True)
        if R<best:
            best=R;np.save(ROOT/'best.npy',np.load(ROOT/row['final_coordinate_file']));run.write(ROOT/'best.json',row);print('CAVITY_BEST',row['grid_id'],best,flush=True)
        run.write(ROOT/'checkpoint.json',dict(completed=len(records),requested=len(executed),grid_budget=64,grid_planned=len(grid),skipped=len(grid)-len(executed),workers_completed=sorted(done),best_radius=best,records=records,seconds=time.monotonic()-tic,scope='Finite paired-transport numerical search; no exact covering or global-optimality claim.'))
        if len(records)%8==0:print('CAVITY_PROGRESS',len(records),len(executed),'best',best,'seconds',time.monotonic()-tic,flush=True)
    for p in processes:p.join()
    run.write(ROOT/'completion.json',dict(completed=len(records),requested=len(executed),best_radius=best,errors=sum('error' in r for r in records),workers_exitcode=[p.exitcode for p in processes],seconds=time.monotonic()-tic))
    print('CAVITY_DONE',len(records),best,'seconds',time.monotonic()-tic,flush=True)

if __name__=="__main__":main()
