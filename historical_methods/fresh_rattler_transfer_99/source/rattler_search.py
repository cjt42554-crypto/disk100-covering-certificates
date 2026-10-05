# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""New99 skeletons: move their inactive-cover disk to actual maximum holes."""
import engine as run
from alignment import compare
import os,sys,json,itertools,time,traceback,multiprocessing as mp,argparse,queue as queue_mod,collections
from pathlib import Path
import numpy as np
from scipy.optimize import nnls,linear_sum_assignment
ROOT=run.ROOT;BASE=ROOT.parent

def pidinfo():return dict(pid=os.getpid(),namespace_pid_status=[s for s in Path('/proc/self/status').read_text().splitlines() if s.startswith('NSpid:')],affinity=sorted(os.sched_getaffinity(0)))
def unordered_hash(C):
    C=np.asarray(C,dtype='<f8');P=C[np.lexsort((C[:,1],C[:,0]))];return run.csha(P)

def prior_hashes():
    ordered={};unordered={};arrays=0;errors=[];excluded={'.git','.aws','.agents','.codex','github_archive_staging','fresh_rattler_transfer_99','__pycache__'}
    frozenmath={'rattler_contact_transfer','best147_boundary_transfer','next146_boundary_transfer','paired_transport_146','topology_patch_hybrid_146','event_geometry_146'}
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

def skeleton_compare(target,source,seeds=144):
    """Explicit99-to99 or99-to100 subset matching, rotation/reflection only."""
    best=None
    for reflection in (False,True):
        V=source.copy()
        if reflection:V[:,1]*=-1
        for initial in np.linspace(0,2*np.pi,seeds,endpoint=False):
            angle=float(initial)
            for it in range(10):
                Q=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]);W=V@Q;rows,cols=linear_sum_assignment(np.sum((target[:,None,:]-W[None,:,:])**2,axis=2));A=target[rows];Z=V[cols];angle=float(np.arctan2(np.sum(A[:,0]*Z[:,1]-A[:,1]*Z[:,0]),np.sum(A*Z)))
            Q=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]);d=np.linalg.norm(target[rows]-V[cols]@Q,axis=1);error=float(np.sum(d*d))
            if best is None or error<best['sum_squared_error']:best=dict(sum_squared_error=error,maximum_binding_center_distance=float(d.max()),rotation_angle=angle,reflection=reflection,target_selected_rows=rows,source_selected_columns=cols,target_unselected_indices=sorted(set(range(len(target)))-set(rows)),coordinate_match_below2e6=bool(d.max()<2e-6))
    best.update(rotation_seeds=seeds,scope='Successful99 coordinate alignment is a geometry witness; failed finite alignment is not a proof of global geometric inequivalence. Radius is not used for matching.');return best

def frozen_checks():
    out={}
    for name in ('rattler_contact_transfer','best147_boundary_transfer','next146_boundary_transfer','paired_transport_146','topology_patch_hybrid_146','event_geometry_146'):
        mf=json.loads((BASE/name/'MANIFEST_SHA256.json').read_text());out[name]=dict(count=len(mf),manifest_sha256=run.sha(BASE/name/'MANIFEST_SHA256.json'),unchanged=all(run.sha(BASE/name/p)==h for p,h in mf.items()))
    return out

def prepare():
    tic=time.monotonic();run.write(ROOT/'status.json',dict(stage='fresh source geometry/99 novelty/deletion and hole inventory',budget=96,**pidinfo()));source_rows=json.loads((ROOT/'inputs/source_inputs.json').read_text());refs=[('old65_known99','inputs/known99_reference.npy'),('old001','inputs/old001_100.npy'),('old003','inputs/old003_100.npy'),('current146','inputs/current146_100.npy')];refdata={};refdiags={}
    for label,path in refs:
        C=np.load(ROOT/path);refdata[label]=C;refdiags[label]=dict(file=path,file_sha256=run.sha(ROOT/path),diagnostic=run.diagnose(C),actual_whole_disk_radius=run.true_cover(C)[0])
    inventory=[];novelty=[];targets=[];grid=[];seen=set();seenu=set();prior,prioru,oldarrays,olderrors=prior_hashes();assert not olderrors;run.write(ROOT/'prior_ordered_hashes.json',prior);run.write(ROOT/'prior_unordered_hashes.json',prioru)
    for row in source_rows:
        i=row['source_index'];free=row['free_index'];C=np.load(ROOT/row['copied_file']);assert run.sha(ROOT/row['copied_file'])==row['input_file_sha256']==run.sha(row['input_original_file']);R,score=run.true_cover(C);g,p,J,keys,sites,top=run.inventory(C);ac=np.flatnonzero(g.max()-g<1e-9);dg=run.diagnose(C);used=dg['binding_indices'];assert len(used)==99 and dg['nonbinding_indices']==[free];M=np.r_[J[ac].T,np.ones((1,len(ac)))];rhs=np.r_[np.zeros(C.size),1.];lam,res=nnls(M,rhs,maxiter=20000);lm={int(a):float(lam[j]) for j,a in enumerate(ac)};D=np.delete(C,free,axis=0);debtR,debtScore=run.true_cover(D);np.save(ROOT/'inputs'/f'source{i:02d}_delete_free{free:03d}.npy',D)
        source=dict(row,radius=R,whole_disk_score=score,diagnostic=dg,parent_skeleton_coordinate_sha256=run.csha(C[used]),parent_skeleton_unordered_coordinate_sha256=unordered_hash(C[used]),parent_active_inventory_sha256=dg['active_signature'],parent_inactive_inventory_sha256=dg['inactive_signature'],without_nonbinding_radius=debtR,nonbinding_deletion_debt=debtR-R,nonbinding_deletion_score=debtScore,NNLS_residual_inf=float(np.max(np.abs(M@lam-rhs))),NNLS_lambda_min=float(lam.min()),NNLS_lambda_max=float(lam.max()),NNLS_lambda_sum=float(lam.sum()),active_keys=[keys[a] for a in ac],NNLS_lambda=lam,statistical_scope='Fresh statistics refer only to this input SHA; nonbinding is not legal deletion or stationarity certification.')
        comparisons=[]
        for label,path in refs:
            K=refdata[label];kd=refdiags[label]['diagnostic'];full=compare(K,C,144);target=K[kd['binding_indices']] if kd['binding_count']==99 else K;sk=skeleton_compare(target,C[used],144);comparisons.append(dict(reference=label,reference_file_sha256=refdiags[label]['file_sha256'],reference_binding_count=kd['binding_count'],full100_geometry_and_inventory_comparison=full,binding99_or99subset_comparison=sk,actual_active_inventory_counts=[dg['active_count'],kd['active_count']],actual_boundary_counts=[dg['boundary_count'],kd['boundary_count']]))
        assert not comparisons[0]['binding99_or99subset_comparison']['coordinate_match_below2e6'];assert not any(x['full100_geometry_and_inventory_comparison']['matches'] for x in comparisons)
        novelty.append(dict(source_index=i,input_file_sha256=row['input_file_sha256'],actual_binding_count=99,nonbinding_index=free,comparisons=comparisons,old65_skeleton_match_found=False,scope='All100 geometry plus actual inventory and separate99 binding comparison, not radius-only basin detection. No finite nonmatch is claimed as a global inequivalence theorem.'))
        holes=[]
        for a in ac:
            assert free not in sites[a];S=C.copy();S[free]=p[a];sr,ss=run.true_cover(S);holes.append(dict(hole_id=int(a),key=keys[a],sites=sites[a],point=p[a],kind='boundary' if keys[a][0] in ('b','f') else 'internal',support_radius=float(np.sqrt(g[a])),stress=lm[int(a)],distance_from_free=float(np.linalg.norm(p[a]-C[free])),start_actual_whole_disk_radius=sr,start_whole_disk_score=ss,actual_cavity_prior_gain=debtR-sr,angle=float(np.arctan2(p[a,1],p[a,0]))))
        boundary=[h for h in holes if h['kind']=='boundary'];internal=[h for h in holes if h['kind']=='internal'];assert len(boundary)==32
        chosen=[];seenholes=set()
        def select(options,count,family):
            added=0
            for h in options:
                if h['hole_id'] in seenholes:continue
                chosen.append(dict(h,selection_family=family));seenholes.add(h['hole_id']);added+=1
                if added==count:break
            assert added==count
        if i==0:
            select(sorted(boundary,key=lambda h:h['actual_cavity_prior_gain'],reverse=True),3,'control_boundary_actual_cavity_gain')
            select(sorted(boundary,key=lambda h:h['stress'],reverse=True),3,'control_boundary_high_stress')
            order=sorted(boundary,key=lambda h:h['stress']);mid=len(order)//2;select(order[mid:]+order[:mid][::-1],3,'control_boundary_middle_stress');select(order,3,'control_boundary_weak_stress')
            select(sorted(internal,key=lambda h:h['actual_cavity_prior_gain'],reverse=True),2,'internal_actual_cavity_gain');select(sorted(internal,key=lambda h:h['stress'],reverse=True),1,'internal_high_stress');select(sorted(internal,key=lambda h:h['stress']),1,'internal_weak_stress')
        else:
            select(boundary,32,'all32_actual_boundary_maxholes')
            select(sorted(internal,key=lambda h:(h['actual_cavity_prior_gain'],-h['distance_from_free']),reverse=True),2,'internal_actual_cavity_gain')
            select(sorted(internal,key=lambda h:h['stress'],reverse=True),2,'internal_high_stress');order=sorted(internal,key=lambda h:h['stress']);mid=len(order)//2;select(order[mid:]+order[:mid][::-1],2,'internal_middle_stress');select(order,2,'internal_weak_stress')
        expected=16 if i==0 else 40;assert len(chosen)==expected
        # Debt repair prioritizes actual attempts; provenance records every chosen hole.
        chosen.sort(key=lambda h:(-h['actual_cavity_prior_gain'],-h['stress'],h['hole_id']))
        for local,h in enumerate(chosen):
            index=(0 if i==0 else 16 if i==1 else 56)+local;worker=index%4;S=C.copy();S[free]=h['point'];cs=run.csha(S);us=unordered_hash(S);reason='exact_prior_initial_sha_duplicate' if cs in prior else 'exact_prior_unordered_geometry_duplicate' if us in prioru else 'exact_new_initial_sha_duplicate' if cs in seen else 'exact_new_unordered_geometry_duplicate' if us in seenu else None
            gr=dict(grid_id=index,source_index=i,free_center_index=free,parent_input_file=row['copied_file'],parent_input_file_sha256=row['input_file_sha256'],parent_coordinate_sha256=run.csha(C),parent_skeleton_coordinate_sha256=source['parent_skeleton_coordinate_sha256'],parent_skeleton_unordered_coordinate_sha256=source['parent_skeleton_unordered_coordinate_sha256'],parent_active_inventory_sha256=dg['active_signature'],parent_inactive_inventory_sha256=dg['inactive_signature'],parent_radius=R,parent_nonbinding_deletion_debt=debtR-R,target=h,target_prior_gain=h['actual_cavity_prior_gain'],target_stress=h['stress'],target_selection_rationale='All actual boundary maxima for priority sources1/2; finite control subset for0. Internal holes cover debt gain and high/middle/weak fresh NNLS stress.',worker=worker,start_coordinate_sha256=cs,start_unordered_geometry_sha256=us,skip=reason is not None,skip_reason=reason,offsets=0,all_center_variables=100)
            if reason is None:
                folder=ROOT/'workers'/f'worker{worker:02d}'
                for d in ('starts','finals','stages','proposals'):(folder/d).mkdir(parents=True,exist_ok=True)
                path=folder/'starts'/f'{index:04d}.npy';assert np.array_equal(S[used],C[used]);np.save(path,S);gr['start_file_sha256']=run.sha(path);seen.add(cs);seenu.add(us)
            grid.append(gr);targets.append(dict(source_index=i,grid_id=index,**h))
        inventory.append(source);run.write(ROOT/'screening'/f'source{i:02d}_all_actual_holes.json',holes);run.write(ROOT/'screening'/f'source{i:02d}_selected_targets.json',chosen);run.write(ROOT/'novelty_checkpoint.json',dict(completed_sources=len(inventory),total_sources=3,novelty=novelty,seconds=time.monotonic()-tic,**pidinfo()));print('NEW99_SOURCE',i,'free',free,'active',len(ac),'R',R,'delete_debt',debtR-R,'old65_binding_match_maxdistance',comparisons[0]['binding99_or99subset_comparison']['maximum_binding_center_distance'],'targets',len(chosen),flush=True)
    frozen=frozen_checks();assert all(v['unchanged'] for v in frozen.values());assert len(grid)==96;checks=dict(stage='prepared',grid_budget=96,planned=len(grid),executable=sum(not r['skip'] for r in grid),skipped=sum(r['skip'] for r in grid),skip_counts=dict(collections.Counter(r['skip_reason'] for r in grid if r['skip'])),source_budgets={0:16,1:40,2:40},source_boundary_targets={0:12,1:32,2:32},source_internal_targets={0:4,1:8,2:8},prior_initial_arrays=oldarrays,prior_initial_unique_ordered=len(prior),prior_initial_unique_unordered=len(prioru),prior_read_errors=olderrors,frozen_trees=frozen,seconds=time.monotonic()-tic,**pidinfo(),scope='New physically changed99 skeletons retain their100th inactive-cover disk. Finite transport inventory; no free deletion, exact covering or global-optimality claim.');run.write(ROOT/'input_checks.json',checks);run.write(ROOT/'source_inventory.json',inventory);run.write(ROOT/'source_novelty.json',novelty);run.write(ROOT/'reference_inventory.json',refdiags);run.write(ROOT/'target_ledger.json',targets);run.write(ROOT/'grid.json',grid);run.write(ROOT/'status.json',checks);print('PREPARED',json.dumps(run.clean(checks)),flush=True);return grid


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
    runtime['worker_pids']=[p.pid for p in processes];run.write(ROOT/'runtime_provenance.json',runtime);print('RATTLER_WORKERS',json.dumps(runtime),flush=True)
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
            best=R;np.save(ROOT/'best.npy',np.load(ROOT/row['final_coordinate_file']));run.write(ROOT/'best.json',row);print('RATTLER_BEST',row['grid_id'],best,flush=True)
        run.write(ROOT/'checkpoint.json',dict(completed=len(records),requested=len(executed),grid_budget=96,grid_planned=len(grid),skipped=len(grid)-len(executed),workers_completed=sorted(done),best_radius=best,records=records,seconds=time.monotonic()-tic,scope='Finite paired-transport numerical search; no exact covering or global-optimality claim.'))
        if len(records)%8==0:print('RATTLER_PROGRESS',len(records),len(executed),'best',best,'seconds',time.monotonic()-tic,flush=True)
    for p in processes:p.join()
    run.write(ROOT/'completion.json',dict(completed=len(records),requested=len(executed),best_radius=best,errors=sum('error' in r for r in records),workers_exitcode=[p.exitcode for p in processes],seconds=time.monotonic()-tic))
    print('RATTLER_DONE',len(records),best,'seconds',time.monotonic()-tic,flush=True)

if __name__=="__main__":main()
