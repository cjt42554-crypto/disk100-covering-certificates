# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Finite structured remove-two/refill-two search from immutable de82 input."""
import engine as run
import os,sys,json,itertools,time,traceback,multiprocessing as mp,argparse,queue as queue_mod
from pathlib import Path
import numpy as np
from scipy.optimize import nnls
ROOT=run.ROOT; BASE=ROOT.parent

def pidinfo():
    return dict(pid=os.getpid(),namespace_pid_status=[s for s in Path('/proc/self/status').read_text().splitlines() if s.startswith('NSpid:')],affinity=sorted(os.sched_getaffinity(0)))

def mapped_inventory(C,labels,lam_by_key=None):
    g,p,J,keys,sites,top=run.inventory(C);out=[]
    for i,(key,st) in enumerate(zip(keys,sites)):
        mk=[key[0]]+[int(labels[s]) for s in key[1:4 if key[0]=='t' else 3 if key[0]=='b' else 2]]
        if key[0]=='b':mk.append(key[-1])
        # Boundary signs remain canonical because retained labels are increasing.
        out.append(dict(hole_id=i,key=key,original_labels_key=mk,sites=st,original_sites=[int(labels[s]) for s in st],point=p[i].tolist(),radius=float(np.sqrt(max(g[i],0))),squared_radius=float(g[i]),kind='boundary' if key[0] in ('b','f') else 'internal',input_stress=float((lam_by_key or {}).get(tuple(mk),0.))))
    return out

def prior_start_hashes():
    """Inspect local prior search traces only; no config/archive or other projects."""
    hashes={};arrays=0;errors=[]
    excluded={'.git','.aws','.agents','.codex','github_archive_staging','paired_transport_146','__pycache__'}
    for folder,dirs,files in os.walk(BASE):
        dirs[:]=[d for d in dirs if d not in excluded]
        for name in files:
            p=Path(folder)/name
            # Save named starts and predictor/event seeds too, including legacy arrays.
            if p.suffix not in ('.npy','.npz'):continue
            if not any(s in (name+' '+str(p.parent.relative_to(BASE))).lower() for s in ('start','seed','predictor','paired_event_search')):continue
            try:
                if p.suffix=='.npy':values=[('',np.load(p,allow_pickle=False))]
                else:
                    with np.load(p,allow_pickle=False) as z:values=[(k,z[k]) for k in z.files if k.lower() in ('c','centers','start','seed','points','x')]
                for key,a in values:
                    if a.shape==(201,):a=a[:-1].reshape(100,2)
                    elif a.shape==(200,):a=a.reshape(100,2)
                    if a.shape==(100,2):
                        hashes.setdefault(run.csha(a),[]).append(str(p.relative_to(BASE))+(':'+key if key else ''));arrays+=1
            except Exception as e:errors.append(dict(path=str(p.relative_to(BASE)),error=str(e)))
    return hashes,arrays,errors

def prepare():
    tic=time.monotonic();run.write(ROOT/'status.json',dict(stage='fresh100 deletion and pair screening',budget=128,**pidinfo()))
    C=np.load(ROOT/'inputs/baseline.npy');R,detail=run.true_cover(C);g,p,J,keys,sites,top=run.inventory(C);ac=np.flatnonzero(g.max()-g<1e-9);dg=run.diagnose(C)
    assert run.sha(ROOT/'inputs/baseline.npy')=='de82a12038c6f9a705b7a985c7e25c53c8e94f5c506946aba9a6c6f03a02852d'
    assert dg['binding_count']==100 and len(ac)==146 and len(top[1])==32
    M=np.r_[J[ac].T,np.ones((1,len(ac)))];rhs=np.r_[np.zeros(C.size),1.];lam,res=nnls(M,rhs,maxiter=10000);lm={tuple(keys[a]):float(lam[j]) for j,a in enumerate(ac)}
    baseline_holes=mapped_inventory(C,np.arange(100),lm);donors=[]
    for i in range(100):
        D=np.delete(C,i,axis=0);rr,dd=run.true_cover(D);incident=[j for j,a in enumerate(ac) if i in sites[a]]
        np.save(ROOT/'inputs'/f'delete_{i:03d}.npy',D)
        donors.append(dict(index=i,deletion_radius=rr,deletion_debt=rr-R,deletion_score=dd,incident_stress=float(lam[incident].sum()),incident_active_count=len(incident),incident_active_keys=[keys[ac[j]] for j in incident],boundary_incident_count=sum(keys[ac[j]][0]=='b' for j in incident),nonbinding=False))
    debt_rank={d['index']:r for r,d in enumerate(sorted(donors,key=lambda d:d['deletion_debt']))};stress_rank={d['index']:r for r,d in enumerate(sorted(donors,key=lambda d:d['incident_stress']))}
    pool=sorted({d['index'] for d in sorted(donors,key=lambda d:d['deletion_debt'])[:18]}|{d['index'] for d in sorted(donors,key=lambda d:d['incident_stress'])[:18]})
    adjacency={tuple(sorted(e)) for t in top[0] for e in itertools.combinations(map(int,t),2)}
    provisional=[]
    for a,b in itertools.combinations(pool,2):
        shared=[j for j,k in enumerate(ac) if a in sites[k] and b in sites[k]];dist=float(np.linalg.norm(C[a]-C[b]));adj=(a,b) in adjacency
        family='adjacent_shared_star' if adj and shared else 'separated_disjoint_stars' if not adj and not shared and dist>.32 else None
        if family is None:continue
        inc=[j for j,k in enumerate(ac) if a in sites[k] or b in sites[k]]
        provisional.append(dict(pair=[a,b],family=family,separation=dist,shared_active_count=len(shared),shared_active_keys=[keys[ac[j]] for j in shared],shared_stress=float(lam[shared].sum()),union_incident_stress=float(lam[inc].sum()),preselection_rank=debt_rank[a]+debt_rank[b]+stress_rank[a]+stress_rank[b],single_debt_sum=donors[a]['deletion_debt']+donors[b]['deletion_debt']))
    screening=[]
    for family in ('adjacent_shared_star','separated_disjoint_stars'):
        opts=sorted((x for x in provisional if x['family']==family),key=lambda x:(x['preselection_rank'],x['union_incident_stress']))[:32]
        for row in opts:
            a,b=row['pair'];labels=np.array([i for i in range(100) if i not in (a,b)]);D=C[labels];rr,dd=run.true_cover(D);holes=mapped_inventory(D,labels,lm);cavity=[h for h in holes if h['squared_radius']>R*R+1e-6]
            row.update(pair_deletion_radius=rr,pair_deletion_debt=rr-R,pair_deletion_score=dd,cavity_hole_count=len(cavity),cavity_internal_count=sum(h['kind']=='internal' for h in cavity),cavity_boundary_count=sum(h['kind']=='boundary' for h in cavity),maximum_cavity_point=dd['maximum_point'])
            k=len(screening);np.save(ROOT/'screening'/f'pair{k:03d}_delete.npy',D);row['screening_index']=k;row['deletion_file']=f'screening/pair{k:03d}_delete.npy';row['deletion_file_sha256']=run.sha(ROOT/row['deletion_file']);run.write(ROOT/'screening'/f'pair{k:03d}_inventory.json',dict(row,retained_labels=labels,actual_holes=holes,cavity_holes=cavity));screening.append(row)
        run.write(ROOT/'screening_checkpoint.json',dict(screened=len(screening),pool=pool,seconds=time.monotonic()-tic,**pidinfo()))
        print('PAIR_SCREENED',family,len(opts),time.monotonic()-tic,flush=True)
    # Four of each topology family: true pair debt first, with a weak-stress and
    # mixed-rank representative; repeated donor participation is capped at two.
    selected=[];usedpairs=set();counts={}
    for family in ('adjacent_shared_star','separated_disjoint_stars'):
        opts=[r for r in screening if r['family']==family]
        ranks={field:{r['screening_index']:j for j,r in enumerate(sorted(opts,key=lambda z:z[field]))} for field in ('pair_deletion_debt','union_incident_stress')}
        for r in opts:r['pair_mixed_rank']=ranks['pair_deletion_debt'][r['screening_index']]+ranks['union_incident_stress'][r['screening_index']]
        for selection,field in [('lowest_actual_pair_debt','pair_deletion_debt'),('second_low_actual_pair_debt','pair_deletion_debt'),('weak_incident_stress','union_incident_stress'),('mixed_actual_debt_stress','pair_mixed_rank')]:
            candidates=sorted(opts,key=lambda z:(z[field],z['pair_deletion_debt']))
            pick=next((r for r in candidates if tuple(r['pair']) not in usedpairs and all(counts.get(i,0)<2 for i in r['pair'])),None)
            if pick is None:pick=next(r for r in candidates if tuple(r['pair']) not in usedpairs)
            selected.append(dict(pick,selection=selection,selection_scope='Finite structured pool; this selection is not a global theorem.'));usedpairs.add(tuple(pick['pair']))
            for i in pick['pair']:counts[i]=counts.get(i,0)+1
    oldhashes,oldarrays,olderrors=prior_start_hashes();run.write(ROOT/'prior_start_hashes.json',oldhashes)
    grid=[];targetsets=[];seen=set();target_rejects=[]
    for pair_id,row in enumerate(selected):
        a,b=row['pair'];labels=np.array([i for i in range(100) if i not in (a,b)]);D=C[labels];holes=mapped_inventory(D,labels,lm);cavity=sorted((h for h in holes if h['squared_radius']>R*R+1e-6),key=lambda h:h['squared_radius'],reverse=True);boundary=[h for h in holes if h['kind']=='boundary'];internal=[h for h in cavity if h['kind']=='internal']
        candidates=[];uniques=[]
        def add(q1,q2,family,rationale):
            P=np.array([q1['point'],q2['point']]);reason=None
            if np.linalg.norm(P[0]-P[1])<.025:reason='two_target_centers_near_duplicate'
            if min(np.max(np.linalg.norm(P-C[[a,b]],axis=1)),np.max(np.linalg.norm(P[::-1]-C[[a,b]],axis=1)))<1e-8:reason='restores_original_two_center_set'
            if any(min(np.max(np.linalg.norm(P-Q,axis=1)),np.max(np.linalg.norm(P[::-1]-Q,axis=1)))<1e-8 for Q in uniques):reason='duplicate_unordered_target_geometry'
            if reason:target_rejects.append(dict(pair_id=pair_id,family=family,reason=reason,targets=[q1,q2]));return
            S=np.r_[D,P];rr,dd=run.true_cover(S);uniques.append(P);candidates.append(dict(targets=[q1,q2],family=family,rationale=rationale,preoptimization_whole_disk_radius=rr,preoptimization_whole_disk_score=dd))
        # A greedy98->99->100 inventory changes the second vacancy after the first fill.
        for q1 in cavity[:3]:
            D99=np.r_[D,np.array(q1['point'])[None,:]];h99=mapped_inventory(D99,np.r_[labels,100],lm);q2=max(h99,key=lambda h:h['squared_radius']);q2=dict(q2,inventory_stage='actual99_after_first_fill',first_point=q1['point']);add(q1,q2,'dynamic_cavity_refill','First actual98 vacancy critical; then actual99 whole-disk maximum after that insertion.')
        # Locally distinct maxima associated with both deleted donor positions.
        local=[]
        for i in (a,b):local.append(sorted(cavity,key=lambda h:np.linalg.norm(np.array(h['point'])-C[i]))[:3])
        for ia,ib in ((0,0),(0,1),(1,0),(1,1),(2,0),(0,2)):
            if len(local[0])>ia and len(local[1])>ib:add(local[0][ia],local[1][ib],'two_local_cavity_critical','Actual98 criticals near the two deleted centers, with alternative local split geometry.')
        # Move pair capacity between a deletion cavity and a real boundary gap.
        q=cavity[0];eligible=[h for h in boundary if np.linalg.norm(np.array(h['point'])-q['point'])>.12];eligible.sort(key=lambda h:(h['input_stress'],h['squared_radius']))
        for name,k in [('weak',0),('middle',len(eligible)//2),('strong',len(eligible)-1)]:
            if eligible:add(q,eligible[k],'cavity_plus_'+name+'_stress_boundary','Largest actual98 cavity plus an actual98 boundary Voronoi maximum selected by fresh input NNLS stress; zero means a newly exposed98 gap.')
        # Two actual boundary gaps with spread; one near the deleted pair.
        boundary_stress=sorted(boundary,key=lambda h:h['input_stress'],reverse=True)
        for q1 in boundary_stress[:3]:
            spread=[h for h in boundary if np.linalg.norm(np.array(h['point'])-q1['point'])>.35]
            if spread:add(q1,max(spread,key=lambda h:(h['input_stress'],h['squared_radius'])),'two_boundary_stress_gaps','Two actual98 boundary maxima with substantial spatial separation and fresh input stress priority.')
        if boundary:
            q1=min(boundary,key=lambda h:min(np.linalg.norm(np.array(h['point'])-C[a]),np.linalg.norm(np.array(h['point'])-C[b])));spread=[h for h in boundary if np.linalg.norm(np.array(h['point'])-q1['point'])>.35]
            if spread:add(q1,max(spread,key=lambda h:h['squared_radius']),'near_pair_plus_boundary_gap','Actual boundary gap near removed pair plus far actual98 boundary gap of largest debt.')
        for q1 in internal[:2]:
            farther=[h for h in internal if np.linalg.norm(np.array(h['point'])-q1['point'])>.06]
            if farther:add(q1,max(farther,key=lambda h:h['squared_radius']),'two_internal_cavity_critical','Two distinct internal critical holes from actual98 cluster vacancy; retains an internal alternative.')
        # Family coverage first, then fill remaining slots by actual start debt.
        chosen=[];used=set()
        buckets=[('dynamic_cavity_refill',1),('two_local_cavity_critical',1),('cavity_plus_strong_stress_boundary',1),('cavity_plus_middle_stress_boundary',1),('cavity_plus_weak_stress_boundary',1),('two_boundary_stress_gaps',1),('near_pair_plus_boundary_gap',1),('two_internal_cavity_critical',1)]
        for fam,num in buckets:
            opts=sorted([(i,z) for i,z in enumerate(candidates) if z['family']==fam],key=lambda z:z[1]['preoptimization_whole_disk_radius'])
            for i,z in opts[:num]:chosen.append(z);used.add(i)
        for i,z in sorted(enumerate(candidates),key=lambda z:z[1]['preoptimization_whole_disk_radius']):
            if len(chosen)>=8:break
            if i not in used:chosen.append(z);used.add(i)
        chosen=chosen[:8];run.write(ROOT/'screening'/f'selected_pair{pair_id:02d}_target_inventory.json',dict(pair=row,all_actual98_holes=holes,all_proposed_targetsets=candidates,chosen_targetsets=chosen))
        for target_id,z in enumerate(chosen):
            targetsets.append(dict(pair_id=pair_id,target_id=target_id,**z))
            for assignment in (0,1):
                index=pair_id*16+target_id*2+assignment;worker=pair_id%4;folder=ROOT/'workers'/f'worker{worker:02d}'
                for d in ('starts','finals','stages','proposals'):(folder/d).mkdir(parents=True,exist_ok=True)
                S=C.copy();targets=z['targets'] if assignment==0 else z['targets'][::-1];S[[a,b]]=[t['point'] for t in targets];cs=run.csha(S);reason='exact_prior_trace_duplicate' if cs in oldhashes else 'exact_new_grid_duplicate' if cs in seen else None
                gr=dict(grid_id=index,pair_id=pair_id,pair=[a,b],pair_family=row['family'],pair_selection=row['selection'],pair_deletion_debt=row['pair_deletion_debt'],pair_union_incident_stress=row['union_incident_stress'],target_id=target_id,assignment=assignment,targets=targets,target_family=z['family'],rationale=z['rationale'],preoptimization_whole_disk_radius=z['preoptimization_whole_disk_radius'],worker=worker,baseline_file_sha256=run.sha(ROOT/'inputs/baseline.npy'),baseline_coordinate_sha256=run.csha(C),start_coordinate_sha256=cs,skip=reason is not None,skip_reason=reason,offsets=0,all_center_variables=100,assignment_scope='Both assignments deliberately included: unordered geometric starts agree but label assignment and start SHA differ.')
                if reason is None:
                    assert np.array_equal(S[labels],C[labels]);np.save(folder/'starts'/f'{index:04d}.npy',S);gr['start_file_sha256']=run.sha(folder/'starts'/f'{index:04d}.npy');seen.add(cs)
                grid.append(gr)
        print('PAIR_TARGETS',pair_id,row['pair'],row['family'],len(chosen),'debt',row['pair_deletion_debt'],flush=True)
    frozen={}
    for name in ('rattler_contact_transfer','best147_boundary_transfer','next146_boundary_transfer'):
        mf=json.loads((BASE/name/'MANIFEST_SHA256.json').read_text());frozen[name]=dict(count=len(mf),manifest_sha256=run.sha(BASE/name/'MANIFEST_SHA256.json'),unchanged=all(run.sha(BASE/name/p)==s for p,s in mf.items()))
    checks=dict(stage='prepared',baseline_radius=R,baseline_whole_disk_score=detail,baseline_diagnostic=dg,input_file_sha256=run.sha(ROOT/'inputs/baseline.npy'),statistical_scope='All fresh NNLS/stress and deletion statistics are de82 baseline only.',fresh_NNLS_residual_inf=float(np.max(np.abs(M@lam-rhs))),lambda_min=float(lam.min()),lambda_max=float(lam.max()),single_deletion_debt_range=[min(d['deletion_debt'] for d in donors),max(d['deletion_debt'] for d in donors)],donor_pool=pool,pair_screened=len(screening),selected_pairs=[r['pair'] for r in selected],prior_start_arrays=oldarrays,prior_unique_hashes=len(oldhashes),prior_hash_read_errors=olderrors,grid_budget=128,grid_planned=len(grid),grid_executable=sum(not r['skip'] for r in grid),grid_skipped=sum(r['skip'] for r in grid),target_rejection_count=len(target_rejects),original_frozen_trees=frozen,seconds=time.monotonic()-tic,**pidinfo(),numerical_only=True)
    run.write(ROOT/'input_checks.json',checks);run.write(ROOT/'donor_inventory.json',donors);run.write(ROOT/'pair_screening.json',screening);run.write(ROOT/'selected_pairs.json',selected);run.write(ROOT/'baseline_inventory.json',baseline_holes);run.write(ROOT/'target_sets.json',targetsets);run.write(ROOT/'target_rejections.json',target_rejects);run.write(ROOT/'grid.json',grid);run.write(ROOT/'status.json',checks);np.savez_compressed(ROOT/'inputs/baseline_geometry.npz',C=C,g=g,p=p,J=J,active=ac,lam=lam,T=top[0],E=top[1],sg=top[2],F=top[3])
    print('PREPARED',json.dumps(run.clean(checks)),flush=True);return grid

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
    runtime['worker_pids']=[p.pid for p in processes];run.write(ROOT/'runtime_provenance.json',runtime);print('PAIR_WORKERS',json.dumps(runtime),flush=True)
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
            best=R;np.save(ROOT/'best.npy',np.load(ROOT/row['final_coordinate_file']));run.write(ROOT/'best.json',row);print('PAIR_BEST',row['grid_id'],best,flush=True)
        run.write(ROOT/'checkpoint.json',dict(completed=len(records),requested=len(executed),grid_budget=128,grid_planned=len(grid),skipped=len(grid)-len(executed),workers_completed=sorted(done),best_radius=best,records=records,seconds=time.monotonic()-tic,scope='Finite paired-transport numerical search; no exact covering or global-optimality claim.'))
        if len(records)%8==0:print('PAIR_PROGRESS',len(records),len(executed),'best',best,'seconds',time.monotonic()-tic,flush=True)
    for p in processes:p.join()
    run.write(ROOT/'completion.json',dict(completed=len(records),requested=len(executed),best_radius=best,errors=sum('error' in r for r in records),workers_exitcode=[p.exitcode for p in processes],seconds=time.monotonic()-tic))
    print('PAIR_DONE',len(records),best,'seconds',time.monotonic()-tic,flush=True)

if __name__=='__main__':main()
