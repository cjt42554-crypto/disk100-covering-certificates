# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""New146 input:12x32 boundary transfers plus15 local internal probes."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import numpy as np,json,hashlib,time,traceback,argparse,multiprocessing as mp,itertools
from scipy.linalg import null_space
from scipy.optimize import nnls,minimize
import geometry as geom
from whole_disk_oracle import true_cover
ROOT=Path(__file__).resolve().parent.parent;BASE=ROOT.parent
REFERENCE=.11544695488558258
THRESHOLD=.1154469548855
MATERIAL=1e-10
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def csha(C):return hashlib.sha256(np.asarray(C,dtype='<f8').tobytes()).hexdigest()
def clean(v):
    if isinstance(v,np.ndarray):return clean(v.tolist())
    if isinstance(v,(np.integer,np.bool_)):return v.item()
    if isinstance(v,(np.floating,float)):return float(v) if np.isfinite(v) else None
    if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    return v
def write(p,value):
    p=Path(p);q=Path(str(p)+'.tmp');q.write_text(json.dumps(clean(value),ensure_ascii=False,indent=2,allow_nan=False)+'\n');q.replace(p)
def inventory(C,top=None):
    if top is None:top=geom.topology(C)
    T,E,sg,F=top;tg,P,tj=geom.triangles(C,T);bg,Q,bj=geom.boundary(C,E,sg);fg,fj=geom.single_boundary(C,F) if len(F) else (np.empty(0),np.empty((0,2)))
    J=np.zeros((len(T)+len(E)+len(F),C.size),dtype=np.result_type(C));tc=(2*T[:,:,None]+np.array([0,1])[None,None,:]).reshape(len(T),6);ec=(2*E[:,:,None]+np.array([0,1])[None,None,:]).reshape(len(E),4)
    J[np.arange(len(T))[:,None],tc]=tj.reshape(len(T),6);J[len(T)+np.arange(len(E))[:,None],ec]=bj.reshape(len(E),4)
    if len(F):J[len(T)+len(E)+np.arange(len(F))[:,None],2*F[:,None]+np.array([0,1])[None,:]]=fj
    keys=[['t',*sorted(map(int,t))] for t in T]+[['b',int(e[0]),int(e[1]),int(s)] for e,s in zip(E,sg)]+[['f',int(i)] for i in F]
    sites=[list(map(int,t)) for t in T]+[list(map(int,e)) for e in E]+[[int(i)] for i in F]
    points=np.r_[P,Q,-C[F]/np.linalg.norm(C[F],axis=1)[:,None]]
    return np.r_[tg,bg,fg],points,J,keys,sites,top
def signature(keys):return hashlib.sha256(json.dumps(sorted(keys),separators=(',',':')).encode()).hexdigest()
def diagnose(C):
    g,p,J,kk,sites,top=inventory(C);active=np.flatnonzero(g.max()-g<1e-9);used=sorted({i for j in active for i in sites[j]});inactive=[i for i in range(len(g)) if i not in set(active)]
    return dict(active_count=len(active),binding_count=len(used),binding_indices=used,nonbinding_indices=sorted(set(range(len(C)))-set(used)),active_keys=[kk[i] for i in active],inactive_keys=[kk[i] for i in inactive],active_signature=signature([kk[i] for i in active]),inactive_signature=signature([kk[i] for i in inactive]),interior_count=len(top[0]),boundary_count=len(top[1]),singleton_count=len(top[3]),active_tolerance_squared=1e-9)
def prepare():
    tic=time.monotonic();C=np.load(ROOT/'inputs/baseline.npy');R,detail=true_cover(C);g,p,J,kk,sites,top=inventory(C);ac=np.flatnonzero(g.max()-g<1e-9);dg=diagnose(C);assert len(ac)==146 and dg['binding_count']==100
    M=np.r_[J[ac].T,np.ones((1,len(ac)))];rhs=np.r_[np.zeros(C.size),1.];lam,res=nnls(M,rhs,maxiter=10000)
    H=np.zeros((C.size,C.size));eps=1e-26
    for j in range(C.size):
        Z=C.astype(complex);Z.ravel()[j]+=1j*eps;AJ=inventory(Z,top)[2];H[:,j]=np.imag(AJ[ac].T@lam)/eps
    w=np.stack((-C[:,1],C[:,0]),axis=1).ravel();Q=null_space(w[None,:]);AQ=J[ac]@Q;ny=Q.shape[1];m=len(ac);K=np.zeros((ny+1+m,ny+1+m));K[:ny,:ny]=Q.T@((H+H.T)/2)@Q;K[:ny,ny+1:]=AQ.T;K[ny+1:,:ny]=AQ;K[ny,ny+1:]=-1;K[ny+1:,ny]=-1
    br=np.flatnonzero(np.array([kk[i][0]=='b' for i in ac]));rhs2=np.zeros((len(K),len(br)));rhs2[ny+1+br,np.arange(len(br))]=1.;sol=np.linalg.lstsq(K,rhs2,rcond=1e-11)[0];motion=(Q@sol[:ny]).reshape(100,2,len(br));response=np.sqrt(np.mean(np.sum(motion*motion,axis=1),axis=1))
    donors=[]
    for i in range(100):
        rr,dd=true_cover(np.delete(C,i,axis=0));incident=[j for j,a in enumerate(ac) if i in sites[a]];incinactive=[j for j in range(len(g)) if j not in set(ac) and i in sites[j]]
        donors.append(dict(index=i,center=C[i].tolist(),without_donor_radius=rr,deletion_debt=rr-R,deletion_whole_disk_score=dd,incident_active_hole_count=len(incident),incident_active_keys=[kk[ac[j]] for j in incident],incident_stress=float(lam[incident].sum()),incident_star_signature=signature([kk[ac[j]] for j in incident]),incident_inactive_hole_count=len(incinactive),incident_inactive_keys=[kk[j] for j in incinactive],maximum_incident_inactive_radius=float(np.sqrt(g[incinactive].max())) if incinactive else None,boundary_stationary_response_rms=float(response[i]),nonbinding=False,scope='Deletion tested by the whole-disk numerical oracle. No center is treated as removable at zero cost.'))
    ranks={field:{r['index']:j for j,r in enumerate(sorted(donors,key=lambda x:x[field]))} for field in ('deletion_debt','incident_stress','boundary_stationary_response_rms')}
    for d in donors:d['mixed_rank_score']=ranks['deletion_debt'][d['index']]+ranks['incident_stress'][d['index']]-.5*ranks['boundary_stationary_response_rms'][d['index']]
    selected=[];seen=set()
    for family,field,count in [('low_deletion_debt','deletion_debt',5),('low_incident_stress','incident_stress',4),('mixed_debt_stress_response','mixed_rank_score',3)]:
        added=0
        for d in sorted(donors,key=lambda x:x[field]):
            if d['index'] in seen:continue
            selected.append(dict(d,selection_family=family,selection_scope='Five low-debt, four low-incident-stress, three mixed-rank distinct indices. This finite heuristic is not a global donor selection theorem.'));seen.add(d['index']);added+=1
            if added==count:break
    assert len(selected)==12
    boundary=[dict(hole_id=j,key=kk[j],point=p[j].tolist(),sites=sites[j],kind='boundary_grid',active=bool(j in set(ac)),support_radius=float(np.sqrt(g[j]))) for j in range(len(top[0]),len(top[0])+len(top[1]))];assert len(boundary)==32 and all(h['active'] for h in boundary)
    prior=BASE/'rattler_contact_transfer';oldhashes=set()
    for sub in ('','max_stress_supplement','boundary_followup'):
        folder=prior/sub
        for path in (folder/'starts').glob('*.npy'):oldhashes.add(csha(np.load(path)))
    prev=BASE/'best147_boundary_transfer'
    for path in (prev/'workers').glob('worker*/starts/*.npy'):oldhashes.add(csha(np.load(path)))
    grid=[];seenstart=set();started=0;internal=[]
    for k,d in enumerate(selected):
        worker=k%4;folder=ROOT/'workers'/f'worker{worker:02d}'
        for name in ('starts','finals','stages','proposals'):(folder/name).mkdir(parents=True,exist_ok=True)
        local=[]
        if k<5:
            options=[j for j in ac if kk[j][0]=='t' and d['index'] not in sites[j]]
            options.sort(key=lambda j:np.sum((p[j]-C[d['index']])**2))
            for j in options[:3]:
                h=dict(hole_id=int(j),key=kk[j],point=p[j].tolist(),sites=sites[j],kind='internal_local_probe',active=True,support_radius=float(np.sqrt(g[j])),distance_from_donor=float(np.linalg.norm(p[j]-C[d['index']])),rationale='One of three nearest actual active internal maxholes outside this low-debt donor star; probes the local cavity allocation without restricting targets to the boundary.')
                local.append(h);internal.append(dict(donor=d['index'],target=h))
        targets=[(k*32+j,h) for j,h in enumerate(boundary)]+[(384+k*3+j,h) for j,h in enumerate(local)]
        for index,h in targets:
            S=C.copy();S[d['index']]=h['point'];cs=csha(S);reason=None
            if d['index'] in h['sites']:reason='target_hole_already_incident_to_donor'
            elif np.linalg.norm(S[d['index']]-C[d['index']])<1e-9:reason='zero_effect_target_near_original_center'
            elif cs in oldhashes:reason='exact_duplicate_of_prior83_or382'
            elif cs in seenstart:reason='exact_duplicate_in_new_grid'
            row=dict(grid_id=index,donor=d['index'],donor_selection=k,donor_family=d['selection_family'],donor_incident_stress=d['incident_stress'],donor_deletion_debt=d['deletion_debt'],donor_response=d['boundary_stationary_response_rms'],target=h,worker=worker,baseline_file_sha256=sha(ROOT/'inputs/baseline.npy'),baseline_coordinate_sha256=csha(C),start_coordinate_sha256=cs,skip=reason is not None,skip_reason=reason,offsets=0,all_center_variables=100)
            if reason is None:
                assert np.array_equal(S[np.arange(100)!=d['index']],C[np.arange(100)!=d['index']]);np.save(folder/'starts'/f'{index:04d}.npy',S);seenstart.add(cs);started+=1
            grid.append(row)
    oldselected=json.loads((prev/'selected_donors.json').read_text());statscope=dict(branch='grid25/de82 new146 numerical baseline only',input_file_sha256=sha(ROOT/'inputs/baseline.npy'),input_coordinate_sha256=csha(C))
    checks=dict(baseline_radius=R,baseline_whole_disk_score=detail,baseline_diagnostic=dg,baseline_file_sha256=sha(ROOT/'inputs/baseline.npy'),statistical_input_scope=statscope,old_tree_manifest_sha256=sha(prior/'MANIFEST_SHA256.json'),old_tree_file_count=429,old_tree_checks_pass=all(sha(prior/p)==s for p,s in json.loads((prior/'MANIFEST_SHA256.json').read_text()).items()),old147_tree_manifest_sha256=sha(prev/'MANIFEST_SHA256.json'),old147_tree_checks_pass=all(sha(prev/p)==s for p,s in json.loads((prev/'MANIFEST_SHA256.json').read_text()).items()),old_start_hash_count=len(oldhashes),old003_selected_donors=[d['index'] for d in oldselected],new_selected_donors=[d['index'] for d in selected],new_only_donor_indices=sorted(set(d['index'] for d in selected)-set(d['index'] for d in oldselected)),independent_NNLS_residual_inf=float(np.max(np.abs(M@lam-rhs))),lambda_sum=float(lam.sum()),lambda_min=float(lam.min()),lambda_max=float(lam.max()),stationary_response_rank=int(np.linalg.matrix_rank(K,tol=1e-10)),stationary_response_dimension=len(K),stationary_response_residual_inf=float(np.max(np.abs(K@sol-rhs2))),H_asymmetry=float(np.max(np.abs(H-H.T))),grid_budget=399,grid_boundary_budget=384,grid_internal_budget=15,grid_executable=started,grid_skipped=len(grid)-started,skip_counts={reason:sum(r['skip_reason']==reason for r in grid) for reason in set(r['skip_reason'] for r in grid if r['skip'])},workers=4,BLAS_threads_per_worker=1,source_hashes={str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'source').glob('*.py'))},seconds=time.monotonic()-tic,scope='Numerical whole-disk donor deletion and finite transfer inventory; no exact certificate or global-optimality claim.')
    write(ROOT/'input_checks.json',checks);write(ROOT/'donor_inventory.json',donors);write(ROOT/'selected_donors.json',selected);write(ROOT/'boundary_inventory.json',boundary);write(ROOT/'internal_inventory.json',internal);write(ROOT/'grid.json',grid);np.savez_compressed(ROOT/'inputs/baseline_geometry.npz',C=C,T=top[0],E=top[1],sg=top[2],F=top[3],g=g,p=p,J=J,active=ac,lam=lam,H=H,Q=Q,KKT=K,stationary_response=sol)
    print('PREPARED',json.dumps(dict(executable=started,skipped=len(grid)-started,donors=[r['index'] for r in selected],new_only_donors=checks['new_only_donor_indices'],debt_range=[min(d['deletion_debt'] for d in donors),max(d['deletion_debt'] for d in donors)],seconds=checks['seconds'])),flush=True)
    return grid
def known99(C):
    K=np.load(ROOT/'inputs/known99_reference.npy');ref=diagnose(K);dd=diagnose(C);used=ref['binding_indices'];ang=np.arctan2(np.sum(C[used,0]*K[used,1]-C[used,1]*K[used,0]),np.sum(C[used]*K[used]));rot=np.array([[np.cos(ang),-np.sin(ang)],[np.sin(ang),np.cos(ang)]]);dist=np.linalg.norm(C-K@rot,axis=1);R,detail=true_cover(C)
    match=dist[used].max()<2e-4 and dd['binding_indices']==used and dd['active_signature']==ref['active_signature'] and abs(R-.11551610349101131)<1e-8
    return dict(match=bool(match),binding_skeleton_max_distance=float(dist[used].max()),all_center_max_distance=float(dist.max()),active_keys_match=dd['active_signature']==ref['active_signature'],inactive_keys_match=dd['inactive_signature']==ref['inactive_signature'],whole_disk_radius=R,free_indices=dd['nonbinding_indices'],scope='Alignment of actual99 binding centers plus full active-key identity and whole-disk oracle; never radius-only detection.')
def optimize(C,folder,row):
    best=C.copy();bestR,initial_detail=true_cover(C);rounds=[];objgrad=np.r_[np.zeros(C.size),1.];index=row['grid_id']
    for outer in range(2):
        radius,top=geom.covering(C);before=geom.top_key(top);model=geom.Model(C,top);x0=np.r_[C.ravel(),radius+1e-10];trust=.18;lower=np.r_[np.maximum(C.ravel()-trust,-1),.09];upper=np.r_[np.minimum(C.ravel()+trust,1),max(.4,radius+.05)];tic=time.monotonic()
        res=minimize(lambda x:x[-1],x0,jac=lambda x:objgrad,method='SLSQP',constraints=[dict(type='ineq',fun=model.con,jac=model.jac)],bounds=list(zip(lower,upper)),options=dict(maxiter=450,ftol=1e-13,disp=False));nextC=res.x[:-1].reshape(-1,2);path=folder/'stages'/f'{index:04d}_{outer:02d}.npy';np.save(path,nextC)
        if not np.isfinite(nextC).all():raise ValueError('nonfinite_coordinates')
        R,detail=true_cover(nextC);rr,nexttop=geom.covering(nextC);changed=before!=geom.top_key(nexttop);mincon=float(model.con(res.x).min());feasible=abs(R-res.x[-1])<1e-10 and mincon>-2e-8 and np.linalg.norm(nextC,axis=1).max()<=1+1e-10
        oldkeys=inventory(C,top)[3];newg,newp,newJ,newkeys,newsites,newtop=inventory(nextC,nexttop);newactive=np.flatnonzero(newg.max()-newg<1e-9);missing=sorted(set(map(tuple,newkeys))-set(map(tuple,oldkeys)));removed=sorted(set(map(tuple,oldkeys))-set(map(tuple,newkeys)));newbinding=sorted(set(tuple(newkeys[j]) for j in newactive)-set(map(tuple,oldkeys)))
        stage=dict(round=outer,coordinate_file=str(path.relative_to(ROOT)),file_sha256=sha(path),coordinate_sha256=csha(nextC),radius=R,whole_disk_score=detail,model_radius=rr,objective=float(res.x[-1]),success=bool(res.success),status=int(res.status),message=str(res.message),iterations=int(res.nit),topology_changed=changed,new_actual_keys=missing,removed_actual_keys=removed,new_actual_binding_keys=newbinding,true_radius_minus_frozen_objective=float(R-res.x[-1]),minimum_frozen_constraint=mincon,feasible=bool(feasible),max_center_norm=float(np.linalg.norm(nextC,axis=1).max()),seconds=time.monotonic()-tic,offsets=0)
        if R<bestR and stage['max_center_norm']<=1+1e-10:bestR=R;best=nextC.copy()
        rounds.append(stage);C=nextC
        if R<THRESHOLD and REFERENCE-R>MATERIAL:
            q=folder/'proposals'/f'grid{index:04d}_round{outer:02d}_R{R:.17f}.npy';np.save(q,C);write(q.with_suffix('.json'),dict(row,stage=stage,proposal_file_sha256=sha(q),proposal_coordinate_sha256=csha(C),scope='Numerical new proposal, exact integer/surd certification required.'));print('POTENTIAL_IMPROVEMENT',row['worker'],index,R,str(q),flush=True)
        if feasible and R>.11551 and diagnose(C)['binding_count']==99:
            kb=known99(C);stage['known99']=kb
            if kb['match']:stage['second_round_skipped_known99_family']=True;break
        if feasible and res.success and not changed:break
    return bestR,best,initial_detail,rounds
def worker_main(worker,rows,queue,cpus):
    if hasattr(os,'sched_setaffinity'):os.sched_setaffinity(0,{cpus[worker%len(cpus)]})
    folder=ROOT/'workers'/f'worker{worker:02d}';records=[];tic=time.monotonic()
    if (folder/'checkpoint.json').exists():records=json.loads((folder/'checkpoint.json').read_text())['records']
    done={r['grid_id']:r for r in records}
    with open(folder/'records.jsonl','a',buffering=1) as fp:
        for row in rows:
            if row['grid_id'] in done:queue.put(dict(type='record',record=done[row['grid_id']]));continue
            tt=time.monotonic();i=row['grid_id']
            try:
                start=np.load(folder/'starts'/f'{i:04d}.npy');R,C,initial,rounds=optimize(start,folder,row);path=folder/'finals'/f'{i:04d}.npy';np.save(path,C);R,detail=true_cover(C);result=dict(row,initial_radius=float(np.sqrt(initial['maximum_squared_distance'])),initial_whole_disk_score=initial,final_radius=R,final_whole_disk_score=detail,final_diagnostic=diagnose(C),final_coordinate_file=str(path.relative_to(ROOT)),final_file_sha256=sha(path),final_coordinate_sha256=csha(C),rounds=rounds,seconds=time.monotonic()-tt,numerical_only=True)
            except Exception as e:result=dict(row,error=str(e),traceback=traceback.format_exc(),seconds=time.monotonic()-tt)
            records.append(result);fp.write(json.dumps(clean(result),ensure_ascii=False,allow_nan=False)+'\n');write(folder/'checkpoint.json',dict(worker=worker,completed=len(records),requested=len(rows),records=records,best_radius=min((r.get('final_radius',1) for r in records),default=1),seconds=time.monotonic()-tic,assigned_cpu=cpus[worker%len(cpus)],BLAS_threads=1));queue.put(dict(type='record',record=result))
    queue.put(dict(type='done',worker=worker))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--prepare-only',action='store_true');args=ap.parse_args();grid=json.loads((ROOT/'grid.json').read_text()) if (ROOT/'grid.json').exists() else prepare()
    if args.prepare_only:return
    executed=[r for r in grid if not r['skip']];cpus=sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else [0,1,2,3];write(ROOT/'runtime_provenance.json',dict(source_sha256={str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'source').glob('*.py'))},worker_cpu_assignments=[cpus[w%len(cpus)] for w in range(4)],workers=4,BLAS_threads_per_worker=1));ctx=mp.get_context('spawn');queue=ctx.Queue();processes=[ctx.Process(target=worker_main,args=(w,[r for r in executed if r['worker']==w],queue,cpus)) for w in range(4)]
    for p in processes:p.start()
    records=[];done=set();tic=time.monotonic();best=REFERENCE
    while len(done)<4:
        item=queue.get()
        if item['type']=='done':done.add(item['worker']);continue
        row=item['record'];records.append(row);R=row.get('final_radius',1)
        if R<THRESHOLD and REFERENCE-R>MATERIAL:
            q=ROOT/'frozen_proposals'/f'grid{row["grid_id"]:04d}_R{R:.17f}.npy';np.save(q,np.load(ROOT/row['final_coordinate_file']));write(q.with_suffix('.json'),dict(row,frozen_file_sha256=sha(q),frozen_coordinate_sha256=csha(np.load(q))));print('FROZEN_PROPOSAL',row['grid_id'],R,str(q),sha(q),flush=True)
        if R<best:
            best=R;src=ROOT/row['final_coordinate_file'];np.save(ROOT/'best.npy',np.load(src));write(ROOT/'best.json',row);print('GRID_BEST',row['grid_id'],best,flush=True)
        write(ROOT/'checkpoint.json',dict(completed=len(records),requested=len(executed),grid_budget=len(grid),skipped=len(grid)-len(executed),workers_completed=sorted(done),best_radius=best,records=records,seconds=time.monotonic()-tic,scope='Finite numerical grid; no exact full-cover or global-optimality claim.'))
        if len(records)%8==0:print('GRID_PROGRESS',len(records),len(executed),'best',best,'seconds',time.monotonic()-tic,flush=True)
    for p in processes:p.join()
    write(ROOT/'completion.json',dict(completed=len(records),requested=len(executed),best_radius=best,errors=sum('error' in r for r in records),workers_exitcode=[p.exitcode for p in processes],seconds=time.monotonic()-tic))
    print('GRID_DONE',len(records),best,'seconds',time.monotonic()-tic,flush=True)
if __name__=='__main__':main()
