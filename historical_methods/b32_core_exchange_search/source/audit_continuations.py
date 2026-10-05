# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Independent full-D and actual graph/skeleton audit of3+10 continuations."""
import engine as own
import os,json,time,collections,multiprocessing as mp,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from independent_oracle import all_triple_cover
ROOT=own.ROOT;BASE=ROOT.parent
# Use the unchanged strict engine and complete-circle validator of branch001.
sys.path.insert(0,str(ROOT/'batches/batch001/source'))
import importlib.util
_spec=importlib.util.spec_from_file_location("general_metric_validation",ROOT/"batches/batch001/source/metric_validation.py")
_general=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_general);validate=_general.validate
from contact_geometry import contact_compare

def init(cpus):os.sched_setaffinity(0,{cpus[(mp.current_process()._identity[0]-1)%len(cpus)]})
def replay(task):
    kind,label,file,stored=task;R,d=all_triple_cover(np.load(ROOT/file));return dict(kind=kind,label=label,file=file,radius=R,stored_radius=stored,error=abs(R-stored),detail=d)

def main():
    tic=time.monotonic();failures=[];startchecks=[];stages=[];tasks=[];rows=[];runtimes=[];graphs=[]
    for branch in('batch001','batch002'):
        p=ROOT/'batches'/branch;cp=json.loads((p/'checkpoint.json').read_text());complete=json.loads((p/'completion.json').read_text());runtime=json.loads((p/'runtime_provenance.json').read_text());checks={q:own.sha(p/q)==h for q,h in runtime['source_sha256'].items()};runtimes.append(dict(branch=branch,source_SHA_checks=checks,completion=complete));assert all(checks.values())
        hm=json.loads((p/'inputs/owner_handoff/MANIFEST_SHA256.json').read_text());assert all(own.sha(p/'inputs/owner_handoff'/q)==h for q,h in hm.items());old=json.loads((p/'prior_ordered_hashes.json').read_text());oldu=json.loads((p/'prior_unordered_hashes.json').read_text());seen=set();seenu=set();sourcepaths=set()
        for r in cp['records']:
            if'error'in r:failures.append(dict(branch=branch,grid_id=r['grid_id'],error=r['error']));continue
            rows.append(dict(branch=branch,record=r));i=r['grid_id'];S=np.load(p/'inputs/owner_handoff'/r['source_catalog_file']);C=np.load(p/'workers'/f'worker{r["worker"]:02d}'/'starts'/f'{i:04d}.npy');cs=own.csha(C);us=own.csha(C[np.lexsort((C[:,1],C[:,0]))]);mpc={int(k):v for k,v in r['target_to_source_catalog_bijection'].items()};outside=sorted(set(range(100))-set(r['modified_target_labels']));v=validate(C,r['target_graph']);sc=dict(branch=branch,grid_id=i,outside_bitwise_retained=np.array_equal(C[outside],S[[mpc[j]for j in outside]]),common_core67_bitwise_retained=np.array_equal(C[[j for j in outside if j<68]],S[[mpc[j]for j in outside if j<68]]),source_file_SHA_match=own.sha(p/'inputs/owner_handoff'/r['source_catalog_file'])==r['source_catalog_file_sha256'],initial_sha_match=cs==r['initial_coordinate_sha256'],initial_unordered_sha_match=us==r['initial_unordered_geometry_sha256'],old_ordered_duplicate=cs in old,old_unordered_duplicate=us in oldu,new_ordered_duplicate=cs in seen,new_unordered_duplicate=us in seenu,actual_start_graph_verified=v['target_metric_realized']);seen.add(cs);seenu.add(us);startchecks.append(sc)
            if not all(sc[k]for k in('outside_bitwise_retained','common_core67_bitwise_retained','source_file_SHA_match','initial_sha_match','initial_unordered_sha_match','actual_start_graph_verified'))or any(sc[k]for k in('old_ordered_duplicate','old_unordered_duplicate','new_ordered_duplicate','new_unordered_duplicate')):failures.append(sc)
            sf=p/'inputs/owner_handoff'/r['source_catalog_file'];sourcepaths.add(str(sf.relative_to(ROOT)));file=str((p/'workers'/f'worker{r["worker"]:02d}'/'starts'/f'{i:04d}.npy').relative_to(ROOT));tasks.append(('start',branch+':'+str(i),file,r['initial_radius']));fp=p/r['final_coordinate_file'];assert own.sha(fp)==r['final_file_sha256'];assert own.csha(np.load(fp))==r['final_coordinate_sha256'];tasks.append(('final',branch+':'+str(i),str(fp.relative_to(ROOT)),r['final_radius']))
            for z in r['rounds']:
                fp=p/z['coordinate_file'];C=np.load(fp);assert own.sha(fp)==z['file_sha256'];assert own.csha(C)==z['coordinate_sha256'];vv=validate(C,r['target_graph']);stages.append(dict(branch=branch,grid_id=i,round=z['round'],actual_graph_match=vv.get('actual_fingerprint')==z['target_and_actual_metric'].get('actual_fingerprint')));assert stages[-1]['actual_graph_match'];tasks.append(('stage',branch+':'+str(i)+':'+str(z['round']),str(fp.relative_to(ROOT)),z['radius']))
        for file in sorted(sourcepaths):tasks.append(('source',branch,file,own.true_cover(np.load(ROOT/file))[0]))
    assert len(rows)==13;cpus=sorted(os.sched_getaffinity(0));results=[]
    with ProcessPoolExecutor(max_workers=4,initializer=init,initargs=(cpus,))as pool:
        for result in pool.map(replay,tasks):
            results.append(result)
            if result['error']>2e-11:failures.append(result)
    refs=[]
    for label,file in [('current146','b32_metric_seed_search/inputs/current146_100.npy'),('old65','b32_metric_seed_search/inputs/known99_reference.npy'),('new14','b32_metric_seed_search/final_handoff/index14_grid0061.npy'),('new15','b32_metric_seed_search/final_handoff/index15_grid0070.npy')]+[(f'event99_source{i}',f'topology_patch_hybrid_146/inputs/aligned_basin{i:02d}.npy')for i in range(3)]:refs.append((label,BASE/file))
    classifications=[];groups=[];deletions=[]
    for item in rows:
        branch=item['branch'];r=item['record'];p=ROOT/'batches'/branch;C=np.load(p/r['final_coordinate_file']);diag=own.diagnose(C);match=None;comparisons=[]
        for label,file in refs+groups:
            cmp=contact_compare(np.load(file),C);comparisons.append(dict(label=label,comparison=cmp))
            if cmp.get('matches'):match=label;break
        if match is None:
            for label,file in refs+groups:
                cmp=contact_compare(np.load(file),C,True)
                if cmp.get('matches'):comparisons.append(dict(label=label,comparison=cmp));match=label;break
        if match is None:
            match='finite_unmatched_'+branch+'_'+str(r['grid_id']);groups.append((match,p/r['final_coordinate_file']))
        classifications.append(dict(branch=branch,grid_id=r['grid_id'],file_sha256=r['final_file_sha256'],radius=r['final_radius'],binding=diag['binding_count'],free_indices=diag['nonbinding_indices'],active_count=diag['active_count'],final_actual_fingerprint=r['final_target_and_actual_metric'].get('actual_fingerprint'),initial_actual_fingerprint=r['optimization_initial_actual_fingerprint'],intended_target=r['intended_target_index'],initial_intended_achieved=r['intended_target_graph_achieved'],label=match,comparisons=comparisons,
            topology_events=sum(z['topology_changed']for z in r['rounds']),new_binding_events=sum(bool(z['new_actual_binding_keys'])for z in r['rounds']),maximum_actual_minus_frozen=max(z['true_radius_minus_frozen_objective']for z in r['rounds']),last_SLSQP_status=r['rounds'][-1]['status'],scope='Geometry+strict actual support inventories, finite failed matching proves no inequivalence or local minimum.'))
        for free in diag['nonbinding_indices']:
            D=np.delete(C,free,axis=0);folder=ROOT/'diagnostics/nonbinding_deletions';folder.mkdir(parents=True,exist_ok=True);f=folder/f'{branch}_grid{r["grid_id"]:04d}_remove{free:03d}.npy';np.save(f,D);R,d=all_triple_cover(D);deletions.append(dict(branch=branch,grid_id=r['grid_id'],free_index=free,source_sha256=r['final_file_sha256'],deleted_file=str(f.relative_to(ROOT)),deleted_file_sha256=own.sha(f),full100_radius=r['final_radius'],deleted99_radius=R,actual_deletion_debt=R-r['final_radius'],independent_detail=d,scope='Diagnostic only; numerical nonbinding is never treated as legal redundancy or permission to remove a covering disk.'))
    material=[x for x in rows if x['record']['final_radius']<own.THRESHOLD and own.REFERENCE-x['record']['final_radius']>own.MATERIAL];obs=dict(paths_with_topology_event=sum(r['topology_events']>0 for r in classifications),paths_with_new_binding=sum(r['new_binding_events']>0 for r in classifications),actual_gap_paths=sum(r['maximum_actual_minus_frozen']>1e-8 for r in classifications));report=dict(status='PASS'if not failures else'FAIL',optimized=13,budget=96,optimizer_exceptions=0,orchestration_failure_count=1,orchestration_failure_recovered_without_optimizer=True,independent_alltriple_calls=len(results)+len(deletions),main_independent_calls=len(results),independent_counts=dict(collections.Counter(r['kind']for r in results)),maximum_score_error=max(r['error']for r in results),start_checks=startchecks,stage_checks=stages,runtimes=runtimes,results=results,classifications=classifications,nonbinding_deletion_diagnostics=deletions,material_improvements=len(material),obstacles=obs,failures=failures,seconds=time.monotonic()-tic,scope='Mathematical array acceptance PASS is separate from the preserved branch001 empty-worker orchestration failure. No exact/global covering, deletion redundancy, or strict local-minimum assertion.');own.write(ROOT/'continuations_replay_report.json',report);print('CONTINUATIONS_AUDIT',json.dumps({k:report[k]for k in('status','optimized','independent_alltriple_calls','maximum_score_error','orchestration_failure_count','material_improvements','obstacles','seconds')}));print('CLASSIFICATIONS',[(r['branch'],r['grid_id'],r['radius'],r['binding'],r['label'])for r in classifications])

if __name__=='__main__':main()
