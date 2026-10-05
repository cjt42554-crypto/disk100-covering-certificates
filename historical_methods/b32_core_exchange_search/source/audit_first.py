# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Independent acceptance and geometry/support comparison of first6 exchanges."""
import engine as run
import novelty
from metric_validation import validate
from independent_oracle import all_triple_cover
from alignment import compare
from pathlib import Path
import json,time,collections
import numpy as np
ROOT=run.ROOT;BASE=ROOT.parent

def main():
    tic=time.monotonic();cp=json.loads((ROOT/'checkpoint.json').read_text());rows=cp['records'];grid=json.loads((ROOT/'grid.json').read_text());complete=json.loads((ROOT/'completion.json').read_text());runtime=json.loads((ROOT/'runtime_provenance.json').read_text());old=json.loads((ROOT/'prior_ordered_hashes.json').read_text());oldu=json.loads((ROOT/'prior_unordered_hashes.json').read_text());results=[];checks=[];failures=[];classifications=[];targetref={14:'b32_metric_seed_search/final_handoff/index14_grid0061.npy',15:'b32_metric_seed_search/final_handoff/index15_grid0070.npy'}
    assert cp['completed']==complete['completed']==6;assert complete['worker_exitcodes']==[0]*4
    runtimechecks={p:run.sha(ROOT/p)==h for p,h in runtime['source_sha256'].items()};assert all(runtimechecks.values())
    def score(kind,index,path,stored):
        R,d=all_triple_cover(np.load(ROOT/path));entry=dict(kind=kind,index=index,file=path,radius=R,stored_radius=stored,error=abs(R-stored),detail=d);results.append(entry)
        if abs(R-stored)>2e-11:failures.append(entry)
    seen=set();seenu=set()
    for row in rows:
        if'error'in row:failures.append(dict(grid_id=row['grid_id'],error=row['error']));continue
        i=row['grid_id'];p=ROOT/'workers'/f'worker{row["worker"]:02d}'/'starts'/f'{i:04d}.npy';C=np.load(p);S=np.load(ROOT/'inputs/first_handoff'/row['source_catalog_file']);mp={int(k):v for k,v in row['target_to_source_catalog_bijection'].items()};outside=sorted(set(range(100))-set(row['modified_target_labels']));cs=run.csha(C);us=novelty.unordered_hash(C);v=validate(C,row['target_graph']);sc=dict(grid_id=i,source_outside_modified_bitwise_preserved=np.array_equal(C[outside],S[[mp[j]for j in outside]]),source_retained_core67_bitwise_preserved=np.array_equal(C[[j for j in outside if j<68]],S[[mp[j]for j in outside if j<68]]),start_file_sha_match=run.sha(p)==row['start_file_sha256'],coordinate_sha_match=cs==row['initial_coordinate_sha256'],unordered_sha_match=us==row['initial_unordered_geometry_sha256'],old_ordered_duplicate=cs in old,old_unordered_duplicate=us in oldu,new_ordered_duplicate=cs in seen,new_unordered_duplicate=us in seenu,actual_target_metric_realized=v['target_metric_realized']);seen.add(cs);seenu.add(us);checks.append(sc)
        if not all(sc[k]for k in('source_outside_modified_bitwise_preserved','source_retained_core67_bitwise_preserved','start_file_sha_match','coordinate_sha_match','unordered_sha_match','actual_target_metric_realized'))or any(sc[k]for k in('old_ordered_duplicate','old_unordered_duplicate','new_ordered_duplicate','new_unordered_duplicate')):failures.append(sc)
        score('start',i,str(p.relative_to(ROOT)),row['initial_radius']);p=ROOT/row['final_coordinate_file'];C=np.load(p);assert run.sha(p)==row['final_file_sha256'];assert run.csha(C)==row['final_coordinate_sha256'];score('final',i,row['final_coordinate_file'],row['final_radius'])
        for z in row['rounds']:
            p=ROOT/z['coordinate_file'];assert run.sha(p)==z['file_sha256'];assert run.csha(np.load(p))==z['coordinate_sha256'];score('stage',f'{i}:{z["round"]}',z['coordinate_file'],z['radius'])
        cmp=compare(np.load(BASE/targetref[row['target_index']]),C,144);classification=dict(grid_id=i,target_index=row['target_index'],input_file_sha256=row['final_file_sha256'],radius=row['final_radius'],diagnostic=run.diagnose(C),full100_reference=targetref[row['target_index']],comparison=cmp,scope='Explicit all100 geometry and strict5e-14 active/inactive inventories. A close geometry with differing nearactive support is not asserted a strict basin identity.');classifications.append(classification)
    source=ROOT/'inputs/first_handoff/inputs/current_catalog10.npy';score('source',10,str(source.relative_to(ROOT)),run.true_cover(np.load(source))[0]);best=min(rows,key=lambda r:r['final_radius']);material=[r for r in rows if r['final_radius']<run.THRESHOLD and run.REFERENCE-r['final_radius']>run.MATERIAL]
    obstacles=dict(paths_with_topology_event=sum(any(z['topology_changed']for z in r['rounds'])for r in rows),paths_with_new_binding=sum(any(z['new_actual_binding_keys']for z in r['rounds'])for r in rows),actual_gap_paths=sum(any(z['true_radius_minus_frozen_objective']>1e-8 for z in r['rounds'])for r in rows),final_not_feasible=sum(not r['rounds'][-1]['feasible']for r in rows))
    report=dict(status='PASS'if not failures else'FAIL',completed=6,budget=96,optimizer_exceptions=sum('error'in r for r in rows),independent_alltriple_calls=len(results),independent_counts=dict(collections.Counter(r['kind']for r in results)),maximum_score_error=max(r['error']for r in results),best_radius=best['final_radius'],best_file=best['final_coordinate_file'],best_file_sha256=best['final_file_sha256'],material_improvements=len(material),final_status_counts=dict(collections.Counter(r['rounds'][-1]['status']for r in rows)),final_binding_counts=dict(collections.Counter(r['final_diagnostic']['binding_count']for r in rows)),start_checks=checks,runtime_source_checks=runtimechecks,geometry_classification=classifications,obstacles=obstacles,independent_results=results,failures=failures,seconds=time.monotonic()-tic,scope='Finite6 new exact-core-exchange metrics followed by original full100 numeric search; no exact upper or per-graph/global optimum claim.')
    run.write(ROOT/'first_batch_replay_report.json',report);print('FIRST6_AUDIT',json.dumps({k:v for k,v in report.items()if k not in('start_checks','runtime_source_checks','geometry_classification','independent_results','failures')}))

if __name__=='__main__':main()
