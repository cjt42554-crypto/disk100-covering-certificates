# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Independent finite numerical replay of completed proxy paths, no optimizer."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,hashlib,itertools,time
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.spatial import Delaunay,ConvexHull
from independent_oracle import all_triple_cover
from whole_disk_oracle import true_cover
import event_jets,graph_tools

ROOT=Path(__file__).resolve().parent.parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ch(C):return hashlib.sha256(np.asarray(C,dtype='<f8').tobytes()).hexdigest()
def write(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def replay(item):
    path,expected=item;C=np.load(ROOT/path);tic=time.monotonic()
    assert sha(ROOT/path)==expected and C.shape==(100,2) and np.isfinite(C).all()
    R,detail=all_triple_cover(C);dynamic,dd=true_cover(C)
    return dict(file=path,file_sha256=expected,coordinate_sha256=ch(C),radius_all_triple=R,radius_dynamic=dynamic,absolute_radius_difference=abs(R-dynamic),passed=abs(R-dynamic)<2e-11,independent_oracle=detail,seconds=time.monotonic()-tic)
def strict_event_geometry(C,ids):
    jet=event_jets.event_jet(C,'incircle',ids);q=jet['center'];r2=jet['radius_squared'];P=C[ids]
    hull=ConvexHull(P);convex=len(hull.vertices)==4
    cycle=P[hull.vertices];insidequad=bool(np.all(np.cross(np.roll(cycle,-1,axis=0)-cycle,q-cycle)>=-2e-12))
    other=[i for i in range(100) if i not in ids];margin=np.sum((C[other]-q)**2,axis=1)-r2
    points=np.sum((P-q)**2,axis=1)-r2
    faces={tuple(sorted(map(int,t))) for t in Delaunay(C).simplices}
    a,b,c,d=ids;before=[tuple(sorted((a,b,c))),tuple(sorted((a,b,d)))];after=[tuple(sorted((a,c,d))),tuple(sorted((b,c,d)))]
    return dict(raw_event_value=float(jet['raw_value']),q=q.tolist(),rho_squared=float(r2),q_norm_squared=float(q@q),q_strict_inside_unit_disk=bool(q@q<1-1e-8),four_points_convex=convex,q_inside_quad_diagnostic=insidequad,minimum_other96_circle_margin=float(margin.min()),nearest_other_site=int(other[int(np.argmin(margin))]),other96_strict_nearest=bool(margin.min()>1e-8),support_radius_squared_spread=float(np.ptp(points)),before_adjacent_faces=[list(t) for t in before],after_adjacent_faces=[list(t) for t in after],before_faces_actual=all(t in faces for t in before),after_faces_actual=all(t in faces for t in after))
def main():
    audit=ROOT/'audit';audit.mkdir(exist_ok=True)
    result=json.loads((ROOT/'PROXY_BATCH_RESULT.json').read_text());assert len(result['paths'])==2
    references={};allrows=[];physical=[];statuses={}
    prior_ord=set(json.loads((ROOT/'prior_ordered_hashes.json').read_text()));prior_unord=set(json.loads((ROOT/'prior_unordered_hashes.json').read_text()))
    for summary in result['paths']:
        rp=ROOT/summary['result_file'];assert sha(rp)==summary['result_sha256'];d=json.loads(rp.read_text());target=d['target'];ids=target['event_labels']
        references[target['source_file']]=target['source_sha256']
        for a in d['attempts']:
            references[a['starts_file']]=a['starts_file_sha256']
            for row in a['rounds']:
                references[row['file']]=row['file_sha256'];statuses[str(row['status'])]=statuses.get(str(row['status']),0)+1
        for row in d['accepted']:references[row['file']]=row['file_sha256']
        selected={r['fraction']:r for r in d['accepted'] if abs(r['fraction']-.9995)<1e-12 or abs(r['fraction']-1)<1e-12 or abs(r['fraction']-1.0005)<1e-12}
        assert len(selected)==3
        gate=[]
        for f,row in sorted(selected.items()):
            C=np.load(ROOT/row['file']);geom=strict_event_geometry(C,ids);ordered=ch(C);unordered=ch(C[np.lexsort((C[:,1],C[:,0]))])
            gate.append(dict(fraction=f,file=row['file'],file_sha256=row['file_sha256'],coordinate_sha256=ordered,old_ordered_duplicate=ordered in prior_ord,old_unordered_duplicate=unordered in prior_unord,geometry=geom))
        pre,tie,post=gate
        valid=all(g['geometry']['q_strict_inside_unit_disk'] and g['geometry']['four_points_convex'] and g['geometry']['other96_strict_nearest'] and not g['old_ordered_duplicate'] and not g['old_unordered_duplicate'] for g in gate)
        valid=valid and pre['geometry']['raw_event_value']>0 and post['geometry']['raw_event_value']<0 and pre['geometry']['before_faces_actual'] and post['geometry']['after_faces_actual']
        physical.append(dict(target_id=target['target_id'],candidate_id=target['candidate_id'],labels=ids,source_sha256=target['source_sha256'],source_file=target['source_file'],frames=gate,numerical_related_physical_event_gate=valid,exact_binary64_chord_certificate_pending=True,implicit_curve_not_certified=True,whole_original_restores=0))
    unique={}
    filemeta=[]
    for path,expected in references.items():
        assert sha(ROOT/path)==expected;C=np.load(ROOT/path);cs=ch(C)
        unique.setdefault(cs,(path,expected));filemeta.append(dict(file=path,file_sha256=expected,coordinate_sha256=cs))
    write(audit/'SAVED_ARRAY_INVENTORY.json',dict(fields=filemeta,total_fields=len(filemeta),unique_coordinate_arrays=len(unique)))
    with ProcessPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(replay,item) for item in unique.values()]
        for f in as_completed(futures):
            allrows.append(f.result());write(audit/'CHECKPOINT.json',dict(completed=len(allrows),total=len(unique),passed=all(r['passed'] for r in allrows)))
    report=dict(status='PASS' if all(r['passed'] for r in allrows) and all(r['numerical_related_physical_event_gate'] for r in physical) else 'FAIL',numerical_only=True,no_new_upper_radius_record=True,all_triple_replays=allrows,total_fields=len(filemeta),unique_arrays=len(unique),physical_events=physical,status_inventory=statuses,proxy_solves=result['total_proxy_solves'],original_restores=0,max_independent_oracle_radius_difference=max(r['absolute_radius_difference'] for r in allrows))
    write(audit/'INDEPENDENT_REPLAY.json',report)
    write(ROOT/'PHYSICAL_EVENT_HANDOFF.json',dict(status=report['status'],review='Finite numerical true-diagram flip gates and all saved fields independently replayed; exact dyadic chord pending',events=physical,independent_report_sha256=sha(audit/'INDEPENDENT_REPLAY.json'),remaining_original_restore_budget=77))
    print(json.dumps(dict(status=report['status'],unique_arrays=len(unique),fields=len(filemeta),maximum_difference=report['max_independent_oracle_radius_difference'],proxy_solves=report['proxy_solves'],physical_events=len(physical)),indent=2),flush=True)
if __name__=='__main__':main()
