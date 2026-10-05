# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Strict inventory supplement and immutable accepted metric-frame handoff."""
import engine as run
from pathlib import Path
import json,hashlib,shutil,collections
import numpy as np
from scipy.optimize import nnls
ROOT=run.ROOT

def main():
    report=json.loads((ROOT/'replay_report.json').read_text());assert report['status']=='PASS'
    cp=json.loads((ROOT/'checkpoint.json').read_text());rows=cp['records'];stats=json.loads((ROOT/'graph_metric_results.json').read_text());catalog=json.loads((ROOT/'inputs/catalog_handoff/catalog.json').read_text());graphs={g['independent_index']:g for g in catalog['graphs']};strict=[];nnlsrows=[]
    for r in rows:
        C=np.load(ROOT/r['final_coordinate_file']);g,q,J,kk,sites,top=run.inventory(C);ac=np.flatnonzero(g.max()-g<=5e-14);used=sorted({i for k in ac for i in sites[k]});near=np.flatnonzero(g.max()-g<1e-9);extra=sorted(set(near)-set(ac))
        strict.append(dict(grid_id=r['grid_id'],independent_index=r['independent_index'],file_sha256=r['final_file_sha256'],strict_active_count=len(ac),strict_binding_count=len(used),strict_free_indices=sorted(set(range(100))-set(used)),legacy_extra_nearactive=[dict(key=kk[k],squared_slack=float(g.max()-g[k]))for k in extra]))
    run.write(ROOT/'screening/strict5e14_final_binding_inventory.json',dict(rows=strict,legacy_squared_tolerance=1e-9,strict_squared_tolerance=5e-14,scope='Numerical binding classification at explicit slack tolerances, not exact contact equality or disk removability.'))
    for r in stats['graphs']:
        C=np.load(ROOT/r['best_file']);g,q,J,kk,sites,top=run.inventory(C);ac=np.flatnonzero(g.max()-g<=5e-14);M=np.r_[J[ac].T,np.ones((1,len(ac)))];rhs=np.r_[np.zeros(200),1];lam,res=nnls(M,rhs,maxiter=20000)
        nnlsrows.append(dict(independent_index=r['independent_index'],input_file=r['best_file'],input_file_sha256=r['best_file_sha256'],strict_active_count=len(ac),strict_binding_count=len({i for a in ac for i in sites[a]}),NNLS_residual_inf=float(np.max(np.abs(M@lam-rhs))),NNLS_lambda_min=float(lam.min()),NNLS_lambda_max=float(lam.max()),active_keys=[kk[a]for a in ac],lambda_values=lam,squared_tolerance=5e-14))
    run.write(ROOT/'screening/strict5e14_graph_best_NNLS.json',dict(rows=nnlsrows,scope='Actual numerical near-maximum inventory at5e-14 squaredslack; no exact equality, stationarity, or strictSOSC assertion.'))
    r14=min((r for r in rows if r['independent_index']==14),key=lambda r:r['final_radius']);C=np.load(ROOT/r14['final_coordinate_file']);legacy=json.loads((ROOT/'screening/graph14_best_diagnostic.json').read_text());g,q,J,keys,sites,top=run.inventory(C);ac=np.flatnonzero(g.max()-g<1e-9);lam=np.array(legacy['lambda_values']);stress=np.zeros(100)
    for k,l in zip(ac,lam):stress[sites[k]]+=l
    variation=[]
    for r in sorted((r for r in rows if r['independent_index']==14),key=lambda r:r['grid_id']):
        A=np.load(ROOT/r['final_coordinate_file']);u,s,vh=np.linalg.svd(A.T@C);Q=u@vh;dist=np.linalg.norm(A@Q-C,axis=1);ix=np.argsort(dist)[::-1]
        variation.append(dict(grid_id=r['grid_id'],file_sha256=r['final_file_sha256'],radius=r['final_radius'],rotation=Q,rotation_determinant=float(np.linalg.det(Q)),all100_max_distance=float(dist.max()),top_moving_centers=[dict(index=int(i),distance=float(dist[i]),baseline_incident_NNLS_stress=float(stress[i]))for i in ix[:8]],active_keys_equal=r['final_diagnostic']['active_signature']==r14['final_diagnostic']['active_signature'],inactive_keys_equal=r['final_diagnostic']['inactive_signature']==r14['final_diagnostic']['inactive_signature']))
    run.write(ROOT/'screening/graph14_metric_variation_diagnostic.json',dict(reference_grid_id=r14['grid_id'],reference_file_sha256=r14['final_file_sha256'],reference_radius=r14['final_radius'],NNLS_incident_stress=stress,lambda_below1e7_rows=[dict(key=keys[ac[i]],weight=float(lam[i]))for i in range(len(lam))if lam[i]<1e-7],variations=variation,scope='Same catalog labels and optimal rotation, all100 retained; nearlyzero numeric stress proves no redundancy, exact zero mode, or common basin.'))
    clarification=dict(legacy_inventory_squared_tolerance=1e-9,strict_inventory_squared_tolerance=5e-14,strict_final_binding_counts=dict(collections.Counter(r['strict_binding_count']for r in strict)),
        graph14_all100_same_catalog_label_rotational_distance_max=max(r['all100_max_distance']for r in variation),
        note='The17 finite unmatched contact-inventory labels are not17 proved physical basins. In particular graph14 has four legacy contact labels while all six all100 metrics align within7.33e-7; explicit near-active slacks account for inventory sensitivity. No free or redundant circle follows from tiny NNLS weights.')
    for name in('RESULT.json','replay_report.json','contact_family_classification.json'):
        d=json.loads((ROOT/name).read_text());d['strict_inventory_and_metric_clarification']=clarification;run.write(ROOT/name,d)
    chosen=[r14,max((r for r in rows if r['independent_index']==14),key=lambda r:next(v['all100_max_distance']for v in variation if v['grid_id']==r['grid_id'])),min((r for r in rows if r['independent_index']==15),key=lambda r:r['final_radius'])]
    target=ROOT/'final_handoff';target.mkdir(exist_ok=True);handoff=[]
    accepted={r['index']:r for r in report['independent_results']if r['kind']=='final'}
    for r in chosen:
        assert accepted[r['grid_id']]['absolute_error']<2e-11;src=ROOT/r['final_coordinate_file'];assert run.sha(src)==r['final_file_sha256'];dst=target/f'index{r["independent_index"]:02d}_grid{r["grid_id"]:04d}.npy';shutil.copy2(src,dst);graph=graphs[r['independent_index']]
        handoff.append(dict(independent_index=r['independent_index'],grid_id=r['grid_id'],file=str(dst.relative_to(ROOT)),file_sha256=run.sha(dst),coordinate_sha256=run.csha(np.load(dst)),radius=r['final_radius'],independent_alltriple_acceptance=accepted[r['grid_id']],
            graph_fingerprint=graph['fingerprint'],catalog_file_sha256=run.sha(ROOT/'inputs/catalog_handoff/catalog.json'),target_labelled_faces_sha256=hashlib.sha256(json.dumps(graph['faces'],separators=(',',':')).encode()).hexdigest(),final_label_to_catalog_label=list(range(100)),final_label_to_canonical_label=graph['vertex_to_canonical'],
            actual_labelled_target_faces_equal=r['final_target_and_actual_metric']['actual_target_faces_equal'],actual_outer_circle_arcs_equal=r['final_target_and_actual_metric']['actual_outer_circle_arc_edges_equal'],raw_seed_file=r['owner_seed_file'],raw_seed_file_sha256=r['owner_file_sha256'],
            strict_diagnostic=next(x for x in nnlsrows if x['independent_index']==r['independent_index'])if r['grid_id']in [r14['grid_id'],chosen[-1]['grid_id']]else None,
            scope='Pinned100-center frame only. Graph14 spread representative is a nearby metric with support-tolerance differences, not a new stable basin. Identity labels verified by actual166face equality.'))
    run.write(target/'HANDOFF.json',dict(status='SEALED',records=handoff,source_handoff_manifest_sha256=run.sha(ROOT/'inputs/catalog_handoff/MANIFEST_SHA256.json'),original_final_files_permanently_immutable=True,scope='Independent numeric alltriple+continuousboundary accepted frames for future exact-integer core-exchange geometry. No new optimizer started, no exact/global covering assertion.'))
    m={str(p.relative_to(target)):run.sha(p)for p in sorted(target.rglob('*'))if p.is_file()and p.name not in('MANIFEST_SHA256.json','READY.json')};run.write(target/'MANIFEST_SHA256.json',m);run.write(target/'READY.json',dict(status='SEALED',files=len(m),manifest_sha256=run.sha(target/'MANIFEST_SHA256.json')))
    print('PINNED_HANDOFF',run.sha(target/'MANIFEST_SHA256.json'),json.dumps([(r['independent_index'],r['grid_id'],r['file_sha256'])for r in handoff]))

if __name__=='__main__':main()
