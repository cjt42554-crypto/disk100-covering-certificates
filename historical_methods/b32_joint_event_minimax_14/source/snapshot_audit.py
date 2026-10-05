# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Numerical all-array and budget replay of the two finished new joint pilots.

No optimizer; independently scores every new coordinate array. Reuses only
the five byte-locked selected-frame receipts from the permanent handoff.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,ast,hashlib,time,collections,traceback
import numpy as np
from independent_oracle import all_triple_cover
from whole_disk_oracle import true_cover
from helpers import inventory,write,sha,csha
from selected_audit import circle
import graph_tools
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'audit'

def main():
    OUT.mkdir(exist_ok=True)
    receipt=json.loads((ROOT/'first_joint_handoff'/'INDEPENDENT_ENDPOINT_AUDIT.json').read_text())
    cache={r['coordinate_sha256']:r for r in receipt['rows']}
    files=sorted(list((ROOT/'workers').rglob('*.npy'))+list((ROOT/'inputs').glob('*.npy')))
    groups=collections.defaultdict(list)
    for p in files:
        C=np.load(p);assert C.shape==(100,2) and np.isfinite(C).all();groups[csha(C)].append(p)
    records=[];new_calls=0
    for index,(cs,paths) in enumerate(groups.items()):
        C=np.load(paths[0]);reuse=cs in cache
        if reuse:
            old=cache[cs];R=old['wholeD_radius_alltriple'];detail=old['independent_alltriple_oracle'];Rdt=old['wholeD_radius_DT_allpair'];dtd=old['independent_boundary_continuous_oracle']
        else:
            R,detail=all_triple_cover(C);Rdt,dtd=true_cover(C);new_calls+=1
        assert abs(R-Rdt)<2e-12
        try:graph=graph_tools.extract(C,False)
        except Exception as exc:graph=dict(status='NUMERICAL_GRAPH_UNCLASSIFIED',exception_type=type(exc).__name__,message=str(exc))
        g,J,keys,_=inventory(C);t=max(g.max(),R*R);slacks=t-g;active=np.flatnonzero(slacks<=5e-14)
        binding=sorted({int(v) for i in active for v in keys[i][1:(4 if keys[i][0]=='t' else 3 if keys[i][0]=='b' else 2)]})
        row=dict(coordinate_sha256=cs,files=[dict(file=str(p.relative_to(ROOT)),file_sha256=sha(p)) for p in paths],wholeD_radius_alltriple=float(R),wholeD_radius_DT_allpair=float(Rdt),independent_alltriple_oracle=detail,independent_boundary_continuous_oracle=dtd,alltriple_reused_from_permanent_selected_handoff=reuse,actual_graph=graph,inventory_keys=keys,inventory_squared_distances=g.tolist(),actual_inventory_vs_fullD_squared_gap=float(R*R-g.max()),active_tolerance_squared=5e-14,active_keys=[keys[i] for i in active],binding_center_indices_at_tolerance=binding,nonbinding_center_indices_at_tolerance=sorted(set(range(100))-set(binding)),nonbinding_does_not_imply_removable=True,numerical_only=True)
        records.append(row);write(OUT/'SNAPSHOT_CHECKPOINT.json',dict(completed_unique=index+1,total_unique=len(groups),new_alltriple_calls=new_calls,last_file=str(paths[0].relative_to(ROOT)),last_radius=R,numerical_only=True))
        if (index+1)%8==0:print(json.dumps(dict(unique=index+1,total=len(groups),new_calls=new_calls)),flush=True)
    bycs={r['coordinate_sha256']:r for r in records};byfile={x['file']:r for r in records for x in r['files']}
    old=set(json.loads((ROOT/'prior_ordered_hashes.json').read_text()));oldu=set(json.loads((ROOT/'prior_unordered_hashes.json').read_text()))
    paths=[];total_solves=0;total_nodes=0;new_predictors=[];source=byfile['inputs/joint02_source.npy']
    source_keys={tuple(k):float(g) for k,g in zip(source['inventory_keys'],source['inventory_squared_distances'])};t0=source['wholeD_radius_alltriple']**2
    incoming=[dict(key=list(k),source_slack=t0-v) for k,v in source_keys.items() if 5e-14<t0-v<=5e-9]
    for name in ['joint02','joint03']:
        result=json.loads((ROOT/'workers'/name/'RESULT.json').read_text());target=result['target'];C0=np.load(ROOT/target['source_file']);assert sha(ROOT/target['source_file'])==target['source_file_sha256'];w=np.column_stack((-C0[:,1],C0[:,0])).ravel();w/=np.linalg.norm(w)
        rounds=[s for a in result['attempts'] for s in a['rounds']]
        assert len(rounds)==result['proxy_solves']<=400 and len(result['accepted'])<=40 and result['failed_nodes']<=40
        total_solves+=len(rounds);total_nodes+=len(result['accepted']);previous=C0;prevkeys=set(map(tuple,source['inventory_keys']));prevgraph=source['actual_graph'].get('fingerprint');node_rows=[]
        for a in result['attempts']:
            p=ROOT/a['start_file'];C=np.load(p);assert sha(p)==a['start_file_sha256'];cs=csha(C);cu=csha(C[np.lexsort((C[:,1],C[:,0]))]);assert cs==a['coordinate_sha256'] and cu==a['unordered_sha256'];dup=cs in old or cu in oldu;assert dup==a['old_exact_duplicate'];new_predictors.append(dict(file=a['start_file'],coordinate_sha256=cs,unordered_sha256=cu,old_exact_duplicate=dup))
            base=C
            for s in a['rounds']:
                p=ROOT/s['file'];C=np.load(p);assert sha(p)==s['file_sha256'] and csha(C)==s['coordinate_sha256'];actual=bycs[csha(C)]['wholeD_radius_alltriple'];assert abs(actual-s['actual_radius'])<2e-12
                S=target['primary_triangle'];q,rho,_=circle(C,S);h=np.array([np.sum((C[i]-q)**2)-rho for i in target['added_sites_in_control_order']]);assert np.max(np.abs(h-np.array(s['raw_event_values'])))<2e-13
                assert np.linalg.norm(C-C0,axis=1).max()<=.2+1e-9 and np.linalg.norm(C-base,axis=1).max()<=.05+1e-9 and abs(w@(C.ravel()-C0.ravel()))<1e-9
                assert actual<=.1165 and s['status']==0 and s['success'];base=C
        for node in result['accepted']:
            p=ROOT/node['file'];C=np.load(p);assert sha(p)==node['file_sha256'] and csha(C)==node['coordinate_sha256'];r=bycs[csha(C)];keys=set(map(tuple,r['inventory_keys']));graph=r['actual_graph'].get('fingerprint');actualmap={tuple(k):float(g) for k,g in zip(r['inventory_keys'],r['inventory_squared_distances'])}
            diagnostic=[]
            for inc in incoming:
                key=tuple(inc['key']);g=actualmap.get(key);diagnostic.append(dict(key=list(key),source_slack=inc['source_slack'],still_actual_inventory=g is not None,actual_slack=None if g is None else r['wholeD_radius_alltriple']**2-g))
            node_rows.append(dict(index=node['index'],phase=node['phase'],goal=node['goal'],file=node['file'],coordinate_sha256=csha(C),actual_R=r['wholeD_radius_alltriple'],actual_graph_fingerprint=graph,restricted_graph_changed_from_previous=graph!=prevgraph,new_actual_inventory_keys=[list(k) for k in sorted(keys-prevkeys)],removed_actual_inventory_keys=[list(k) for k in sorted(prevkeys-keys)],source_incoming_contact_tracking=diagnostic,actual_active_count=len(r['active_keys']),actual_binding_count=len(r['binding_center_indices_at_tolerance']),max_motion_from_previous_accepted=float(np.linalg.norm(C-previous,axis=1).max()),no_exact_event_certificate=True))
            previous=C;prevkeys=keys;prevgraph=graph
        paths.append(dict(name=name,result_file_sha256=sha(ROOT/'workers'/name/'RESULT.json'),proxy_solves=len(rounds),accepted_nodes=len(result['accepted']),status_histogram=dict(collections.Counter(str(r['status']) for r in rounds)),failed_nodes=result['failed_nodes'],original_restores=result['original_restores'],termination=result['termination'],incoming_source_contacts=incoming,accepted_inventory_monitor=node_rows))
    assert total_solves==31 and total_nodes==31 and len(new_predictors)==31 and not any(r['old_exact_duplicate'] for r in new_predictors)
    assert len({r['coordinate_sha256'] for r in new_predictors})==31 and len({r['unordered_sha256'] for r in new_predictors})==31
    sourcecode=[]
    for p in sorted((ROOT/'source').glob('*.py')):ast.parse(p.read_text());sourcecode.append(dict(file=str(p.relative_to(ROOT)),sha256=sha(p),ast_pass=True))
    audit=dict(status='PASS_ALL_SAVED_CENTER_FIELDS_NUMERICAL_AND_RUNTIME',coordinate_fields=len(files),unique_arrays=len(records),independent_alltriple_calls_new=new_calls,permanent_selected_receipts_reused=len(cache),rows=records,path_reconciliation=paths,predictor_novelty=new_predictors,source_AST=sourcecode,total_proxy_solves=31,total_original_restores=0,all_triples_per_unique=161700,all_pairs_per_unique=4950,sampled_coverage_claim=False,no_strict_cover_or_event_claim=True,no_global_claim=True)
    write(OUT/'ALL_SNAPSHOTS_AUDIT.json',audit)
    print(json.dumps(dict(status=audit['status'],coordinate_fields=len(files),unique_arrays=len(records),new_alltriple_calls=new_calls,proxy_solves=31,original_restores=0)),flush=True)
if __name__=='__main__':main()
