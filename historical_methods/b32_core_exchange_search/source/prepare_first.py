# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Verify6 immutable real core-site/wall exchanges before original optimization."""
import engine as run
import novelty
from metric_validation import validate
from pathlib import Path
import json,copy,collections
import numpy as np
ROOT=run.ROOT;BASE=ROOT.parent;HANDOFF=ROOT/'inputs/first_handoff'

def main():
    ready=json.loads((HANDOFF/'READY.json').read_text());assert ready['status']=='READY_VERIFIED';manifest=json.loads((HANDOFF/'MANIFEST_SHA256.json').read_text());assert all(run.sha(HANDOFF/p)==h for p,h in manifest.items())
    h=json.loads((HANDOFF/'HANDOFF.json').read_text());catalog=json.loads((BASE/'b32_metric_seed_search/inputs/catalog_handoff/catalog.json').read_text());assert run.sha(BASE/'b32_metric_seed_search/inputs/catalog_handoff/catalog.json')==h['catalog_sha256'];graphs={r['independent_index']:r for r in catalog['graphs']};old=json.loads((ROOT/'prior_ordered_hashes.json').read_text());oldu=json.loads((ROOT/'prior_unordered_hashes.json').read_text());seen=set();seenu=set();grid=[];original=np.load(ROOT/'inputs/current146_100.npy');assert run.sha(ROOT/'inputs/current146_100.npy')==h['current_baseline_original_sha256']
    for i,r in enumerate(h['records']):
        C=np.load(HANDOFF/r['coordinate_relative_file']);assert run.sha(HANDOFF/r['coordinate_relative_file'])==r['file_sha256'];assert run.csha(C)==r['coordinate_sha256'];us=novelty.unordered_hash(C);assert us==r['unordered_geometry_sha256'];S=np.load(HANDOFF/r['source_catalog10_coordinate_relative_file']);assert run.sha(HANDOFF/r['source_catalog10_coordinate_relative_file'])==r['source_catalog10_coordinate_file_sha256'];mp={int(k):v for k,v in r['target_to_source_catalog_bijection'].items()};mo={int(k):v for k,v in r['target_to_source_original_label'].items()};assert set(mp)==set(mo)==set(range(100));assert set(mp.values())==set(mo.values())==set(range(100));changed=set(r['target_modified_labels']);outside=sorted(set(range(100))-changed);assert np.array_equal(C[outside],S[[mp[j]for j in outside]]);assert np.array_equal(S[[mp[j]for j in range(100)]],original[[mo[j]for j in range(100)]])
        actualchanged=sorted(j for j in range(100)if not np.array_equal(C[j],S[mp[j]]));assert actualchanged==sorted(changed);assert len(outside)==r['retained_real_centers_bitwise'];assert len(set(range(68))-changed)==67
        graph=copy.deepcopy(graphs[r['target_graph_index']]);orientation=r['expected_metric_orientation'];assert orientation in(-1,1)
        if orientation==-1:graph['faces']=[[a,c,b]for a,b,c in graph['faces']]
        v=validate(C,graph);v.pop('actual_graph',None);R,d=run.true_cover(C);cs=run.csha(C);reason=None
        if cs in old:reason='old_ordered_coordinates'
        elif us in oldu:reason='old_unordered_geometry'
        elif cs in seen or us in seenu:reason='new_exact_geometry_duplicate'
        elif not v['target_metric_realized']:reason='actual_target_metric_support_failure'
        assert abs(R-r['initial_whole_disk_radius'])<2e-11
        row=dict(grid_id=i,worker=i%4,owner_record=r,owner_coordinate_relative_file=r['coordinate_relative_file'],owner_file_sha256=r['file_sha256'],source_catalog_file=r['source_catalog10_coordinate_relative_file'],source_catalog_file_sha256=r['source_catalog10_coordinate_file_sha256'],source_original_file_sha256=h['current_baseline_original_sha256'],
            target_index=r['target_graph_index'],target_graph=graph,target_fingerprint=graph['fingerprint'],metric_orientation=orientation,original_catalog_faces_were_reversed=orientation==-1,
            target_to_source_catalog_bijection=mp,target_to_source_original_label=mo,modified_target_labels=actualchanged,retained_real_centers_bitwise=len(outside),retained_common_core_centers_bitwise=67,
            old_ordered_duplicate=cs in old,old_unordered_duplicate=us in oldu,initial_coordinate_sha256=cs,initial_unordered_geometry_sha256=us,initial_radius=R,initial_whole_disk_score=d,initial_target_and_actual_metric=v,
            owner_integer_target_face_strict_audit=r['integer_target_face_audit'],owner_whole_disk_independent_score=r['initial_independent_whole_disk_score'],skip=reason is not None,skip_reason=reason,offsets=0,all_center_variables=100)
        if reason is None:
            folder=ROOT/'workers'/f'worker{row["worker"]:02d}'
            for name in('starts','finals','stages','proposals'):(folder/name).mkdir(parents=True,exist_ok=True)
            p=folder/'starts'/f'{i:04d}.npy';np.save(p,C);row['start_file_sha256']=run.sha(p);seen.add(cs);seenu.add(us)
        grid.append(row)
    run.write(ROOT/'grid.json',grid);checks=dict(status='READY_VERIFIED',budget=96,owner_starts=len(grid),executable=sum(not r['skip']for r in grid),skip_reasons=dict(collections.Counter(r['skip_reason']for r in grid if r['skip'])),source_handoff_manifest_sha256=run.sha(HANDOFF/'MANIFEST_SHA256.json'),target_counts=dict(collections.Counter(r['target_index']for r in grid)),retained_real_centers=[r['retained_real_centers_bitwise']for r in grid],old_ordered_unique=len(old),old_unordered_unique=len(oldu),scope='Six exact integer-motif-supported new metrics; all100 nooffset original optimizer follows, no harmonic replay or radius certificate.');run.write(ROOT/'input_checks.json',checks);run.write(ROOT/'preparation_status.json',checks);print('PREPARED_EXCHANGE',json.dumps(run.clean(checks)))

if __name__=='__main__':main()
