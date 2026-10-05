# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Immutable review/source gates and complete local trace hash supplement."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import ast,json
import numpy as np
from proxy_driver import ROOT,sha,csha,write
import event_jets

def main():
    review=ROOT.parent/'b32_event_targeting_14_review';ready=json.loads((review/'READY.json').read_text());assert ready['status']=='PERMANENTLY_SEALED_SOURCE_JET_REVIEW'
    assert sha(review/'MANIFEST_SHA256.json')==ready['manifest_sha256']
    m=json.loads((review/'MANIFEST_SHA256.json').read_text())
    for f,h in m.items():assert sha(review/f)==h
    targets=json.loads((review/'QUALIFIED_TARGETS.json').read_text())['targets'];assert len(targets)==2
    jetssource=review/'inputs/owner_first_targets/TARGET_JETS.npz';jetslocal=ROOT/'inputs/TARGET_JETS.npz';jetslocal.write_bytes(jetssource.read_bytes())
    pinned=ROOT/'inputs/pinned14.npy';pinned.write_bytes((review/'inputs/pinned14.npy').read_bytes())
    prepared=[]
    for row in targets:
        source=review/row['source_file'];assert sha(source)==row['source_sha256'];local=ROOT/'inputs'/f'{row["target_id"]}_source.npy';local.write_bytes(source.read_bytes());C=np.load(local);jet=event_jets.event_jet(C,row['event_kind'],row['event_labels'])
        assert abs(jet['raw_value']-row['source_raw_event_value'])<1e-14 and abs(jet['gradient_norm']-row['source_raw_event_gradient_norm'])<1e-13
        prepared.append(dict(row,source_file=str(local.relative_to(ROOT)),pinned_reference_file='inputs/pinned14.npy',ranking_direction_file='inputs/TARGET_JETS.npz',ranking_direction_file_sha256=sha(jetslocal),ranking_direction_index=row['owner_target_id']))
    old=ROOT.parent/'b32_core_exchange_search';ordered=json.loads((old/'batches/batch002/prior_ordered_hashes.json').read_text());unordered=json.loads((old/'batches/batch002/prior_unordered_hashes.json').read_text());tracechecks=[]
    for dirname in ['b32_core_exchange_search','b32_weak_contact_continuation','b32_weak_contact_event_extension','b32_event_targeting_14_review']:
        p=ROOT.parent/dirname;manifest=json.loads((p/'MANIFEST_SHA256.json').read_text());fields=0
        for relative,h in manifest.items():
            f=p/relative
            if f.suffix not in ('.npy','.npz'):continue
            assert sha(f)==h
            if f.suffix=='.npy':arrays=[('$',np.load(f,allow_pickle=False))]
            else:
                with np.load(f,allow_pickle=False) as z:arrays=[(k,z[k]) for k in z.files]
            for field,C in arrays:
                if C.shape!=(100,2):continue
                # Direction arrays need no covering claim; retaining their hashes
                # is harmless conservative duplicate rejection for predictors.
                ch=csha(C);uh=csha(C[np.lexsort((C[:,1],C[:,0]))]);ordered.setdefault(ch,[]).append(str(f)+':'+field);unordered.setdefault(uh,[]).append(str(f)+':'+field);fields+=1
        tracechecks.append(dict(directory=dirname,manifest_sha256=sha(p/'MANIFEST_SHA256.json'),retained100x2_fields=fields))
    write(ROOT/'prior_ordered_hashes.json',ordered);write(ROOT/'prior_unordered_hashes.json',unordered)
    source_manifest={str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'source').glob('*')) if p.is_file()}
    for p in (ROOT/'source').glob('*.py'):ast.parse(p.read_text())
    write(ROOT/'EXECUTION_SOURCE_SNAPSHOT.json',dict(source_files=source_manifest,review_manifest_sha256=ready['manifest_sha256'],retained_ordered_hashes=len(ordered),retained_unordered_hashes=len(unordered),trace_checks=tracechecks,worker_limit=2,proxy_launch_authorization='Root explicit: after reviewPASS directly first2pilot; original zero-event restores remain gated.'))
    write(ROOT/'PREPARED_TARGETS.json',dict(status='READY_AFTER_INDEPENDENT_TARGET_REVIEW',review_manifest_sha256=ready['manifest_sha256'],targets=prepared,source_snapshot_sha256=sha(ROOT/'EXECUTION_SOURCE_SNAPSHOT.json')))
    write(ROOT/'NOT_LAUNCHED.json',dict(status='SOURCE_REVIEW_GATES_PASSED_READY_TO_LAUNCH_FIRST2',proxy_solves=0,new_unconstrained_cover_solves=0,maximum_initial_pilots=2))
    print('TWO_PILOTS_SOURCE_REVIEW_GATES_PASS',ready['manifest_sha256'],len(ordered),len(unordered),flush=True)
if __name__=='__main__':main()
