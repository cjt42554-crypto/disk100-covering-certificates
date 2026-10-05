# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Verify prior receipts and independently audit only new copied center arrays."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,numpy as np
from concurrent.futures import ProcessPoolExecutor,as_completed
from proxy_driver import sha,csha,write
from independent_oracle import all_triple_cover
from whole_disk_oracle import true_cover
ROOT=Path(__file__).resolve().parent.parent
def run(item):
    p,field,expected=item;z=np.load(ROOT/p);C=z if field=='$' else z[field];assert csha(C)==expected;R,d=all_triple_cover(C);r,t=true_cover(C)
    return dict(file=p,field=field,coordinate_sha256=expected,file_sha256=sha(ROOT/p),all_triple_radius=R,dynamic_radius=r,absolute_difference=abs(R-r),passed=abs(R-r)<2e-11,independent_oracle=d)
def main():
    audit=ROOT/'final_audit';audit.mkdir(exist_ok=True);cache={};receipts=[]
    for p,key in [('audit/INDEPENDENT_REPLAY.json','all_triple_replays'),('audit_original/INDEPENDENT_RESTORE_AUDIT.json','all_triple_replays')]:
        d=json.load(open(ROOT/p));assert d['status']=='PASS';receipts.append(dict(file=p,file_sha256=sha(ROOT/p)))
        for r in d[key]:
            assert r['passed'] and sha(ROOT/r['file'])==r['file_sha256'];C=np.load(ROOT/r['file']);assert csha(C)==r['coordinate_sha256'];cache[r['coordinate_sha256']]=r
    fields=[];unique={}
    for p in sorted(ROOT.rglob('*')):
        if p.suffix not in ('.npy','.npz'):continue
        z=np.load(p);arrays=[('$',z)] if p.suffix=='.npy' else [(k,z[k]) for k in z.files]
        for key,C in arrays:
            if C.shape!=(100,2):continue
            assert np.isfinite(C).all();r=dict(file=str(p.relative_to(ROOT)),field=key,file_sha256=sha(p),coordinate_sha256=csha(C));fields.append(r);unique.setdefault(r['coordinate_sha256'],r)
    new=[(r['file'],r['field'],h) for h,r in unique.items() if h not in cache];records=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for f in as_completed([pool.submit(run,item) for item in new]):records.append(f.result());write(audit/'CHECKPOINT.json',dict(new_completed=len(records),new_total=len(new),prior_receipt_unique_arrays=len(cache)))
    allpassed=all(x['passed'] for x in records);write(audit/'SAVED_ARRAY_ACCEPTANCE.json',dict(status='PASS' if allpassed else 'FAIL',mathematical_scope='All stored100x2 center fields, including copied finite probes. Prior independent receipts verified by file+coordinateSHA; only newly copied arrays replayed, no repeated optimization or trajectory.',center_fields=len(fields),unique_arrays=len(unique),reused_verified_receipt_unique_arrays=len(set(unique)&set(cache)),new_independent_replays=len(records),new_replay_records=records,prior_receipts=receipts,center_field_inventory=fields,maximum_new_radius_oracle_difference=max((r['absolute_difference'] for r in records),default=0.),numerical_only=True))
    print(json.dumps(dict(status='PASS' if allpassed else 'FAIL',fields=len(fields),unique=len(unique),reuse=len(set(unique)&set(cache)),new=len(records)),indent=2))
if __name__=='__main__':main()
