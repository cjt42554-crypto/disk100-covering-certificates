#!/usr/bin/env python3
"""Exact centerwise transport of the already proved closed-D covering."""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime, timezone
import hashlib
import json
import sys
import time

assert __debug__ and sys.flags.no_site and sys.dont_write_bytecode
sys.path.append('/usr/lib/python3/dist-packages')
import numpy as np
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
assert not (ROOT/'EXECUTION_STARTED.json').exists()
source=BASE/'radial_guarded_exact_cover'
assert sha(source/'MANIFEST_SHA256.json')=='e7ada40fe44ed923ef2d088fceb7ff72d7e05b217717d59cb3df2957d7f924ba'
for n,h in read(source/'MANIFEST_SHA256.json').items():assert sha(source/n)==h,n
q=read(source/'QUALIFIED_SOURCE_RECEIPT.json')
assert q['status']=='PASS' and q['exact_closed_unit_disk_cover'] is True
old=source/'inputs/centers.npy'
new=BASE/'radial_guarded_original_restore/runs/joint03_guarded_final/final.npy'
producer=BASE/'radial_guarded_original_restore/runs/joint03_guarded_final/RESULT.json'
assert sha(old)=='f4a4a88461396d45a96877f9a17c090299b8a82479fe598ebb6e15a7edaae565'
assert sha(new)=='ea5ebe01e8aea5eeb86d83bd72a9ca86071e0ef1ac2109b3501d0d3cf505c398'
assert sha(producer)=='37952b5e0915418e89cc3654cd4c32152d64b92123c3752a9c6fc6e715a5bca4'
assert read(producer)['final_file_sha256']==sha(new)
for p in [old,new]:
    a=np.load(p,allow_pickle=False)
    assert a.shape==(100,2) and a.dtype==np.dtype('<f8') and np.isfinite(a).all()
with (ROOT/'EXECUTION_STARTED.json').open('x') as f:f.write(json.dumps({'utc':datetime.now(timezone.utc).isoformat(),'proof_attempt':1,'source_manifest_sha256':sha(source/'MANIFEST_SHA256.json')})+'\n')
tic=time.monotonic()
A=np.load(old,allow_pickle=False);B=np.load(new,allow_pickle=False)
delta=F(1,100000000)
rows=[]
for i in range(100):
    differences=[F.from_float(float(B[i,j]))-F.from_float(float(A[i,j])) for j in range(2)]
    d2=sum(x*x for x in differences)
    assert d2<=delta*delta,i
    rows.append({'center':i,'difference_xy':[[x.numerator,x.denominator] for x in differences],
                 'distance_squared':[d2.numerator,d2.denominator], 'distance_squared_le_delta_squared':True})
cap=F(*q['exact_cover_radius_fraction'])+delta
assert cap<F('0.1165')
result={'status':'PASS_EXACT_CENTERWISE_TRANSPORT_CLOSED_D_COVER',
    'schema':'EXACT_BINARY64_SOURCE_COVER_BY_CENTERWISE_TRANSPORT_V1',
    'old_source_npy_sha256':sha(old),'new_source_npy_sha256':sha(new),
    'old_source_manifest_sha256':sha(source/'MANIFEST_SHA256.json'),
    'old_qualified_receipt_sha256':sha(source/'QUALIFIED_SOURCE_RECEIPT.json'),
    'producer_result_sha256':sha(producer),'centers':100,'exact_displacement_checks':100,
    'delta_fraction':[delta.numerator,delta.denominator],
    'exact_closed_D_cover_radius_fraction':[cap.numerator,cap.denominator],
    'exact_closed_D_cover_radius_decimal':'.116404861624',
    'exact_closed_D_cover_radius_lt_0p1165':True,
    'theorem':'For each x in closed D choose a source center a_i with |x-a_i|<=R; the corresponding new center b_i satisfies |a_i-b_i|<=delta by exact dyadic squared distance, so |x-b_i|<=R+delta.',
    'rows':rows,'seconds':time.monotonic()-tic,
    'new_global_candidate_oracle_calls':0,'optimizer_calls':0,'ordinary_calls':0,
    'old_cover_verifier_replays':0,'record_upper_improved':False,'global_optimum_proved':False}
with (ROOT/'RESULT.json').open('x') as f:f.write(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='rows'},sort_keys=True))
