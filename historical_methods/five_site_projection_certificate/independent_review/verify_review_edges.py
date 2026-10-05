# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Independent arithmetic/input-binding checks; no full-cover acceptor execution."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,importlib.util,json

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
SOURCE=ROOT/'verify_projection.py';EXPECTED='63a9c0f0379796c941dd068a7aa2c714a335143910b9c44300110c8322bdf6b4'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(SOURCE)==EXPECTED
spec=importlib.util.spec_from_file_location('projection_reviewed',SOURCE)
projection=importlib.util.module_from_spec(spec);spec.loader.exec_module(projection)
original_spec=projection.importlib.util.spec_from_file_location
cover_calls=[];mode='forbid';reference=None
class StopBeforeCover(Exception):pass
def hooked_spec(name,location,*args,**kw):
    spec=original_spec(name,location,*args,**kw)
    if name=='cover_acceptor':
        original_exec=spec.loader.exec_module
        def exec_module(module):
            original_exec(module)
            def cover_mock(rows,radius):
                cover_calls.append(dict(rows=len(rows),radius=radius,mode=mode))
                if mode=='reuse-existing-for-IO-only':
                    assert radius=='2' and rows==reference
                    return json.loads((ROOT/'regressions/positive_certificate.json').read_text())['source_closed_disk_cover']
                raise StopBeforeCover('Full cover execution deliberately intercepted')
            module.verify=cover_mock
        spec.loader.exec_module=exec_module
    return spec
projection.importlib.util.spec_from_file_location=hooked_spec
tests=[]
def rejected(label,callback,substring):
    before=len(cover_calls)
    try:callback()
    except ValueError as e:assert substring in str(e);tests.append(dict(case=label,status='PASS',expected_rejection=str(e)))
    else:raise AssertionError('Missing expected rejection: '+label)
    assert len(cover_calls)==before,'Acceptor invoked before malformed input rejection'

# The prior parser flattened 100 triples into 150 unrelated pairs including (0,0).
# The intended first-two-coordinate sites all miss x=(0,0) by distance at least 3.
rows=[['3','0','0'],['0','3','0'],['-3','0','0'],['0','-3','0'],['1.8','-2.4','0']]+[['4','4','0']]*95
malformed=HERE/'malformed_100_triples.json';malformed.write_text(json.dumps(rows)+'\n')
points=[tuple(F(x) for x in row[:2]) for row in rows]
q,r2=projection.circle(points,[0,1,2]);assert q==(0,0) and r2==9
vectors={i:(points[i][0],points[i][1]) for i in range(5)}
projection.ray_order(vectors)
assert all(projection.norm(points[i])>r2 for i in range(5,100))
flat=[F(x) for row in rows for x in row];paired=list(zip(flat[::2],flat[1::2]))
assert len(paired)==150 and (F(0),F(0)) in paired
assert min(projection.norm(p) for p in points)==9
rejected('malformed_three_coordinate_JSON_before_acceptor',lambda:projection.certify(malformed,[0,1,2],[3,4],'1'),'exactly two coordinates')
counterexample=dict(status='VALID_ARITHMETIC_COUNTEREXAMPLE_TO_PRIOR_DIMENSION_BINDING',
 prior_source_sha256='5d9cb7906d916a43f1808aa4ec1a0067a41f8f377174be44421d880d57c8f528',
 original_rows=100,original_coordinate_columns=3,intended_primary_radius_squared=[9,1],
 intended_minimum_origin_distance_squared=[9,1],claimed_radius=[1,1],flattened_acceptor_centers=150,
 flattened_acceptor_contains_origin=True,proof='The flattened set covers the unit disk by its origin-centered radius1 disk; intended100 planar centers miss the origin. Added centers are already on circle so displacement bound D=0.',
 full_cover_acceptor_executions=0,final_guard_rejected_before_cover_call=True)
(HERE/'dimension_counterexample.json').write_text(json.dumps(counterexample,indent=2)+'\n')

positive=ROOT/'regressions/positive.csv'
rejected('duplicate_selected_sites',lambda:projection.certify(positive,[0,1,2],[2,4],'2'),'distinct')
rejected('out_of_range_selector',lambda:projection.certify(positive,[0,1,2],[3,100],'2'),'outside range')
rejected('nearer_site_original_negative_receipt_replayed_geometry_only',lambda:projection.certify(ROOT/'regressions/nearer_blocker.csv',[0,1,2],[3,4],'2'),'blocks')
class ChangingInput:
    count=0
    def __str__(self):return str(positive)
    def read_bytes(self):
        self.count+=1;return positive.read_bytes()+(b'changed' if self.count>1 else b'')
rejected('source_bytes_change_during_parser_read',lambda:projection.certify(ChangingInput(),[0,1,2],[3,4],'2'),'changed during read')

for x,bits in [(F(0),0),(F(4),0),(F(1,9),100),(F(2),100),(F(3,7),3),(F(10**40+1,10**25+3),11)]:
    lo,hi=projection.sqrt_interval(x,bits);assert lo>=0 and lo*lo<=x<=hi*hi
tests.append(dict(case='six_sqrt_interval_zero_square_nonsquare_large_rational_edges',status='PASS'))
R=F(1,3);d=F(2,7);rl,rh=projection.sqrt_interval(R*R);dl,dh=projection.sqrt_interval(d*d)
D=max(rh-dl,dh-rl,F(0));assert abs(R-d)<=D
assert all(abs(x-y)<=D for x in (rl,rh) for y in (dl,dh))
tests.append(dict(case='radial_displacement_rectangle_bound',status='PASS'))
shift=(F(7,5),F(-2,7));tri=[(shift[0]+1,shift[1]),(shift[0],shift[1]+1),(shift[0]-1,shift[1])]
center,square=projection.circle(tri,[0,1,2]);assert center==shift and square==1
tests.append(dict(case='translated_exact_primary_circle_formula',status='PASS'))
rejected('same_ray_projection_collision',lambda:projection.ray_order({0:(F(1),F(0)),1:(F(2),F(0))}),'same circle ray')
rejected('center_on_polygon_boundary',lambda:projection.ray_order({0:(F(1),F(0)),1:(F(3,5),F(4,5)),2:(F(0),F(1)),3:(F(-3,5),F(4,5)),4:(F(-1),F(0))}),'not strictly inside')

# Final valid two-column geometry reaches the cover gate, which stays intercepted.
try:projection.certify(positive,[0,1,2],[3,4],'2')
except StopBeforeCover:tests.append(dict(case='valid_original_positive_geometry_reaches_cover_gate_without_executing_it',status='PASS'))
else:raise AssertionError('Full cover gate interception failed')

# Only test raw-byte and normalized-reader bindings: reuse the existing actual
# source-cover result on identical decoded rows, never execute another acceptor.
mode='reuse-existing-for-IO-only'
reader_spec=original_spec('read_only_cover_parser',projection.VERIFIER)
reader=projection.importlib.util.module_from_spec(reader_spec);reader_spec.loader.exec_module(reader)
reference,reader_sha=reader.read_centers(str(positive))
crlf=HERE/'positive_crlf_IO_only.csv';crlf.write_bytes(positive.read_bytes().replace(b'\n',b'\r\n'))
result=projection.certify(crlf,[0,1,2],[3,4],'2')
assert result['source_csv_sha256']==sha(crlf)
assert result['source_cover_reader_normalized_text_sha256']==reader_sha==sha(positive)
assert result['source_csv_sha256']!=reader_sha
assert result['maximum_center_displacement_upper']==[1,2] and result['cover_radius']==[5,2]
tests.append(dict(case='CRLF_source_byte_hash_and_reader_normalized_hash_distinct',status='PASS',
    source_cover_result_reused_only_for_IO_test=True,new_cover_acceptor_executions=0))
projection.importlib.util.spec_from_file_location=original_spec
assert sha(SOURCE)==EXPECTED
report=dict(status='PASS',reviewed_source_sha256=EXPECTED,arithmetic_and_input_tests=tests,
    tests_passed=len(tests),actual_full_cover_acceptor_executions=0,
    intercepted_cover_gate_calls=cover_calls,
    existing_positive_and_nearer_blocker_receipts_reused=True,
    new_optimizer_or_full_disk_or_old_math_runs=0,
    source_original_unchanged=True,
    dimension_counterexample='dimension_counterexample.json',
    CRLF_binding_check=dict(raw_byte_sha256=sha(crlf),normalized_text_sha256=reader_sha,
        source_cover_return_value_is_mocked_reuse_for_IO_only=True,no_new_mathematical_cover_acceptance_claimed=True))
(HERE/'EDGE_CHECK_REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(dict(status='PASS',tests=len(tests),source_sha256=EXPECTED,full_cover_acceptor_executions=0)))
