# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Independent all-triple replay and immutable numerical evidence seal."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,hashlib,time,ast
import numpy as np
from independent_oracle import all_triple_cover
from whole_disk_oracle import true_cover
import graph_tools,branches
from diagnose import ROOT,write,sha
from continue_cones import csha

def main():
    start=time.monotonic();stages=[]
    for fn in sorted((ROOT/'stages').glob('curve*/RESULT.json')):
        r=json.loads(fn.read_text());stages.extend(r['stages'])
    for fn in sorted((ROOT/'all_contact_extension/stages').glob('curve*/RESULT.json')):
        r=json.loads(fn.read_text());stages.extend(r['stages'])
    expected={r['coordinate_file']:r for r in stages};configs=[];types=[];directions=[]
    for f in sorted(ROOT.rglob('*.npy')):
        C=np.load(f,allow_pickle=False);relative=str(f.relative_to(ROOT))
        if C.shape!=(100,2):continue
        if relative.startswith('matrices/'):
            directions.append(dict(file=relative,file_sha256=sha(f),array_coordinate_sha256=csha(C),role='direction_not_a_covering_configuration'));continue
        if relative.startswith('stages/') or relative.startswith('all_contact_extension/stages/'):
            assert relative in expected,relative
        configs.append((relative,'$',C,sha(f)))
    for f in sorted(ROOT.rglob('*.npz')):
        with np.load(f,allow_pickle=False) as z:
            for key in z.files:
                C=z[key]
                if C.shape!=(100,2):continue
                relative=str(f.relative_to(ROOT))
                if key=='C':configs.append((relative,key,C.copy(),sha(f)))
                else:directions.append(dict(file=relative,field=key,file_sha256=sha(f),array_coordinate_sha256=csha(C),role='direction_not_a_covering_configuration'))
    seen={};audited=[];errors=[]
    for i,(name,field,C,filehash) in enumerate(configs):
        ch=csha(C)
        if ch not in seen:
            R,detail=all_triple_cover(C);dynamic,dyn_detail=true_cover(C)
            seen[ch]=dict(independent_radius=R,dynamic_radius=dynamic,error=abs(R-dynamic),detail=detail,passed=abs(R-dynamic)<2e-11)
        result=seen[ch];row=dict(file=name,field=field,file_sha256=filehash,coordinate_sha256=ch,**result)
        if name in expected:
            prior=expected[name];assert prior['file_sha256']==filehash and prior['coordinate_sha256']==ch
            row['stored_radius']=prior['actual_radius'];row['stored_radius_error']=abs(result['independent_radius']-prior['actual_radius'])
            gg=graph_tools.extract(C,True);row['restricted_graph_fingerprint']=gg['fingerprint'];row['stored_graph_matches']=gg['fingerprint']==prior['actual_graph']['fingerprint']
            row['stored_model_gap_replayed']=result['independent_radius']-prior['model_radius']
            row['passed']=bool(row['passed'] and row['stored_radius_error']<2e-11 and row['stored_graph_matches'])
        audited.append(row)
        if not row['passed']:errors.append(row)
        if (i+1)%20==0:write(ROOT/'audit/CHECKPOINT.json',dict(completed=len(audited),total=len(configs),unique_oracle_calls=len(seen),failures=errors));print('AUDIT',len(audited),'of',len(configs),'unique',len(seen),flush=True)
    # Verify every immutable external source referenced by the own source bindings.
    sourcechecks=[]
    binding=json.loads((ROOT/'SOURCE_BINDING.json').read_text())
    for row in binding:
        actual=sha(ROOT.parent/row['source']);sourcechecks.append(dict(path=row['source'],expected=row['sha256'],actual=actual,passed=actual==row['sha256']))
    for relative in ['inputs/owner_modes']:
        p=ROOT/relative;m=json.loads((p/'MANIFEST_SHA256.json').read_text())
        for name,h in m.items():sourcechecks.append(dict(path=str((p/name).relative_to(ROOT)),expected=h,actual=sha(p/name),passed=sha(p/name)==h))
    olddirs=['rattler_contact_transfer','best147_boundary_transfer','next146_boundary_transfer','paired_transport_146','topology_patch_hybrid_146','fresh_rattler_transfer_99','cavity_repair_transfer_99','b32_metric_seed_search','b32_core_exchange_search','event_geometry_146','fresh_rattler_geometry_bridge','fresh99_second_order_diagnostic']
    frozenchecks=[]
    for name in olddirs:
        p=ROOT.parent/name;m=json.loads((p/'MANIFEST_SHA256.json').read_text());bad=[]
        for relative,h in m.items():
            if not (p/relative).is_file() or sha(p/relative)!=h:bad.append(relative)
        frozenchecks.append(dict(path=name,manifest_sha256=sha(p/'MANIFEST_SHA256.json'),payload_files=len(m),passed=not bad,failures=bad))
    static=[]
    for fn in sorted((ROOT/'source').glob('*.py')):
        ast.parse(fn.read_text());static.append(dict(file=str(fn.relative_to(ROOT)),sha256=sha(fn),parse_passed=True))
    audit=dict(status='PASS' if not errors and all(x['passed'] for x in sourcechecks+frozenchecks) else 'FAIL',config_fields=len(configs),independent_unique_alltriple_calls=len(seen),stage_fields=len(stages),maximum_oracle_error=max(x['error'] for x in seen.values()),audit_rows=audited,direction_type_ledger=directions,errors=errors,source_checks=sourcechecks,old_frozen_directory_checks=frozenchecks,source_AST_checks=static,elapsed_seconds=time.monotonic()-start,numerical_only=True,continuous_boundary_method='All4950 pair bisector/unit-circle crossings (both signs), all100 antipodes and axes, actual nearest100 squared distance; no sampled acceptance.')
    write(ROOT/'audit/INDEPENDENT_ACCEPTANCE.json',audit);assert audit['status']=='PASS'
    diag=json.loads((ROOT/'DIAGNOSTIC.json').read_text());inc=json.loads((ROOT/'INCOMING_CONTACT_REVIEW.json').read_text());first=json.loads((ROOT/'CONTINUATION_RESULT.json').read_text());ext=json.loads((ROOT/'all_contact_extension/CONTINUATION_RESULT.json').read_text())
    final=dict(status='NUMERICAL_EVIDENCE_PASS_NO_NEW_UPPER',input_sha256=diag['input_sha256'],source_radius=diag['whole_radius'],current_strict_upper_reference=.1154469548855829,cover_optimizer_runs=0,previous_structure_optimizer_runs=19,raw_pose_rank_sensitive=True,rank_review=inc['rank_review'],tiny_positive_lambda_count=diag['weak_count'],weak_lambda_threshold_heuristic=1e-7,near_inactive_count=inc['near_inactive_count'],complete_strong_tangent_dimension=inc['strong_complete_dimension'],complete_strong_Hmin=inc['strong_complete_Hmin'],FD_h1e6_error_complete=inc['rank_review'][1]['rank_tolerance_sweep'][2]['FD_h1e6_operator_error'],initial146_only_cone_trials=12,incoming22_cone_trials=12,incoming22_failed_trial_count=sum(not x['success'] for x in inc['finite_cone_trials']),incoming22_three_distinct_rays=True,mirror_equivalence_max_center_error=json.loads((ROOT/'CURVE01_GEOMETRY_COMPARISON.json').read_text())['maxdistance'],initial_continuation_curves=3,physical_mirror_deduplicated_extensions=2,initial_stages=sum(len(r['stages']) for r in first['curves']),extended_stages=sum(len(r['stages']) for r in ext['curves']),total_actual_stage_count=len(stages),pure_geometry_Newton_failures=sum(len(r['failures']) for r in first['curves']+ext['curves']),restricted_or_circle_events=sum(x['restricted_DT_changed'] or x['circle_inventory_changed'] for x in stages),extension_endpoint_records=[dict(curve=r['curve'],termination=r['termination'],phase=r['stages'][-1]['phase'],actual_radius=r['stages'][-1]['actual_radius'],file=r['stages'][-1]['coordinate_file'],file_sha256=r['stages'][-1]['file_sha256'],graph_fingerprint=r['stages'][-1]['actual_graph']['fingerprint'],actual_minus_model_radius=r['stages'][-1]['actual_minus_model_radius']) for r in ext['curves']],no_optimizer_qualified_topology_event=True,no_new_U=True,independent_audit=dict(status=audit['status'],unique_calls=len(seen),config_fields=len(configs),maximum_error=audit['maximum_oracle_error'],report_sha256=sha(ROOT/'audit/INDEPENDENT_ACCEPTANCE.json')),scope='Finite numerical weak/incoming-contact curvature and continuation diagnosis. All-positive exact criticality is distinguished from the tiny-lambda and incoming relaxations. No exact-zero stress, strict SOSC, optimal family, exhaustion or global claim. Raw near-null dimension is explicitly retained. All100 centers kept; no disk deletion, no covering optimizer and no restarts.')
    write(ROOT/'RESULT.json',final)
    (ROOT/'README.txt').write_text('Numerical B32 index14 weak-contact diagnosis and two physically distinct repaired continuations.\nInput: pinned14 grid61 SHA78d5a0ae08b07010da7c52e9634d1a14ecdcf12cf0e04ff5d1895895231452c0.\nOriginal t=R^2 analytic stress Hessian agrees with three finite-difference steps. Raw spatial rank is tolerance sensitive: complete augmented critical/gauge tangents have54/all-positive and68/strong dimensions, not53/67; added near-null curvature is positive numerically.\n14 tiny positive multipliers remain positive finite data, not exact zeros. Initial146-only cone directions are not true cover tangents because eight near-inactive real holes can take over linearly. Incoming22-inequality conservative model produces12 finite cost trials (one infeasible status8 preserved); three distinct signed directions include two reflection-equivalent curves.\n121 first path stages stop at actual old-hole takeover. Two physically distinct paths extend from saved pre-takeover arrays with all remaining old Voronoi branch inequalities included.48 new stages keep actual wholeD equal to repaired model up to finite Rcap .11565. No restricted-DT/circle event occurred; disk-exterior full-hull events are recorded separately.\nNo whole original objective optimizer ran, no previous optimizer or saved stage reran, no new U and no qualified event handoff. All169 actual stages plus unique baseline/source configurations independently accepted by all161700 triples and continuous all-pair boundary candidate oracle.\nReplay in a new directory: source/diagnose.py, source/incoming_diagnose.py, source/continue_cones.py, source/extend_all_contacts.py, source/audit_and_seal.py in that order. Current directory permanently frozen; do not replay into it.\nFinite failures and positive curvature are numerical evidence, not a strict SOSC certificate or global obstruction.\n')
    # Manifest excludes itself and the READY receipt, which refers to its hash.
    m={str(f.relative_to(ROOT)):sha(f) for f in sorted(ROOT.rglob('*')) if f.is_file() and f.name not in ('MANIFEST_SHA256.json','READY.json')}
    write(ROOT/'MANIFEST_SHA256.json',m)
    for name,h in m.items():assert sha(ROOT/name)==h,name
    write(ROOT/'READY.json',dict(status='PERMANENTLY_SEALED_NUMERICAL_EVIDENCE',payload_files=len(m),manifest_sha256=sha(ROOT/'MANIFEST_SHA256.json'),audit_report_sha256=sha(ROOT/'audit/INDEPENDENT_ACCEPTANCE.json'),new_upper=False,whole_optimizer_runs=0))
    print('SEALED',len(m),sha(ROOT/'MANIFEST_SHA256.json'),'unique_audits',len(seen),'stages',len(stages),'error',audit['maximum_oracle_error'],flush=True)
if __name__=='__main__':main()
