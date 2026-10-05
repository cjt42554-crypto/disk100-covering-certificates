# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Exact exterior event classification and immutable finite-extension receipt."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,ast,importlib.util
from fractions import Fraction
import numpy as np
from diagnose import ROOT,sha,write

def main():
    helperdir=ROOT/'inputs/exact_helper';sys.path.insert(0,str(helperdir))
    spec=importlib.util.spec_from_file_location('sealed_exact_segment',helperdir/'exact_segment_events.py');h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    records=json.loads((ROOT/'exact_chords/RESULT.json').read_text())['records'];classified=[]
    for rec in records:
        if not rec.get('roots'):continue
        report=json.loads((ROOT/rec['output']).read_text());A=np.load(ROOT/rec['A']);B=np.load(ROOT/rec['B']);Q=int(report['common_denominator'])
        ints=[]
        for C in (A,B):
            rows=[]
            for row in C:
                rows.append([int(n)*(Q//int(d)) for n,d in (float(x).as_integer_ratio() for x in row)])
            ints.append(rows)
        P=[tuple([ints[0][j][k],ints[1][j][k]-ints[0][j][k]] for k in (0,1)) for j in range(100)]
        for root in report['all_predicate_roots']:
            labels=root['labels'];assert len(labels)==3 and root['kind'] in ('TRIANGLE_ORIENTATION_ZERO','FULL_HULL_COLLINEARITY')
            a,b,c=(P[j] for j in labels);u=h.psub(b,a);v=h.psub(c,a);D=h.mul([2],h.cross(u,v));uu=h.dot(u,u);vv=h.dot(v,v)
            pr=(h.sub(h.mul(uu,v[1]),h.mul(vv,u[1])),h.sub(h.mul(u[0],vv),h.mul(v[0],uu)))
            N=(h.add(h.mul(a[0],D),pr[0]),h.add(h.mul(a[1],D),pr[1]));side=h.primitive(h.sub(h.mul([Q*Q],h.mul(D,D)),h.dot(N,N)))
            sequence=h.sturm(side);rootcount=h.variation(sequence,Fraction(0))-h.variation(sequence,Fraction(1));negative0=h.signat(side,Fraction(0))<0;negative1=h.signat(side,Fraction(1))<0
            assert rootcount==0 and negative0 and negative1
            classified.append(dict(curve=rec['curve'],segment=rec['segment'],predicate_kind=root['kind'],labels=labels,root_isolation=root['root_isolation'],common_denominator=str(Q),associated_circumcenter_unit_side_integer_polynomial=list(map(str,side)),unit_side_definition='Q^2*D^2-|N|^2; negative means circumcenter outside closed unit disk whenever orientation is nonzero',unit_side_roots_in_closed_chord=0,unit_side_strict_negative_at_both_endpoints=True,unit_side_strict_negative_through_chord_by_Sturm=True,event_classification='exterior collinearity; no finite circumcenter in the unit disk',qualified_new_restricted_event=False,scope='This proof classifies the associated three-site locus on this straight chord; it does not certify a Newton arc, nor later-cell predicates after the first root.'))
    write(ROOT/'exact_chords/EXTERIOR_EVENT_CLASSIFICATION.json',dict(events=classified,all_events_exterior=True,distinct_chords=len({(v['curve'],v['segment']) for v in classified}),event_loci=len(classified)))
    final_integer=[]
    for curve in (0,2):
        r=json.loads((ROOT/f'continuation/stages/curve{curve:02d}/RESULT.json').read_text());last=r['stages'][-1];C=np.load(ROOT/last['coordinate_file']);verified=h.verify(C);assert verified['status']=='PASS';final_integer.append(dict(curve=curve,file=last['coordinate_file'],sha256=last['file_sha256'],integer_graph=verified,above_cap_guard=True))
    write(ROOT/'exact_chords/FINAL_GUARD_INTEGER_GRAPHS.json',final_integer)
    audit=json.loads((ROOT/'audit/INDEPENDENT_ACCEPTANCE.json').read_text());assert audit['status']=='PASS' and not audit['exact_chord_worker_errors']
    parents=json.loads((ROOT/'INPUT_HANDOFF.json').read_text());runs=[json.loads((ROOT/f'continuation/stages/curve{curve:02d}/RESULT.json').read_text()) for curve in (0,2)]
    endpoints=[];checkpoints=[];duplicate_errors=[]
    for run in runs:
        admissible=[s for s in run['stages'] if s['actual_radius']<=.1160];last=admissible[-1]
        endpoints.append(dict(curve=run['curve'],phase=last['phase'],actual_radius=last['actual_radius'],file=last['coordinate_file'],file_sha256=last['file_sha256'],restricted_graph_fingerprint=last['actual_graph']['fingerprint'],actual_minus_model_radius=last['actual_minus_model_radius'],max_center_norm=last['max_center_norm'],relevant_event=False,qualified_cover_optimizer_seed=False))
        for cap in (.11575,.11585,.1160):
            above=next(s for s in run['stages'] if s['actual_radius']>cap);below=[s for s in run['stages'] if s['actual_radius']<=cap]
            checkpoints.append(dict(curve=run['curve'],cap=cap,last_saved_at_or_below=below[-1] if below else None,first_saved_crossing_guard=above,observed_restricted_or_circle_event=False))
        for row in run['stages']:
            if any(row['old_duplicates'].values()):duplicate_errors.append(dict(file=row['coordinate_file'],duplicates=row['old_duplicates']))
    assert not duplicate_errors
    write(ROOT/'RADIUS_CAP_CHECKPOINTS.json',checkpoints)
    prior=json.loads((ROOT/'inputs/old_six_physical_event_routes.json').read_text());prior_compare=dict(prior_metadata_sha256=sha(ROOT/'inputs/old_six_physical_event_routes.json'),old_six_event_types=[r['event_identity']['event_type'] for r in prior],old_event_support_cardinalities=[len(r['event_identity']['quad_original_label_set']) for r in prior],current_event_types=sorted({r['predicate_kind'] for r in classified}),current_event_support_cardinalities=sorted({len(r['labels']) for r in classified}),matched_old_convex_four_site_incircle_events=0,reason='Current recorded exact events are exterior three-site orientation/hull collinearities; all six legacy events are four-site empty-circle diagonal exchanges. Source row labels are not equated across configurations. No new restricted event or optimizer handoff is claimed.')
    write(ROOT/'OLD_SIX_PHYSICAL_EVENT_COMPARISON.json',prior_compare)
    frozen=[]
    for name in ['b32_core_exchange_search','b32_weak_contact_continuation','b32_metric_seed_search','event_geometry_146']:
        p=ROOT.parent/name;m=json.loads((p/'MANIFEST_SHA256.json').read_text());bad=[f for f,hsh in m.items() if not (p/f).is_file() or sha(p/f)!=hsh]
        frozen.append(dict(directory=name,manifest_sha256=sha(p/'MANIFEST_SHA256.json'),payloads=len(m),passed=not bad,failed_paths=bad));assert not bad
    write(ROOT/'FROZEN_SOURCE_RECHECK.json',frozen)
    for path in (ROOT/'source').glob('*.py'):ast.parse(path.read_text())
    result=dict(status='ACCEPTED_FINITE_LOCAL_EXTENSION_NO_RELEVANT_EVENT',source_parent_manifest_sha256=parents[0]['source_parent_manifest_sha256'],source_radius=.11551473632645183,current_strict_upper_reference=.1154469548855829,independent_physical_curves=2,old169_stages_recomputed=0,old19_optimizers_repeated=0,new_cover_optimizers=0,original_structure_optimization_budget_used=19,remaining_authorized_whole_optimization_budget=77,new_saved_stage_count=audit['new_stage_count'],saved_frames_at_or_below_cap=sum(r['within_requested_radius_cap'] for r in audit['independent_alltriple_records']),guard_frames_above_cap=sum(not r['within_requested_radius_cap'] for r in audit['independent_alltriple_records']),independent_alltriple_calls=audit['new_alltriple_calls'],maximum_oracle_error=audit['maximum_oracle_error'],exact_chords=len(records),exact_chord_predicates=audit['exact_chord_predicates'],exact_no_initial_tie_chords=sum(r['status']=='PASS_EXACT_CHORD_FIRST_EVENT' for r in records),exact_chords_without_DT_circle_hull_roots=sum(r.get('roots')==0 for r in records),exterior_event_chords=len({(v['curve'],v['segment']) for v in classified}),exterior_event_loci=len(classified),observed_restricted_or_circle_events=sum(s['restricted_DT_changed'] or s['circle_inventory_changed'] for r in runs for s in r['stages']),endpoints_at_or_below_cap=endpoints,Newton_corrector_failures=sum(len(r['failures']) for r in runs),old_exact_ordered_unordered_duplicate_errors=duplicate_errors,no_new_U=True,no_qualified_optimizer_restart=True,audit_sha256=sha(ROOT/'audit/INDEPENDENT_ACCEPTANCE.json'),exact_exterior_classification_sha256=sha(ROOT/'exact_chords/EXTERIOR_EVENT_CLASSIFICATION.json'),scope='Finite continuation of two fixed mechanisms to actualRcap .116. Numerical whole-domain scores and exact first-event certificates for each straight endpoint chord are distinct. Entire implicit Newton arcs, post-first-root cells, global impossibility and optimality are not certified. All100 centers retained; new original objective optimization requires a genuinely new related event, and none was observed.')
    write(ROOT/'RESULT.json',result)
    (ROOT/'README.txt').write_text('Two original incoming-cone14 mechanisms continued in a new directory. Pre-cap parent arrays and their already accepted next guard arrays were SHA copied; those former phases were reused, not recomputed. First new phases .060/.088.\n81 new saved stages:79 at/below .116 and2 cap guards above. All81 independently accepted by161700triples+4950pair unit-circle intersections+100antipodes/axes, actualnearest100 distance; no sampled acceptance. Max numeric oracle mismatch is in RESULT.\n81 exact integer/Sturm chord first-event audits use the sealed owner helper;77 chords have no DT/circle/hull roots. Four chords have exterior three-site orientation/hull collinearities only. All six associated triangle circle-side polynomials stay strictly negative with no roots on their chords, independently proving these loci are exterior. Events are separated from six old four-site empty-circle diagonal routes.\nExact chord results do not certify the nonlinear Newton arcs or later cell evolution after a chord first event. All saved restricted graphs remain independentindex14. No qualified relevant event, no newU, no new covering optimizer, no random modes, no old169 orold19 restarts. This is finite local failure evidence.\nReplay into a fresh copy: source/continue_events.py, source/audit_new_extension.py, source/classify_and_seal.py. The preceding frozen directories and old inventory are read-only dependencies. Do not run into this permanently frozen directory.\n')
    manifest={str(f.relative_to(ROOT)):sha(f) for f in sorted(ROOT.rglob('*')) if f.is_file() and f not in (ROOT/'MANIFEST_SHA256.json',ROOT/'READY.json')}
    write(ROOT/'MANIFEST_SHA256.json',manifest)
    for f,hsh in manifest.items():assert sha(ROOT/f)==hsh
    write(ROOT/'READY.json',dict(status='PERMANENTLY_SEALED_FINITE_EXTENSION',payload_files=len(manifest),manifest_sha256=sha(ROOT/'MANIFEST_SHA256.json'),independent_audit_sha256=sha(ROOT/'audit/INDEPENDENT_ACCEPTANCE.json'),cover_optimizer_runs=0,new_upper=False))
    print('SEALED',len(manifest),sha(ROOT/'MANIFEST_SHA256.json'),'stages',audit['new_stage_count'],'error',audit['maximum_oracle_error'],flush=True)
if __name__=='__main__':main()
