# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Audit only new stages; reuse all accepted parent frames without rerunning them."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,hashlib,time,traceback,importlib.util
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
import branches,graph_tools
from independent_oracle import all_triple_cover
from whole_disk_oracle import true_cover
from diagnose import ROOT,write,sha

def chord_worker(job):
    helperdir=ROOT/'inputs/exact_helper';sys.path.insert(0,str(helperdir))
    spec=importlib.util.spec_from_file_location('owned_exact_segment',helperdir/'exact_segment_events.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    out=ROOT/job['output'];a=ROOT/job['A'];b=ROOT/job['B']
    if out.exists():
        r=json.loads(out.read_text());assert r['endpoint_A_file_sha256']==sha(a) and r['endpoint_B_file_sha256']==sha(b)
        return dict(job,reused_existing_exact_probe=True,status=r['status'],result_sha256=sha(out),first_event_kinds=[x['kind'] for x in r['first_event_cluster']],roots=r['predicate_roots_in_segment'],predicates=r['polynomial_predicates_checked'])
    try:
        r=module.audit(np.load(a),np.load(b));r.update(endpoint_A_file=str(a),endpoint_B_file=str(b),endpoint_A_file_sha256=sha(a),endpoint_B_file_sha256=sha(b))
        write(out,r)
        return dict(job,reused_existing_exact_probe=False,status=r['status'],result_sha256=sha(out),first_event_kinds=[x['kind'] for x in r['first_event_cluster']],roots=r['predicate_roots_in_segment'],predicates=r['polynomial_predicates_checked'])
    except Exception as e:
        r=dict(job,status='EXCEPTION',type=type(e).__name__,message=str(e),traceback=traceback.format_exc());write(out,r);return r

def main():
    tic=time.monotonic();parents=json.loads((ROOT/'INPUT_HANDOFF.json').read_text());runs=[];jobs=[];newstages=[]
    for curve in (0,2):
        r=json.loads((ROOT/f'continuation/stages/curve{curve:02d}/RESULT.json').read_text());runs.append(r)
        p=next(v for v in parents if v['curve']==curve);A=p['reused_stage_local_file']
        for stage in r['stages']:
            B=stage['coordinate_file'];output=f'exact_chords/curve{curve:02d}_segment{stage["stage"]:03d}.json'
            jobs.append(dict(curve=curve,segment=stage['stage'],A=A,B=B,output=output,A_sha256=sha(ROOT/A),B_sha256=sha(ROOT/B)))
            A=B;newstages.append(stage)
    write(ROOT/'exact_chords/JOB_LEDGER.json',jobs)
    exact_results=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(chord_worker,j):j for j in jobs}
        for f in as_completed(futures):
            exact_results.append(f.result());write(ROOT/'exact_chords/CHECKPOINT.json',dict(completed=len(exact_results),total=len(jobs),records=exact_results))
            if len(exact_results)%20==0:print('EXACT_CHORDS',len(exact_results),'of',len(jobs),flush=True)
    exact_results.sort(key=lambda v:(v['curve'],v['segment']));write(ROOT/'exact_chords/RESULT.json',dict(records=exact_results,worker_limit=4,BLAS_threads_per_worker=1,scope='Exact first event for each straight binary64 endpoint chord only. Any first event changes the applicable cell; later old-cell roots are not interpreted as new-cell events. Newton arcs between endpoints are not certified.'))
    alltriple=[];errors=[]
    for i,stage in enumerate(newstages):
        p=ROOT/stage['coordinate_file'];C=np.load(p);assert sha(p)==stage['file_sha256']
        R,details=all_triple_cover(C);dynamic,dyninfo=true_cover(C);graph=graph_tools.extract(C,True)
        row=dict(file=stage['coordinate_file'],file_sha256=sha(p),independent_radius=R,stored_radius=stage['actual_radius'],dynamic_radius=dynamic,error=max(abs(R-dynamic),abs(R-stage['actual_radius'])),details=details,actual_graph_fingerprint=graph['fingerprint'],within_requested_radius_cap=R<=.116,passed=abs(R-stage['actual_radius'])<2e-11 and abs(R-dynamic)<2e-11 and graph['fingerprint']==stage['actual_graph']['fingerprint'])
        alltriple.append(row)
        if not row['passed']:errors.append(row)
        if (i+1)%20==0:write(ROOT/'audit/CHECKPOINT.json',dict(completed=len(alltriple),total=len(newstages),failures=errors));print('ALLTRIPLE',len(alltriple),'of',len(newstages),flush=True)
    sourcechecks=[]
    for row in json.loads((ROOT/'SOURCE_BINDING.json').read_text()):
        actual=sha(ROOT.parent/row['source']);copy=sha(ROOT/row['copy']);sourcechecks.append(dict(row,current_source_sha256=actual,current_copy_sha256=copy,passed=actual==copy==row['sha256']))
    acceptance=dict(status='PASS' if not errors and all(v['passed'] for v in sourcechecks) else 'FAIL',new_stage_count=len(newstages),new_alltriple_calls=len(alltriple),no_prior_stage_oracle_repeated=True,maximum_oracle_error=max(v['error'] for v in alltriple),independent_alltriple_records=alltriple,source_checks=sourcechecks,errors=errors,exact_chord_worker_errors=[r for r in exact_results if r['status']=='EXCEPTION'],exact_chord_count=len(exact_results),exact_chord_predicates=sum(v.get('predicates',0) for v in exact_results),exact_chord_first_event_kinds={k:sum(k in v.get('first_event_kinds',[]) for v in exact_results) for k in sorted({k for v in exact_results for k in v.get('first_event_kinds',[])})},elapsed_seconds=time.monotonic()-tic,numerical_cover_acceptance_only=True,exact_chords_not_a_Newton_arc_certificate=True)
    write(ROOT/'audit/INDEPENDENT_ACCEPTANCE.json',acceptance)
    print('DONE',len(newstages),'error',acceptance['maximum_oracle_error'],'exact_errors',len(acceptance['exact_chord_worker_errors']),'kinds',acceptance['exact_chord_first_event_kinds'],flush=True)
if __name__=='__main__':main()
