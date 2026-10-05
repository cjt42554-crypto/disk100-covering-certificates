# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Replay saved numerical geometry without repeating any optimization path."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,hashlib,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import forced_paths as f
ROOT=f.ROOT


def main():
    tic=time.monotonic();checks=json.loads((ROOT/'input_checks.json').read_text());geom=np.load(ROOT/'full48_geometry.npz');C0=geom['C'];Q=geom['Q'];w=geom['w'];t0=float(geom['t0'])
    frozen=json.loads((ROOT.parent/'event_release_diagnostic/frozen_linear_manifest.json').read_text())
    source_status={name:f.sha(ROOT.parent/'event_release_diagnostic'/name)==digest for name,digest in frozen.items()}
    copy_status={name:f.sha(ROOT/name)==digest for name,digest in checks['input_files'].items()}
    events=json.loads((ROOT/'inputs/events.json').read_text());reports=[];candidate_files=[];failures=[];scored=0;maxscore=0.
    for folder in sorted(ROOT.glob('path*_*')):
        if not(folder/'result.json').exists():continue
        result=json.loads((folder/'result.json').read_text());event=result['event'];supp=[];margins=[]
        for row in result['steps']:
            snap=folder/row['snapshot'];arr=np.load(snap);C=arr['C'];R,detail=f.true_cover(C);scored+=1;error=abs(R-row['true_radius']);maxscore=max(maxscore,error)
            ep=f.exact_predicates(C,event);em=float(f.incircle(C,event));relative=em/event['base_value_corrected'];inventorym,inventory=f.event_inventory(C,events);margins.append(inventorym)
            top=f.geometry.topology(C);g,A,kk=f.all_holes(C,top);ac=np.flatnonzero(g.max()-g<8e-10);grad=f.event_jet(C,event)[1]/event['base_value_corrected'];Beq=np.r_[np.c_[A[ac]@Q,-np.ones(len(ac))],np.r_[grad@Q,0.][None,:]];sv=np.linalg.svd(Beq,compute_uv=False);rank=int(np.linalg.matrix_rank(Beq,tol=1e-9))
            cr=row.get('stationary_diagnostic',{});zerostress=cr.get('zero_stress_keys',[])
            entry=dict(snapshot=row['snapshot'],phase=row['phase'],fraction=row['event_fraction'],coordinate_sha256=f.center_sha(C),file_sha256=f.sha(snap),independent_true_radius=R,stored_score_error=error,event_relative_residual=abs(relative-(1-row['event_fraction'])),exact_binary64_predicates=ep,complete_event_inventory=inventory,epigraph_event_equality_rows=len(Beq),epigraph_event_equality_rank=rank,epigraph_event_rank_deficiency=len(Beq)-rank,epigraph_event_smallest_singular_value=float(sv[-1]),zero_stress_keys=zerostress,curvature_scope='Stored positive curvature is on the active equality subspace; zero multiplier inequalities additionally require critical-cone analysis before a local-minimum certificate',zero_stress_requires_critical_cone_analysis=bool(zerostress))
            supp.append(entry)
            if error>2e-12:failures.append(dict(path=folder.name,snapshot=row['snapshot'],reason='score_mismatch',error=error))
            if ep['incircle_sign']!=row['exact_predicates']['incircle_sign']:failures.append(dict(path=folder.name,snapshot=row['snapshot'],reason='exact_sign_mismatch'))
            if row['phase']!='failed_attempt':
                if entry['event_relative_residual']>2e-9 or R>.13+1e-12 or np.linalg.norm(C-C0,axis=1).max()>.3+1e-10 or np.linalg.norm(C,axis=1).max()>1+1e-9:failures.append(dict(path=folder.name,snapshot=row['snapshot'],reason='accepted_state_budget_or_event_failure'))
        f.write_json(folder/'supplementary_checks.json',supp)
        np.savez_compressed(folder/'all_baseline_event_margins_replayed.npz',event_ids=np.array([e['id'] for e in events]),margins=np.array(margins),snapshots=np.array([r['snapshot'] for r in result['steps']]))
        candidate_status=[]
        for restored in result['restorations']:
            for stage in restored['rounds']:
                path=folder/stage['coordinate_file'];C=np.load(path);R,detail=f.true_cover(C);scored+=1;err=abs(R-stage['independent_true_radius']);maxscore=max(maxscore,err)
                candidate_files.append(dict(path=str(path.relative_to(ROOT)),file_sha256=f.sha(path),coordinate_sha256=f.center_sha(C),independent_true_radius=R,event_sign=f.exact_predicates(C,event)['incircle_sign'],stage=stage['stage']))
                if err>2e-12:failures.append(dict(path=str(path),reason='restoration_stage_score_mismatch',error=err))
            path=folder/restored['coordinate_file'];C=np.load(path);R,detail=f.true_cover(C);scored+=1;err=abs(R-restored['true_radius']);maxscore=max(maxscore,err)
            status=dict(path=str(path.relative_to(ROOT)),file_sha256=f.sha(path),coordinate_sha256=f.center_sha(C),independent_true_radius=R,event_sign=f.exact_predicates(C,event)['incircle_sign'],final_event_constraint_present=restored['final_event_constraint_present'],fully_unconstrained_stage_valid=restored['fully_unconstrained_stage_valid'],stage='final_all_actual_holes_zero_offsets')
            candidate_files.append(status);candidate_status.append(status)
            if err>2e-12 or status['final_event_constraint_present'] or not status['fully_unconstrained_stage_valid']:failures.append(dict(path=str(path),reason='final_restore_validation_failure',error=err))
        er=result.get('event_validation');states=[r for r in result['steps'] if r['phase']=='event'];evcheck=None
        if states:
            Ce=np.load(folder/states[-1]['snapshot'])['C'];evcheck=f.validate_event(Ce,event)
            if not evcheck['valid']:failures.append(dict(path=folder.name,reason='event_geometry_validation_failure'))
        accepted=[r for r in result['steps'] if r['phase']!='failed_attempt'];small=[]
        for row in accepted:
            x=row['event_fraction']
            if 0<x<=.003:
                observed=row['true_radius']**2-t0;pred=result['initial_prediction']['predicted_squared_barrier']*x*x
                small.append(dict(fraction=x,observed_squared_radius_cost=observed,quadratic_model_squared_radius_cost=pred,ratio=observed/pred if pred else None,scope='Numerical small-step consistency, no Taylor remainder bound'))
        reports.append(dict(path=folder.name,event=event['id'],accepted_nodes=len(accepted),failed_attempt_nodes=len(result['steps'])-len(accepted),last_accepted_fraction=accepted[-1]['event_fraction'],stop=result['stop'],event_valid=evcheck and evcheck['valid'],event_inside_norm=evcheck and evcheck.get('common_circle_center_norm'),event_binary64_exact_zero=evcheck and evcheck['exact_predicates']['incircle_exactly_zero'],exact_binary64_sign_bracket=er and er.get('exact_binary64_sign_bracket'),sampled_peak_true_radius=max(r['true_radius'] for r in accepted),restorations=candidate_status,small_step_model_checks=small,min_positive_equality_curvature=min([min(r['stationary_diagnostic']['tangent_H_eigenvalues']) for r in result['steps'] if r.get('stationary_diagnostic',{}).get('tangent_H_eigenvalues')],default=None),zero_stress_nodes=sum(bool(r.get('stationary_diagnostic',{}).get('zero_stress_keys')) for r in result['steps']),scope='Saved-state numerical validation only; sampled peak not an interval supremum certificate'))
    if not all(source_status.values()) or not all(copy_status.values()):failures.append(dict(reason='frozen_or_copied_input_hash_mismatch'))
    audit=dict(status='PASS' if not failures else 'FAIL',source_frozen_hashes=source_status,copied_input_hashes=copy_status,optimization_paths_repeated=0,completed_paths=len(reports),independent_scored_arrays=scored,maximum_independent_score_difference=maxscore,validated_numerical_events=sum(bool(r['event_valid']) for r in reports),paths=reports,candidate_files=candidate_files,failures=failures,seconds=time.monotonic()-tic,one_CPU_BLAS_thread=True,source_driver_sha256=f.sha(ROOT/'forced_paths.py'),guarded_core_sha256=f.sha(ROOT.parent/'pro_core_audit/guarded_core.py'),scope='Complete saved geometry replay, exact rational signs of saved doubles, numerical inside/nearest/circle checks; no exact event or full-segment or global optimality claim',provenance_note='First three paths were computed before adding full baseline boundary-event inventory and epigraph-rank fields. This audit adds those diagnostics on their preserved snapshots without repeating optimization. The finite optimization formulas are unchanged.')
    f.write_json(ROOT/'audit_report.json',audit);f.write_json(ROOT/'candidate_hashes.json',candidate_files)
    manifest={str(p.relative_to(ROOT)):f.sha(p) for p in sorted(ROOT.rglob('*')) if p.is_file() and p.name!='MANIFEST_SHA256.json'}
    f.write_json(ROOT/'MANIFEST_SHA256.json',manifest)
    print(json.dumps({k:audit[k] for k in ['status','completed_paths','independent_scored_arrays','maximum_independent_score_difference','validated_numerical_events','optimization_paths_repeated','seconds']}),flush=True)


if __name__=='__main__':main()
