# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Replay saved physical-circle states with the complete all100 oracle."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import physical_paths as p
ROOT=p.ROOT;f=p.f

def main():
    tic=time.monotonic();m=np.load(ROOT/'inputs/matrices.npz');C0=m['C'];Q=m['Q'];w=m['w'];t0=float(m['t']);checks=json.loads((ROOT/'input_checks.json').read_text());reports=[];failures=[];scored=0;maxscore=0.;candidate_hashes=[]
    inputchecks={name:p.sha(ROOT/name)==digest for name,digest in checks['input_sha256'].items()};frozen={}
    for dirname in ('event_forced_paths','paired_event_search'):
        tree=ROOT.parent/dirname;manifest=json.loads((tree/'MANIFEST_SHA256.json').read_text());bad=[name for name,digest in manifest.items() if p.sha(tree/name)!=digest];frozen[dirname]=dict(unchanged=not bad,manifest_sha256=p.sha(tree/'MANIFEST_SHA256.json'),files=len(manifest),mismatches=bad)
    for resultpath in sorted(ROOT.glob('path*/result.json')):
        result=json.loads(resultpath.read_text());case=result['case'];folder=resultpath.parent;nodes=[];event=None
        for row in result['steps']:
            path=folder/row['snapshot'];a=np.load(path);C=a['C'];x=a['x'];q=a['q'];rho=float(a['rho_squared']);R,detail=f.true_cover(C);scored+=1;err=abs(R-row['true_radius']);maxscore=max(maxscore,err);prob=p.PhysicalProblem(C0,Q,w,t0,case,row['s'],C0,trust=.4);prob.core.top=p.geometry.topology(C);g,A,kk=f.all_holes(C,prob.core.top);physical=float(prob.con(x)[len(g)+300:].min());eq=float(np.max(np.abs(prob.eq(x))));cycle=[float(p.area(C,t)/a0) for t,a0 in zip(case['cycle_triangles'],case['cycle_initial_areas'])]
            node=dict(snapshot=row['snapshot'],phase=row['phase'],s=row['s'],true_radius=R,stored_score_error=err,physical_eq_scaled_residual=eq,physical_inequality_minimum=physical,common_area_ratio=float(p.area(C,case['common'])/case['common_initial_area']),minimum_cycle_area_ratio=min(cycle),q_norm=float(np.linalg.norm(q)),rho_squared=rho,physical_J_rank=int(np.linalg.matrix_rank(prob.eqjac(x),tol=1e-9)),exact_binary64_predicates=[f.exact_predicates(C,e) for e in case['events']],coordinate_sha256=f.center_sha(C),file_sha256=p.sha(path),maximum_center_norm=float(np.linalg.norm(C,axis=1).max()),maximum_baseline_center_shift=float(np.linalg.norm(C-C0,axis=1).max()))
            nodes.append(node)
            if err>2e-12 or node['coordinate_sha256']!=row['coordinate_sha256'] or node['file_sha256']!=row['file_sha256']:failures.append(dict(path=str(path),reason='score_or_snapshot_hash_mismatch'))
            if [z['incircle_sign'] for z in node['exact_binary64_predicates']]!=[z['incircle_sign'] for z in row['exact_binary64_predicates']]:failures.append(dict(path=str(path),reason='exact_saved_sign_mismatch'))
            if row['phase']!='failed_attempt' and (eq>3e-9 or physical<-2e-8 or R>.125+1e-10 or node['q_norm']>.999+1e-9 or node['maximum_center_norm']>1+1e-9 or node['maximum_baseline_center_shift']>.3+1e-9 or node['physical_J_rank']!=5):failures.append(dict(path=str(path),reason='accepted_physical_state_budget_or_rank_failure'))
            if row['phase']=='five_site_circle':
                event=p.validate(C,q,rho,case)
                if not event['valid']:failures.append(dict(path=str(path),reason='physical_five_site_event_invalid',validation=event))
        before=[r for r in nodes if r['phase']=='corrected'];after=[r for r in nodes if r['phase']=='after_circle'];brackets=None
        if before and after:
            brackets=[[before[-1]['exact_binary64_predicates'][i]['incircle_sign'],after[-1]['exact_binary64_predicates'][i]['incircle_sign']] for i in range(2)]
            if brackets!=[[1,-1],[1,-1]]:failures.append(dict(path=str(folder),reason='physical_event_sign_bracket_failure',brackets=brackets))
        restored=[]
        for rr in result['restorations']:
            for stage in rr['rounds']:
                path=folder/stage['coordinate_file'];C=np.load(path);R,detail=f.true_cover(C);scored+=1;err=abs(R-stage['true_radius']);maxscore=max(maxscore,err)
                candidate_hashes.append(dict(path=str(path.relative_to(ROOT)),stage=stage['stage'],true_radius=R,coordinate_sha256=f.center_sha(C),file_sha256=p.sha(path),known_family=stage['known_family']))
                if err>2e-12:failures.append(dict(path=str(path),reason='restoration_stage_score_mismatch'))
            path=folder/rr['coordinate_file'];C=np.load(path);R,detail=f.true_cover(C);scored+=1;err=abs(R-rr['true_radius']);maxscore=max(maxscore,err);entry=dict(path=str(path.relative_to(ROOT)),true_radius=R,coordinate_sha256=f.center_sha(C),file_sha256=p.sha(path),fully_unconstrained_valid=rr['fully_unconstrained_valid'],physical_constraints=False if not rr.get('skipped') else None,offsets=0,exact_binary64_event_signs=[f.exact_predicates(C,e)['incircle_sign'] for e in case['events']],actual_all100_oracle=True);restored.append(entry);candidate_hashes.append(entry)
            if err>2e-12:failures.append(dict(path=str(path),reason='restored_candidate_score_mismatch'))
        report=dict(path_index=case['path_index'],ids=case['ids'],union=case['union'],initial_hull_vertices=case['cycle'],initial_gaps=case['gaps'],local_prediction_radius_barrier=case['predicted_radius_barrier'],accepted_nodes=result['accepted_nodes'],failed_attempt_nodes=result['failed_attempt_nodes'],last_s=result['last_s'],valid_physical_five_site_circle=bool(event and event['valid']),both_events_crossed=bool(event and event['valid'] and after),event_validation=event,exact_binary64_sign_brackets=brackets,sampled_peak_true_radius=result['sampled_peak_true_radius'],minimum_accepted_common_area_ratio=min(n['common_area_ratio'] for n in nodes if n['phase']!='failed_attempt'),minimum_accepted_cycle_area_ratio=min(n['minimum_cycle_area_ratio'] for n in nodes if n['phase']!='failed_attempt'),stop=result['stop'],restorations=restored,nodes=nodes,scope='Numerical physical geometry and saved covering nodes; no exact event/path/global-optimum certification')
        reports.append(report);p.write(folder/'independent_replay.json',report)
    if not all(inputchecks.values()) or not all(t['unchanged'] for t in frozen.values()):failures.append(dict(reason='frozen_or_copied_input_hash_mismatch'))
    audit=dict(status='PASS' if not failures else 'FAIL',completed_paths=len(reports),accepted_nodes=sum(r['accepted_nodes'] for r in reports),failed_attempt_nodes=sum(r['failed_attempt_nodes'] for r in reports),valid_physical_five_site_circles=sum(r['valid_physical_five_site_circle'] for r in reports),crossed_physical_circles=sum(r['both_events_crossed'] for r in reports),independently_scored_arrays=scored,maximum_score_difference=maxscore,copied_input_hashes=inputchecks,frozen_trees=frozen,optimization_paths_repeated=0,actual_coverage_center_count=100,source_sha256=p.sha(ROOT/'physical_paths.py'),paths=reports,failures=failures,seconds=time.monotonic()-tic,scope='Exact rational signs of saved doubles and complete numerical all100 covering replay only. Auxiliary circle constraints are not the coverage oracle. No certification improvement/global optimality claim')
    p.write(ROOT/'audit_report.json',audit);p.write(ROOT/'candidate_hashes.json',candidate_hashes);print(json.dumps({k:audit[k] for k in ['status','completed_paths','accepted_nodes','failed_attempt_nodes','valid_physical_five_site_circles','crossed_physical_circles','independently_scored_arrays','maximum_score_difference','seconds']}),flush=True)

if __name__=='__main__':main()
