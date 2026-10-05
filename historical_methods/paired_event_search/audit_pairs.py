# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Independent numerical replay of the twelve saved coupled event paths."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,time,itertools
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import paired_paths as p
f=p.solo;ROOT=p.ROOT


def main():
    tic=time.monotonic();checks=json.loads((ROOT/'input_checks.json').read_text());C0=np.load(ROOT/'inputs/baseline_corrected.npy');reports=[];candidate_hashes=[];failures=[];scored=0;scoreerror=0.;accepted=0;attemptfail=0
    copychecks={name:p.sha(ROOT/name)==digest for name,digest in checks['input_sha256'].items()}
    frozenroot=ROOT.parent/'event_forced_paths';manifest=json.loads((frozenroot/'MANIFEST_SHA256.json').read_text());frozenchecks={name:p.sha(frozenroot/name)==digest for name,digest in manifest.items()}
    for resultpath in sorted(ROOT.glob('pair*/result.json')):
        result=json.loads(resultpath.read_text());folder=resultpath.parent;events=result['pair']['events'];rows=result['steps'];nodes=[]
        for row in rows:
            path=folder/row['snapshot'];C=np.load(path)['C'];R,detail=f.true_cover(C);scored+=1;err=abs(R-row['true_radius']);scoreerror=max(scoreerror,err)
            values=np.array([f.incircle(C,e)/e['base_value_corrected'] for e in events]);res=float(np.max(np.abs(values-(1-row['fraction']))));predicates=[f.exact_predicates(C,e) for e in events]
            union=result['pair']['union_sites'];ds=np.linalg.norm(C[union,None,:]-C[None,union,:],axis=2);ds[np.diag_indices(len(union))]=np.inf;at=np.unravel_index(int(np.argmin(ds)),ds.shape)
            node=dict(snapshot=row['snapshot'],phase=row['phase'],fraction=row['fraction'],independent_true_radius=R,stored_score_error=err,event_relative_residual=res,exact_binary64_predicates=predicates,minimum_cavity_pair_distance=float(ds[at]),closest_cavity_sites=[union[at[0]],union[at[1]]],coordinate_sha256=f.center_sha(C),file_sha256=p.sha(path),actual_center_norm=float(np.linalg.norm(C,axis=1).max()),maximum_baseline_shift=float(np.linalg.norm(C-C0,axis=1).max()))
            nodes.append(node)
            if err>2e-12 or node['file_sha256']!=row['file_sha256'] or node['coordinate_sha256']!=row['coordinate_sha256']:failures.append(dict(path=str(path),reason='score_or_hash_mismatch'))
            if [x['incircle_sign'] for x in predicates]!=[x['incircle_sign'] for x in row['exact_binary64_predicates']]:failures.append(dict(path=str(path),reason='saved_exact_binary64_sign_mismatch'))
            if row['phase']=='failed_attempt':attemptfail+=1
            else:
                accepted+=1
                if R>.125+1e-10 or res>3e-9 or node['actual_center_norm']>1+1e-9 or node['maximum_baseline_shift']>.3+1e-9:failures.append(dict(path=str(path),reason='accepted_geometry_or_event_budget_failure'))
        numeric_roots=[row for row in rows if row['phase']=='double_event'];validation=None
        if numeric_roots:
            Ce=np.load(folder/numeric_roots[-1]['snapshot'])['C'];validation=p.pair_validation(Ce,events)
            if validation['valid']!=result['event_validation']['valid']:failures.append(dict(path=str(folder),reason='event_classification_mismatch'))
        before=[row for row in rows if row['phase']=='corrected'];after=[row for row in rows if row['phase']=='after_both_events'];brackets=None
        if before and after:
            brackets=[[before[-1]['exact_binary64_predicates'][i]['incircle_sign'],after[-1]['exact_binary64_predicates'][i]['incircle_sign']] for i in range(2)]
            if brackets!=[[1,-1],[1,-1]]:failures.append(dict(path=str(folder),reason='crossing_sign_bracket_failure',brackets=brackets))
        restores=[]
        for restored in result['restorations']:
            for stage in restored['rounds']:
                path=folder/stage['coordinate_file'];C=np.load(path);R,detail=f.true_cover(C);scored+=1;err=abs(R-stage['true_radius']);scoreerror=max(scoreerror,err)
                revised=p.known_basin(C)
                candidate_hashes.append(dict(path=str(path.relative_to(ROOT)),file_sha256=p.sha(path),coordinate_sha256=f.center_sha(C),true_radius=R,stage=stage['stage'],known_single_basin_legacy=stage['known_single_basin'],known_binding_skeleton_replayed=revised,known_basin_legacy_false_negative=bool(revised['match'] and not stage['known_single_basin']['match'])))
                if err>2e-12:failures.append(dict(path=str(path),reason='restoration_stage_score_mismatch'))
            path=folder/restored['coordinate_file'];C=np.load(path);R,detail=f.true_cover(C);scored+=1;err=abs(R-restored['true_radius']);scoreerror=max(scoreerror,err)
            entry=dict(path=str(path.relative_to(ROOT)),file_sha256=p.sha(path),coordinate_sha256=f.center_sha(C),independent_true_radius=R,fully_unconstrained_stage_valid=restored['fully_unconstrained_stage_valid'],final_event_constraints_present=restored.get('final_event_constraints_present'),known_single_basin_legacy=restored.get('known_single_basin'),known_binding_skeleton_replayed=p.known_basin(C),exact_binary64_event_signs=[f.exact_predicates(C,e)['incircle_sign'] for e in events]);restores.append(entry);candidate_hashes.append(entry)
            if err>2e-12:failures.append(dict(path=str(path),reason='final_restoration_score_mismatch'))
        strictcurves=[row['stationary_diagnostic']['tangent_H_eigenvalues'] for row in rows if row.get('stationary_diagnostic',{}).get('tangent_H_eigenvalues')]
        knownstops=[stage for rr in result['restorations'] for stage in rr['rounds'] if stage.get('remaining_optimization_skipped_known_single_basin')]
        report=dict(pair_index=result['pair']['pair_index'],ids=result['pair']['ids'],family=result['pair']['family'],normalized_cross_coupling=result['pair']['normalized_cross_coupling'],gram_eigenvalues=result['pair']['gram_eigenvalues'],predicted_radius_barrier=result['pair']['predicted_radius_barrier'],accepted_nodes=result['accepted_nodes'],failed_attempt_nodes=result['failed_attempt_nodes'],last_fraction=result['last_accepted_fraction'],stop=result['stop'],validated_numeric_event_pair=bool(validation and validation['valid']),validated_five_site_cavity=bool(validation and validation.get('five_site_cavity_cocircular_valid')),crossed_both_events=bool(validation and validation['valid'] and after),exact_binary64_sign_brackets=brackets,event_validation=validation,sampled_peak_true_radius=result['sampled_peak_true_radius'],restorations=restores,minimum_saved_equality_tangent_curvature=min([min(ev) for ev in strictcurves],default=None),known_basin_optimization_stops=len(knownstops),raw_root_minimum_cavity_site_distance=next((n['minimum_cavity_pair_distance'] for n in nodes if n['phase']=='double_event'),None),nodes=nodes,scope='Numerical root acceptance and sampled cover; no interval/path/global-optimum certificate')
        missed=[]
        for rr in result['restorations']:
            first=None;after_count=0
            for stage in rr['rounds']:
                revised=p.known_basin(np.load(folder/stage['coordinate_file']))
                if revised['match'] and not stage['known_single_basin']['match']:
                    if first is None:first=stage['round']
                    elif stage['round']>first:after_count+=1
            if first is not None:missed.append(dict(first_matching_round=first,extra_solver_rounds_after_first_match=after_count,unbound_center=65,reason='Legacy alignment fitted all 100 centers; unused center 65 drifts while the other 99 binding centers match the known skeleton'))
        report['known_basin_legacy_detection_misses']=missed
        reports.append(report);p.write(folder/'independent_replay.json',report)
    if not all(copychecks.values()) or not all(frozenchecks.values()):failures.append(dict(reason='copied_or_frozen_input_hash_mismatch'))
    audit=dict(status='PASS' if not failures else 'FAIL',completed_pairs=len(reports),accepted_nodes=accepted,failed_attempt_nodes=attemptfail,independently_scored_arrays=scored,maximum_score_difference=scoreerror,validated_numeric_event_pairs=sum(r['validated_numeric_event_pair'] for r in reports),validated_five_site_cavities=sum(r['validated_five_site_cavity'] for r in reports),crossed_both_event_pairs=sum(r['crossed_both_events'] for r in reports),invalid_primitive_zero_pairs=sum(not r['validated_numeric_event_pair'] for r in reports),known_binding_skeleton_detection_missed_pairs=sum(bool(r['known_basin_legacy_detection_misses']) for r in reports),extra_solver_rounds_after_known_skeleton=sum(x['extra_solver_rounds_after_first_match'] for r in reports for x in r['known_basin_legacy_detection_misses']),known_detection_fix_scope='Future detector fits only the 99 binding-center skeleton and reports free center 65. Original logs are preserved; no optimization was repeated for this fix.',copied_input_hashes=copychecks,frozen_event_forced_tree_unchanged=all(frozenchecks.values()),frozen_tree_file_count=len(frozenchecks),frozen_manifest_sha256=p.sha(frozenroot/'MANIFEST_SHA256.json'),optimization_paths_repeated=0,source_sha256=p.sha(ROOT/'paired_paths.py'),helper_sha256=p.sha(frozenroot/'forced_paths.py'),paths=reports,failures=failures,seconds=time.monotonic()-tic,scope='Complete saved numerical replay; exact rational signs only certify saved binary64 predicates. Simultaneous raw primitive zero is explicitly insufficient for actual five-site events; no global optimum claim')
    p.write(ROOT/'audit_report.json',audit);p.write(ROOT/'candidate_hashes.json',candidate_hashes)
    print(json.dumps({k:audit[k] for k in ['status','completed_pairs','accepted_nodes','failed_attempt_nodes','independently_scored_arrays','maximum_score_difference','validated_numeric_event_pairs','validated_five_site_cavities','crossed_both_event_pairs','invalid_primitive_zero_pairs','frozen_event_forced_tree_unchanged','seconds']}),flush=True)


if __name__=='__main__':main()
