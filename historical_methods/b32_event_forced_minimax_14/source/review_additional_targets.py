# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Limited independent additional-pool jets/qualification review, no optimizer."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.spatial import ConvexHull,Delaunay
import event_jets,graph_tools
from whole_disk_oracle import true_cover
from proxy_driver import write,sha,csha
ROOT=Path(__file__).resolve().parent.parent
def physical(C,ids):
    e=event_jets.event_jet(C,'incircle',ids);q=e['center'];r2=e['radius_squared'];hull=ConvexHull(C[ids]);P=C[np.array(ids)[hull.vertices]];inside=bool(np.all(np.cross(np.roll(P,-1,axis=0)-P,q-P)>=-1e-12));others=[i for i in range(100) if i not in ids]
    return dict(q=q.tolist(),q_norm_squared=float(q@q),q_inside_unit_disk=bool(q@q<1-1e-9),quad_convex=len(hull.vertices)==4,q_inside_quad_diagnostic=inside,minimum_other96_circle_power=float(np.min(np.sum((C[others]-q)**2,axis=1)-r2)),raw_event_value=float(e['raw_value']))
def main():
    base=ROOT/'additional_review';D=base/'inputs';h=json.load(open(D/'HANDOFF.json'));m=json.load(open(D/'MANIFEST_SHA256.json'));assert all(sha(D/p)==digest for p,digest in m.items());Z=np.load(D/'TARGET_JETS.npz');out=[]
    for r in h['records']:
        ids=list(map(int,r['ordered_triangle']))+[int(r['fourth_site'])];p=D/r['source_coordinate_relative_file'];assert sha(p)==r['source_coordinate_file_sha256'];C=np.load(p);e=event_jets.event_jet(C,'incircle',ids);j=r['jet_gradient_index'];gradient=Z['raw_gradient'][j];H=Z['raw_hessian'][j];e0=r['raw_source_scale'];cols=np.array([(2*i+k) for i in ids for k in (0,1)]);errors=[]
        for step in (1e-5,3e-6,1e-6):
            fd=[];fH=[]
            for col in cols:
                A=C.ravel().copy();B=A.copy();A[col]+=step;B[col]-=step;u=event_jets.event_jet(A.reshape(100,2),'incircle',ids);v=event_jets.event_jet(B.reshape(100,2),'incircle',ids);fd.append((u['raw_value']-v['raw_value'])/(2*step));fH.append((u['raw_jacobian']-v['raw_jacobian'])/(2*step))
            fH=np.array(fH).T;errors.append(dict(h=step,gradient_absolute_max_error=float(np.max(np.abs(np.array(fd)-e['raw_jacobian'][cols]))),hessian_absolute_max_error=float(np.max(np.abs(fH-H[:,cols])))))
        T=C.copy();q=e['center'];v=T[ids[3]]-q;T[ids[3]]=q+v*np.sqrt(e['radius_squared']/(v@v));tie=physical(T,ids);source=physical(C,ids);dep=event_jets.differential_dependency(T,'incircle',ids);V=Z['ranking_jet_direction'][j];dtzero=float(e['raw_jacobian']@V/e0);rotation=np.column_stack((-C[:,1],C[:,0])).ravel();R,details=true_cover(C)
        passed=abs(e['raw_value']-e0)<1e-12 and np.max(np.abs(gradient-e['raw_jacobian']))<1e-12 and min(x['gradient_absolute_max_error'] for x in errors)<1e-8 and min(x['hessian_absolute_max_error'] for x in errors)<1e-7 and abs(dtzero+1)<1e-7 and tie['q_inside_unit_disk'] and tie['quad_convex'] and tie['minimum_other96_circle_power']>1e-8 and abs(tie['raw_event_value'])<1e-12
        out.append(dict(target_id=r['target_id'],candidate_id=r['candidate_id'],labels=ids,source_file=r['source_coordinate_relative_file'],source_file_sha256=sha(p),source_coordinate_sha256=csha(C),source_is_stationary_reference=r['source_is_stationary_reference'],source_radius_actual=R,source_radius_owner=r['source_radius'],source_actual_graph=graph_tools.extract(C,True),raw_clearance=e['raw_value'],raw_source_gradient_norm=e['gradient_norm'],owner_raw_gradient_max_error=float(np.max(np.abs(gradient-e['raw_jacobian']))),owner_normalized_gradient_max_error=float(np.max(np.abs(Z['normalized_gradient'][j]-e['raw_jacobian']/e0))),rotation_derivative=float(e['raw_jacobian']@rotation),ranking_direction_normalized_event_derivative=dtzero,derivative_checks=errors,source_geometry=source,ideal_floating_tie_geometry=tie,ideal_floating_tie_not_cover_restart=True,event_dependency_at_ideal_tie=dep,reference_KKT_Hessian_transferred_to_cap=False,review_passed=bool(passed),eligible_as_source_event_target_only=bool(passed),additional_proxy_launch=False))
    report=dict(status='PASS' if all(x['review_passed'] for x in out) else 'FAIL',owner_handoff_manifest_sha256=sha(D/'MANIFEST_SHA256.json'),owner_targets=4,reviewed_targets=out,whole_optimizer_calls=0,finite_diagnostic_steps=(1e-5,3e-6,1e-6),actual_solved_event_endpoints=0,new_strict_certificate=False,scope='Sign/scaling/FD/source geometry and ideal floating tie feasibility only. No cap stationarity, no event path, no cover restart qualified.')
    write(base/'ADDITIONAL_REVIEW.json',report);print(json.dumps(dict(status=report['status'],targets=[dict(id=x['target_id'],passed=x['review_passed'],minFD=min(y['gradient_absolute_max_error'] for y in x['derivative_checks']),minH=min(y['hessian_absolute_max_error'] for y in x['derivative_checks']),sourceR=x['source_radius_actual']) for x in out]),indent=2))
if __name__=='__main__':main()
