# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Common-chart derivatives, gauge and nonzero post cones. Zero cover calls."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as F
import json,itertools,numpy as np
from exact_cones import circle,positions
from review_joint import ROOT,sha,dec,enc,write,event_jet,finite_check,row_key
def main():
    I=ROOT/'inputs';B=I/'path_spec';A=I/'nonzero_alternatives';M=json.load(open(A/'MANIFEST_SHA256.json'));assert all(sha(A/p)==s for p,s in M.items());alternates={r['joint_id']:r for r in json.load(open(A/'NONZERO_PARAMETER_ALTERNATIVES.json'))['records']};spec=json.load(open(B/'PATH_SPEC.json'));Z=np.load(B/'SOURCE_EVENT_JETS.npz');joint={r['joint_id']:r for r in json.load(open(I/'JOINT_FIVE_SITE_CONSTRUCTIONS.json'))['records']};out=[];contracts=[]
    for row in spec['records']:
        jid=row['joint_id'];jr=joint[jid];p=B/row['source_coordinate_relative_file'];assert sha(p)==row['source_file_sha256'];C=np.load(p);shared=row['shared_primary_triangle'];added=row['added_sites'];source=[tuple(F(*float(x).as_integer_ratio()) for x in P) for P in C];q,r2=circle([source[i] for i in shared]);e0=[sum((source[i][j]-q[j])**2 for j in (0,1))-r2 for i in added];field,P,k=positions(source,shared,added,q,r2,{i:F(0) for i in added});jets=[]
        for j,(site,d) in enumerate(zip(added,row['source_event_definitions'])):
            e,J=event_jet(C,shared,site);ownerj=d['source_jet_index'];Hs=Z['raw_hessian'][ownerj];cols=(2*np.array(shared+[site])[:,None]+[0,1]).ravel();errors=[]
            for step in (1e-5,3e-6,1e-6):
                fd=[];fh=[]
                for col in cols:
                    L=C.ravel().copy();R=L.copy();L[col]+=step;R[col]-=step;u,v=event_jet(L.reshape(100,2),shared,site),event_jet(R.reshape(100,2),shared,site);fd.append((u[0]-v[0])/(2*step));fh.append((u[1]-v[1])/(2*step))
                errors.append(dict(step=step,gradient_max_error=float(np.max(np.abs(np.array(fd)-J[cols]))),hessian_max_error=float(np.max(np.abs(np.array(fh).T-Hs[:,cols])))))
            assert abs(e-d['driver_source_raw_scale'])<1e-14 and np.max(np.abs(J-Z['raw_gradient'][ownerj]))<1e-12 and min(r['gradient_max_error'] for r in errors)<1e-8 and min(r['hessian_max_error'] for r in errors)<1e-7
            assert dec(d['exact_source_event_scale'])==e0[j];assert abs(np.linalg.norm(J)-d['driver_source_gradient_norm'])<1e-12
            jets.append(dict(added_site=site,raw_source_value=e,exact_raw_source_value=enc(e0[j]),source_gradient_norm=float(np.linalg.norm(J)),raw_owner_gradient_max_error=float(np.max(np.abs(J-Z['raw_gradient'][ownerj]))),finite_derivative_checks=errors,rotation_derivative=float(J@np.column_stack((-C[:,1],C[:,0])).ravel()),full201_event_row_t_column=0))
        # Prove requested nonzero direction in the exact tie's local cone.
        alt=alternates[jid];eta=[F(x) for x in alt['preferred_nonzero_direction']];faces=alt['faces'];signs=[]
        for f in faces:
            for v in set(jr['union_five_sites'])-set(f):
                coeff=field.height_det_coeff(P,list(f)+[v],added);signs.append(field.sign(field.add(field.scale(coeff[0],eta[0]),field.scale(coeff[1],eta[1]))))
        assert len(signs)==6 and min(signs)>0
        adapted=dict(branch_id=alt['branch_id'],faces=faces,rational_raw_event_direction=alt['preferred_nonzero_direction'],finite_raw_event_epsilon=alt['epsilon']);finite=finite_check(source,shared,added,q,r2,e0,adapted)
        eN=dec(alt['epsilon']);goals=[float(eN*eta[j])/jets[j]['raw_source_value'] for j in (0,1)];assert np.max(np.abs(np.array(goals)-np.array(alt['driver_normalized_power_goals'])))<1e-12
        G=np.array([event_jet(C,shared,i)[1] for i in added]);rotation=np.column_stack((-C[:,1],C[:,0])).ravel();gauge=rotation/np.linalg.norm(rotation);aug=np.vstack((np.column_stack((G,np.zeros(2))),np.r_[gauge,0]));s=np.linalg.svd(aug,compute_uv=False);assert s[-1]>1e-3
        out.append(dict(joint_id=jid,source_sha256=row['source_file_sha256'],source_is_pinned14=True,cap_KKT_Hessian_not_used=True,jets=jets,augmented_two_events_plus_one_rotation_gauge_singular_values=s.tolist(),augmented_rank=3,linear_equality_null_dimension=198,local_hole_event_dependence_from_review_joint=True,nonzero_alternative_exact_raw_cone_signs=signs,nonzero_alternative_finite_local_checks=finite,wholeD_calls=0,optimizer_calls=0))
        contracts.append(dict(joint_id=jid,source_file=str(p.relative_to(ROOT)),source_file_sha256=sha(p),primary_triangle=shared,added_sites_in_control_order=added,raw_source_e0=[x['raw_source_value'] for x in jets],fixed_source_gradient_norm=[x['source_gradient_norm'] for x in jets],event_raw_formula='h_i(C)=||C_i-O_primary(C)||^2-r_primary(C)^2',event_equality_solver_row='[h_i(C)-target_raw_i]/fixed_source_gradient_norm_i; derivative [grad_h_i,0]/same fixed norm',prejoint_target_raw='(1-f)*raw_source_e0, each component; f from0 to1',postjoint_raw_eta=list(map(int,eta)),postjoint_epsilon=enc(eN),postjoint_target_raw=[enc(eN*x) for x in eta],postjoint_normalized_goals=goals,postjoint_f_values=[1-x for x in goals],desired_post_diagonals=alt['diagonals'],desired_post_faces=faces,old_diagonals_retained=0,exact_local_rational_control_direction_review_pass=True,original_cover_restart_qualified=False,required_dynamic_endpoint_gates=['Nondegenerate primary triangle','Five convex points with compatible cyclic order','Actual event q strictly inside unit disk','Other95 strict nearest at joint circle','Actual selected3 faces orientation/insideD/other97 strict empty','Complete actual restrictedDT/circlearc and actual wholeD oracle','Fresh recomputation of control-plane fan when implicit minimax frame differs from this radial existence frame'],do_not_treat_saved_constant_fan_as_global_path_certificate=True))
    write(ROOT/'JETS_NONZERO_ALTERNATIVE_REVIEW.json',dict(status='PASS',paths=out,new_exact_other97_checks=582,wholeD_calls=0,optimizer_calls=0));write(ROOT/'QUALIFIED_CONSTRAINTS.json',dict(status='READY_SOURCE_AND_LOCAL_FIVE_CONE_CONSTRAINTS',paths=contracts,proxy_launch_authorized_in_this_review=False,wholeD_calls=0,optimizer_calls=0,scope='Gate for parent-authorized finite joint minimax pilot only. No measured radius or original-problem restart claim.'))
    print(json.dumps(dict(status='PASS',joint_ids=[2,3],alternative_exact_checks=582,wholeD_calls=0,optimizer_calls=0),indent=2))
if __name__=='__main__':main()
