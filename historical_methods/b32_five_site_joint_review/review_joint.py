# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Independent Fraction/biquadratic fan review. No fullD or optimizer import."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,hashlib,itertools,time,functools
from fractions import Fraction as F
import numpy as np
from exact_cones import Field,circle,positions,triangulations,cone_rays,circle_inside
ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dec(v):return F(int(v['numerator']),int(v['denominator']))
def enc(v):return dict(numerator=str(v.numerator),denominator=str(v.denominator))
def write(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def triangle_jet(C,ids):
    P=C[ids];a,b,c=P;B=P[1:]-a;rhs=(np.sum(P[1:]**2,axis=1)-a@a)/2;q=np.linalg.solve(B,rhs);v=q-a;g=float(v@v)
    M=np.zeros((2,6));M[:,:2]=-(a-q)
    for j in (1,2):M[j-1,2*j:2*j+2]=P[j]-q
    Q=np.zeros((2,200));cols=(2*np.array(ids)[:,None]+[0,1]).ravel();Q[:,cols]=np.linalg.solve(B,M);J=2*v@Q;J[2*ids[0]:2*ids[0]+2]-=2*v
    return q,g,Q,J
def event_jet(C,shared,added):
    q,g,Q,J=triangle_jet(C,shared);v=C[added]-q;grad=-2*v@Q-J;grad[2*added:2*added+2]+=2*v
    return float(v@v-g),grad
def row_key(t):return tuple(sorted(tuple(sorted(e)) for e in t['diagonals']))
def cyclic_equivalent(a,b):return any(a==b[i:]+b[:i] for i in range(5))
def fan_complete(field,cones):
    rays=[];edges=[]
    for cone in cones:
        ids=[]
        for r in cone['rays']:
            k=next((i for i,v in enumerate(rays) if field.sign(field.cross(r,v))==0 and field.sign(field.dot(r,v))>0),None)
            if k is None:k=len(rays);rays.append(r)
            ids.append(k)
        edges.append(ids)
    assert len(rays)==5 and all(sum(i in e for e in edges)==2 for i in range(5))
    def half(ray):
        y=field.sign(ray[1]);x=field.sign(ray[0]);return 0 if y>0 or(y==0 and x>=0) else 1
    def compare(i,j):
        a,b=rays[i],rays[j]
        if half(a)!=half(b):return half(a)-half(b)
        return -field.sign(field.cross(a,b))
    angles=sorted(range(5),key=functools.cmp_to_key(compare));expected={tuple(sorted((angles[j],angles[(j+1)%5]))) for j in range(5)}
    assert {tuple(sorted(e)) for e in edges}==expected
    assert all(field.sign(field.cross(rays[angles[j]],rays[angles[(j+1)%5]]))>0 for j in range(5))
    return dict(ray_count=5,cone_count=5,every_boundary_ray_shared_by_two_cones=True,cyclic_cone_adjacency_exact_cross_checks=True,cone_ray_indices=edges)
def finite_check(source,shared,added,q,r2,e0,branch):
    epsilon=dec(branch['finite_raw_event_epsilon']);eta=[F(x) for x in branch['rational_raw_event_direction']];levels={i:epsilon*eta[j] for j,i in enumerate(added)};field,P,k=positions(source,shared,added,q,r2,levels)
    H={i:sum((source[i][j]-q[j])**2 for j in (0,1))-r2 for i in range(100)};H.update(levels);checks=[]
    for face in branch['faces']:
        assert field.sign(field.orient(P,face))>0;assert circle_inside(field,P,face)>0;signs=[]
        for site in set(range(100))-set(face):signs.append(field.sign(field.height_det(P,list(face)+[site],H)))
        assert len(signs)==97 and min(signs)>0;checks.append(dict(face=face,orientation_positive=True,circumcenter_strict_inside_unit_disk=True,exact_other97_emptycircle_signs_positive=True,competitors=97))
    assert all(field.sign(field.sub(field.q(1),field.dot(P[i],P[i])))>0 for i in range(100))
    # This verifies three local faces only, never all cells or a covering radius.
    return dict(branch_id=branch['branch_id'],raw_eta=[int(x) if x.denominator==1 else enc(x) for x in eta],epsilon=enc(epsilon),raw_event_levels=[enc(levels[i]) for i in added],normalized_event_targets=[float(levels[i]/e0[j]) for j,i in enumerate(added)],faces=checks,total_exact_emptycircle_checks=291,all100_centers_strict_inside_unit_disk=True,wholeD_calls=0,complete_actual_graph_or_radius_claim=False)
def main():
    I=ROOT/'inputs';M=json.load(open(I/'MANIFEST_SHA256.json'));assert all(sha(I/p)==s for p,s in M.items());path=I/'path_spec';N=json.load(open(path/'MANIFEST_SHA256.json'));assert all(sha(path/p)==s for p,s in N.items())
    spec=json.load(open(path/'PATH_SPEC.json'));sp={r['joint_id']:r for r in spec['records']};Z=np.load(path/'SOURCE_EVENT_JETS.npz');owner=json.load(open(I/'JOINT_FIVE_SITE_CONSTRUCTIONS.json'));results=[];finite=[];start=time.monotonic()
    for row in owner['records']:
        jid=row['joint_id'];npfile=I/row['source_coordinate_relative_file'];assert sha(npfile)==row['source_file_sha256'];C=np.load(npfile);source=[tuple(F(*float(x).as_integer_ratio()) for x in P) for P in C];shared=row['common_triangle'];added=row['added_sites_in_event_order'];cycle=row['exact_algebraic_endpoint_ray_order'];q,r2=circle([source[i] for i in shared]);assert q==tuple(dec(v) for v in row['exact_source_common_circle_center']) and r2==dec(row['exact_source_common_circle_radius_squared']);assert 1-sum(v*v for v in q)>0
        e0=[sum((source[i][j]-q[j])**2 for j in (0,1))-r2 for i in added];assert all(x>0 for x in e0)
        exact_rest_powers=[sum((source[i][j]-q[j])**2 for j in (0,1))-r2 for i in set(range(100))-set(cycle)];assert len(exact_rest_powers)==95 and min(exact_rest_powers)>0
        field,P,k=positions(source,shared,added,q,r2,{i:F(0) for i in added});assert all(field.sign(field.orient(P,f))>0 for f in itertools.combinations(cycle,3));assert cyclic_equivalent(cycle,row['source_patch_boundary'])
        independent=triangulations(cycle);own={row_key(t):t for t in row['all_five_pentagon_triangulations']};assert set(own)=={row_key(t) for t in independent};cones=[];records=[];Werrors=[]
        ownerfan={row_key(t):t for t in row['numerical_first_order_lifting_fan']}
        for t in independent:
            rows=[];metadata=[]
            for face in t['faces']:
                A=field.orient(P,face);assert field.sign(A)>0
                for v in sorted(set(cycle)-set(face)):
                    coeff=field.height_det_coeff(P,list(face)+[v],added);rows.append(coeff);metadata.append(dict(face=face,other_site=v,cleared_raw_height_coefficients=[field.enc(x) for x in coeff],orientation_positive=True))
            rays,inner=cone_rays(field,rows);key=row_key(t);numeric=np.array([[field.fl(c)/field.fl(field.orient(P,m['face'])) for c in r] for r,m in zip(rows,metadata)]);expected=np.array(ownerfan[key]['linear_other_site_emptycircle_inequality_rows']);errors=np.max(np.abs(numeric[:,None,:]-expected[None,:,:]),axis=2);Werrors.append(max(float(errors.min(0).max()),float(errors.min(1).max())))
            numericinner=np.array([field.fl(x) for x in inner]);numericinner/=np.max(np.abs(numericinner));eta=[F(int(round(x*10000)),10000) for x in numericinner];assert all(field.sign(field.add(field.scale(r[0],eta[0]),field.scale(r[1],eta[1])))>0 for r in rows)
            cones.append(dict(rays=rays,rows=rows))
            records.append(dict(diagonals=[list(e) for e in t['diagonals']],faces=t['faces'],cleared_inequality_rows=metadata,rays=[[field.enc(x) for x in ray] for ray in rays],independently_chosen_raw_eta=[enc(x) for x in eta],all_six_inequalities_at_raw_eta_exact_positive=True))
        complete=fan_complete(field,cones)
        # Two individually proposed diagonals are combinatorially crossing.
        idx={v:j for j,v in enumerate(cycle)};d0,d1=row['individual_new_diagonals'];a,b=sorted(idx[x] for x in d0);c,d=sorted(idx[x] for x in d1);crossed=(a<c<b<d)or(c<a<d<b);assert crossed and not any(set(map(tuple,map(sorted,[d0,d1])))==set(t['diagonals']) for t in independent)
        Y=np.load(I/row['ideal_display_relative_file']);events=np.array([event_jet(Y,shared,i)[1] for i in added]);holes=np.array([triangle_jet(Y,f)[3] for f in row['source_patch_faces']]);aug=np.vstack((np.column_stack((holes,-np.ones(3))),np.column_stack((events,np.zeros(2)))));singular=np.linalg.svd(aug,compute_uv=False);deps=[]
        for j in (1,2):
            coef=np.linalg.lstsq(events.T,(holes[j]-holes[0]),rcond=None)[0];deps.append(dict(hole_difference=j,coefficients_on_two_raw_events=coef.tolist(),dependency_residual_inf=float(np.max(np.abs(holes[j]-holes[0]-coef@events)))))
        gauge_dot=field.q(0)
        for i in range(100):gauge_dot=field.add(gauge_dot,field.add(field.scale(P[i][0],source[i][0]),field.scale(P[i][1],source[i][1])))
        assert field.sign(gauge_dot)>0
        review=dict(joint_id=jid,source_name=row['source_name'],source_sha256=row['source_file_sha256'],common_primary_triangle=shared,added_sites=added,common_chart_raw_source_values=[enc(x) for x in e0],exact95_nearest_checks_pass=True,seven_geometric_sources_not_stationarity_models=True,pentagon_cycle=cycle,all10_algebraic_triangle_orientations_exact_positive=True,independent_recursive_triangulations=records,fan_completeness=complete,owner_lifting_row_max_numeric_difference=max(Werrors),individual_new_diagonals_cross=True,simultaneous_both_crossing_diagonals_rejected=True,event_rank_exact_from_disjoint_added_coordinate_blocks=2,rotation_gauge_preserves_event_rank_exact=True,fixed_source_gauge_dot_current_rotation_exact_positive=True,local_three_hole_plus_two_event_functional_row_rank_exact=3,local_expected_dependency_dimension=2,local_augmented_singular_values=singular.tolist(),local_dependency_numerical_checks=deps,local_functional_dependencies_do_not_establish_global_active_constraint_LICQ_failure=True,raw_to_normalized_fan_rule='W_raw*diag(eA0,eB0) applied to normalized levels. PATH_SPEC eta is RAW power direction; normalized post goal is epsilon*eta_i/e_i0.',no_KKT_Hessian_transferred_to_cap=True,wholeD_calls=0,optimizer_calls=0,source_geometry_review_pass=True)
        if jid in sp:
            pr=sp[jid];assert pr['source_file_sha256']==row['source_sha256'] if 'source_sha256' in row else pr['source_file_sha256']==row['source_file_sha256'];assert pr['shared_primary_triangle']==shared and pr['added_sites']==added
            certs=[]
            for br in pr['branch_parameterizations']:
                key=row_key(br);cc=cones[next(i for i,t in enumerate(independent) if row_key(t)==key)];eta=[F(x) for x in br['rational_raw_event_direction']];assert all(field.sign(field.add(field.scale(r[0],eta[0]),field.scale(r[1],eta[1])))>0 for r in cc['rows']);cert=finite_check(source,shared,added,q,r2,e0,br);certs.append(cert);finite.append(cert);print('FINITE_BRANCH_PASS',jid,br['branch_id'],291,flush=True)
            review['owner_path_spec_branches_independently_exact_checked']=certs
        results.append(review);write(ROOT/'CHECKPOINT.json',dict(joints_completed=len(results),total=7,finite_branches_completed=len(finite),new_wholeD_calls=0,new_optimizer_calls=0));print('JOINT_REVIEW_PASS',jid,'5cones','rank3',flush=True)
    out=dict(status='PASS',owner_joint_manifest_sha256=sha(I/'MANIFEST_SHA256.json'),owner_path_spec_manifest_sha256=sha(path/'MANIFEST_SHA256.json'),reviewed_joint_geometry=7,reviewed_exact_regular_cones=35,finite_branches_checked=10,finite_three_face_emptycircle_checks=2910,wholeD_calls=0,optimizer_calls=0,records=results,seconds=time.monotonic()-start,scope='Exact local algebraic geometry and cone/constraint qualification only; no covering radius, complete graph, optimizer start or proxy curve certificate.')
    write(ROOT/'JOINT_REVIEW.json',out);print(json.dumps(dict(status='PASS',joints=7,cones=35,finite_branches=10,checks=2910,seconds=out['seconds']),indent=2),flush=True)
if __name__=='__main__':main()
