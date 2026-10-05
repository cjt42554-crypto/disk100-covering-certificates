# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Independent numerical selected-frame acceptance and exact-byte handoff.

No optimizer. Exact dyadic fan predicates below concern a lifted 5-point
configuration only; they are not an exact five-circle or whole-cover proof.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal
import json,hashlib,shutil,itertools
import numpy as np
from scipy.spatial import ConvexHull,Delaunay
from independent_oracle import all_triple_cover
from whole_disk_oracle import true_cover
import graph_tools

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'first_joint_handoff'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def csha(C):return hashlib.sha256(np.asarray(C,dtype='<f8').tobytes()).hexdigest()
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def det(M):
    if len(M)==1:return M[0][0]
    return sum(((-1)**j)*M[0][j]*det([r[:j]+r[j+1:] for r in M[1:]]) for j in range(len(M)))
def circle(C,ids):
    P=C[ids];u=P[1]-P[0];v=P[2]-P[0]
    z=np.cross(u,v);uu=float(u@u);vv=float(v@v)
    q=P[0]+np.array([uu*v[1]-vv*u[1],u[0]*vv-v[0]*uu])/(2*z)
    return q,float(np.sum((q-P[0])**2)),float(z)
def exact_fan(C,t):
    """Cofactors of det[x,y,height,1]; -det/orient = lower-plane slack."""
    added=t['added_sites_in_control_order'];ids=t['primary_triangle']+added
    rows=[]
    for face in t['desired_post_faces']:
        for v in sorted(set(ids)-set(face)):
            quad=face+[v];M=[[F.from_float(float(C[i,0])),F.from_float(float(C[i,1])),F(0),F(1)] for i in quad]
            A=det([[r[0],r[1],r[3]] for r in M[:3]])
            assert A>0
            coeff=[]
            for a in added:
                index=quad.index(a) if a in quad else None
                cofactor=F(0) if index is None else (-1)**(index+2)*det([[x for j,x in enumerate(r) if j!=2] for k,r in enumerate(M) if k!=index])
                coeff.append(-cofactor/A)
            eta=t['postjoint_raw_eta'];margin=sum(x*y for x,y in zip(coeff,eta))
            rows.append(dict(face=face,other_site=v,coefficients=[str(x) for x in coeff],margin_for_source_raw_eta=str(margin),numeric_margin=float(margin)))
    candidates=[]
    for a,b in itertools.product(range(-12,0),repeat=2):
        if np.gcd(abs(a),abs(b))!=1:continue
        values=[sum(F(x)*v for x,v in zip(row['coefficients'],[a,b])) for row in rows]
        if min(values)>0:candidates.append(dict(raw_eta=[a,b],minimum_exact_margin=str(min(values)),normalized_margin=float(min(values))/np.hypot(a,b)))
    candidates.sort(key=lambda r:-r['normalized_margin'])
    return dict(method='Exact Fraction cofactor lower-lift inequalities at the saved dyadic positions; not a circle certificate',rows=rows,source_raw_eta=t['postjoint_raw_eta'],source_raw_eta_inside_actual_preferred_cone=all(F(r['margin_for_source_raw_eta'])>0 for r in rows),small_negative_rational_alternatives=candidates[:5],alternatives_not_optimized_or_launch_authorized=True)
def gate(C,t,kind):
    S=t['primary_triangle'];added=t['added_sites_in_control_order'];ids=S+added
    q,rho,A=circle(C,S);h=[float(np.sum((C[i]-q)**2)-rho) for i in added]
    hull=ConvexHull(C[ids]);convex=len(hull.vertices)==5
    remaining=[i for i in range(100) if i not in ids]
    nearest=float(np.min(np.sum((C[remaining]-q)**2,axis=1)-rho))
    actual_faces={tuple(sorted(map(int,f))) for f in Delaunay(C).simplices}
    checks=[]
    for face in t['desired_post_faces']:
        qq,rr,orient=circle(C,face);others=[i for i in range(100) if i not in face]
        power=float(np.min(np.sum((C[others]-qq)**2,axis=1)-rr))
        present=tuple(sorted(face)) in actual_faces
        checks.append(dict(face=face,orientation=orient,q_squared_norm=float(qq@qq),other97_minimum_power=power,actual_DT_face_present=present,strict_gate=bool(orient>1e-10 and qq@qq<1-1e-8 and power>1e-10 and present)))
    geometry=bool(A>1e-10 and convex and q@q<1-1e-8 and nearest>1e-8 and max(map(abs,h))<1e-8)
    fan=exact_fan(C,t)
    return dict(kind=kind,primary_orientation=A,primary_q=q.tolist(),primary_radius_squared=rho,primary_q_squared_norm=float(q@q),raw_powers=h,convex5=convex,cyclic_labels=[ids[i] for i in hull.vertices],other95_minimum_power=nearest,near_joint_geometry_numerical_gate=geometry,preferred_source_eta_fan_gate=fan['source_raw_eta_inside_actual_preferred_cone'],preferred_post_faces=checks,preferred_post_three_strict_faces_gate=all(r['strict_gate'] for r in checks),actual_fan=fan,exact_joint_or_cover_certificate=False)
def main():
    assert not OUT.exists(),'Never rewrite an existing handoff'
    OUT.mkdir();(OUT/'coordinates').mkdir();(OUT/'source').mkdir()
    for name in ['selected_audit.py','independent_oracle.py','whole_disk_oracle.py','graph_tools.py']:
        shutil.copyfile(ROOT/'source'/name,OUT/'source'/name)
    specs=[]
    for joint in ['joint02','joint03']:
        result=json.loads((ROOT/'workers'/joint/'RESULT.json').read_text())
        for node in result['accepted']:
            if (node['phase']=='prejoint' and node['goal']>=.9999) or (node['phase']=='postjoint' and node['goal']==2**-18):
                specs.append((joint,result,node))
    rows=[]
    for joint,result,node in specs:
        src=ROOT/node['file'];assert sha(src)==node['file_sha256'];C=np.load(src);assert C.shape==(100,2) and np.isfinite(C).all()
        label=joint+('_near1' if node['phase']=='prejoint' and node['goal']==1 else '_pre' if node['phase']=='prejoint' else '_post')
        dest=OUT/'coordinates'/(label+'.npy');shutil.copyfile(src,dest);assert sha(dest)==sha(src)
        csv=OUT/'coordinates'/(label+'.csv')
        texts=[[format(Decimal.from_float(float(x)),'f') for x in p] for p in C]
        assert all(F(s)==F.from_float(float(x)) for p,text in zip(C,texts) for x,s in zip(p,text))
        csv.write_text('\n'.join(','.join(p) for p in texts)+'\n')
        R,detail=all_triple_cover(C);Rdt,dt=true_cover(C);passed=abs(R-Rdt)<2e-12 and abs(R-node['actual_radius'])<2e-12
        geom=gate(C,result['target'],label);graph=graph_tools.extract(C,True)
        # Deliberately a loose decimal radius supplied for root's exact test.
        # An independently numerical radius is not promoted to a strict upper.
        proposed=format(Decimal.from_float(float(max(R,Rdt)))+Decimal('0.000000000002'),'f')
        row=dict(name=label,primary=result['target']['primary_triangle'],added=result['target']['added_sites_in_control_order'],source_initial_file_sha256=result['target']['source_file_sha256'],original_file=node['file'],original_file_sha256=sha(src),coordinate_file=str(dest.relative_to(OUT)),coordinate_file_sha256=sha(dest),coordinate_sha256=csha(C),csv_file=str(csv.relative_to(OUT)),csv_file_sha256=sha(csv),csv_every_coordinate_exactly_equals_saved_binary64_rational=True,phase=node['phase'],phase_control=node['goal'],proxy_status=node['status'],raw_control=node['raw_goal'],wholeD_radius_alltriple=R,wholeD_radius_DT_allpair=Rdt,wholeD_agreement_pass=passed,independent_alltriple_oracle=detail,independent_boundary_continuous_oracle=dt,geometry_gate=geom,actual_graph=graph,proposed_rational_radius_for_root_exact_test=proposed,exact_coverage_certificate=False,original_restart_qualified=False)
        rows.append(row);print(json.dumps(dict(name=label,R=R,pass_numerical=passed,joint_geometry=geom['near_joint_geometry_numerical_gate'],preferred_eta=geom['preferred_source_eta_fan_gate'],preferred_post=geom['preferred_post_three_strict_faces_gate'])),flush=True)
    assert len(rows)==5 and all(r['wholeD_agreement_pass'] for r in rows)
    report=dict(status='PASS_NUMERICAL_SELECTED_FRAMES',rows=rows,unique_arrays=5,alltriple_calls=5,proxy_solves=0,original_restores=0,no_exact_joint_or_cover_certificate=True,scope='Independent selected endpoints, all161700 triples+all4950 pairs+antipodes, and an exact dyadic control fan diagnosis only')
    write(OUT/'INDEPENDENT_ENDPOINT_AUDIT.json',report)
    handoff=dict(status='READY_IMMUTABLE_SELECTED_ENDPOINTS',rows=[{k:r[k] for k in ['name','primary','added','coordinate_file','coordinate_file_sha256','coordinate_sha256','csv_file','csv_file_sha256','original_file','original_file_sha256','source_initial_file_sha256','phase','phase_control','wholeD_radius_alltriple','proposed_rational_radius_for_root_exact_test']} for r in rows],joint02='Preferred post reached with actual strict three-face numerical gate. Near1 is numerical tie only.',joint03='Near1 actual five-circle geometry gate passes. Original source raw eta fails its preferred cone in actual geometry; no post solve attempted.',proxy_total_original_phase=31,original_restores=0,exact_projection_delegated_to_root=True,no_whole_path_certificate=True)
    write(OUT/'HANDOFF.json',handoff)
    manifest={str(p.relative_to(OUT)):sha(p) for p in sorted(OUT.rglob('*')) if p.is_file()}
    write(OUT/'MANIFEST_SHA256.json',manifest);assert all(sha(OUT/p)==s for p,s in manifest.items())
    write(OUT/'READY.json',dict(status='READY_PERMANENT_IMMUTABLE',payload_files=len(manifest),manifest_sha256=sha(OUT/'MANIFEST_SHA256.json'),all_payload_sha_pass=True,numerical_only=True,no_original_restore_launched=True))
    print(json.dumps(dict(ready=str(OUT/'READY.json'),manifest_sha256=sha(OUT/'MANIFEST_SHA256.json'),files=len(manifest))),flush=True)
if __name__=='__main__':main()
