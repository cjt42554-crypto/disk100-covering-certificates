# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Five-cone finite event parameter gates using exact rational interval arithmetic.
No covering radius, QP, nonlinear optimizer or prior path replay is called.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[k]='1'
import sys,json,hashlib,shutil,itertools,math
from pathlib import Path
from fractions import Fraction as F
import numpy as np
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent.parent;OWNER=ROOT.parent
sys.path.insert(0,str(OWNER/'frozen_research_support_v2/source'))
from event_jets import jets,event

def enc(x):return dict(numerator=str(x.numerator),denominator=str(x.denominator))
def de(x):return F(int(x['numerator']),int(x['denominator']))
def E(z):return dict(lower=enc(z[0]),upper=enc(z[1]))
def fr(x):return F(*float(x).as_integer_ratio())
def exact(x):return (x,x)
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def neg(a):return (-a[1],-a[0])
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    v=[x*y for x in a for y in b];return min(v),max(v)
def square(a):return (F(0) if a[0]<=0<=a[1] else min(a[0]*a[0],a[1]*a[1]),max(a[0]*a[0],a[1]*a[1]))
def div(a,b):
    assert not b[0]<=0<=b[1];return mul(a,(1/b[1],1/b[0]))
def dot(a,b):return add(mul(a[0],b[0]),mul(a[1],b[1]))
def norm(a):return add(square(a[0]),square(a[1]))
def vsub(a,b):return [sub(x,y) for x,y in zip(a,b)]
def cross(a,b):return sub(mul(a[0],b[1]),mul(a[1],b[0]))
def orient(a,b,c):return cross(vsub(b,a),vsub(c,a))
def sqrt_interval(k,bits=100):
    assert k>0;M=1<<bits;n=math.isqrt((k.numerator*M*M)//k.denominator);lo=F(n,M);hi=F(n+1,M);assert lo*lo<=k<hi*hi;return (lo,hi)
def circle(P):
    a,b,c=P;u=vsub(b,a);v=vsub(c,a);D=mul(exact(F(2)),cross(u,v));assert D[0]>0;U=norm(u);V=norm(v);z=[div(sub(mul(U,v[1]),mul(V,u[1])),D),div(sub(mul(u[0],V),mul(v[0],U)),D)];q=[add(x,y) for x,y in zip(a,z)];return q,norm(z)
def construct(P,q,r2,added,h):
    X=[[exact(x) for x in p] for p in P];ks={}
    for d in added:
        v=[P[d][j]-q[j] for j in range(2)];n=sum(x*x for x in v);k=(r2+h[d])/n;assert k>0;s=sqrt_interval(k);ks[d]=dict(scale_squared=enc(k),positive_square_root_interval=E(s),sqrt_interval_squares_proved=True)
        X[d]=[add(exact(q[j]),mul(s,exact(v[j]))) for j in range(2)]
    return X,ks

def circuit_coefficients(X,tri,j,added):
    quad=tri+[j];out=[]
    for d in added:
        if d not in quad:out.append(exact(F(0)));continue
        k=quad.index(d);rest=[X[v] for n,v in enumerate(quad) if n!=k];a=orient(*rest);out.append(a if (k+3)%2==0 else neg(a))
    return out

def certify_branch(C,P,q,r2,added,triangles,eta,eps):
    X,roots=construct(P,q,r2,added,{v:eps*eta[k] for k,v in enumerate(added)});union=set(v for t in triangles for v in t);orientations=[]
    for t in itertools.combinations(list(union),3):
        a=orient(*(X[v] for v in t));orientations.append(E(a))
    cert=[]
    for f in triangles:
        o=orient(*(X[v] for v in f));assert o[0]>0;qf,R2=circle([X[v] for v in f]);inside=sub(exact(F(1)),norm(qf));assert inside[0]>0
        pp=[]
        for j in range(100):
            if j in f:continue
            power=sub(norm(vsub(X[j],qf)),R2)
            if power[0]<=0:return None
            pp.append(dict(other_site=j,power_interval=E(power)))
        cert.append(dict(face=f,orientation_interval=E(o),circumcenter_unitdisk_inside_interval=E(inside),other97_emptycircle_intervals=pp,all97_strict_positive=True))
    for j,x in enumerate(X):
        if sub(exact(F(1)),norm(x))[0]<=0:return None
    Y=np.array([[float((a+b)/2) for a,b in row] for row in X]);return X,Y,roots,cert

def integer_partial(C,faces):
    rs=[[float(x).as_integer_ratio() for x in row] for row in C];Q=max(d for row in rs for _,d in row);P=[tuple(n*(Q//d) for n,d in row) for row in rs];out=[]
    D2=lambda a,b:sum(x*y for x,y in zip(a,b))
    cr=lambda a,b:a[0]*b[1]-a[1]*b[0]
    for f in faces:
        a,b,c=[P[v] for v in f];u=(b[0]-a[0],b[1]-a[1]);v=(c[0]-a[0],c[1]-a[1]);D=2*cr(u,v);assert D>0;U=D2(u,u);V=D2(v,v);w=(U*v[1]-V*u[1],u[0]*V-v[0]*U);N=(a[0]*D+w[0],a[1]*D+w[1]);R2=D2(w,w);powers=[]
        for j,p in enumerate(P):
            if j in f:continue
            z=(p[0]*D-N[0],p[1]*D-N[1]);power=D2(z,z)-R2;assert power>0;powers.append(power)
        inside=(Q*D)**2-D2(N,N);assert inside>0;out.append(dict(face=f,orientation_integer=str(D//2),circumcenter_inside_integer=str(inside),strict_other97_integer_tests=97,minimum_emptycircle_integer=str(min(powers))))
    assert all(D2(p,p)<Q*Q for p in P);return dict(common_denominator=str(Q),all100_centers_strict_insideD=True,faces=out,complete_actual_graph_claim=False)

def main():
    assert not (ROOT/'PATH_SPEC.json').exists()
    src=OWNER/'five_site_joint_event';D=json.load(open(src/'JOINT_FIVE_SITE_CONSTRUCTIONS.json'))['records'];jobs=[D[2],D[3]];(ROOT/'inputs').mkdir(exist_ok=True);(ROOT/'finite_references').mkdir(exist_ok=True);records=[];gs=[];hs=[];source_scales=[]
    for sid,joint in enumerate(jobs):
        f=src/joint['source_coordinate_relative_file'];assert hashlib.sha256(f.read_bytes()).hexdigest()==joint['source_file_sha256'];inp=ROOT/'inputs'/f.name
        if not inp.exists():shutil.copy2(f,inp)
        assert hashlib.sha256(inp.read_bytes()).hexdigest()==joint['source_file_sha256']
        C=np.load(f);P=[tuple(fr(x) for x in c) for c in C];q=tuple(de(x) for x in joint['exact_source_common_circle_center']);r2=de(joint['exact_source_common_circle_radius_squared']);added=joint['added_sites_in_event_order'];tri=joint['common_triangle'];union=joint['union_five_sites'];ideal,idealroots=construct(P,q,r2,added,{v:F(0) for v in added})
        exact_orient=[]
        for t in itertools.combinations(joint['exact_algebraic_endpoint_ray_order'],3):
            o=orient(*(ideal[v] for v in t));assert o[0]>0;exact_orient.append(dict(face=list(t),exact_positive_orientation_interval=E(o)))
        events=[]
        for d in added:
            e,g,H,c=jets(C,'dt_flip',tri+[d]);HH=np.zeros((200,200));HH[np.ix_(c,c)]=H;ei=de(joint['exact_source_event_levels'][str(d)]);gn=float(np.linalg.norm(g));events.append(dict(ordered_common_triangle=tri,fourth_site=d,raw_definition='||c_d-O_common_triangle(C)||^2-r_common_triangle(C)^2',exact_source_event_scale=enc(ei),driver_source_raw_scale=e,driver_source_gradient_norm=gn,exact_vs_binary64_scale_absolute_difference=abs(e-float(ei)),normalized_definition='e_d(C)/driver_source_raw_scale',source_normalized_value=1.,source_jet_index=len(gs)));gs.append(g);hs.append(HH);source_scales.append(e)
        branchrows=[]
        old={tuple(sorted(d)) for d in joint['old_diagonals']}
        for bid,branch in enumerate(joint['all_five_pentagon_triangulations']):
            coeff=[]
            for face in branch['faces']:
                for d in set(union)-set(face):coeff.append(dict(face=face,other_site=d,coeff=circuit_coefficients(ideal,face,d,added)))
            eta=None;proof=None
            for maxabs in range(1,9):
                for w in itertools.product(range(-maxabs,maxabs+1),repeat=2):
                    if max(abs(w[0]),abs(w[1]))!=maxabs or not any(w):continue
                    signs=[add(mul(a['coeff'][0],exact(F(w[0]))),mul(a['coeff'][1],exact(F(w[1])))) for a in coeff]
                    if all(s[0]>0 for s in signs):eta=w;proof=signs;break
                if eta is not None:break
            assert eta is not None
            eps=F(1,1<<18);finite=None
            for shrink in range(12):
                finite=certify_branch(C,P,q,r2,added,branch['faces'],eta,eps)
                if finite is not None:break
                eps/=2
            assert finite is not None
            X,Y,roots,facecert=finite;partial=integer_partial(Y,branch['faces']);p=ROOT/'finite_references'/f'joint{joint["joint_id"]:02d}_branch{bid:02d}_UNMEASURED_RADIUS.npy';np.save(p,Y)
            levels=[eps*F(v) for v in eta];goals=[float(levels[k])/events[k]['driver_source_raw_scale'] for k in range(2)];raws=[event(Y,'dt_flip',tri+[v]) for v in added];oldkept=len(old&{tuple(sorted(d)) for d in branch['diagonals']})
            branchrows.append(dict(branch_id=bid,diagonals=branch['diagonals'],faces=branch['faces'],source_old_diagonals_retained=oldkept,two_non_crossing_new_edges=oldkept==0,rational_raw_event_direction=list(eta),exact_fixed_direction_circuit_tests=[dict(face=a['face'],other_site=a['other_site'],coefficient_on_raw_event_a=E(a['coeff'][0]),coefficient_on_raw_event_b=E(a['coeff'][1]),strict_direction_interval=E(s)) for a,s in zip(coeff,proof)],finite_raw_event_epsilon=enc(eps),finite_raw_event_levels=[enc(x) for x in levels],finite_driver_normalized_power_goals=goals,finite_driver_gradient_norm_scaled_power_goals=[float(levels[k])/events[k]['driver_source_gradient_norm'] for k in range(2)],finite_homotopy_f_values=[1-g for g in goals],algebraic_added_scale_intervals=roots,finite_local_restricted_DT_exact_interval_faces=facecert,finite_binary64_reference_relative_file=str(p.relative_to(ROOT)),finite_binary64_reference_file_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),finite_reference_raw_event_values=raws,finite_reference_raw_event_absolute_errors=[abs(raws[k]-float(levels[k])) for k in range(2)],finite_binary64_three_face_integer_check=partial,whole_domain_radius=None,qualified_cover_optimizer_restart=False,complete_actual_restricted_graph_or_circlearcs_verified=False))
            print('joint',joint['joint_id'],'branch',bid,'eta',eta,'eps',float(eps),'oldkept',oldkept,'rawerr',[abs(raws[k]-float(levels[k])) for k in range(2)],flush=True)
        records.append(dict(joint_id=joint['joint_id'],source_file_sha256=joint['source_file_sha256'],source_coordinate_relative_file='inputs/'+f.name,source_is_pinned14=True,shared_primary_triangle=tri,added_sites=added,component_candidate_ids=joint['component_candidate_ids'],candidate504_relation=('Common-triangle E72 has same zero-circle locus as original504 raw E29 but a different positive source scale; use this new shared-primary definition' if joint['joint_id']==3 else None),exact_algebraic_joint_circle_roots=idealroots,exact_all10_triangle_positive_orientations=exact_orient,source_event_definitions=events,branch_parameterizations=branchrows,
            continuation_spec=dict(prejoint_normalized_powers='(1-f,1-f), 0<=f<=1',postjoint_raw_powers='epsilon*(eta_a,eta_b); chosen fixed rational eta strictly inside one certified regular cone',postjoint_normalized_powers='epsilon*eta_i/driver_source_raw_scale_i',variable_dimension=201,spatial_dimension=200,rotation_quotient_dimension=199,epigraph='t=R^2',gauge='one spatial rotation gauge only; a common origin rotation preserves the local gates',original_cover_constraints='all actual holes/incoming; permit changing active set',required_fresh_finite_gates=['primary orientation nonzero','convex ordered five-site patch','selected3 triangle orientations positive','selected3 circumcenters strictly insideD','all97 other sites strictly outside each selected3 circle','fullactual restrictedDT/circlearc classification before any newgraph or cover restart claim'],constant_local_radial_family='c_d(h)=q+sqrt((r2+h_d)/||c_d(source)-q||^2)*(c_d(source)-q); shared3 and other95 fixed, only a geometric existence family',implicit_minimax_path_first_event_certified=False)))
    np.savez(ROOT/'SOURCE_EVENT_JETS.npz',raw_gradient=np.array(gs),raw_hessian=np.array(hs),driver_raw_source_scale=np.array(source_scales))
    out=dict(status='EXACT_FIVE_CONE_FINITE_LOCAL_GATES',records=records,joint_count=2,branch_count=10,exact_finite_face_checks=30,exact_finite_emptycircle_otherpoint_checks=2910,binary64_partial_emptycircle_integer_checks=2910,source_geometric_jets=4,source_scale_normalization_fixed=True,new_wholeD_calls=0,QP_calls=0,optimizer_calls=0,complete_actual_graph_or_cover_restart_claim=False)
    (ROOT/'PATH_SPEC.json').write_text(json.dumps(out,indent=2)+'\n');(ROOT/'README.txt').write_text('This is an exact interval / integer local event-parameter specification, not a covering radius evaluation or optimizer restart. Only joint2(80/81) andjoint3(40/504 locus) are extended. All five rational regular-cone interior directions are certified per joint using100-bit rational square-root enclosures. Every finite epsilon endpoint has its three local triangle orientations, disk-interior circumcenters, and all97 competing-center circle powers certified strictly; the rounded binary64 reference arrays also pass exact integer selected-three-face tests. Full boundary arcs, the rest of the actual restrictedDT, and the whole-domain covering radius are intentionally unknown. Use common-primary power functions and the explicit fixed source/gradient scales in PATH_SPEC. Three-hole plus two-event LICQ is not asserted. All original near-active/tiny-positive multiplier constraints belong to the future full201/rotation-quotient200 epigraph solve, not these local geometry certificates.\n')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();files={str(p.relative_to(ROOT)):sha(p) for p in sorted(ROOT.rglob('*')) if p.is_file()};(ROOT/'MANIFEST_SHA256.json').write_text(json.dumps(files,indent=2)+'\n');ready=dict(status='PERMANENTLY_SEALED_EXACT_FIVE_CONE_PATH_SPEC',payload_files=len(files),manifest_sha256=sha(ROOT/'MANIFEST_SHA256.json'),joints=2,branches=10,wholeD_calls=0,QP_calls=0,optimizer_calls=0,qualified_cover_restarts=0);(ROOT/'READY.json').write_text(json.dumps(ready,indent=2)+'\n')
    for p in ROOT.rglob('*'):
        if p.is_file():p.chmod(0o444)
    print(json.dumps(ready,indent=2))
if __name__=='__main__':main()
