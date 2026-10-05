# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Pure exact-circle / pentagon combinatorics. No wholeD, QP or optimizer."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[k]='1'
import sys,json,hashlib,itertools,functools,shutil
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import numpy as np
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent.parent;BASE=ROOT.parent

def fr(x):return F(*float(x).as_integer_ratio())
def enc(x):return dict(numerator=str(x.numerator),denominator=str(x.denominator))
def vec(v):return [enc(x) for x in v]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def sign(x):return (x>0)-(x<0)
def circle(P):
    a,b,c=P;u=sub(b,a);v=sub(c,a);D=2*cross(u,v);assert D
    U=dot(u,u);V=dot(v,v);w=((U*v[1]-V*u[1])/D,(u[0]*V-v[0]*U)/D);q=tuple(x+y for x,y in zip(a,w));return q,dot(w,w)
def radical_sign(x,y,k):
    assert k>0
    if not x:return sign(y)
    if not y:return sign(x)
    if sign(x)==sign(y):return sign(x)
    return sign(x)*sign(x*x-y*y*k)
def angle_order(vs):
    def half(v):return 0 if v[1]>0 or(v[1]==0 and v[0]>=0) else 1
    def cmp(a,b):
        if half(vs[a])!=half(vs[b]):return half(vs[a])-half(vs[b])
        z=cross(vs[a],vs[b]);return -sign(z) if z else 0
    return sorted(vs,key=functools.cmp_to_key(cmp))
def cyclic_same(a,b):return any(a==b[i:]+b[:i] for i in range(len(b)))
def patch_boundary(faces):
    directed={}
    for f in faces:
        for j in range(3):
            k=(f[j],f[(j+1)%3]);assert k not in directed;directed[k]=f[(j+2)%3]
    boundary=[k for k in directed if k[::-1] not in directed];nxt=dict(boundary);assert len(boundary)==5 and len(nxt)==5
    start=min(nxt);v=start;out=[]
    while v not in out:out.append(v);v=nxt[v]
    assert v==start and len(out)==5;return out

def pentagon_triangulations(cycle):
    n=5;bedges={tuple(sorted((cycle[j],cycle[(j+1)%n]))) for j in range(n)};idx={v:i for i,v in enumerate(cycle)}
    diags=[tuple(sorted(e)) for e in itertools.combinations(cycle,2) if tuple(sorted(e)) not in bedges]
    def crossing(e,f):
        if set(e)&set(f):return False
        a,b=sorted(idx[x] for x in e);c,d=sorted(idx[x] for x in f);return (a<c<b<d)or(c<a<d<b)
    out=[]
    for pair in itertools.combinations(diags,2):
        if crossing(*pair):continue
        edges=bedges|set(pair);faces=[list(t) for t in itertools.combinations(cycle,3) if all(tuple(sorted(e)) in edges for e in itertools.combinations(t,2))];assert len(faces)==3
        out.append(dict(diagonals=[list(d) for d in pair],faces=faces))
    assert len(out)==5;return out,crossing

def fan(tris,C,union,added):
    heights=np.zeros((5,2));uids={v:i for i,v in enumerate(union)}
    for j,v in enumerate(added):heights[uids[v],j]=1
    result=[]
    for t in tris:
        W=[]
        for face in t['faces']:
            ii=[uids[x] for x in face];M=np.vstack((C[face].T,np.ones(3)))
            for v in set(union)-set(face):
                bary=np.linalg.solve(M,np.r_[C[v],1]);W.append(heights[uids[v]]-bary@heights[ii])
        W=np.array(W);angles=np.linspace(0,2*np.pi,1440,endpoint=False);v=np.column_stack((np.cos(angles),np.sin(angles)));margin=(v@W.T).min(1);ok=np.flatnonzero(margin>1e-8);best=int(np.argmax(margin));result.append(dict(diagonals=t['diagonals'],linear_other_site_emptycircle_inequality_rows=W.tolist(),sampled_strict_sector_angles=[float(angles[k]) for k in ok],best_sample_angle=float(angles[best]),best_sample_minimum_margin=float(margin[best]),status='NUMERICAL_FIRST_ORDER_REGULAR_LIFTING_FAN_ONLY'))
    return result

def main():
    data=json.load(open(BASE/'frozen_research_support_v2/candidate_targets.json'))['records'];jobs=[(0,80,k) for k in (23,79,81)]+[(0,40,k) for k in (3,41,42)]+[(1,505,504)]
    (ROOT/'inputs').mkdir();(ROOT/'ideal_display').mkdir();records=[]
    for jid,(sid,a,b) in enumerate(jobs):
        ra,rb=data[a],data[b];assert ra['source_id']==rb['source_id']==sid
        C=np.load(ra['source_file']);sourcep=ROOT/'inputs'/Path(ra['source_file']).name
        if not sourcep.exists():shutil.copy2(ra['source_file'],sourcep)
        P=[tuple(fr(x) for x in c) for c in C];S=set(ra['labels'])&set(rb['labels']);U=set(ra['labels'])|set(rb['labels']);assert len(S)==3 and len(U)==5
        shared=next(f for f in ra['adjacent_faces'] if set(f)==S);assert any(set(f)==S for f in rb['adjacent_faces']);neighbors=[next(f for f in r['adjacent_faces'] if set(f)!=S) for r in (ra,rb)];faces=[shared]+neighbors;boundary=patch_boundary(faces);q,r2=circle([P[i] for i in shared]);assert r2>0 and 1-dot(q,q)>0
        added=[next(iter(set(r['labels'])-S)) for r in (ra,rb)];directions={i:sub(P[i],q) for i in U};order=angle_order(directions);same=cyclic_same(order,boundary)
        distinct=all(cross(directions[i],directions[j])!=0 or dot(directions[i],directions[j])<0 for i,j in itertools.combinations(U,2));assert distinct
        rest=sorted(set(range(100))-U);powers={i:dot(sub(P[i],q),sub(P[i],q))-r2 for i in rest};assert min(powers.values())>0
        e0={i:dot(directions[i],directions[i])-r2 for i in added};assert all(v>0 for v in e0.values());Y=C.copy();endpoint=[]
        for i in added:
            k=r2/dot(directions[i],directions[i]);A=1-dot(q,q)-r2;B=-2*dot(q,directions[i]);sg=radical_sign(A,B,k);Y[i]=np.array([float(x) for x in q])+np.sqrt(float(k))*np.array([float(x) for x in directions[i]]);endpoint.append(dict(site=i,scale_squared=enc(k),positive_scale_branch=True,exact_unitdisk_inside_sign=sg,numerical_center_norm=float(np.linalg.norm(Y[i]))))
        display=ROOT/'ideal_display'/f'joint{jid:02d}_FLOATING_DISPLAY_NOT_SEED.npy';np.save(display,Y)
        qcircle_inside=1-dot(q,q)-r2;whole_circle_inside=qcircle_inside>0 and qcircle_inside*qcircle_inside>4*dot(q,q)*r2
        tris,crossing=pentagon_triangulations(order);newdiags=[ra['opposites'],rb['opposites']];old=[ra['edge'],rb['edge']];crossnew=crossing(*newdiags)
        if same:assert any({tuple(sorted(d)) for d in t['diagonals']}=={tuple(sorted(d)) for d in old} for t in tris)
        qinsideall=all(cross(directions[order[i]],directions[order[(i+1)%5]])>0 for i in range(5))
        quad_grad_coeffs=[]
        # Exact determinant coefficients at the algebraic display are numerical;
        # the structural formula itself is exact in a shared triangle chart.
        for quad in itertools.combinations(order,4):
            M=np.column_stack((np.ones(4),Y[list(quad)]));cs=[]
            for d in added:
                if d not in quad:cs.append(0.);continue
                k=list(quad).index(d);cs.append(float(((-1)**(k+3))*np.linalg.det(np.delete(M,k,axis=0))))
            quad_grad_coeffs.append(dict(quad=list(quad),determinant_last_height_column_coefficients_on_event_pair=cs))
        row=dict(joint_id=jid,source_name=ra['source_name'],source_file_sha256=ra['source_file_sha256'],source_coordinate_relative_file=str(sourcep.relative_to(ROOT)),component_candidate_ids=[a,b],component_event_quads=[ra['labels'],rb['labels']],common_triangle=shared,union_five_sites=sorted(U),added_sites_in_event_order=added,source_patch_faces=faces,source_patch_boundary=boundary,exact_algebraic_endpoint_ray_order=order,source_patch_boundary_matches_exact_final_ray_order=same,
            exact_source_common_circle_center=vec(q),exact_source_common_circle_radius_squared=enc(r2),exact_common_circle_insideD_margin=enc(1-dot(q,q)),exact_remaining95_nearestcircle_checks=95,remaining95_power_signs_all_positive=True,exact_min_remaining95_circle_power=enc(min(powers.values())),exact_source_event_levels={str(i):enc(v) for i,v in e0.items()},
            endpoint_algebraic_centers=endpoint,all_five_algebraic_endpoint_centers_exact_insideD=all(e['exact_unitdisk_inside_sign']>0 for e in endpoint),whole_common_circle_exact_insideD=whole_circle_inside,all_five_circle_points_exact_distinct=True,five_point_convexity_exact_from_distinct_common_circle=True,common_circle_center_inside_five_convex_hull=qinsideall,quad_hull_containment_used_as_gate=False,
            ideal_display_relative_file=str(display.relative_to(ROOT)),floating_display_not_exact_binary64_joint_tie=True,qualified_cover_restart=False,
            joint_geometric_event_construction_qualified=bool(same and all(e['exact_unitdisk_inside_sign']>0 for e in endpoint)),old_diagonals=old,individual_new_diagonals=newdiags,individual_new_diagonals_cross_on_convex_five_cycle=crossnew,all_five_pentagon_triangulations=tris,numerical_first_order_lifting_fan=fan(tris,Y,order,added),four_subset_determinant_jets=quad_grad_coeffs,
            independent_event_rank_exact=2,three_epigraph_holes_plus_two_events_rank_exact=3,three_holes_plus_two_events_gradient_dependency_dimension_exact=2,scope='Exact algebraic empty nearest five-circle construction and combinatorics. No current binary64 event tie, full implicit path first-event, whole-domain radius, optimizer or new graph claim.')
        records.append(row);print(jid,ra['source_name'],a,b,'sites',sorted(U),'qualified',row['joint_geometric_event_construction_qualified'],'sourcecyclematch',same,'inside',[e['exact_unitdisk_inside_sign'] for e in endpoint],'newdiags_cross',crossnew,flush=True)
    out=dict(status='COMPLETE_PURE_GEOMETRY',records=records,construction_count=len(records),qualified_constructions=sum(r['joint_geometric_event_construction_qualified'] for r in records),new_wholeD_calls=0,QP_calls=0,optimizer_calls=0,selected_event_pilot_repeated=False)
    (ROOT/'JOINT_FIVE_SITE_CONSTRUCTIONS.json').write_text(json.dumps(out,indent=2)+'\n')
    (ROOT/'DERIVATION.txt').write_text('''For shared noncollinear source triangle S and added sites a,b, let q,r2 be the exact rational circumcircle of S from binary64 source coordinates. Put v_i=c_i-q, e_i0=||v_i||2-r2>0. The algebraic radial family c_i(f)=q+sqrt((r2+(1-f)e_i0)/||v_i||2)*v_i satisfies e_i(C(f))/e_i0=1-f exactly. The three sites S and every other95site remain fixed. At f=1 the five sites are exactly cocircular and all other95 powers at q remain strictly positive. q is strictly insideD. For f in[0,1), q remains nearest to the three shared sites; this does not certify absence of other actual graph or covering events elsewhere along the family. InsideD of the two algebraic endpoints is checked by the exact sign of A+B*sqrt(k); the full radial segments then stay inside the convex disk. The displayed npys are rounded visual references and have no exact tie or radius claim.

Five sites on a positive-radius circle with distinct rays form a strictly convex pentagon. Their cyclic order is established by exact rational ray comparisons. Source patch boundary is checked against that order. The five possible pentagon triangulations are enumerated combinatorially. In particular the two individually flipped replacement diagonals often cross, so two single-flip graph targets cannot simply be combined. Joint event constraints require fresh actual DT post-side qualification.

For the two common-triangle circlepower events, the a-coordinate columns of grad e_a are2(c_a-q), while grad e_b has zeros there, and conversely at b. Both are nonzero, proving event row rank2 even after quotienting the rotation gauge. Three adjacent epigraph-hole rows minus the two event rows have two exact dependencies at the common circle, by the generic co-circle radius-gradient identity. One hole row is independent of the events because its t-column is-1. Thus the five local rows have rank3, not LICQ.

Subtract the affine common-triangle circle equation from the last height column of every four-site determinant det[1,x,y,x2+y2]. The three shared-site heights vanish; the two remaining heights are exactly e_a,e_b. Each of the five four-subset determinants is therefore A(C)*e_a+B(C)*e_b, and its differential at the joint tie lies in their rank2 span. The saved coefficient arrays and regular-lifting cones use rounded algebraic display positions and are numerical first-order diagnostics. They do not certify a path, exact sector endpoints or a whole-domain covering radius.
''')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();files={str(p.relative_to(ROOT)):sha(p) for p in sorted(ROOT.rglob('*')) if p.is_file()};(ROOT/'MANIFEST_SHA256.json').write_text(json.dumps(files,indent=2)+'\n');ready=dict(status='PERMANENTLY_SEALED_PURE_FIVE_SITE_GEOMETRY',payload_files=len(files),manifest_sha256=sha(ROOT/'MANIFEST_SHA256.json'),constructions=len(records),qualified_geometric_constructions=out['qualified_constructions'],new_wholeD_calls=0,QP_calls=0,optimizer_calls=0,exact_current_binary64_tie_claim=False);(ROOT/'READY.json').write_text(json.dumps(ready,indent=2)+'\n')
    for p in ROOT.rglob('*'):
        if p.is_file():p.chmod(0o444)
    print(json.dumps(ready,indent=2))
if __name__=='__main__':main()
