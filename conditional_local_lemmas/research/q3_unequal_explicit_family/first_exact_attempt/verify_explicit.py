#!/usr/bin/env python3
"""One exact algebraic unequal fan; no optimizer, parameter scan or cover search."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,sys

if not __debug__:raise RuntimeError('Assertions must remain enabled')
HERE=Path(__file__).resolve().parent
BACKEND=HERE.parent/'coupled_edges/verify_support.py'
PIN='90fd6e3426d8c113bacedd55f9a13f5e47981e9a48d21bfef48183ec4602efec'
CFILE=HERE.parent/'quadratic_lambda14_conjugate/certificate.json'
CPIN='f325d42bfe651599881b7dc652567729a1586389ee385c38060840bc41696bc2'
if hashlib.sha256(BACKEND.read_bytes()).hexdigest()!=PIN:raise ValueError('primitive changed')
if hashlib.sha256(CFILE.read_bytes()).hexdigest()!=CPIN:raise ValueError('conjugate input changed')
sys.path.insert(0,str(BACKEND.parent))
from verify_support import I,PI,sinq,cosq,asin

R=Q(11321,100000);LAM=Q(7,50);T=Q(1071883,20000000)
A=Q(359,500);B=Q(1691,2500);ETA=Q(1,10**7)
HALF_TAN=Q(461,10000)
TARGET=Q(19931551,10**9)

def need(v,m):
    if not v:raise ValueError(m)
def encq(q):return [q.numerator,q.denominator]
def unpack(x):return I(Q(*x[0]),Q(*x[1]))
def add(a,b):return tuple(a[j]+b[j] for j in (0,1))
def sub(a,b):return tuple(a[j]-b[j] for j in (0,1))
def dot(a,b):return a[0]*b[0]+a[1]*b[1]
def norm2(a):return dot(a,a)
def dist2(a,b):return norm2(sub(a,b))
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def reflect(a):return a[0],-a[1]
def ep(a):return [z.enc() for z in a]

def cubic_cosine():
    # cos(3u)=-cos(2a), u=(pi-2a)/3 in (0,pi/3).
    f=lambda z:4*z*z*z-3*z+(2*A*A-1)
    lo=Q(85,100);hi=Q(87,100)
    need(f(lo)<0<f(hi) and 12*lo*lo-3>0,'cubic root isolation')
    while hi-lo>Q(1,10**28):
        md=(lo+hi)/2
        if f(md)<0:lo=md
        else:hi=md
    return I(lo,hi),[encq(lo),encq(hi)]

def verify():
    before={str(p):p.read_bytes() for p in (BACKEND,CFILE)}
    V=1-R*R;n=I(V).sqrt()
    sa=I(1-A*A).sqrt();sb=I(1-B*B).sqrt()
    a=PI/2-asin(I(A));b=PI/2-asin(I(B))
    k,kroot=cubic_cosine();su=(1-k*k).sqrt()
    rho=(1-R*R*(2-B*B)).sqrt()-R*(1+B)-ETA
    p=(rho,I(0))
    x2=2*A*A-1;y2=2*A*sa
    c0=(rho+R*(x2+B),-R*(y2+sb))
    c1=(rho+R*(1+B),-R*sb)
    c2=reflect(c1);c3=reflect(c0)
    d0=(rho+2*R*k*(1-2*k*k),4*R*k*k*su)
    d1=(rho-2*R*k,I(0));d2=reflect(d0)
    centers=[p,c0,c1,c2,c3,d0,d1,d2]
    # Exact algebraic pi/16 rotation (half-angle definitions).
    c8=((1+I(2).sqrt()/2)/2).sqrt()
    s8=((1-I(2).sqrt()/2)/2).sqrt()
    cs=((1+c8)/2).sqrt();ss=((1-c8)/2).sqrt()
    direction=(I(1),I(0));rings=[]
    for j in range(32):
        if j not in (0,1,31):
            q=(n*direction[0],n*direction[1])
            rings.append((j,q));centers.append(q)
        direction=(cs*direction[0]-ss*direction[1],ss*direction[0]+cs*direction[1])
    need(len(centers)==37 and len(rings)==29,'one interior degree-seven star plus 33 boundary sites')

    v=sub(c1,c0);l2=norm2(v);delta=(norm2(c1)-norm2(c0))/2
    disc=l2-delta*delta
    need(l2.a>0 and disc.a>0,'canonical ordinary pair root exists')
    zz=disc.sqrt()
    q0=((delta*v[0]+zz*v[1])/l2,(delta*v[1]-zz*v[0])/l2)
    q1=(I(1),I(0));q2=reflect(q0)
    need(q0[0].a>0 and q0[1].b<0,'negative-side canonical root')
    other=((delta*v[0]-zz*v[1])/l2,(delta*v[1]+zz*v[0])/l2)
    need(other[0].b<0 and other[1].a>0,'other pair root is in far upper-left quadrant')
    sinrot=2*HALF_TAN/(1+HALF_TAN*HALF_TAN)
    cosrot=(1-HALF_TAN*HALF_TAN)/(1+HALF_TAN*HALF_TAN)
    need(cosrot*cosrot+sinrot*sinrot==1,'rational unit rotation')
    need(sinrot>sinq(Q(23,250)).b and cosrot>cosq(Q(1,10)).b,
         'terminal owned segment is strictly longer than .092 and shorter than .1')
    left=(cosrot*q0[0]+sinrot*q0[1],-sinrot*q0[0]+cosrot*q0[1])
    right=reflect(left)
    need(-q0[1].b>sinq(Q(23,250)).b and -q0[1].a<sinq(Q(1,5)).a,
         'middle assigned arcs strictly between .092 and .2')

    cap_slacks=[];norm_low=[]
    for c in (c0,c1,c2,c3):
        cap_slacks.append(I(V)-norm2(c));norm_low.append(norm2(c)-(1-R)*(1-R))
    need(all(z.a>0 for z in cap_slacks+norm_low),'all four reflected boundary norms are strict')
    need(1-rho.b>R,'p can never cover a unit-circle point')
    interior_slacks=[I((1-R)*(1-R))-norm2(z) for z in (p,d0,d1,d2)]
    need(all(z.a>0 for z in interior_slacks),'four interior sites cannot be closest on the covered circle')
    common=[]
    for qi,(ci,cj) in zip((q0,q1,q2),((c0,c1),(c1,c2),(c2,c3))):
        sl=[I(R*R)-dist2(qi,z) for z in (ci,cj)]
        need(all(z.a>0 for z in sl),'canonical transition common coverage strict')
        common.append([z.enc() for z in sl])
    terminal=[I(R*R)-dist2(left,c0),I(R*R)-dist2(right,c3)]
    need(all(z.a>Q(1,100000) for z in terminal),'true terminal coverage on strictly longer arcs')

    # No ring center can compete anywhere in [arg(left),arg(right)]:
    # nearest angular ring directions are +/-pi/8 and distance there
    # is minimized at the corresponding endpoint.
    firstpos=next(q for j,q in rings if j==2)
    firstneg=next(q for j,q in rings if j==30)
    ring_protected=[dist2(left,firstneg)-R*R,dist2(right,firstpos)-R*R]
    need(all(z.a>0 for z in ring_protected),'ring disks stay outside the entire protected fan interval')
    left_ring_edge=(n*c8+R*s8,R*c8-n*s8)
    right_ring_edge=reflect(left_ring_edge)
    overlap=[I(R*R)-dist2(left_ring_edge,c0),I(R*R)-dist2(right_ring_edge,c3)]
    need(all(z.a>0 for z in overlap),'boundary-ring footprints overlap the outer fan disks')
    sin32=((1-cs)/2).sqrt()
    need(sin32.b<R,'pi/16 ring spacings are shorter than twice their radius-R footprint half-angle')

    # Non-adjacent ownership cuts for the four real assigned arcs.
    # q dot(c2-c0)-delta02 increases for negative angles since both
    # difference coordinates are positive, and is negative at q0.
    v02=sub(c2,c0);delta02=(norm2(c2)-norm2(c0))/2
    guard02=dot(q0,v02)-delta02
    need(v02[0].a>0 and v02[1].a>0 and guard02.b<0,'c2 never cuts the negative outer assigned arc')
    v13=sub(c3,c1);delta13=(norm2(c3)-norm2(c1))/2
    guard13=v13[0]*cosq(Q(1,5))-delta13
    need(v13[0].b<0 and v13[1].a>0 and guard13.b<0,
         'c3 never cuts the negative middle assigned arc')
    # c0/c3 and c1/c2 comparisons are exact vertical equal-norm pairs.
    need(c0[0].a>0 and c0[1].b<0 and c1[0].a>0 and c1[1].b<0,
         'all fan arcs lie in the concave positive-dot coverage window')
    dot_window=min(c0[0].a,c1[0].a)*cosq(Q(3,10)).a-max(-c0[1].a,-c1[1].a)*sinq(Q(3,10)).b
    need(dot_window>0,'short-arc coverage from endpoint coverage is valid')

    Os=[(rho+R*B,-R*sb),(rho+R,I(0)),(rho+R*B,R*sb),
        (rho+R*x2,R*y2),(rho-R*k,R*su),(rho-R*k,-R*su),(rho+R*x2,-R*y2)]
    cycle=[1,2,3,4,5,6,7]
    faces=[];min_empty=None
    for j,O in enumerate(Os):
        f=[0,cycle[j],cycle[(j+1)%7]]
        area=cross(sub(centers[f[1]],p),sub(centers[f[2]],p))/2
        need(area.a>0,'seven triangle orientations strict')
        empty=[]
        for h,z in enumerate(centers):
            if h in f:continue
            margin=dist2(z,O)-R*R
            need(margin.a>0,'all nonincident sites outside each exact circumcircle')
            empty.append(margin.a)
        m=min(empty);min_empty=m if min_empty is None else min(min_empty,m)
        faces.append(dict(vertices=f,circumcenter=ep(O),exact_circumradius=encq(R),area=area.enc(),empty_margin_lower=encq(m)))

    area_sum=sum((unpack(f['area']) for f in faces[:3]),I(0))/R/R
    S=area_sum+2*sa*A+sb*B-4*LAM*(2*(1-A*A)+(1-B*B))+T*(2*a+b)
    C=unpack(json.loads(CFILE.read_text())['conjugate'])
    gap=3*C-S
    need(Q(1,100)<gap.a<=gap.b<Q(19,1000)<TARGET,'same-score local loss strictly below requested .019931551')
    # Compute the EXISTING middle assigned-lens slack, without spending it.
    sh=((1-q0[0])/2).sqrt();ch=((1+q0[0])/2).sqrt()
    hmax=(R*R-sh*sh).sqrt()
    radial=c1[0]*ch-c1[1]*sh
    h=ch-radial
    P=sh*(hmax-h)
    need(h.a>0 and P.a>0,'exact inward middle assigned lens has positive old H-area slack')
    existing_middle_payment=2*P/R/R
    need(existing_middle_payment.a>Q(4,5),'large existing middle-cell lens loss is present in this witness')
    need(all(q.read_bytes()==raw for q,raw in before.items()),'frozen primitive/conjugate unchanged')
    return dict(status='PASS_ONE_ALGEBRAIC_UNEQUAL_Q3_D7_LOCAL_OBSTRUCTION',
                R=encq(R),lambda_fraction=encq(LAM),t=encq(T),cos_a=encq(A),cos_b=encq(B),radial_eta=encq(ETA),
                algebraic_inner_cosine_root=kroot,inner_cosine_cubic=['4','0','-3',str(2*A*A-1)],
                p=ep(p),four_boundary_centers=[ep(z) for z in (c0,c1,c2,c3)],
                three_interior_neighbors=[ep(z) for z in (d0,d1,d2)],
                ring_centers=[{'angle_index_pi_over_16':j,'center':ep(z)} for j,z in rings],
                boundary_center_count=33,interior_site_count=4,total_site_count=37,
                actual_apex_angles=[a.enc(),b.enc(),a.enc()],canonical_switches=[ep(z) for z in (q0,q1,q2)],
                other_first_pair_root=ep(other),strict_terminal_point=ep(left),
                terminal_rotation_half_tangent=encq(HALF_TAN),boundary_norm_cap_slacks=[z.enc() for z in cap_slacks],
                boundary_norm_lower_slacks=[z.enc() for z in norm_low],interior_circle_exclusion_slacks=[z.enc() for z in interior_slacks],
                common_switch_cover_slacks=common,terminal_cover_slacks=[z.enc() for z in terminal],
                ring_protected_interval_exclusion_slacks=[z.enc() for z in ring_protected],
                ring_fan_footprint_overlap_slacks=[z.enc() for z in overlap],
                nearest_guard02=guard02.enc(),nearest_guard13=guard13.enc(),
                seven_exact_empty_circumcircle_faces=faces,min_empty_squared_margin=encq(min_empty),
                whole_score=S.enc(),single_face_conjugate=C.enc(),same_score_loss=gap.enc(),
                target_fee=encq(TARGET),target_failure_margin=(TARGET-gap).enc(),
                existing_middle_cell_lens_slack_normalized=existing_middle_payment.enc(),
                primitive_sha256=PIN,accepted_conjugate_sha256=CPIN,
                optimizer_calls=0,graph_catalogues_run=0,old_batches_run=0,
                scope='Explicit 37-site diagram with 33 covered boundary-cell arcs, an interior degree-seven q3 apex and seven empty circumradius-R star faces. S loss<.019. The entire unit circle is covered. The whole closed unit disk, a 100-site optimum, all other face circumradii, the global profile and outer-Q condition are NOT certified. This blocks the requested standalone local S fee from these local premises; existing middle-cell lens slack is computed once and not appended as a new fee.')

if __name__=='__main__':print(json.dumps(verify(),indent=2))
