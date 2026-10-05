#!/usr/bin/env python3
"""Exact full-apex conjugate for the fixed new lambda=.14 support row."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,sys

if not __debug__:
    raise RuntimeError('Assertions must remain enabled')
HERE=Path(__file__).resolve().parent
PRIMITIVE=HERE.parent/'coupled_edges/verify_support.py'
PIN='90fd6e3426d8c113bacedd55f9a13f5e47981e9a48d21bfef48183ec4602efec'
if hashlib.sha256(PRIMITIVE.read_bytes()).hexdigest()!=PIN:
    raise ValueError('Frozen arithmetic primitive changed')
sys.path.insert(0,str(PRIMITIVE.parent))
from verify_support import I,PI,sinq,cosq

LAM=Q(7,50)
T=Q(1071883,20000000)

def need(v,m):
    if not v:raise ValueError(m)

def trig(a,fn):
    a=I(a)
    md=(a.a+a.b)/2
    radius=(a.b-a.a)/2
    c=fn(md)
    # sin and cos are globally 1-Lipschitz. The primitive's Taylor
    # remainder bounds are valid at these rational arguments.
    return I(c.a-radius,c.b+radius)

def score(a,t,obtuse=False):
    a=I(a);s=trig(a,sinq);c=trig(a,cosq)
    if obtuse:return s-4*LAM*s*s+t*(PI-a)
    return s*(1+2*c)-4*LAM*s*s+t*a

def derivative(a,t,obtuse=False):
    s,c=sinq(a),cosq(a)
    if obtuse:return c*(1-8*LAM*s)-t
    return c+2*cosq(2*a)-4*LAM*sinq(2*a)+t

def maximum(t,obtuse=False):
    lo=Q(0);hi=PI.a/3
    need(derivative(lo,t,obtuse).a>0,'positive left derivative')
    need(derivative(hi,t,obtuse).b<0,'negative right derivative')
    iterations=0
    while hi-lo>Q(1,10**25):
        md=(lo+hi)/2
        d=derivative(md,t,obtuse)
        if d.a>0:lo=md
        elif d.b<0:hi=md
        else:raise ValueError('derivative sign unresolved before fixed tolerance')
        iterations+=1
    peak=I(lo,hi)
    lower=score((lo+hi)/2,t,obtuse).a
    upper=score(peak,t,obtuse).b
    value=I(lower,upper)
    need(value.b-value.a<Q(1,10**23),'support enclosure width')
    return dict(peak_angle=peak.enc(),value=value.enc(),bisection_iterations=iterations)

def unpack(x):return I(Q(*x[0]),Q(*x[1]))

def verify():
    before=PRIMITIVE.read_bytes()
    need(LAM==Q(7,50) and T==Q(1071883,20000000),'fixed parameters')
    need(Q(17,10)**2<3 and Q(14,25)**2<Q(1,2),
         'rational tail and curvature guards')
    need(Q(3,125)<T<Q(1,2),'both complete branch tails decrease')
    need(-2*Q(17,10)+4*LAM<0,'acute second-derivative guard')
    acute=maximum(T)
    obtuse=maximum(T,True)
    acute0=maximum(Q(0))
    C=unpack(acute['value']);Co=unpack(obtuse['value'])
    C0=unpack(acute0['value']);Co0=I(Q(25,56))
    need(C.a>Co.b and C0.a>Co0.b,'strict branch dominance')
    need(PRIMITIVE.read_bytes()==before,'frozen primitive unchanged')
    return dict(status='PASS_FIXED_LAMBDA14_FULL_APEX_CONJUGATE',
                lambda_fraction=[LAM.numerator,LAM.denominator],
                tilt_fraction=[T.numerator,T.denominator],
                acute_branch=acute,obtuse_branch=obtuse,
                conjugate=C.enc(),acute_obtuse_margin=(C-Co).enc(),
                zero_tilt_acute=acute0,zero_tilt_obtuse=Co0.enc(),
                zero_tilt_conjugate=C0.enc(),
                conjugate_difference=(C-C0).enc(),
                primitive_sha256=PIN,old_batches_run=0,optimizers_run=0,
                scope='Full-apex single-face conjugate at lambda=7/50 and the fixed '
                      'positive tilt, plus zero tilt. No joint fan fee, arc bootstrap, '
                      'global radius exclusion or optimality certificate.')

if __name__=='__main__':
    print(json.dumps(verify(),indent=2))
