#!/usr/bin/env python3
"""New Fraction-only rectangle guards; no frozen mathematical imports."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def pair(x):
    return F(x[0], x[1])
def interval(x):
    return pair(x[0]), pair(x[1])
def pack(x):
    if isinstance(x, F):
        return [x.numerator, x.denominator]
    if isinstance(x, dict):
        return {k:pack(v) for k,v in x.items()}
    if isinstance(x, (list, tuple)):
        return [pack(v) for v in x]
    return x

def verify():
    if not __debug__:
        raise RuntimeError('Run without -O: required exact assertions must remain active.')
    pins=json.loads((ROOT/'DEPENDENCY_PINS.json').read_text())
    before={name:sha(ROOT/name) for name in pins}
    assert before==pins
    c=json.loads((ROOT/'dependencies/actual/certificate.json').read_text())
    radius=json.loads((ROOT/'dependencies/radius/endpoint_arithmetic.json').read_text())
    assert c['status']=='PASS_SAME_QUADRATIC_ACTUAL_SIGNED_ANCHOR_FAN_FEE_1_40'
    assert pair(c['lambda'])==F(7,50)
    assert pair(c['tilt'])==F(1071883,20000000)
    assert pair(c['same_score_fee'])==F(1,40)
    L,H=F(7,10),F(24,25)
    oldfee=F(1,40)
    newfee=F(16027,800000)
    eps=F(1,200)
    perturb=5*eps
    lamlo,lamhi=F(27,200),F(29,200)
    t0=F(1071883,20000000)
    tlo,thi=t0-eps,t0+eps
    peaklo,peakhi=interval(c['peak_angle'])
    assert F(4,5)<peaklo<=peakhi<F(17,20)
    C0lo,C0hi=interval(c['accepted_conjugate'])
    assert C0lo>F(147,100)
    curvlo,curvhi=interval(c['curvature_strict_lower_enclosure'])
    assert curvlo>4
    marginL=interval(c['localization_margins'][0])[0]
    marginH=interval(c['localization_margins'][1])[0]
    assert marginL>F(137,10000)
    assert marginH>F(14,1000)
    new_gap_min=oldfee+min(marginL,marginH)-perturb*(H-L)
    assert new_gap_min>newfee
    dL=4*(F(4,5)-L)-perturb
    dH=4*(H-F(17,20))-perturb
    assert dL>0 and dH>0
    concave_guard=F(7,10)+4*F(4,5)-4*lamhi
    assert concave_guard>0
    tail_derivative=F(1,2)-1+thi
    assert tail_derivative<0
    Cref_lower=C0lo-4*eps-eps*H
    assert Cref_lower>F(7,5)
    obtuse_upper=1+thi*F(22,7)
    degenerate_upper=F(1,2)+thi*F(22,7)
    assert F(7,5)-obtuse_upper>newfee
    assert Cref_lower-degenerate_upper>newfee
    assert degenerate_upper<1
    assert 2*H+F(1,10)<3
    for enclosed in c['cut_ratio_endpoint_enclosures']:
        assert interval(enclosed)[0]>F(81,50)
    cut_slope=F(1,2)*F(81,50)**2-thi
    assert cut_slope>1
    assert F(1,10)>newfee
    assert pair(c['paid_gamma_slope_strict_lower'])-eps>1
    sqrt_lower=F(19,100)
    sqrt_upper=F(3873,20000)
    root_square=F(3,80)
    assert sqrt_lower**2<root_square<sqrt_upper**2
    k=1-eps
    assert k==F(199,200)
    assert perturb==F(1,40)
    secant_lower=k-F(3,80)/(2*sqrt_lower)
    assert secant_lower>0
    fee_from_guard=k*oldfee-perturb*sqrt_upper
    assert fee_from_guard==newfee
    assert newfee-F(1,50)==F(27,800000)>0
    Rmax=F(949600,4987489)
    assert F(231,2000)<Rmax
    # Endpoint arithmetic is accepted dependency evidence, not recomputed by
    # executing its source; its exact fraction is pinned along with its proof.
    assert radius['status']=='PASS_NEW_EXACT_RADIUS_ENDPOINT_ARITHMETIC'
    assert pair(radius['lambda'])==F(7,50)
    assert pair(radius['tilt'])==t0
    assert pair(radius['same_local_normalized_fee'])==oldfee
    assert pair(radius['radius_interval']['lower'])==0
    assert radius['radius_interval']['lower_included'] is False
    assert pair(radius['radius_interval']['upper'])==Rmax
    assert radius['radius_interval']['upper_included'] is True
    after={name:sha(ROOT/name) for name in pins}
    assert after==before
    return pack({'status':'PASS_UNIFORM_PARAMETER_ACTUAL_A_LOCAL_FEE',
        'lambda_interval':[lamlo,lamhi], 'tilt_interval':[tlo,thi],
        'radius_upper':F(231,2000), 'fee':newfee,
        'excess_over_1_50':newfee-F(1,50),
        'new_endpoint_gap_lower':new_gap_min,
        'new_G_endpoint_derivative_margins':[dL,dH],
        'restricted_concavity_lower':concave_guard,
        'decreasing_tail_derivative_upper':tail_derivative,
        'acute_support_lower':Cref_lower,
        'obtuse_support_upper':obtuse_upper,
        'convex_cut_loss_slope_lower':cut_slope,
        'degenerate_score_upper':degenerate_upper,
        'sensitivity':'Dnew >= (199/200) Dold - (1/40) sqrt(3 Dold/2)',
        'monotone_secant_lower':secant_lower,
        'square_root_bounds':[sqrt_lower,sqrt_upper],
        'dependency_count':len(pins),'dependency_hashes_before_after_match':True,
        'frozen_checker_calls':0,'frozen_math_source_imports':0,
        'optimizer_calls':0,'geometry_search_calls':0,
        'scope':'Actual signed A common-apex three-face local fee; same-parameter full-apex support. No independent J fee, full-cover row or global radius claim.'})

if __name__=='__main__':
    print(json.dumps(verify(),indent=2,sort_keys=True))
