#!/usr/bin/env python3
"""New broad-corollary rational guards only; accepted geometry external."""
from fractions import Fraction as F
import json
def verify():
    if not __debug__:
        raise RuntimeError('Run ordinary Python without -O.')
    eps=F(7,500)
    t0=F(1071883,20000000)
    tmin,tmax=t0-eps,t0+eps
    lammin,lammax=F(7,50)-eps,F(7,50)+eps
    fee=F(41,2000)
    assert lammin==F(63,500) and lammax==F(77,500)
    assert tmin==F(791883,20000000) and tmax==F(1351883,20000000)
    assert F(387,10000)-5*eps*F(13,50)==fee
    Cref=F(147,100)-4*eps-eps*F(24,25)
    assert Cref==F(17507,12500)>F(7,5)
    assert F(7,5)-(1+tmax*F(22,7))>fee
    assert Cref-(F(1,2)+tmax*F(22,7))>fee
    assert F(7,10)+4*F(4,5)-4*lammax>0
    assert tmax<F(1,2)
    assert F(1,2)*F(81,50)**2-tmax>1
    assert F(2,5)-5*eps>0 and F(11,25)-5*eps>0
    oldfee=F(19,500)
    assert F(23,100)**2<F(3,2)*oldfee<F(191,800)**2
    f_lower=(1-eps)*oldfee-5*eps*F(191,800)
    assert f_lower==F(41511,2000000)>fee
    assert 1-eps-F(21,200)/F(46,100)>0
    curvature=4-8*eps
    assert curvature==F(486,125)
    assert 1/(2*curvature)==F(125,972)
    assert 5*eps/curvature==F(7,388)
    def p(x): return [x.numerator,x.denominator]
    return {'status':'PASS_NEW_BROAD_FRACTION_GUARDS_ONLY',
        'lambda_interval':[p(lammin),p(lammax)],
        'tilt_interval':[p(tmin),p(tmax)],'fee':p(fee),
        'sensitivity_threshold_strict_lower':p(f_lower),
        'radius_upper':[794400,5378209],
        'cap_coefficient':p(1/(2*curvature)),
        'peak_displacement_upper':p(5*eps/curvature),
        'old_checker_calls':0,'old_math_imports':0,
        'scope':'Broad actual signed-A local corollary. Accepted .038 theorem/radius and full analytic proof required. No whole-cover radius exclusion.'}
if __name__=='__main__':
    print(json.dumps(verify(),sort_keys=True,indent=2))
