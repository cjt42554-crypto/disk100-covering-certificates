#!/usr/bin/env python3
"""One new arithmetic-only corollary guard; no geometry source import."""
from fractions import Fraction as F
import json

def verify():
    if not __debug__:
        raise RuntimeError('Run ordinary Python without -O.')
    oldfee=F(19,500)
    fee=F(25473,800000)
    upper=F(191,800)
    lower=F(23,100)
    assert lower**2<F(3,2)*oldfee<upper**2
    assert F(199,200)*oldfee-F(1,40)*upper==fee
    assert F(199,200)-F(3,80)/(2*lower)>0
    assert F(161,5000)-fee==F(287,800000)>0
    assert F(15109287,70000000)>fee
    assert F(1,10)>fee
    assert F(7,5)-F(47890713,70000000)>fee
    K=F(993,1250)
    S=F(1103,1000)**2+F(51,25)**2
    Rmax=K/S
    assert S==F(5378209,1000000)
    assert Rmax==F(794400,5378209)
    assert F(231,2000)<Rmax
    def p(x):
        return [x.numerator,x.denominator]
    return {'status':'PASS_NEW_FRACTION_COROLLARY_ARITHMETIC_ONLY',
        'fee':p(fee),'old_accepted_fee_required':p(oldfee),
        'square_root_bracket':[p(lower),p(upper)],
        'acute_localization_gap_remaining':p(F(161,5000)-fee),
        'radius_upper':p(Rmax),'old_checker_calls':0,
        'old_math_source_imports':0,'optimizer_or_search_calls':0,
        'scope':'Scalar algebra for a local actual-A corollary; full .038 theorem, radius and parameter proofs/reviews required separately.'}

if __name__=='__main__':
    print(json.dumps(verify(),sort_keys=True,indent=2))
