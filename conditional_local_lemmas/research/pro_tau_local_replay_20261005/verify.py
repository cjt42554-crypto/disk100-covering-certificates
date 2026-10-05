#!/usr/bin/env python3
"""Exact scalar guards for the two-actual-anchor tau lemma.

Only Python's standard library is used. There is no optimizer, sampled
configuration, old three-face checker, or global-cover enumeration here.

All interval endpoints are integer multiples of 10**(-70). Trigonometric
point evaluations use a finite Taylor polynomial and the Lagrange remainder.
The geometric identities and the continuous-domain reductions are proved
in PROOF.md; this program verifies their explicit numerical guards.
"""
from __future__ import annotations
from fractions import Fraction as Q
from math import factorial, isqrt
from pathlib import Path
import json

SCALE = 10**70

def ceildiv(a: int, b: int) -> int:
    if b <= 0:
        raise ValueError('ceildiv requires a positive denominator')
    return -((-a)//b)

class I:
    __slots__ = ('lo', 'hi')
    def __init__(self, lo: int, hi: int):
        if lo > hi:
            raise ValueError('reversed interval')
        self.lo, self.hi = lo, hi
    @classmethod
    def of(cls, x: Q | int | 'I') -> 'I':
        if isinstance(x, I):
            return x
        x = Q(x)
        return cls((x.numerator*SCALE)//x.denominator,
                   ceildiv(x.numerator*SCALE, x.denominator))
    @classmethod
    def hull(cls, a: Q, b: Q) -> 'I':
        aa, bb = cls.of(a), cls.of(b)
        return cls(aa.lo, bb.hi)
    def __add__(self, other):
        other=I.of(other); return I(self.lo+other.lo, self.hi+other.hi)
    __radd__=__add__
    def __neg__(self): return I(-self.hi,-self.lo)
    def __sub__(self, other): return self + (-I.of(other))
    def __rsub__(self, other): return I.of(other) + (-self)
    def __mul__(self, other):
        other=I.of(other)
        vals=[self.lo*other.lo,self.lo*other.hi,self.hi*other.lo,self.hi*other.hi]
        return I(min(vals)//SCALE,ceildiv(max(vals),SCALE))
    __rmul__=__mul__
    def reciprocal(self):
        if self.lo <= 0 <= self.hi:
            raise ValueError('division by an interval containing zero')
        if self.hi < 0:
            return -((-self).reciprocal())
        return I(SCALE*SCALE//self.hi,ceildiv(SCALE*SCALE,self.lo))
    def __truediv__(self, other): return self * I.of(other).reciprocal()
    def __rtruediv__(self, other): return I.of(other) * self.reciprocal()
    def __pow__(self, n: int):
        if not isinstance(n,int) or n < 0:
            raise ValueError('only nonnegative integer powers are supported')
        out=I.of(1); base=self
        while n:
            if n & 1: out=out*base
            n//=2
            if n: base=base*base
        return out
    def sqrt(self):
        if self.lo < 0: raise ValueError('negative square-root lower endpoint')
        low=isqrt(self.lo*SCALE)
        x=self.hi*SCALE; high=isqrt(x)
        if high*high != x: high+=1
        return I(low,high)
    def abs(self):
        if self.lo>=0: return self
        if self.hi<=0: return -self
        return I(0,max(-self.lo,self.hi))
    def loq(self): return Q(self.lo,SCALE)
    def hiq(self): return Q(self.hi,SCALE)
    def exact(self):
        return {'lower':str(self.loq()), 'upper':str(self.hiq())}

def imin(*xs):
    xs=[I.of(x) for x in xs]
    return I(min(x.lo for x in xs),min(x.hi for x in xs))
def imax(*xs):
    xs=[I.of(x) for x in xs]
    return I(max(x.lo for x in xs),max(x.hi for x in xs))

# Taylor degree 79 for sine, 78 for cosine, both expanded at zero.
# The derivatives of sine/cosine have absolute value <=1 everywhere.
def sin_point(x: Q) -> I:
    x=Q(x)
    if abs(x)>3: raise ValueError('point outside certified evaluation range')
    val=Q(0)
    for k in range(40):
        val+=(-1)**k*x**(2*k+1)/factorial(2*k+1)
    err=abs(x)**80/factorial(80)
    return I.hull(val-err,val+err)
def cos_point(x: Q) -> I:
    x=Q(x)
    if abs(x)>3: raise ValueError('point outside certified evaluation range')
    val=Q(0)
    for k in range(40):
        val+=(-1)**k*x**(2*k)/factorial(2*k)
    err=abs(x)**79/factorial(79)
    return I.hull(val-err,val+err)

R=Q(5661,50000)
LAM=Q(7,50)
T0=Q(1071883,20000000)
D_MAX=Q(9,400)
BUDGET_SKEW=Q(5,2)*Q(14951,15000)**2
KAPPA=Q(41,200)
REM=Q(13,20)
TAU_COST=Q(11,50)
CHI_BAR=Q(261,500)
CHI_SLOPE=Q(71,100)
FEE=Q(178523206195158271774904671901,
      8011730250000000000000000000000)
SL,SU=Q(753,1000),Q(113,125)
WL,WU=-Q(12,125),Q(813,8000)
ALPHA_LO,ALPHA_HI=Q(2067,2000),Q(11661,8000)
H_LO,H_HI=Q(2199,8000),Q(137,250)
A_LO,A_HI=Q(411,2000),Q(5029,8000)
DELTA_MAX=Q(3,20)

GUARDS=[]

def guard(name: str, gap, margin: Q=Q(0)):
    gap=I.of(gap); margin=Q(margin)
    if not gap.loq()>margin:
        raise AssertionError(f'{name}: guard failed: {gap.exact()} not > {margin}')
    GUARDS.append({'name':name,'gap':gap.exact(),'required_strict_margin':str(margin)})
    return gap

def Ffun(x: Q):
    s=sin_point(x); c=cos_point(x)
    return s*(1+2*c)-4*LAM*s*s+T0*x

def Fprime(x: Q):
    return cos_point(x)+2*cos_point(2*x)-4*LAM*sin_point(2*x)+T0

def rho_product(alpha: Q):
    b=sin_point(alpha); u=(1-R*R*b*b).sqrt()
    rho=u-R*cos_point(alpha)
    return rho, rho*b

def run():
    GUARDS.clear()
    # Self-contained support/localization guards. No three-face gap is used.
    guard('support: second apex branch < 6/5', Q(6,5)-(1+T0*Q(22,7)))
    guard('support: F(5/6)>7/5', Ffun(Q(5,6))-Q(7,5))
    guard('peak lies to the right of .828', Fprime(Q(207,250)))
    guard('peak lies to the left of .829', -Fprime(Q(829,1000)))
    guard('theta below .7 loses more than .035',
          Ffun(Q(5,6))-Ffun(Q(7,10)),Q(7,200))
    guard('theta above .96 loses more than .035',
          Ffun(Q(5,6))-Ffun(Q(24,25)),Q(7,200))
    guard('apex loss globally increases on theta window',
          1-cos_point(Q(7,10))-T0)
    guard('apex loss derivative >1 up to z=.1',
          Q(1,2)*Q(3871,2500)**2-T0-1)
    guard('old F curvature >4 on [.7,.96]',
          sin_point(Q(7,10))+4*sin_point(Q(48,25))
          +Q(28,25)*cos_point(Q(48,25))-4)
    # Budget -> the uniform parameter box.
    guard('Dmax smaller than localization threshold',Q(7,200)-D_MAX)
    guard('individual theta lower .72 from budget',
          (Q(207,250)-Q(18,25))**2-D_MAX/2)
    guard('individual theta upper .94 from budget',
          (Q(47,50)-Q(829,1000))**2-D_MAX/2)
    guard('sin(.72)^2 >13/30',sin_point(Q(18,25))**2-Q(13,30))
    guard('sin(.94)<.81',Q(81,100)-sin_point(Q(47,50)))
    guard('sigma coefficient >=5/4',
          2*sin_point(Q(18,25))**2/sin_point(Q(297,400))-Q(5,4))
    guard('|d_i|<.14 from sigma',sin_point(Q(7,50))**2-Q(9,500))
    guard('|eta|<.096 from sigma',Q(12,125)**2
          -Q(2,5)*Q(15000,14951)**2*D_MAX)
    # The exact arithmetic endpoints of the box.
    assert ALPHA_LO == Q(3,2)*SL+WL
    assert ALPHA_HI == Q(3,2)*SU+WU
    assert H_LO == SL/2-WU and H_HI == SU/2-WL
    assert A_LO == SL/2+WL-DELTA_MAX/2
    assert A_HI == SU/2+WU+DELTA_MAX/2
    guard('all alpha,h,A are in first quadrant',Q(3,2)-ALPHA_HI)
    guard('B*=sin(alpha)>6/7',sin_point(ALPHA_LO)-Q(6,7))
    guard('cos(alpha)>.1',cos_point(ALPHA_HI)-Q(1,10))
    guard('sin(A_i)>.2025',sin_point(A_LO)-Q(81,400))
    guard('u>.99 along homotopy',1-R*R-Q(99,100)**2)
    guard('B(t)>.83 along homotopy',
          Q(6,7)-Q(15,13)*D_MAX-Q(83,100))
    guard('negative translation branch excluded',
          Q(1,10)*Q(83,100)-Q(15,13)*D_MAX)
    neg_x=R*Q(99,100)-R*R-R*Q(7,40)
    guard('signed tangential directions for linear support',
          neg_x/Q(5)-Q(23,20)*R*R)
    # Eight interval slabs, not samples. The delta dependence has already
    # been minimized analytically as A+Q*cos(delta)+P*sin(delta).
    slab_lowers=[]
    for j in range(8):
        sl=SL+(SU-SL)*Q(j,8)
        su=SL+(SU-SL)*Q(j+1,8)
        amin=Q(3,2)*sl+WL; amax=Q(3,2)*su+WU
        hmin=sl/2-WU; hmax=su/2-WL
        rho_min,prod_min=rho_product(amin)
        rho_max,prod_max=rho_product(amax)
        A0=KAPPA/2-R*prod_max
        Qlo=-KAPPA/2*cos_point(2*sl)+R*R*prod_min*cos_point(hmax)
        sin2lo=imin(sin_point(2*sl),sin_point(2*su))
        Plo=KAPPA/2*sin2lo-R*rho_max*cos_point(hmin)
        Phi=KAPPA/2-R*rho_min*cos_point(hmax)
        Pabs=imax(Plo.abs(),Phi.abs())
        # cos(delta) is positive. min(Qlo,Qlo*cos(delta_max))
        # is valid regardless of the sign of Qlo.
        lower=A0+imin(Qlo,Qlo*cos_point(DELTA_MAX)) \
                  -Pabs*sin_point(DELTA_MAX)
        guard(f'linear tau coefficient slab {j+1}/8: kappa_i<41/200',
              lower,Q(1,4000))
        slab_lowers.append({'s_lower':str(sl),'s_upper':str(su),
                            'gap':lower.exact()})
    # Uniform Taylor-remainder estimate for f(b)=sqrt(1-R²|b|²) Jcw b/|b|.
    guard('2/sqrt(3)<7/6 after squaring',Q(7,6)**2-Q(4,3))
    M1=Q(100,83)
    guard('radial first derivative <M1',M1-R*R*Q(100,99))
    M2bound=Q(7,6)*Q(100,99)*Q(100,83)**2+R*R*Q(100,99)**3
    guard('directional Hessian norm <7/4',Q(7,4)-M2bound)
    remainder=(R*R+R*Q(7,4)+2*R*M1)*Q(15,13)**2
    guard('mean-norm quadratic remainder <13/20 T²',REM-remainder)
    guard('total tau cost <11/50 T',TAU_COST-KAPPA-REM*D_MAX)
    # Signed halfarc enclosure, optional in the fee proof.
    guard('halfarc quadratic coefficient <81/100',
          Q(81,100)-Q(15,13)**2/(2*Q(83,100)))
    # A scalar closing argument, not a configuration optimization.
    qmax=R*cos_point(ALPHA_LO)/(1-R*R*sin_point(ALPHA_LO)**2).sqrt()
    guard('chi s derivative <3/5',Q(3,5)-(Q(1,2)+Q(3,2)*qmax))
    guard('baseline chi <.522',
          sin_point(Q(43,400))-R*sin_point(Q(2487,2000)))
    guard('Cauchy coefficient <.71',CHI_SLOPE**2-Q(9,100)-1/BUDGET_SKEW)
    guard('scalar norm envelope is increasing in used shape budget',
          R*CHI_SLOPE*sin_point(CHI_BAR)/Q(3,20)-TAU_COST,Q(1,25))
    final_angle=CHI_BAR+CHI_SLOPE*Q(3,20)
    assert final_angle==Q(1257,2000)
    mean_upper=1+R*R-2*R*cos_point(final_angle)
    guard('final mean outer squared norm <83/100',Q(83,100)-mean_upper,Q(3,10000))
    inherited_rho=Q(911088875,10**9)
    guard('given highM norm floor squared >83/100',inherited_rho**2-Q(83,100))
    # Direct arc-to-norm calculation is reported separately, not needed above.
    ah=Q(147,2000)
    arc_radicand=R*R-sin_point(ah)**2
    guard('full arc .147 has a real radial inner endpoint',arc_radicand)
    arc_rho=cos_point(ah)-arc_radicand.sqrt()
    guard('full assigned arc>.147 implies norm>.911124459',
          arc_rho-Q(911124459,10**9))
    guard('new fee exceeds requested fee by >1/5000',D_MAX-FEE,Q(1,5000))
    summary={
        'R':str(R),'lambda':str(LAM),'t':str(T0),
        'new_two_face_gap':str(D_MAX),
        'requested_fee':str(FEE),
        'fee_margin_exact':str(D_MAX-FEE),
        'outer_squared_norm_floor':str(Q(83,100)),
        'mean_norm_upper_under_D_le_gap':mean_upper.exact(),
        'arc_implied_norm_floor':arc_rho.exact(),
        'linear_tau_coefficient':str(KAPPA),
        'quadratic_tau_remainder':str(REM),
        'uniform_total_tau_coefficient':str(TAU_COST),
        'sigma_eta_coefficient':str(BUDGET_SKEW),
        'linear_support_slabs':slab_lowers,
        'guards':GUARDS,
        'guard_count':len(GUARDS),
        'scope':'Two actual adjacent triangles only; not a global 100-disk certificate.'
    }
    return summary

if __name__=='__main__':
    data=run()
    path=Path(__file__).with_name('verification.json')
    path.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    for item in GUARDS:
        print('PASS:',item['name'])
    print(f'PASS: {len(GUARDS)} exact guards, including 8 full interval slabs.')
    print('THEOREM: if |c1|^2<=1-R^2 and |c0|^2,|c2|^2>=83/100,')
    print('         then 2*C-S2 > 9/400 for the stated actual-anchor geometry.')
    print('         No middle-halfarc variance is charged.')
    print('NOT A GLOBAL r_D(100) CERTIFICATE.')
