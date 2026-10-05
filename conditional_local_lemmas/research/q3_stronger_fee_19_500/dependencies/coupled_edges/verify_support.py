#!/usr/bin/env python3
"""Exact rational arithmetic for a conditional, nonsymmetric boundary lemma.
Does NOT alone prove a universal lower bound: the m=100, B=33..35,
and no-ear hypotheses must be established separately at the same radius.
"""
from fractions import Fraction as Q
from math import factorial,isqrt
from pathlib import Path
import json
DIGITS=35
SCALE=10**DIGITS
R=Q(9041,80000)
K=Q(6,25)
T1=Q(3,100)

def outward(a,b):
 return Q(a.numerator*SCALE//a.denominator,SCALE),Q(-((-b.numerator*SCALE)//b.denominator),SCALE)
class I:
 def __init__(self,a,b=None):
  if isinstance(a,I):self.a,self.b=a.a,a.b;return
  a=Q(a);b=a if b is None else Q(b)
  assert a<=b
  self.a,self.b=outward(a,b)
 def __add__(x,y):
  y=I(y);return I(x.a+y.a,x.b+y.b)
 __radd__=__add__
 def __neg__(x):return I(-x.b,-x.a)
 def __sub__(x,y):return x+-I(y)
 def __rsub__(x,y):return I(y)+-x
 def __mul__(x,y):
  y=I(y);v=[x.a*y.a,x.a*y.b,x.b*y.a,x.b*y.b];return I(min(v),max(v))
 __rmul__=__mul__
 def __truediv__(x,y):
  y=I(y);assert y.a*y.b>0
  return x*I(1/y.b,1/y.a)
 def __rtruediv__(x,y):return I(y)/x
 def sqrt(x):
  assert x.a>=0
  la=isqrt(x.a.numerator*SCALE*SCALE//x.a.denominator)
  ub=isqrt(x.b.numerator*SCALE*SCALE//x.b.denominator)+1
  a,b=Q(la,SCALE),Q(ub,SCALE)
  assert a*a<=x.a and b*b>x.b
  return I(a,b)
 def enc(x):return [[x.a.numerator,x.a.denominator],[x.b.numerator,x.b.denominator]]

def atan(q,N=60):
 x=Q(1,q);v=sum(((-1)**j*x**(2*j+1)/Q(2*j+1) for j in range(N)),Q(0))
 e=(-1)**N*x**(2*N+1)/Q(2*N+1)
 return I(min(v,v+e),max(v,v+e))
PI=16*atan(5)-4*atan(239)
assert Q(3)<PI.a<PI.b<Q(22,7)
def sinq(x):
 v=sum(((-1)**j*x**(2*j+1)/factorial(2*j+1) for j in range(30)),Q(0))
 e=abs(x)**61/factorial(61)
 return I(v-e,v+e)
def cosq(x):
 v=sum(((-1)**j*x**(2*j)/factorial(2*j) for j in range(30)),Q(0))
 e=abs(x)**60/factorial(60)
 return I(v-e,v+e)
def sin(x):
 x=I(x);assert 0<=x.a<=x.b<PI.a/2
 return I(sinq(x.a).a,sinq(x.b).b)
def cos(x):
 x=I(x);assert 0<=x.a<=x.b<PI.a/2
 return I(cosq(x.b).a,cosq(x.a).b)
def asinq(x):
 assert 0<=x<1
 lo,hi=Q(0),PI.b/2
 for _ in range(125):
  md=(lo+hi)/2
  sm=sinq(md)
  if sm.b<x:lo=md
  elif sm.a>x:hi=md
  else:
   # At an unresolved tiny midpoint, use derivative cos>=sqrt(1-x^2)/2
   # only after bounding md inside a fixed interval and derivative away from zero.
   if hi-lo<Q(1,10**30):break
   raise ValueError('insufficient sine precision during arcsin enclosure')
  if hi-lo<Q(1,10**30):break
 return I(lo,hi)
def asin(x):
 x=I(x);return I(asinq(x.a).a,asinq(x.b).b)

def verify():
 assert 0<R<Q(1,5) and K==Q(6,25)
 # The analytically proved concavity numerator estimate for the bad-arc penalty.
 penalty_second_numerator=-Q(266,625)+Q(18744,390625)
 assert penalty_second_numerator==-Q(147506,390625)<0
 y=(-(1-K)+I((1-K)*(1-K)+32).sqrt())/8
 xs=(1-y*y).sqrt()
 intercept=xs*((1-K)+2*y)
 # G' = K: 4 y^2 +(1-K)y-2=0, and positive y is unique.
 assert y.a>0 and y.b<1 and xs.a>0 and xs.b<1
 tc=asin(I(R)*I(K).sqrt())
 assert tc.a>T1
 s1=sin(T1);c1=cos(T1);h1=(R*R-s1*s1).sqrt();alpha1=asin(s1/R)
 rows=[]
 for B in (33,34,35):
  t0=PI/B;s=sin(t0);c=cos(t0);h=(R*R-s*s).sqrt()
  sineta=s*(c-h)/R
  coseta=(h*c+s*s)/R
  alpha0=asin(s/R)
  b=R*R*K*coseta
  assert b.b<R*R
  H0=s*h+t0-s*c
  F0=H0+b*(alpha0-t0)
  D0=c*(R*R-2*s*s)/h+2*s*s+b*(c/h-1)
  F1=s1*h1+T1-s1*c1+b*(alpha1-T1)
  D1=c1*(R*R-2*s1*s1)/h1+2*s1*s1+b*(c1/h1-1)
  u=s1/R;v=(1-u*u).sqrt()
  P1=2*R*R*v*(K-u*u)
  PD1=2*R*u*c1*(3*u*u-K-2)/v
  # Difference of affine upper support for F+penalty at t1,
  # from upper support for F at t0. Check both bad-arc endpoints.
  a=F0-t0*D0-F1-P1+T1*(D1+PD1)
  z=D0-D1-PD1
  margins=[a,a+tc*z]
  assert all(m.a>0 for m in margins),(B,[m.enc() for m in margins])
  assert tc.b<t0.a
  # At the mean angle the bad-arc penalty is zero.
  assert s.a*s.a>R*R*K
  area_upper=(198-2*B)*3*I(3).sqrt()/4*R*R+B*R*R*(intercept+K*sineta)+B*H0
  deficit=PI-area_upper
  assert deficit.a>0,(B,deficit.enc())
  rows.append(dict(B=B,cos_eta_enclosure=coseta.enc(),small_arc_support_margins=[m.enc() for m in margins],area_deficit=deficit.enc()))
 return dict(status='PASS',radius=[R.numerator,R.denominator],k=[K.numerator,K.denominator],t1=[T1.numerator,T1.denominator],scope='Conditional geometric boundary lemma: m=100, no ears, B=33,34,35. These hypotheses need independent global proof at this radius.',pi_enclosure=PI.enc(),support_x=xs.enc(),support_intercept=intercept.enc(),penalty_concavity_numerator_upper=[penalty_second_numerator.numerator,penalty_second_numerator.denominator],cases=rows)
if __name__=='__main__':
 result=verify();p=Path(__file__).with_name('support_certificate.json');p.write_text(json.dumps(result,indent=2)+'\n')
 print('PASS: rational scalar support checks for r=9041/80000, k=6/25, B=33,34,35.')
 for row in result['cases']:
  m=Q(*row['area_deficit'][0]);print('B',row['B'],'certified deficit >',m)
 print('CONDITIONAL: m=100, B=33..35 and no ears must be certified separately at this radius.')
