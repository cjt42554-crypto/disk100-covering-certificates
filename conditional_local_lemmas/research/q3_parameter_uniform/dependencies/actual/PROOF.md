# A same-quadratic actual signed-anchor three-face loss of 1/40

Fix

```
R=11321/100000, lambda=7/50, t=1071883/20000000.
```

Let three positive-area triangles `(p,c_i,c_(i+1))`, `i=0,1,2`, be
consecutive in one consistently oriented fan with their actual common
apex `p`. Each triangle has circumradius at most `R`. Choose the orientation
so their successive apex increments `gamma_i` are positive and their
oriented areas are positive. For each base there is an actual point `q_i`
with

```
|q_i|=1, |q_i-c_i|<=R, |q_i-c_(i+1)|<=R.
```

Define its signed opposite-side anchor and normalized half-base by

```
A_i=-det(c_(i+1)-c_i,q_i-c_i)/2,
x_i=|c_(i+1)-c_i|/(2R).
```

An anchor on the apex side has negative `A_i`; that case is included.
Let `C` be the full-apex independent support for the **same** `lambda,t`
proved in `../quadratic_lambda14_conjugate/PROOF.md`. Then

```
sum_i [area_i/R^2+A_i/R^2-4lambda*x_i^2+t*gamma_i]
                 < 3C-1/40.                            (1)
```

Here `C` is approximately `1.473500630352807`; the actual loss `1/40=.025`
exceeds `.019931551`. The conclusion covers an expanded common-unitpoint
domain: no center norm bound, reflected configuration, assigned-arc
length, canonical tie, nearest-site rule, or bound on `p` is assumed.
It is a local score theorem. Its global substitution requires the
separate signed area identity and eligible disjoint fan incidence.

## 1. Full angle/skew/anchor domain, before any localization

For one triangle put `theta=asin x`, `y_x=cos theta` and let its two
base angles be `(pi-gamma)/2 +/- d`. The sine rule gives

```
theta<=gamma<=pi-theta,
area/R^2=x^2*cot(gamma/2)-sigma,
sigma=2*x^2*sin(d)^2/sin(gamma)>=0,
|d|<(pi-gamma)/2<pi/2.                              (2)
```

Indeed the circumradius is `Rx/sin gamma`, while the area equals
`2R^2 x^2 sin(baseangle_1)sin(baseangle_2)/sin gamma`. In base-midpoint
coordinates divided by `R`, write the endpoints as `(-x,0),(x,0)` and
the anchor as `(s,h)`, with positive height away from the apex. The two
disk inequalities are equivalent to

```
s^2+h^2+2x|s|<=y_x^2.
tau=x*y_x-A/R^2=x*(y_x-h)>=0.                      (3)
```

Thus its exact actual score is `F(x,gamma)-sigma-tau`, where

```
F(x,gamma)=x^2*cot(gamma/2)+x*sqrt(1-x^2)
                         -4lambda*x^2+t*gamma.
F_gamma_gamma=(x^2/2)*csc(gamma/2)^2*cot(gamma/2)>0.
```

The full convex endpoint reduction gives

```
G(theta)=sin theta+sin(2theta)-4lambda*sin(theta)^2+t*theta,
Bo(theta)=sin theta-4lambda*sin(theta)^2+t*(pi-theta).
```

The accepted source proves that both branch maxima occur on `[0,pi/3]`,
where each branch is strictly concave, and that each branch decreases
on `[pi/3,pi/2]`. Its exact certificate gives

```
C=max G, Co=max Bo, C-Co>9/10,
theta*=argmax G in (.7,.96), .676<cos(theta*)<.677.   (4)
```

This is the new full-apex conjugate at `lambda=.14`; the old implementation
restricted to `lambda<=1/8` is not used. In particular, this proof does
not assert that `G` is concave on all of `(0,pi/2)`. The conjugate source,
proof, and saved certificate are pinned and read without rerunning its
checker.

For `theta<pi/2`, `F` lies below the convex interpolation between `G`
and `Bo`. Hence every actual score is at most `C`. The endpoint
`theta=pi/2` has score at most `Co`, and the `x=0` limiting scores are
at most `t*pi<=Co`. Degenerate closures therefore cannot occur in the
high-score contradiction below. All obtuse/skew cases are included
before localization; acuteness is a consequence of the paid deficit.

## 2. One deficit pays angle, apex slack, skew and actual-anchor loss

Suppose (1) fails. Let `score_i` denote its summands, and put

```
d_i=C-score_i>=0, D=sum d_i<=D0=1/40,
L=7/10, H=24/25, zcut=1/10.
```

Since `score_i>=C-D0>Co`, convex endpoint interpolation forces
`G(theta_i)>=score_i` and therefore `C-G(theta_i)<=d_i`. The exact
guards `C-G(L)>D0`, `C-G(H)>D0`, the unique peak, and the global branch
monotonicity from Section 1 imply `L<theta_i<H`. On this window,

```
-G''=sin theta+4sin(2theta)+8lambda*cos(2theta)
 >=sin L+4min(sin(2L),sin(2H))+8lambda*cos(2H)>4.    (5)
```

This uses sine concavity on `[2L,2H]` and cosine monotonicity there.

For `gamma` between `theta` and `theta+zcut`,

```
-F_gamma >= (1/2)*r(theta)^2-t,
r(theta)=sin theta/sin((theta+zcut)/2).
```

The sign of `r'` is the sign of
`N=cos(theta)sin(u)-(1/2)sin(theta)cos(u)`, with
`u=(theta+zcut)/2`. Differentiation gives
`N'=-(3/4)sin(theta)sin(u)<0`. Thus `r` has no interior minimum.
Exact endpoint guards give `r(L),r(H)>81/50`, and consequently

```
-F_gamma > (1/2)*(81/50)^2-t >1.                   (6)
```

The cut lies inside the allowed interval, since `2H+zcut<pi`.
At the cut `F<=G-zcut<=C-zcut<C-D0`. On the whole remaining interval
`[theta+zcut,pi-theta]`, convexity bounds `F` by the larger of its cut
value and its obtuse endpoint, both strictly below `C-D0`. This
contradicts the high actual score if `gamma` reaches that cut. Thus
`z_i=gamma_i-theta_i<zcut`. Integrating (5) from `theta*` and (6)
over each actual apex interval yields the **single simultaneous budget**

```
2*sum e_i^2+sum z_i+sum sigma_i+sum tau_i <=D<=1/40,
e_i=theta_i-theta*.                                (7)
```

In particular each actual `gamma_i<=theta_i+D0<.985<pi/2`. No
additional independent loss is charged for a variable already in (7).

## 3. Exact tips and linear actual-anchor error, for either signed side

Let `a_i` be the outer lens tip at its base, namely `(0,y_x)` in (3).
Write `Delta=y_x-h=tau/x`. Equation (3) implies

```
2x|s|<=2y_x*Delta-Delta^2<=2y_x*Delta,
|q_i-a_i|/R<=Delta*sqrt(1+y_x^2/x^2)=tau_i/x_i^2.    (8)
```

This is linear in the actual signed-anchor slack; negative anchor
height simply enlarges `tau`. There is no absolute-anchor replacement.

Let `m_i` be the bisector direction of the two actual rays from `p`
and `r_i=x_i/sin gamma_i` its normalized circumradius. The outer tip
has the exact complex expression

```
(a_i-p)/R=e^(im_i)[U_i*cos d_i+i*V_i*sin d_i],
U_i=r_i(1+cos gamma_i)+cos theta_i,
V_i=r_i(1-cos gamma_i)-cos theta_i.                 (9)
```

For completeness, the sine rule gives the normalized base midpoint
relative to its ray bisector as
`r_i[(1+cos gamma_i)cos d_i+i(1-cos gamma_i)sin d_i]`.
Its actual outward normal is `e^(im_i-id_i)`, so adding the outer-tip
normal contribution `cos(theta_i)e^(im_i-id_i)` proves (9). Reversing
the skew sign convention conjugates the skew term and changes no norm
estimate. The shared rays give the exact lifted phase identities

```
m_0-m_1=-(gamma_0+gamma_1)/2,
m_2-m_1= (gamma_1+gamma_2)/2.                       (10)
```

The following rational bounds hold on the localized theta window:

```
S0=.64<sin theta<S1=.82,
C1=.573<cos theta<C0=.765,
cot theta<6/5, cot(theta/2)<11/4,
Qmin=2.14< Q_i=1+2cos theta_i <Qmax=2.54.
```

Since `0<=1-r_i<=cot(theta_i)z_i` and
`U_i=sin theta_i*cot(gamma_i/2)+cos theta_i`,

```
0<U_i<=Q_i,
|U_i-Q_i|<=(11/4)z_i.                             (11)
```

The derivative bound in (11) follows from
`|U_gamma|=(sin theta/2)csc(gamma/2)^2<=cot(theta/2)`.
The exact identity

```
V_i=(1-2cos theta_i)-(1-r_i)(1-cos theta_i)
                       +r_i(cos theta_i-cos gamma_i)
```

and `z_i<=D0` give

```
Vlo=1-2C0-(6/5)(1-C1)D0 >-.55,
Vhi=1-2C1+D0 <.55.
```

Set `cs=625/512=1/(2S0^2)`, `ct=625/256=1/S0^2`. By (2),
`sin(d_i)^2<=cs*sigma_i`, and `1-cos d_i<=sin(d_i)^2`. Thus the
norm squared of the complex skew correction from `U_i` to the tip is
at most

```
U_i^2 sin(d_i)^4+V_i^2 sin(d_i)^2
 <=cs[.55^2+Qmax^2*cs*D0]*sigma_i
 <.79^2*sigma_i.                                  (12)
```

The additional anchor error from (8) is at most `ct*tau_i`.

## 4. A full-box acute-tip Jacobian, with both endpoint imaginary signs

Rotate coordinates so the actual middle bisector `m_1` is zero.
For three independent theta values in the box `[L,H]^3`, define

```
a=(theta_0+theta_1)/2, b=(theta_1+theta_2)/2,
Bacute=(Q_0 e^(-ia)+Q_2 e^(ib))/2-Q_1,
bacute=(Q_2 e^(ib)-Q_0 e^(-ia))/2.
```

The endpoint derivative columns for `Bacute` are
`e^(-ia)[-sin theta_0-iQ_0/4]` and
`e^(ib)[-sin theta_2+iQ_2/4]`. Its middle column is
`2sin theta_1+(i/4)(Q_2e^(ib)-Q_0e^(-ia))`.
The middle column for `bacute` is
`(i/4)(Q_2e^(ib)+Q_0e^(-ia))`; its endpoint columns have the same
absolute component bounds as those of `Bacute`.

Use these exact rational component bounds:

```
jt=S1*C0+Qmax*S1/4,
ju=max(S1^2-Qmin*C1/4, Qmax*C0/4-S0^2),
jv=2S1-Qmin*S0/2,
jw=(Qmax*C0-Qmin*C1)/4,
jvb=(Qmax*S1-Qmin*S0)/4,
jwb=Qmax*C0/2.
```

The endpoint real magnitude is at most `jt`. Its imaginary component
is `sin(theta_0)sin(a)-Q_0cos(a)/4`, or its sign reversal. Its upper
bound is the first entry of `ju` and its negative lower magnitude is
the second entry. Both are included: the imaginary derivative can
change sign in this new box, so the older positive-sign simplification
is not used. Exact guards show the first entry dominates the second.

The middle real component of `Bacute` is positive because
`2S0-Qmax*S1/2>0`, and its upper bound is `jv`. The remaining middle
component bounds follow by taking differences or sums of two positive
products in their rational intervals. Hence the whole-box squared
Frobenius bounds are

```
JB2=2(jt^2+ju^2)+jv^2+jw^2,
Jb2=2(jt^2+ju^2)+jvb^2+jwb^2.                     (13)
```

The straight segment from `(theta*,theta*,theta*)` to the actual
theta triple stays in the box. Integrating its Jacobian and using the
operator norm bound by the Frobenius norm gives changes at most
`sqrt(JB2)*sqrt(sum e_i^2)` and `sqrt(Jb2)*sqrt(sum e_i^2)`.
This is a complete localized-box map bound, not an optimizer Hessian
approximation or a claimed pointwise global quadratic inequality.

## 5. Shared-budget optimization of all translation-free errors

The actual combinations are

```
B=((q_0+q_2)/2-q_1)/R,
b=(q_2-q_0)/(2R).
```

Their ideal values at `theta_i=theta*` are
`B*=(-A,0)`, `b*=(0,J)`, with

```
y=cos theta*, Q=1+2y,
A=Q(1-y), J=Q*sqrt(1-y^2),
.76<A<.763, 1.73<J<1.74.                          (14)
```

The common translation `p` cancels exactly. Comparing actual tips
to acute tips at the same actual phases costs at most `(11/4)z_i`
by (11); changing the outer phases costs at most `Qmax` times their
phase errors `(z_0+z_1)/2`, `(z_1+z_2)/2`. After the half-endpoint
weights in `B,b`, safe aggregated z bounds are

```
LB*sum z, LB=11/4+Qmax/2=201/50=4.02,
Lb*sum z, Lb=(11/4+Qmax)/2=529/200=2.645.
```

The corresponding anchor errors are at most `ct*sum tau` and
`(ct/2)*sum tau`; these coefficients are smaller than `LB,Lb`.
The skew weights are `(1/2,1,1/2)` for `B` and `(1/2,0,1/2)` for
`b`. Their squared sums are `3/2` and `1/2`. Cauchy with the angular
payment coefficient two in (7) gives

```
alphaB2=JB2/2+(3/2)*.79^2,
alphab2=Jb2/2+(1/2)*.79^2,
W=sum z+sum tau,
|B-B*|<=sqrt(alphaB2)*sqrt(D-W)+LB*W,
|b-b*|<=sqrt(alphab2)*sqrt(D-W)+Lb*W.               (15)
```

Each term uses exactly its part of (7), once. The exact guards give
`4LB^2D0<alphaB2` and `4Lb^2D0<alphab2`. The derivative of
`alpha*sqrt(D-W)+Llinear*W` is therefore negative throughout
`0<=W<D`; equivalently use the square-root tangent at `W=0`.
Both maxima in (15) occur at zero linear spending. Consequently

```
|B-B*|<=sqrt(alphaB2*D0)<.27,
|b-b*|<=sqrt(alphab2*D0)<.24.                       (16)
```

This optimizes errors against one paid deficit rather than adding
four independent worst-case errors. No BC lens loss has been used.

## 6. The three actual unit anchors contradict that paid fit

For any three unit points their actual `B,b` obey

```
2|B_perp_to_b|<=R(|B|^2+|b|^2).                   (17)
```

To prove this without an external curvature theorem, let `a=(q_0+q_2)/2`
and physical half-chord `h=(q_2-q_0)/2`. Unit endpoint equations give
`a dot h=0` and `|a|^2+|h|^2=1`. For `h!=0`, choose a unit vector `e`
normal to it, and write `a=alpha e`, `q_1 dot e=beta`. Both scalars
lie in `[-1,1]`. The identity

```
(1-alpha*beta)^2-(alpha-beta)^2
                        =(1-alpha^2)(1-beta^2)>=0
```

proves `|alpha-beta|<=1-alpha*beta`. Also
`|a-q_1|^2+|h|^2=2(1-alpha*beta)`. Since `a-q_1=RB`, `h=Rb`,
this is (17). Our estimate (16) ensures `b!=0`; repeated anchors
elsewhere in the original domain have not been excluded by fiat.

From (14)--(16), `|b_x|/|b|<.24/(1.73-.24)`. The exact square guard
`1-[.24/(1.73-.24)]^2>.98^2` implies that the ideal horizontal axis
projects onto the actual `b`-perpendicular by more than `.98`. Hence

```
2|B_perp_to_b|>2(.76*.98-.27)=1187/1250=.9496,
R(|B|^2+|b|^2)<R[(.763+.27)^2+(1.74+.24)^2].       (18)
```

The left lower bound minus the right upper bound is strictly greater
than `19/50=.38`, verified by rational arithmetic. Equations (17) and
(18) contradict each other. This excludes the closed endpoint `D=D0`
as well, proving (1). Limiting degenerate faces cannot rescue a
high-score triple: Section 1 already excludes them before divisions
or the localized fit, and all subsequent closed error bounds are strict
in the final contradiction.

## 7. Exact verification, scope and the once-only global interface

`verify_actual_quadratic.py` checks all new localization, trig ratio,
curvature, both-sign Jacobian, skew, shared-budget, ideal-size and final
unit-compatibility guards with outward rational intervals. It reads
only the already accepted full-apex conjugate certificate and imports
only the pinned interval/Taylor primitive definitions. All imported
source bytes are checked before and after. Its single author run
passes; no old checker batch, optimizer or geometric search is run.

The score in (1) includes each marked triangle and its actual signed
anchor once, with exactly the latest quadratic coefficient and apex
tilt. In the signed shoelace identity

```
center_polygon_area+sum actual_signed_anchor+sum actual_BC_lens=pi,
```

one may replace these three marked/actual-anchor score terms by (1).
The BC lens term is retained once with its existing support and loss;
the slack `tau` here belongs to the old **edge anchor** allowance, not
that lens term. There is no conversion from the false standalone J
fee, no auxiliary division by an anchor-to-lens bridge constant, and
no separately added raw triangle fee. Several local payments require
disjoint charged marked-face/edge-anchor sets. Other unmarked angle
credits, BI credits and centered credits also retain their original
single allocation.

This theorem neither verifies the complete parent global budget nor
all profiles, nonordinary cases or topology. The applicable signed
shoelace reduction, disjoint fan incidence, coefficient matching and
remaining global rows are separate obligations. It does not give a
100-circle full-unit-disk construction, global radius exclusion, or
the exact optimal polynomial. The root-owned signed budget interface
in `../q3_actual_anchor_support/GLOBAL_BUDGET_INTERFACE.md` states the
same once-only principle for the earlier linear objective; the present
theorem changes the objective explicitly rather than transferring that
older numerical fee.
