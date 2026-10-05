# A paid mean-norm envelope and the exact remaining tau interface

This is a new analytic addendum, separate from the frozen zero-tau
proof and `WORKING_INTERFACE.md`. It proves the shape-envelope
calculation below. The final positive-tau conclusion is **conditional**
on a uniform directional derivative bound that is not proved here.
No mathematical program, imported source, checker or parameter scan
was executed to produce this addendum.

Use the definitions and fixed value `R=5661/50000` of the frozen
zero-tau proof. Write `W=tau0+tau1`, `B=D-W`, and suppose
`D<=223/10000`. The accepted actual per-face deficit budget gives

    4E^2+Delta^2+2 sum x_i sin(d_i)^2 <= B.       (1)

Removing paid apex slack in the simultaneous manner of Section 3
of the frozen proof gives the upper norm-driving angle

    chi0(s,eta)=s/2+eta+asin(R sin(3s/2-eta)).

The original shape ranges remain `.75<s<.904`, `|Delta|<.15`.
The sum identity and `x_i>.64` also give
`2.56 sin^2 eta<=D`, hence `|eta|<.095`: indeed
`sin(.095)>.094` and `2.56(.094)^2>.0223`.
Negative eta may be replaced by zero for an upper norm bound and
a lower payment, because `chi0` is increasing in eta on this range.

## 1. Weighted skew minimization loses only quadratically in Delta

Put `delta=Delta/2`, `X=x0+x1`, `Y=x0-x1`.
For every common skew `m`, without any symmetry assumption,

    2 sum x_i sin(d_i)^2
      = X-[X cos(2eta) cos(2m)-Y sin(2eta) sin(2m)]
      >= X-sqrt(X^2 cos^2(2eta)+Y^2 sin^2(2eta))
      >= (8x0x1/X) sin^2 eta cos^2 eta.          (2)

The last inequality uses `sqrt(1-u)<=1-u/2` for `0<=u<=1`,
with `u=4x0x1 sin^2(2eta)/X^2`.
The exact weight is

    8x0x1/X
       =4 sin s cos delta [1-cot^2 s tan^2 delta].

Because `|delta|<.075`, `cos delta>.997` and
`|tan delta|<.076`. Also `cot s<1.1`: use
`sin(.75)>.68` and the fourth-order upper polynomial
`cos(.75)<.74`. Consequently

    cos delta [1-cot^2 s tan^2 delta]
      >(.997)[1-(1.1)^2(.076)^2]>.99.

Moreover `cos^2 eta>1-(.095)^2>.99`.
Thus (2), with `sin s>.68`, implies the convenient non-strict
single-budget bound

    4E^2+Delta^2+(1666/625) sin^2 eta <= B,      (3)

where `1666/625=2.6656=(4)(.68)(.98)`.
This avoids the linear `|Delta|` coefficient of the simpler
minimum-norm relaxation in the frozen proof.

## 2. Mean norm uses the same angular payment

Set `eta+=max(eta,0)`. For `0<=eta+<.095`,
`eta+<=1.002 sin eta+`, by the lower sine polynomial.
On the connecting rectangle from `(theta*,0)` to `(s,eta+)`,
the frozen derivative bounds remain

    chi0_s<.6,  0<chi0_eta<1,
    chi0(theta*,0)<.523.

Cauchy applied only to (3) therefore gives

    chi0(s,eta+)
      <=.523+sqrt[(.6^2/4+1.002^2/(1666/625))(B-Delta^2)]
      < .523+sqrt[(47/100)(B-Delta^2)].          (4)

For the last rational guard, `.09+1.004004/2.6656<.47` follows
from `(.38)(2.6656)=1.012928>1.004004`.

At the positive unit-tip branch, the exact shape mean is

    M(0)=1+R^2-2R cos chi cos Delta,

and actual paid apex slack gives `chi<=chi0(s,eta+)`.
Since `chi0>.375` and `chi0<.662`, and
`t0=(3/2)Delta^2<.034`,

    cos chi0 cos Delta >= cos(chi0+(3/2)Delta^2). (5)

To check (5), if Delta is nonzero then
`sin t0>.999 t0`, `tan chi0>chi0>.375`, and

    cos chi0-cos(chi0+t0)
       >=sin chi0 sin t0
       >cos chi0(.375)(.999)(3/2)Delta^2
       >cos chi0 Delta^2/2
       >=cos chi0(1-cos Delta).

At Delta zero (5) is equality; no strict zero-product inequality
is required.

The function
`sqrt[(47/100)(B-v)]+(3/2)v`, `0<=v<=B<=.0223`, is maximized
at `v=0`. Indeed the square-root loss per v is at least
`(47/100)/(2sqrt[(47/100)B])>(3/2)`, since
`sqrt[(47/100)B]<.103<(47/100)/3`.
Combining this fact, (4), and (5) proves

    M(0) <= N(B),
    N(B)=1+R^2-2R cos(.523+sqrt[(47/100)B]),     (6)

with endpoint `N(0)=1+R^2-2R cos(.523)`.
This is a mean-square bound, not a minimum-norm bound, and the
Delta payment has not been spent a second time.

## 3. A coefficient of .30 suffices for all positive tau

For `0<B<=.0223`, put `y=sqrt[(47/100)B]<.103`.
Then

    N'(B)=R(47/100) sin(.523+y)/y > .30.        (7)

The function `sin(.523+y)/y` decreases on `0<y<=.103`:
its numerator derivative is `y cos(.523+y)-sin(.523+y)<0`,
because `y<.103` while `sin(.523)>.49`.
Thus its value is at least `sin(.626)/.103`.
The cubic lower sine polynomial gives `sin(.626)>.585`, and

    R(.47)(.585)=.031129839 >(.30)(.103)=.0309.

At the maximum budget, `.523+sqrt[(.47)(.0223)]<.626`.
The sixth-order lower cosine polynomial gives `cos(.626)>.81`,
hence

    N(.0223)<1+R^2-2R(.81)=.8294023684.

The actual full-arc condition `.147` requires each outer norm
greater than `.911088875`, whose square exceeds `.829921`.

Consequently, if a correct whole-domain bound along the exact
fixed-fan actual-anchor interpolation were to prove

    M'(t)<=.30 W  for every 0<=t<=1,            (8)

then (6)-(7) would imply

    M(1)<=N(D-W)+.30W<=N(D)<=N(.0223)

and exclude both strong outer arcs at every actual tau.
The limiting case `D-W=0` follows by continuity of N;
W zero needs no derivative assertion. A coefficient `.22` would
be sufficient but is not necessary.

**Condition (8) is the remaining unproved step.** This note neither
asserts (8), nor a general positive-tau fee, nor a global lower
radius bound. The actual fan path and its branch proof have their
own working interface and still need independent checking.

## 4. Useful exact directional cancellation, pending sign and bound

In the phase convention of the frozen zero-tau proof, choose
directed base tangents `e_i` and normals `n_i` toward the apex.
The exact phase formulas give

    n0-cot(theta0)e0=-i b0/x0,
    n1+cot(theta1)e1=+i b1/x1.

Thus if the path gradient signs satisfy `g0 dot e0<=0` and
`g1 dot e1>=0`, the lens derivative coefficients reduce exactly to

    ((-i b0) dot g0)/x0^2,
    ((+i b1) dot g1)/x1^2.                     (9)

The required signs must be proved over the whole path before
using (9) as the selected support branches. They are not assumed
as part of (6).
At the unit-tip endpoint, algebraic substitution gives the first
coefficient as

    R(u-R cos alpha)/x0^2
       [sin alpha(1-R cos h cos Delta/u)+cos h sin Delta],

and the second has x1 in its denominator and minus the final
`cos h sin Delta` term. These identities exhibit the cancellation
between the tangent allowance and the corresponding base half-sine.
They do not yet give a uniform bound for (8).

## 5. Execution boundary

Author mathematical program executions: 0. Mathematical source
imports/replays/old checkers/optimization/search: 0. Writing and
hashing this separately versioned paper note are filesystem
bookkeeping only. Earlier independent zero-tau review is not a
review of this new addendum or its conditional positive-tau bridge.
