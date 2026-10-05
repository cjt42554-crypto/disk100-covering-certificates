# Uniform-radius extension of the actual signed-anchor quadratic fan fee

This is a local geometric theorem. It extends the radius parameter of the
frozen `q3_actual_quadratic_support` proof; it does not establish a global
covering bound, a boundary support, eligible fan incidence, or an optimum.
The original sources and certificates are unchanged. No frozen mathematical
module or verifier was executed for this addendum.

## 1. Exact statement

Fix

```
lambda = 7/50,  t = 1071883/20000000,  D0 = 1/40,
Rmax = 949600/4987489.
```

For every real `R` with `0 < R <= Rmax`, let the three positive-area
triangles `(p,c_i,c_(i+1))`, `i=0,1,2`, be consecutive in one consistently
oriented fan, with their actual common apex `p` and circumradii at most
`R`. Choose the orientation so their actual apex increments `gamma_i`
and oriented areas are positive. For each base suppose that there is an
actual point `q_i` such that

```
|q_i| = 1,  |q_i-c_i| <= R,  |q_i-c_(i+1)| <= R.
```

Put

```
x_i = |c_(i+1)-c_i|/(2R),
A_i = -det(c_(i+1)-c_i,q_i-c_i)/2,
score_i = area_i/R^2 + A_i/R^2 - 4lambda*x_i^2 + t*gamma_i.
```

Define the exact scalar support, independently of `R`, by

```
G(theta) = sin(theta)+sin(2theta)-4lambda*sin(theta)^2+t*theta,
Bo(theta) = sin(theta)-4lambda*sin(theta)^2+t*(pi-theta),
C = max_{0 <= theta <= pi/2} max(G(theta),Bo(theta)),
Co = max_{0 <= theta <= pi/2} Bo(theta).
```

Then

```
sum_{i=0}^2 score_i < 3C - 1/40.                     (U1)
```

In particular the same fee holds throughout `0 < R <= 19/100`. The upper
endpoint `Rmax` is included. This is the largest radius endpoint supplied
by the *particular rational envelopes* in the frozen proof's equation
(18), not a sharp threshold for the local theorem.

The anchor is signed: an anchor on the apex side has negative `A_i`.
There is no hypothesis on center norms or on `p`, reflection, assigned
arc lengths, ownership rules, or distinctness of all anchors in the
initial domain. The ordinary-face restrictions of a whole-cover
substitution remain external.

## 2. Normalize the geometry without normalizing the unit circle away

For each nondegenerate triangle put `theta_i=asin(x_i)`. The sine rule,
using the circumradius bound, gives

```
0 < x_i <= 1,  theta_i <= gamma_i <= pi-theta_i.
```

If the two base angles are `(pi-gamma_i)/2 +/- delta_i`, the exact area is

```
area_i/R^2 = x_i^2*cot(gamma_i/2) - sigma_i,
sigma_i = 2*x_i^2*sin(delta_i)^2/sin(gamma_i) >= 0.
```

In base-midpoint coordinates divided by `R`, the endpoints are
`(-x_i,0),(x_i,0)` and the actual anchor is `(s_i,h_i)`, where positive
height points away from the apex. Its two coverage inequalities give

```
s_i^2+h_i^2+2x_i*|s_i| <= 1-x_i^2,
tau_i = x_i*cos(theta_i)-A_i/R^2
      = x_i*(cos(theta_i)-h_i) >= 0.
```

Thus every exact normalized score is

```
F(theta_i,gamma_i)-sigma_i-tau_i,
F(theta,gamma) = sin(theta)^2*cot(gamma/2)
             +sin(theta)*cos(theta)-4lambda*sin(theta)^2+t*gamma.
```

All these formulas and inequalities are independent of `R`. They require
only `R>0`. The unit-circle condition `|q_i|=1` is deliberately retained in
physical coordinates; it enters the argument in Section 4.

The second derivative of `F` in `gamma` is positive. Consequently its
full feasible-interval maximum is an endpoint value, `G` or `Bo`, and
every score is at most `C`. This treats initially obtuse or skew triangles
before the following localization.

## 3. The frozen paid-fit lemma has no radius hypothesis

Here are the ingredients of the paid-fit argument, including its exact
constants. They explain why the existing accepted scalar facts apply to
the radius-parameterized domain, rather than merely asserting a transfer
of fixed-radius certificate metadata.

The accepted scalar conjugate proves that both branches decrease on
`[pi/3,pi/2]`, are strictly concave on `[0,pi/3]`, and that

```
C-Co > 9/10,  theta* = argmax G in (7/10,24/25),
169/250 < cos(theta*) < 677/1000.
```

Its definition and certificate contain no radius. This concavity statement
is restricted to `[0,pi/3]`; no global concavity of `G` is used.

Suppose (U1) fails and write `d_i=C-score_i >= 0` and
`D=sum d_i <= D0`. Each score is at least `C-D0>Co`. The full convex
endpoint reduction forces `G(theta_i)>=score_i`. The accepted strict
guards `C-G(7/10)>D0` and `C-G(24/25)>D0`, together with the scalar
branch monotonicity, localize every theta to

```
L=7/10 < theta_i < H=24/25.
```

On this interval `-G''>4`. For `zcut=1/10`, the ratio
`sin(theta)/sin((theta+zcut)/2)` has no interior minimum and both endpoint
ratios exceed `81/50`. Indeed its derivative numerator has derivative
`-(3/4)sin(theta)sin((theta+zcut)/2)<0`. Hence the loss slope
`-F_gamma>(1/2)(81/50)^2-t>1` up to the cut. The cut is within the
feasible interval because `2H+zcut<pi`. Its value and the obtuse endpoint
are below `C-D0`; convexity excludes the entire remaining gamma interval.
With `z_i=gamma_i-theta_i`, integration gives the one simultaneous budget

```
2 sum(theta_i-theta*)^2 + sum z_i + sum sigma_i + sum tau_i
                                                    <= D <= 1/40. (U2)
```

In particular `gamma_i < H+D0 < pi/2`. None of this contains `R`.

The exact normalized outer lens tip `a_i` has bisector-phase expression

```
(a_i-p)/R = exp(i*m_i)[U_i*cos(delta_i)+i*V_i*sin(delta_i)],
rho_i = sin(theta_i)/sin(gamma_i),
U_i = rho_i*(1+cos(gamma_i))+cos(theta_i),
V_i = rho_i*(1-cos(gamma_i))-cos(theta_i).
```

The actual anchor-to-tip error satisfies
`|q_i-a_i|/R <= tau_i/x_i^2`, including either sign of `h_i`.
The actual shared fan has the exact lifted phase relations

```
m_0-m_1 = -(gamma_0+gamma_1)/2,
m_2-m_1 =  (gamma_1+gamma_2)/2.
```

The localized rational guards are

```
S0=16/25 < sin(theta_i) < 41/50=S1,
C1=573/1000 < cos(theta_i) < 153/200=C0,
Qmin=107/50 < Q_i=1+2cos(theta_i) < 127/50=Qmax.
```

They imply `|U_i-Q_i| <= (11/4)z_i`, `|V_i|<11/20`, and the squared
skew correction is at most `(79/100)^2 sigma_i`, using the accepted
strict coefficient guard. At `sigma_i=0` the displacement is zero, so
no strict product inequality is asserted. The anchor coefficient is
`ct=625/256=1/S0^2`. All corrections here are in units of `R`.

After rotating so `m_1=0`, define the two actual translation-free vectors

```
B = ((q_0+q_2)/2-q_1)/R,
b = (q_2-q_0)/(2R).
```

The corresponding acute maps on the entire box `[L,H]^3` use
`a=(theta_0+theta_1)/2`, `beta=(theta_1+theta_2)/2` and are

```
Bacute = (Q_0 exp(-i*a)+Q_2 exp(i*beta))/2-Q_1,
bacute = (Q_2 exp(i*beta)-Q_0 exp(-i*a))/2.
```

The endpoint imaginary Jacobian bound retains both possible signs. The
fixed rational component bounds and their squared Frobenius consequences
are

```
jt=S1*C0+Qmax*S1/4,
ju=max(S1^2-Qmin*C1/4,Qmax*C0/4-S0^2),
jv=2S1-Qmin*S0/2,      jw=(Qmax*C0-Qmin*C1)/4,
jvb=(Qmax*S1-Qmin*S0)/4, jwb=Qmax*C0/2,
JB2=2(jt^2+ju^2)+jv^2+jw^2=76960399529/20000000000,
Jb2=2(jt^2+ju^2)+jvb^2+jwb^2=77583868411/20000000000.
```

The straight segment from the ideal peak triple to the actual triple
stays in this box. Jacobian integration and Cauchy with (U2), including
the skew weights `(1/2,1,1/2)` and `(1/2,0,1/2)`, give

```
alphaB2=JB2/2+(3/2)*(79/100)^2=114406399529/40000000000,
alphab2=Jb2/2+(1/2)*(79/100)^2=90065868411/40000000000,
LB=201/50,  Lb=529/200,  W=sum z_i+sum tau_i,
|B-B*| <= sqrt(alphaB2)*sqrt(D-W)+LB*W,
|b-b*| <= sqrt(alphab2)*sqrt(D-W)+Lb*W.
```

The accepted rational guards `4LB^2 D0<alphaB2` and
`4Lb^2 D0<alphab2` make these maxima occur at `W=0`; `D=0` is immediate.
The strict squared error guards imply

```
|B-B*| < EB=27/100,  |b-b*| < Eb=6/25.               (U3)
```

The common apex cancels from these actual vectors. At the peak, writing
`y=cos(theta*)`, their ideal values are

```
B*=(-A,0),  b*=(0,J),
A=(1+2y)(1-y),  J=(1+2y)*sqrt(1-y^2),
Amin=19/25 < A < 763/1000=Amax,
Jmin=173/100 < J < 87/50=Jmax.                      (U4)
```

Every scalar and normalized geometric estimate in this section is the
same accepted proof with arbitrary positive `R`; there is no change of
coefficients, no separate budget for any of its four losses, and no
charge to a BC lens term. The fixed-radius author's final `.38` margin
is not used in the extension.

## 4. Unit anchors supply the only radius-dependent restriction

For completeness, the necessary three-unit-point condition is proved
here. Put `a=(q_0+q_2)/2` and physical half-chord `h=(q_2-q_0)/2`.
The two endpoint norm equations give

```
a dot h=0,  |a|^2+|h|^2=1.
```

By (U3)--(U4), `|b|>Jmin-Eb=149/100>0`; hence `h=Rb` is nonzero.
Choose a unit normal `e` to it, write `a=alpha e` and
`q_1 dot e=beta`, where `alpha,beta` lie in `[-1,1]`. The elementary
identity

```
(1-alpha*beta)^2-(alpha-beta)^2
                             =(1-alpha^2)(1-beta^2)>=0
```

and `1-alpha*beta>=0` imply
`|alpha-beta|<=1-alpha*beta`. Also
`|a-q_1|^2+|h|^2=2(1-alpha*beta)`. Since `a-q_1=RB` and `h=Rb`,
division by positive `R` proves

```
2 |B_perp_to_b| <= R (|B|^2+|b|^2).                 (U5)
```

This is where the physical unit-circle condition enters. From (U3), the
horizontal component of `b` has magnitude less than `Eb`, while
`|b|>Jmin-Eb`. The fixed square guard

```
1-[Eb/(Jmin-Eb)]^2 > (49/50)^2
```

therefore makes the horizontal axis project onto the normal to `b` by
more than `49/50`. Applying the reverse triangle inequality to `B*`
and (U3) gives the strict lower bound

```
2 |B_perp_to_b| > K=2[(19/25)(49/50)-27/100]=1187/1250.
```

The strict upper norm bounds from (U3)--(U4) give

```
|B|^2+|b|^2 < S=(763/1000+27/100)^2+(87/50+6/25)^2
               =(1033/1000)^2+(99/50)^2
               =4987489/1000000.
```

For `0<R<=K/S=949600/4987489`, the hypothesized high-score triple would
therefore have to satisfy

```
K < 2|B_perp_to_b| <= R(|B|^2+|b|^2) < RS <= K,
```

a contradiction. Both fit inequalities remain strict when `D=D0`, so
the conclusion is strict and the radius endpoint is closed even though
`K-Rmax*S=0`. No positive numerical margin at that endpoint is claimed.

The requested simpler interval follows from the exact comparison

```
K-(19/100)*S = 197709/100000000 > 0,
Rmax-19/100 = 197709/498748900 > 0.
```

At larger radii this particular pair of coarse envelopes does not
contradict (U5); that failure says nothing about validity of the fee there.

## 5. Degeneracies and the meaning of a uniform radius

The positive-area argument excludes `theta=pi/2` and the zero-base limit
before any localized division: their scores are at most `Co`, below the
high-score threshold. In the fit, `x_i>S0`, `sin(gamma_i)>0`, and
`|b|>149/100`, so all divisions are justified. Repeated anchors were
allowed initially; `q_0=q_2` would contradict the derived lower bound on
`|b|`, while other repetitions are covered by (U5).

For every fixed positive `R`, a zero-area limiting face, including a
collapsed side different from the base, has

```
limsup score <= 1/2+t*pi
             < 1/2+t*(22/7)=46790713/70000000 < 7/10.
```

Indeed coverage always gives `A_i/R^2<=x_i sqrt(1-x_i^2)<=1/2`, the
quadratic term is nonpositive, and limiting apex angles are in `[0,pi]`.
The accepted `C>147/100` gives `C-D0>289/200>7/10`. The other two scores
are at most `C`, so a high-score triple cannot contain such a zero-area
face. A limiting triple with all positive areas has defined angles and
satisfies the above proof. This closure bound also has no radius
restriction apart from positivity.

The theorem quantifies over each positive radius using that same radius
in the circumradius bound, endpoint coverage, base normalization, area
normalization and unit condition. It is not a claim that the score stays
constant when a fixed configuration is relabeled with another radius.
`R=0` is excluded; no value is assigned to expressions divided by zero.
Physical zero area in a sequence with `R` tending to zero is not by itself
zero *normalized* area, so that circumstance must not be misclassified
as the fixed-positive-radius degeneracy argument above.

## 6. Dependency audit and exact computational scope

The frozen checker declares `R=11321/100000`. Static inspection of its
entire body shows that this variable is used only in
`rhs=R*((Amax+EB)^2+(Jmax+Eb)^2)` and radius metadata. Its final test
`margin>19/50` is tied to that original radius; the extension replaces
only that stronger reporting guard by the strict geometric contradiction
`0<R<=K/S`. Every earlier localization, ratio, curvature, sine/cosine,
tip, two-sign Jacobian, paid-budget, ideal-size and axis guard is radius
independent at the stated fixed `lambda,t`.

The accepted conjugate checker and its certificate contain no radius
parameter. Its proof first states the support at arbitrary positive `R`
and then removes it by the normalized variables. The imported primitive
file declares its own unrelated `R=9041/80000` for its separate boundary
`verify()` routine. Inspection of the selected definitions `I`, `PI`,
`sinq`, and `cosq`, and their helper dependencies, shows that they do not
reference that radius: only rational arithmetic, fixed precision,
factorials, integer square roots and rational Taylor/Machin constants
are involved. Neither that boundary routine nor the conjugate checker
was called here.

The original author certificate and independent full-source/fresh-replay
reviews remain certificates of their declared fixed-radius computation.
They are used here as frozen evidence for the unchanged normalized fit,
not rewritten as range certificates. `SOURCE_PINS.json` identifies the
exact original materials read. The separate new
`radius_endpoint_fraction.py` uses only the standard-library `Fraction`
and exact hard-coded rational envelopes to calculate `S`, `K`, `Rmax`
and the `19/100` comparison; it performs no trigonometry, imported source
execution, repeated guard campaign, optimizer, scan, or geometry search.

Whole-cover use still requires the same signed area identity, matching
quadratic and angle coefficients, actual eligible ordinary faces and
anchors, and disjoint charged face/anchor sets. This addendum extends no
external boundary constant, profile row, topology or nongeneric
whole-cover reduction. The local fee remains `1/40`; in an already valid
whole-cover identity its physical area payment would be `R^2/40`.
