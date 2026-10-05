# Uniform-radius extension of the actual signed-anchor fee 19/500

This is a local geometric parameter extension of the separately frozen
and independently accepted `q3_stronger_fee_19_500` theorem. Its complete
independent analytic review and one fresh replay are recorded in
`independent_audit/INDEPENDENT_REVIEW.md`, with SHA-256
`9a7b5e434a3b4a7f6989fa8d054afe9bd2554f20d98041cdec7a3ce01480aee1`.
The exact accepted normalized facts remain the premise of this analytic
extension. Their acceptance is separately bound in `PREMISE_STATUS.json`;
no fixed-radius author's certificate is relabeled as a range certificate.

No stronger or earlier verifier, frozen mathematical module, trigonometric
guard, optimizer, scan, or geometry search was executed for this addendum.
The earlier frozen `q3_radius_uniform` directory is preserved unchanged.

## 1. Exact parameterized local statement

Fix

```
lambda=7/50,  t=1071883/20000000,  D0=19/500,
Rmax=794400/5378209.
```

For any real `0<R<=Rmax`, take three positive-area triangles
`(p,c_i,c_(i+1))`, `i=0,1,2`, consecutive in one consistently oriented
fan with their actual common apex `p`. Their actual apex increments
`gamma_i` and oriented areas are positive, and each circumradius is at
most `R`. On each base take an actual anchor `q_i` satisfying

```
|q_i|=1,  |q_i-c_i|<=R,  |q_i-c_(i+1)|<=R.
```

Define

```
x_i=|c_(i+1)-c_i|/(2R),
A_i=-det(c_(i+1)-c_i,q_i-c_i)/2,
score_i=area_i/R^2+A_i/R^2-4lambda*x_i^2+t*gamma_i,
G(theta)=sin(theta)+sin(2theta)-4lambda*sin(theta)^2+t*theta,
Bo(theta)=sin(theta)-4lambda*sin(theta)^2+t*(pi-theta),
C=max_{0<=theta<=pi/2} max(G(theta),Bo(theta)),
Co=max_{0<=theta<=pi/2} Bo(theta).
```

Using the independently accepted stronger normalized facts identified
in Section 2, the conclusion is

```
sum_{i=0}^2 score_i < 3C-19/500.                    (S1)
```

In particular it holds throughout `0<R<=231/2000=.1155`. The endpoint
`794400/5378209` is included. This is the endpoint supplied by the
*particular rational envelopes* in the stronger proof's equation (18);
it is not claimed to be a sharp validity threshold for the local fee.

The anchors are signed, including negative `A_i` on the apex side. No
center norm, reflection, assigned arc length, ownership rule, bound on
`p`, or global covering hypothesis is imposed. Initially repeated
anchors are allowed; the high-score fit excludes the problematic
endpoint repetition rather than assuming it away.

The coefficient pair is fixed exactly here. A neighborhood or rectangle
in `lambda,t` requires its own uniform support/guard proof. The fee
`19/500` exceeds `3/100` by `1/125`; this arithmetic does not itself
establish that separate coefficient-rectangle theorem.

## 2. The stronger normalized premises and their radius audit

For a positive face put `theta=asin(x)` and let its base angles be
`(pi-gamma)/2 +/- delta`. The sine rule and circumradius condition give

```
0<x<=1,  theta<=gamma<=pi-theta,
area/R^2=x^2*cot(gamma/2)-sigma,
sigma=2*x^2*sin(delta)^2/sin(gamma)>=0.
```

In normalized base-midpoint coordinates, write its endpoints as
`(-x,0),(x,0)` and its actual anchor as `(s,h)`, where height is positive
away from the apex. Coverage gives

```
s^2+h^2+2x|s|<=1-x^2,
tau=x*cos(theta)-A/R^2=x*(cos(theta)-h)>=0.
```

The signed identity `A/R^2=xh` is valid on either side. The actual score
is exactly `F-sigma-tau`, with

```
F=x^2*cot(gamma/2)+x*cos(theta)-4lambda*x^2+t*gamma.
```

Its positive gamma second derivative makes its complete maximum an
endpoint value `G` or `Bo`. Thus every face score is at most `C`, with
no initial acute or isosceles restriction. All formulas so far use only
positive scaling by `R`; every nonconstant expression is dimensionless.

The accepted scalar conjugate at the fixed `lambda,t` has no radius
parameter. It gives `C-Co>9/10`, the unique acute peak `theta*` in
`(7/10,24/25)`, and

```
169/250 < y=cos(theta*) < 677/1000.
```

It uses strict concavity on `[0,pi/3]` and decreasing tails on
`[pi/3,pi/2]`, never concavity on the entire interval `[0,pi/2]`.

If (S1) fails, each nonnegative deficit `d_i=C-score_i` is at most
`D0`; write `D=sum d_i<=D0`. The stronger theorem's two newly checked
strict guards `C-G(7/10)>D0` and `C-G(24/25)>D0` localize all theta
values to `(L,H)=(7/10,24/25)`. Its full gamma cut at `zcut=1/10`,
endpoint ratio guards `sin(theta)/sin((theta+zcut)/2)>81/50`, and
whole-window curvature `-G''>4` imply the single paid budget

```
2 sum e_i^2+sum z_i+sum sigma_i+sum tau_i <=D<=19/500,
e_i=theta_i-theta*,  z_i=gamma_i-theta_i>=0.          (S2)
```

The cut proof excludes every larger gamma by full-interval convexity,
including initially obtuse faces. Its new explicit acute-apex guard is

```
gamma_i < H+D0=499/500 < pi/2.
```

These facts depend on `D0,lambda,t`, not `R`. The stronger source's
explicit zero-area guard is also dimensionless (Section 4).

The normalized exact outer tip `a_i` is

```
(a_i-p)/R=exp(i*m_i)[U_i*cos(delta_i)+i*V_i*sin(delta_i)],
rho_i=sin(theta_i)/sin(gamma_i),
U_i=rho_i*(1+cos(gamma_i))+cos(theta_i),
V_i=rho_i*(1-cos(gamma_i))-cos(theta_i),
m_0-m_1=-(gamma_0+gamma_1)/2,
m_2-m_1= (gamma_1+gamma_2)/2.
```

The midpoint follows from the ray-length sine rule, and the outward
normal is `exp(i*m_i-i*delta_i)`; this fixes both signs in the tip
formula. The actual anchor-to-tip error is
`|q_i-a_i|/R<=tau_i/x_i^2`, including `h_i<0`. The tip formulas and
shared lifted phase identities have no absolute position restriction.

The stronger source keeps the angle-box guards

```
16/25<sin(theta_i)<41/50,
573/1000<cos(theta_i)<153/200,
107/50<Q_i=1+2cos(theta_i)<127/50,
|U_i-Q_i|<=(11/4)z_i,
-686839/1250000<=V_i<=-27/250,
```

where the last interval is exactly its analytic `Vlo,Vhi` bound and
already implies `|V_i|<11/20`. Its new squared skew coefficient guard is

```
skew_sq=3081055/4194304 < (43/50)^2.
```

Accordingly the squared skew displacement is at most
`(43/50)^2 sigma_i`. At zero sigma the displacement is zero; a strict
product inequality is not asserted there. This is the needed reading
of the stronger source's displayed equation (12) at zero skew, without
changing that frozen source. The anchor coefficient
`ct=625/256` is unchanged. Both corrections are measured in units of
the same `R`, not in absolute lengths.

After rotating so `m_1=0`, put

```
B=((q_0+q_2)/2-q_1)/R,
b=(q_2-q_0)/(2R).
```

The common apex cancels. Integrating the complete acute-map Jacobians
on the convex box `[L,H]^3`, retaining both imaginary derivative signs,
uses the unchanged squared Frobenius bounds

```
JB2=76960399529/20000000000,
Jb2=77583868411/20000000000.
```

Cauchy against the *one* budget (S2), with skew weight-square sums
`3/2` and `1/2`, gives the new coefficients and linear allowances

```
alphaB2=JB2/2+(3/2)*(43/50)^2=121336399529/40000000000,
alphab2=Jb2/2+(1/2)*(43/50)^2=92375868411/40000000000,
LB=201/50, Lb=529/200, W=sum z_i+sum tau_i,
|B-B*|<=sqrt(alphaB2)*sqrt(D-W)+LB*W,
|b-b*|<=sqrt(alphab2)*sqrt(D-W)+Lb*W.
```

The stronger rational guards `4LB^2 D0<alphaB2` and
`4Lb^2 D0<alphab2` make both maxima occur at `W=0`; the zero-budget
case is immediate. Its new strict squared-fit margins prove

```
|B-B*|<EB=17/50,  |b-b*|<Eb=3/10.                  (S3)
```

No angle, skew, apex slack or actual anchor slack is paid twice.
No BC lens loss, BI credit or other global term is used. The ideal
vectors and size bounds are unchanged:

```
B*=(-A,0), b*=(0,J),
A=(1+2y)(1-y), J=(1+2y)*sqrt(1-y^2),
Amin=19/25<A<763/1000=Amax,
Jmin=173/100<J<87/50=Jmax.                           (S4)
```

This establishes the precise stronger normalized premises used below.
No one of the expressions in this section requires the original
fixed radius `.11321`.

## 3. The unit-circle restriction and exact closed endpoint

For any three unit anchors, let physical endpoint mean and half-chord be
`a=(q_0+q_2)/2` and `h=(q_2-q_0)/2`. Then

```
a dot h=0, |a|^2+|h|^2=1.
```

In the hypothesized high-score domain, (S3)--(S4) give
`|b|>Jmin-Eb=143/100>0`; since `R>0`, `h=Rb` is nonzero. Choose a
unit normal `e` to it and write `a=alpha e`, `q_1 dot e=beta`.
Both scalars are in `[-1,1]` and satisfy

```
(1-alpha*beta)^2-(alpha-beta)^2=(1-alpha^2)(1-beta^2)>=0.
```

Because `1-alpha*beta>=0`, this implies
`|alpha-beta|<=1-alpha*beta`. Also
`|a-q_1|^2+|h|^2=2(1-alpha*beta)`. Since `a-q_1=RB` and `h=Rb`,
division by positive `R` proves

```
2|B_perp_to_b| <= R(|B|^2+|b|^2).                   (S5)
```

This is the only nontrivial radius-dependent condition. The stronger
axis guard is

```
1-[Eb/(Jmin-Eb)]^2 > (97/100)^2.
```

Together with (S3) it makes the horizontal axis project onto the
normal to actual `b` by more than `97/100`. The reverse triangle
inequality applied to ideal `B*` yields

```
2|B_perp_to_b| > K=2[(19/25)(97/100)-17/50]
                  =993/1250.
```

The ideal-size and error bounds are all strict, so

```
|B|^2+|b|^2 < S=(763/1000+17/50)^2+(87/50+3/10)^2
               =(1103/1000)^2+(51/25)^2
               =5378209/1000000.
```

For `0<R<=K/S=794400/5378209`, the high-score assumption therefore
forces the impossible chain

```
K < 2|B_perp_to_b| <= R(|B|^2+|b|^2) < R*S <= K.
```

Both strict fits hold at `D=D0`. The endpoint is included although
`K-Rmax*S=0`; no positive coarse-envelope margin at that endpoint
is claimed. In particular the stronger fixed-radius author's reporting
margin `>9/50=.18` is not transferred to the whole interval.

The simpler interval through `.1155` follows exactly from

```
K-(231/2000)*S = 346433721/2000000000 > 0,
Rmax-231/2000 = 346433721/10756418000 > 0.
```

The original radius `11321/100000` is below `231/2000`, so is inside
this extension. At larger radii the particular envelopes no longer
force (S5) to fail; that does not disprove the local fee there.

## 4. Positive-radius closure and all divisions

Initially `theta=pi/2` has support at most `Co`; the zero-base limiting
scores are at most `t*pi<=Co`. Both are below `C-D0`. Localization
therefore occurs before divisions by `x` or any localized sine. In the
high-score fit, `x>16/25`, `sin(gamma)>0`, `cos(delta)>0`, and
`|b|>143/100`; every denominator and normal used above is justified.
Repeated endpoint anchors `q_0=q_2` are initially allowed and then
excluded by the derived bound on `b`. The three-unit-point proof
otherwise permits repeated anchors.

For every fixed positive `R`, a zero-area limiting face has

```
limsup score <= 1/2+t*pi
             <1/2+t*(22/7)=46790713/70000000<7/10.
```

Coverage gives `A/R^2<=x sqrt(1-x^2)<=1/2`, the quadratic term is
nonpositive, and a subsequence has limiting apex angle in `[0,pi]`.
Since the accepted `C>147/100` gives
`C-D0>179/125>7/10`, and the other two scores are at most `C`, every
zero-area limiting face is excluded from a high-score triple. This
includes collapse of a side other than the base. Limits with three
positive areas retain defined parameters and obey the above proof.
The same argument applies to limits with a positive limiting radius
and vanishing normalized area.

The theorem quantifies at every positive `R`, always using that same
radius in coverage, circumradius, areas and base normalization. The
normalized score is not asserted to be unchanged when a fixed physical
configuration is relabeled with a different radius. `R=0` is excluded.
Physical area tending to zero while `R` tends to zero does not by itself
mean normalized area tends to zero; no such limit is silently assigned
the fixed-positive-radius zero-area score bound.

## 5. Certificate dependency and computation audit

Static inspection of the entire stronger verifier shows that its fixed
`R=11321/100000` is used only in the final
`rhs=R*((Amax+EB)^2+(Jmax+Eb)^2)` and returned radius metadata.
Its `margin>9/50` guard is an original fixed-radius reporting claim,
replaced analytically here by `0<R<=K/S` and the strict envelopes.
Every preceding guard is independent of `R`: scalar support, full
localization, apex cut, acute gamma, zero-area bound, trig window,
tip `U,V`, both-sign Jacobian, shared spending, strict fit, ideal sizes,
and axis projection. The new radius proof does not repeat those guards.

The bundled scalar conjugate proof defines its normalized support at
arbitrary positive `R`; its checker and certificate have no radius
parameter. All four bundled dependency hashes are the same already
accepted hashes used by the earlier theorem. The bundled primitive
declares an unrelated boundary radius `9041/80000`, but `I`, `PI`,
`sinq`, `cosq` and their helpers never reference it. They use only
fixed-precision rational arithmetic, factorials, integer square roots,
and rational Taylor/Machin bounds. That radius belongs to the separate
primitive boundary `verify()` routine; it supplies no radius restriction
for the imported definitions. No such source was imported or executed
for this addendum.

`SOURCE_PINS.json` binds the exact stronger sources and dependency bytes.
The stronger author's certificate remains fixed-radius metadata. The
separate original independent replay receipt is evidence for those
normalized facts, not an execution made here. Final independent static
acceptance of that theorem is bound in `PREMISE_STATUS.json`.
The radius extension has its own independent static child review.

The new `radius_endpoint_fraction.py` uses only standard-library Fraction
arithmetic on the quoted rational envelopes to calculate `K,S,Rmax`,
the `.1155` comparisons and the fee surplus over `.03`. It does no
trigonometry, source import, old or new proof-verifier call, scalar-guard
replay, scan, optimizer, or geometry search. Its single new execution
receipt is `endpoint_arithmetic.json`.

Whole-cover use still needs a valid signed area identity, matching
coefficients, actual eligible ordinary fans and unit anchors, and
disjoint charged marked-face/anchor sets. This addendum changes no
external boundary constant, profile row, topology, incidence theorem,
or nongeneric whole-cover reduction. The normalized local payment is
`19/500`; its corresponding physical area payment in an already valid
identity is `(19/500)R^2`. It supplies no global radius exclusion or
exact optimum.
