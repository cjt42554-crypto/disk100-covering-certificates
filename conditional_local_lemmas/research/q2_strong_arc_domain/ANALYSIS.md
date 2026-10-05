# Static progress for the stronger two-face assigned-arc domain

Fix the exact parameters

```
R=5661/50000, lambda=7/50, t=1071883/20000000,
delta=A/B,
A=178523206195158271774904671901,
B=8011730250000000000000000000000,
Knew=(37/20)*(33/32)/R^2.
```

The current parent interface uses the actual middle assigned **HALF**
arc `a_mid`, the joint deficit

```
D_joint=2C-S0-S1+Knew*(a_mid-pi/33)^2,
```

and three consecutive actual assigned full arcs greater than `147/1000`.
Its reflected norm restriction is `|c_j|<=sqrt(1-R^2)`. These are
conditional local premises here. This document proves neither their
global applicability nor the external `Q<.35754` or complete-cover
reduction. It contains no new mathematical execution.

The result established below is a strict exclusion of the **entire
exact-peak mirrored, circumradius-R, outward-tip family** at
`D_joint<=delta`. It is stronger than excluding just the previously
certified `sin(epsilon)=1/12` point. It is not a proof that this family
minimizes the general two-face deficit. The general-fan reduction in
Section 4 remains an exact interface, rather than a closed fee theorem.

## 1. Necessary norm bound from an actual covered arc

Let a radius-R disk with center norm `rho` cover a connected assigned
unit-circle arc of full angle `alpha>147/1000`. A disk with `R<1`
meets the unit circle in one short cap, so `beta=alpha/2` lies in
`(0,asin R]`. Let `e` be the unit radial midpoint of the assigned arc.
Adding the two endpoint squared-distance inequalities gives

```
rho^2+1-2*cos(beta)*(c dot e)<=R^2.
```

Since `c dot e<=rho` and `cos(beta)>0`, this implies

```
rho >= f(beta)=cos(beta)-sqrt(R^2-sin(beta)^2).
```

This argument does not require that the assigned arc be the disk's
complete footprint, radial, or centered on `c`. On the interior of
the feasible interval,

```
f'(beta)=sin(beta)*(cos(beta)/sqrt(R^2-sin(beta)^2)-1)>0,
```

because `cos(beta)^2-(R^2-sin(beta)^2)=1-R^2>0`.
Continuity gives the corresponding strict comparison when an actual
arc reaches the maximum feasible endpoint. Put `z=147/2000=.0735`.
The alternating Taylor bound and elementary cosine bound give

```
sin z>z-z^3/6>.0734,
cos z>=1-z^2/2=.997298875.
R^2=.0128187684,
R^2-(.0734)^2=.0074312084<(.08621)^2=.0074321641.
```

Every finite decimal here and below denotes its exact rational value.
Consequently

```
rho>f(z)>.997298875-.08621=.911088875>.911.       (1)
```

The original parent norm-cut note proves the same necessary bound.
The derivation is repeated here to distinguish an actual assigned arc
from the triangle's normalized base half-length. No assignment rule,
degree, or whole-cover conclusion is inferred from (1).

## 2. Exact-peak mirrored family and its complete loss budget

Let `theta*` be the already accepted unique acute maximizer of

```
G(theta)=sin(theta)+sin(2theta)-4lambda*sin(theta)^2+t*theta,
C=G(theta*).
```

Its accepted enclosure is `.676<y=cos(theta*)<.677`. Write

```
x=sin(theta*)=sqrt(1-y^2),
c=cos(theta*/2)=sqrt((1+y)/2),
s=sin(theta*/2)=sqrt((1-y)/2),
Q=1+2y, V=1-2y<0.
```

Let the admissible base skews be `-epsilon,+epsilon`, with
`|epsilon|<(pi-theta*)/2<pi/2`; in particular the following cosine
is positive. For arbitrary such epsilon, set

```
b=sin(epsilon), a=cos(epsilon),
k=ca-sb, l=sa+cb,
T=Qca-Vsb, H=Qsa+Vcb, M=Vca-Qsb,
u=sqrt(1-R^2 H^2).
```

When the resulting real coordinates and positive side lengths are
defined, the exact family is

```
p=(u-RT,0),
c0=(u-Rk,-2R(ca+sb)x),
c1=(u+RM,0),
c2=(c0_x,-c0_y),
q0=(u,-RH), q1=(u,RH).
```

Both triangles have apex angle `theta*` and circumradius exactly R.
The identities

```
k^2+l^2=1, M^2+H^2=1,
q0-c0=R(k,l), q0-c1=R(-M,-H)
```

show that the two unit anchors are the actual outward lens tips,
with signed-anchor slack zero. Each marked actual score is
`C-2xb^2`, hence its pair deficit is exactly

```
D_actual=2C-S0-S1=4xb^2.                          (2)
```

This is the already independently accepted family's symbolic
identity, specialized to no fixed epsilon. It includes positive,
zero and negative epsilon. Translation onto the other unit-anchor
branch is addressed in Section 3, rather than silently omitted.
Adding the parent's nonnegative middle variance can only increase
the deficit, so `D_joint<=delta` implies `4xb^2<=delta`.

The following rational bounds follow by squaring the accepted y
enclosure and the half-angle identities:

```
.73<x<.74,
.915<c<.916,
.4<s<.403,
2.352<Q<2.354,
|V|<.354.
```

The exact parent target satisfies `delta<23/1000`, as is verified by
the finite integer comparison `1000A<23B`. Thus (2) gives

```
4*.73*b^2<.023,
.023<4*.73*(.09)^2=.023652,
|b|<.09,
a=sqrt(1-b^2)>.995,                              (3)
```

where `.995^2=.990025<1-(.09)^2=.9919`.
This conclusion spends the actual skew deficit once; it does not
spend the joint variance a second time.

## 3. Uniform rational norm exclusion for that whole family

Using `c0=q0-R(k,l)`, its exact outer norm is

```
rho_outer^2=1+R^2-2R*u*k+2R^2*H*l.               (4)
```

First suppose `b>=0`. From (3) and the preceding scalar intervals,

```
k>.915*.995-.403*.09=.874155,
0<l<.403+.916*.09=.48544<.486,
0<H<=Qs<2.354*.403=.948662<.949.
```

For the asserted positive H, a direct lower bound is
`H>2.352*.4*.995-.354*.916*.09>.9`.
Since `R<.114`,

```
R^2 H^2<.013*.901=.011713<1-.994^2=.011964,
u>.994.
```

Equation (4), at the exact R, consequently gives

```
rho_outer^2
 <1+R^2-2R*(.994)*(.874155)+2R^2*(.949)*(.486)
 =.8278871630468752
 <.8281=(.910)^2.                               (5)
```

All constants in this displayed comparison are rational products;
there is no new trigonometric evaluation or radical guard.
The zero-skew case is included: strict scalar interval bounds remain
strict even when b=0.

Now suppose `b<=0`. We instead have

```
k>.915*.995=.910425,
0<l<.403,
0<H<2.354*.403+.354*.916*.09=.97784576<.978.
```

For example `l>.4*.995-.916*.09=.31556>0` and
`H>=Qsa>2.352*.4*.995>0`. Since

```
R^2 H^2<.013*.957=.012441<1-.993^2=.013951,
u>.993,
```

equation (4) now gives the even stronger bound

```
rho_outer^2
 <1+R^2-2R*(.993)*(.910425)+2R^2*(.978)*(.403)
 =.8182098527881312
 <(.910)^2.                                     (6)
```

Reflection gives the same outer norm for c2. Thus the whole positive
translation branch has both outer norms `<.910`, contradicting (1)
if either outer actual assigned arc is greater than .147.

For completeness, the other allowed unit-anchor translation replaces
the common first coordinate of q0,q1 by `-u`, while retaining their
same relative-tip geometry. Its outer norm equation becomes

```
rho_outer,negative^2=1+R^2+2R*u*k+2R^2*H*l>1.
```

The bounds above prove `u,k,H,l>0` throughout the target budget.
That branch is incompatible with the parent's reflected norm cap
`sqrt(1-R^2)<1`. A global rotation of the complete diagram changes
no norm, score, or arc length. The two branches therefore exhaust
the two-unit-anchor translations of this symmetric tip family.

**Restricted-family theorem.** At the exact coefficients and R,
no member of the exact-peak, mirrored-skew, circumradius-R,
actual-outward-tip family can obey both `D_joint<=delta`, the parent
reflected norm cap, and an outer assigned full arc `>.147`.
This conclusion uses only nonnegativity of the joint variance; it
does not need to compute or allocate that variance locally.

## 4. Exact general-fan reduction without a symmetry assumption

For a general positive consecutive two-face fan, rotate the actual
shared middle ray onto the positive horizontal axis. Write

```
theta_i=asin x_i, gamma_i=theta_i+z_i,
r_i=x_i/sin gamma_i,
d_i=the signed base-angle skew,
Wj=(c_j-p)/R, v_i=(q_i-p)/R.
```

These d_i are geometric skews, not the face deficits. With the skew
convention of the accepted actual-tip formula, the exact ray vectors
and shared-side condition are

```
W0=2r0*cos(gamma0/2+d0)*exp(-i gamma0),
W1=2r0*cos(gamma0/2-d0)
  =2r1*cos(gamma1/2+d1)>0,
W2=2r1*cos(gamma1/2-d1)*exp(+i gamma1).            (7)
```

In particular the two skews and angles cannot be varied as unrelated
triangles with two independent middle side lengths.
Let

```
U_i=r_i*(1+cos gamma_i)+cos theta_i,
V_i=r_i*(1-cos gamma_i)-cos theta_i.
```

The exact outward lens tips relative to p are

```
h0=exp(-i gamma0/2)*(U0*cos d0+i V0*sin d0),
h1=exp(+i gamma1/2)*(U1*cos d1+i V1*sin d1),
v_i=h_i+eta_i,
|eta_i|<=tau_i/x_i^2.                            (8)
```

The error bound is the accepted signed-anchor lens inequality and
holds for both signed sides. The actual eta_i remain subject to
the original two endpoint disk inequalities; their directions are
not freely assigned by equality in this norm bound.

More explicitly, let N_i be the actual unit outward base normal
and T_i=i*N_i its directed unit base tangent. In the corresponding
normalized midpoint coordinates the actual anchor has tangent
coordinate v_i,tan and height `cos(theta_i)-Delta_i`, where

```
Delta_i=tau_i/x_i>=0,
eta_i=v_i,tan*T_i-Delta_i*N_i,
v_i,tan^2+2*x_i*|v_i,tan|
                  <=2*cos(theta_i)*Delta_i-Delta_i^2.  (8a)
```

Equation (8a) retains the direction and inward normal sign of the
actual-anchor correction, including negative-height anchors. It
is stronger information than treating eta_i as an arbitrary disk
of radius tau_i/x_i^2.

Define the actual combinations

```
m=(v0+v1)/2, b=(v1-v0)/2,
E=(W0+W2)/2-m, F=(W2-W0)/2,
a=(q0+q1)/2=p+Rm.
```

Relative to the two exact tips, both the mean and half-chord
corrections have the half-endpoint weight:

```
|m-m_tip|, |b-b_tip|
 <=(tau0/x0^2+tau1/x1^2)/2
 <=(625/512)*(tau0+tau1)
```

on the accepted localized window `sin(theta_i)>.64`.
E changes by minus the mean correction, while F is unchanged.
These are two bounds on the same errors, not two available losses.

The two **actual unit-anchor** equations give exactly

```
a dot b=0,
|a|^2+R^2|b|^2=1.                               (9)
```

If `b!=0`, choose a unit normal n to b; then
`a=chi*u*n`, where `u=sqrt(1-R^2|b|^2)` and `chi=+1` or `-1`.
The parent actual middle arc is nonzero, so b!=0 whenever q0,q1
are its genuine distinct endpoints. The precise average outer norm
identity, with the common translation eliminated, is

```
(|c0|^2+|c2|^2)/2
 =u^2+2*chi*R*u*(n dot E)+R^2*(|E|^2+|F|^2).    (10)
```

To verify it, write `c0=a+R(E-F)`, `c2=a+R(E+F)` and average their
squared norms. All cross terms involving F cancel. No inequality or
symmetrization is used. The two strong outer arcs imply by (1) the
necessary strict inequality

```
u^2+2*chi*R*u*(n dot E)+R^2*(|E|^2+|F|^2)
                         >(.911088875)^2.       (11)
```

Thus a sufficiently sharp upper bound on this **average** would
exclude the general pair; it is unnecessary first to show that the
smallest outer norm is maximized by a mirror family.

At exact peak, zero skew and acute endpoints,

```
m*=(Qc,0), b*=(0,Qs),
E*=(-c,0), F*=(0,2cx),
u*=sqrt(1-R^2 Q^2 s^2).
```

Equation (10) recovers the known outer norm formula, so it correctly
captures the cancellation hidden by separate bounds on p and each
outer endpoint.

When q0,q1 are the actual successive middle-cell endpoints in a
connected short assigned arc, its actual HALF angle obeys

```
sin(a_mid)=R|b|, a_mid=asin(R|b|).                (12)
```

Without that endpoint assignment, (12) must not be substituted for
an arbitrary chosen unit-anchor pair. With the assignment, the joint
variance is a function of this same |b|, rather than an independent
angle variable.

## 5. One budget and the remaining closure obligation

For `D_actual=2C-S0-S1<=delta<.023<1/40`, the independently accepted
single-face localization arguments apply face by face and sum over
two faces. They depend on these coefficients and dimensionless
geometry, not on the old proof's radius label. All degenerate,
obtuse and opposite-signed cases are included before localization.
The result is

```
.7<theta_i<.96,
e_i=theta_i-theta*, z_i=gamma_i-theta_i>=0,
sigma_i=2*x_i^2*sin(d_i)^2/sin(gamma_i)>=0,
tau_i=x_i*sqrt(1-x_i^2)-A_i/R^2>=0,
2*(e0^2+e1^2)+z0+z1+sigma0+sigma1+tau0+tau1
                                             <=D_actual.
```

The correctly combined joint budget is consequently

```
2*(e0^2+e1^2)+z0+z1+sigma0+sigma1+tau0+tau1
 +Knew*(asin(R|b|)-pi/33)^2 <=D_joint<=delta.      (13)
```

It is valid only with the actual middle endpoint condition of (12)
and a once-only parent allocation of this variance. There is no
additional local assigned-cap lens loss in (13). In particular tau
already paid in (8) may not be charged again as BC lens loss, and
the variance in (13) cannot also be spent as a free second budget.

Equations (7)--(13), plus the true anchor lens constraints, are the
remaining exact local interface. To prove a general fee above delta,
one must establish that the left side of (11) is at most its right
side throughout this constrained budget and both translation branches
allowed by the reflected norm cap. That inequality is **not proved
here**. Neither the exact-peak mirror exclusion nor the single
counterexample certificate establishes extremality of a mirror family.
The weaker generic Lipschitz estimates from the old three-face method
do not by themselves yield such a sharp outer-norm estimate.

The useful resolved steps are: a strict entire-family exclusion,
an exact shared-side equation, an exact two-unit-anchor translation
elimination, an average-norm target that avoids a symmetry theorem,
and an actual-middle-arc variance expressed through the same chord.
The unresolved step is a sufficiently sharp general constrained
average-norm bound; global incidence, Q and whole-cover reduction
remain separate obligations even if that step is later proved.

## 6. Execution and evidence scope

This directory is a new static addendum. No mathematical code was
executed, no frozen module imported, no old/new checker rerun, and no
optimization or parameter scan was performed. The rational products
in (5) and (6) are displayed so they are independently reviewable.
Metadata hashing and file reads/writes do not count as mathematical
executions. The frozen q2 counterexample and both earlier q3 radius
addenda were neither altered nor replayed.

The only mathematical premise imported by reference is the already
accepted scalar support/peak and its dimensionless single-budget
localization. Source pins record the exact saved inputs used for
this read. They do not imply a fresh replay of those inputs.
