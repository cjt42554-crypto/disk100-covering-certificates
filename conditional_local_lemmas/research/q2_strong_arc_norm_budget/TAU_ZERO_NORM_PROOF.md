# General two-face strong-arc exclusion when actual anchor slack is zero

Fix `R=5661/50000`, `lambda=7/50`,
`t=1071883/20000000`, and the accepted complete scalar support
`C=C_lambda(t)`. Let `(p,c0,c1),(p,c1,c2)` be consecutive consistently
oriented positive-area faces, each of circumradius at most R.
The actual common unit anchors satisfy

    |q_i|=1, |q_i-c_i|<=R, |q_i-c_(i+1)|<=R.

Their actual signed scores and normalized bases are

    x_i=|c_(i+1)-c_i|/(2R),
    A_i=-det(c_(i+1)-c_i,q_i-c_i)/2,
    S_i=area_i/R^2+A_i/R^2-4lambda*x_i^2+t gamma_i,
    D=2C-S0-S1.

Assume the reflected norm cap on the common middle center
`|c1|<=sqrt(1-R^2)`. The outer centers must each cover a connected
actual unit-circle assigned arc of full length greater than `.147`.
Assume additionally, and essentially for this lemma,

    tau_i=x_i sqrt(1-x_i^2)-A_i/R^2=0, i=0,1.    (T0)

Then

    S0+S1 < 2C-223/10000,
    equivalently D>223/10000=.0223.               (1)

The bound `.0223` is strictly greater than the parent's required
`.02228272802808835...`. Adding the nonnegative middle half-arc
variance with `K=(37/20)(33/32)/R^2` only increases the deficit.
This lemma is general over angles, actual apex slack and skew; it
does not assume an isosceles or equal-angle fan. It is, however,
only a **zero actual-anchor-slack** result. The extension to arbitrary
positive tau is an outstanding separate interface in Section 7.

No old checker, conjugate checker, optimizer, parameter enumeration
or new mathematical program was executed to produce this note.
The argument uses accepted per-face constants plus elementary
explicit rational Taylor inequalities. It awaits static independent
review; it is not presented as a machine-replayed new certificate.

## 1. One high-score deficit budget before angle localization

Suppose (1) fails and put `D0=223/10000`, so `0<=D<=D0<1/40`.
Every face has `S_i<=C` by the accepted full-apex scalar support.
Thus every individual deficit is at most D0. The accepted per-face
endpoint localization and apex cut, whose original valid budget
was `1/40`, give

    theta_i=asin x_i in (.7,.96),
    gamma_i=theta_i+z_i, z_i>=0,
    sigma_i=2sin(theta_i)^2 sin(d_i)^2/sin(gamma_i),
    2 sum(theta_i-theta*)^2+sum z_i+sum sigma_i<=D.  (B)

Here the actual base-angle skews are chosen so the common middle
ray length obeys

    |c1-p|/R
       =2r0 cos(gamma0/2+d0)
       =2r1 cos(gamma1/2-d1),
    r_i=x_i/sin gamma_i<=1.

The scalar peak has `.828<theta*<.829`. Its saved certificate and
the accepted per-face proof are bundled under `accepted_inputs/`
without replay. The curvature and slope used for (B) are the same
accepted `-G''>4` and `-F_gamma>1`, not a global concavity premise.

The accepted sine lower bound `sin theta_i>.64` gives

    sin^2 d_i <=(625/512)sigma_i
                <=(625/512)D0<.028.

In particular `|d_i|<.18`, since `sin(.18)>.17` and `.17^2>.028`.
Define

    s=(theta0+theta1)/2,
    E=s-theta*, Delta=theta0-theta1,
    Z=(z0+z1)/2,
    eta=(d0-d1)/2, m=(d0+d1)/2.

Then `2 sum(theta_i-theta*)^2=4E^2+Delta^2`.
Because `sqrt(D0)/2<.075`,

    .75<s<.904,
    |Delta|<.15, |eta|<.18, 0<=Z<=.01115.

Every zero-area limiting face is excluded before these divisions:
its limiting score is at most `1/2+t*pi<.7<C-D0`. This supplies
the same closure as in the accepted proof without any generic
assumption on a complete covering.

## 2. Exact coupled center-norm formula at actual tau=0

Condition (T0) forces q_i to be the actual outward lens tip:
in normalized base coordinates the covered-lens constraint is
`s_base^2+h_base^2+2x|s_base|<=1-x^2`; equality
`A/R^2=x sqrt(1-x^2)` forces `h_base=sqrt(1-x^2)` and `s_base=0`.

Rotate the actual common ray `p c1` onto the positive real axis.
The unit directions from c1 to the two actual tips are

    b0=exp i(d0-gamma0/2-theta0),
    b1=exp i(d1+gamma1/2+theta1).

The unit directions from the outer centers to the same tips are

    a0=exp i(d0-gamma0/2+theta0),
    a1=exp i(d1+gamma1/2-theta1).

These follow directly from midpoint plus outward base normal;
they are not formal independent-J anchors. Put

    alpha=3s/2+Z/2-eta,
    h=s/2+eta-Z/2,
    alpha+h=2s,
    u=sqrt(1-R^2 sin^2 alpha),
    beta=asin(R sin alpha),
    chi=h+beta.

Let phi be the mean phase of b0,b1. Their mean and half-difference
are `cos(alpha) exp(i phi)` and `i sin(alpha) exp(i phi)`.
Unit-anchor equations force

    (q0+q1)/2=+-u exp(i phi),
    (q1-q0)/2=iR sin(alpha) exp(i phi).

The common middle center is consequently

    c1=(+-u-R cos(alpha)) exp(i phi).

The above paid ranges give `0<alpha<1.55<pi/2`, so `cos alpha>0`.
The negative sign would imply
`|c1|=u+R cos alpha>u>=sqrt(1-R^2)`, contrary to its norm cap.
Thus the positive unit-anchor branch is forced without a bound on p.

Subtracting the outer unit offset vectors gives the exact identities

    |c0|^2=1+R^2-2R cos(chi+Delta),
    |c2|^2=1+R^2-2R cos(chi-Delta).              (N)

In particular the requested mean-square identity is

    (|c0|^2+|c2|^2)/2
       =1+R^2-2R cos chi cos Delta.

The minimum is sharper. The paid ranges give
`0<chi-|Delta|<=chi+|Delta|<pi/2`; for example
`h>.375-.18-.01115/2>.189>|Delta|`, while
`chi+|Delta|<.452+.18+.114+.15<.9`.
Hence

    min(|c0|^2,|c2|^2)
       =1+R^2-2R cos(chi-|Delta|).              (M)

Angle asymmetry helps the minimum bound. It must not be charged
as an independent positive norm perturbation. The exact common
middle ray relation was retained in the construction; its radius
ratios disappear from (N) because every actual tip is at distance
exactly R from c1.

## 3. Actual apex slack can be removed without creating skew credit

First,

    d chi/dZ=-1/2+R cos(alpha)/(2u)<0,

since `u>=sqrt(1-R^2)>.993` and `R<.114`.
Thus setting Z=0 can only increase the norm-driving angle chi.
This step alone would not justify changing its skew payment.
The following simultaneous budget keeps that payment valid.

Since `sin(theta_i+z_i)<=sin theta_i+z_i`,

    sigma_i
      >=2sin(theta_i)^2 sin(d_i)^2/(sin theta_i+z_i)
      >=2(sin theta_i-z_i)sin(d_i)^2.

Consequently (B) gives

    D>=4E^2+Delta^2+2 sum sin(theta_i)sin(d_i)^2
                +sum z_i[1-2sin(d_i)^2].

The coefficients in brackets are positive, because
`sin^2 d_i<.028`. Dropping only these nonnegative paid z terms yields

    4E^2+Delta^2+2 sum sin(theta_i)sin(d_i)^2<=D. (B0)

The remaining shape can therefore be bounded by

    chi0(s,eta)=s/2+eta+asin(R sin(3s/2-eta)).

This is an upper relaxation of the actual fan, not an asserted
realizable new fan or an unproved symmetry reduction.

## 4. Removing angle asymmetry pays for the apparent cheaper skew

Let `v=|Delta|`. The sine Lipschitz bound gives
`sin theta_i>=sin s-v/2`, and

    sin^2(m+eta)+sin^2(m-eta)
       =1-cos(2m)cos(2eta)>=2sin^2 eta,

because `|eta|<.18` gives `cos(2eta)>0`. Thus (B0) implies

    4E^2+v^2+(4sin s-2v)sin^2 eta<=D.           (C)

The function chi0 is increasing in eta throughout the paid domain,
because its derivative is `1-R cos(3s/2-eta)/u>0`.
If eta<0, replace eta by zero for an upper norm bound and a lower
payment. Otherwise let eta'>=0 be defined by

    sin eta'=sqrt(1-v/(2sin s)) sin eta.

The factor is real since `v<.15` and `sin s>.68`.
It ensures

    4E^2+4sin s sin^2 eta'<=D.                  (C')

Moreover eta'<=eta<.18, and

    eta-eta'
      <=[sin eta-sin eta']/cos eta
      <=v tan eta/(2sin s)<v/6.

The last inequality uses `tan(.18)<1/5` and `sin s>.68`.
Along this eta interval, `0<d chi0/d eta<1`. Therefore

    chi0(s,eta)-v <=chi0(s,eta')-5v/6
                         <=chi0(s,eta').       (S)

Combining (M), removal of Z and (S), the outer minimum is bounded
by the scalar angle chi0(s,eta') with the single payment (C').
This explicitly accounts for the benefit of unequal angles rather
than assuming equality of theta, gamma, skew, radii or radial sides.

## 5. The resulting one-budget norm bound is below .911

Since `sin s>sin(.75)>.68`, (C') implies

    4E^2+(68/25)sin^2 eta'<=D<=223/10000.

In particular eta'<.092: the elementary lower Taylor bound gives
`sin(.092)>.091`, and `(68/25)*(.091)^2>.0223`.
For `0<=eta'<.092`,

    eta'<=1.002 sin eta'.

This follows from `sin z>=z(1-z^2/6)` and
`1/[1-(.092)^2/6]<1.002`; its zero endpoint is non-strict.

On the whole rectangle `.75<=s<=.904`, `0<=eta'<=.092`, the
angle `alpha0=3s/2-eta'` lies between 1 and pi/2. Hence

    0<d chi0/d eta'<1,
    d chi0/ds
       =1/2+3R cos(alpha0)/(2sqrt(1-R^2 sin^2 alpha0))
       <1/2+3(.114)(.55)/(2(.993))<.6.

Here `cos(1)<.55`, and these derivative bounds hold along the
entire segment from `(theta*,0)` to `(s,eta')`.

The ideal starting angle is less than .523. Indeed
`theta*<.829`, `cos theta*>.676`, so

    (1+2cos theta*)sin(theta*/2)<2.354*.403=.948662,
    R*.948662<.1075,
    sin(.108)>=.108-(.108)^3/6>.1075.

Thus

    chi0(theta*,0)=theta*/2+asin(R(1+2cos theta*)sin(theta*/2))
                   <.4145+.108<.523.

Cauchy against the same paid budget now gives

    chi0(s,eta')
       <.523+.6|E|+1.002sin eta'
       <=.523+sqrt[(.6^2/4+1.002^2/(68/25))D]
       =.523+sqrt[(312201/680000)D]
       <.523+.102=.625.

The last guard is the exact rational inequality

    (312201/680000)(223/10000)
       =69620823/6800000000 <(.102)^2.

Finally the alternating cosine lower polynomial gives

    cos(.625)
      >=1-(.625)^2/2+(.625)^4/24-(.625)^6/720>.81.

By (M) and the relaxed angle bound,

    min(|c0|^2,|c2|^2)
       <1+R^2-2R(.81)
       =.8294023684 <(.911)^2=.829921.          (U)

Thus at least one actual outer center has norm below .911.

## 6. Both actual full assigned arcs require norm above .911

Let beta=.147/2=.0735. A norm-rho center covering an actual
connected arc of full length at least .147 must satisfy

    rho>=cos beta-sqrt(R^2-sin^2 beta).

This follows from the cosine-law footprint condition and is
independent of the fan's base half-sines. The direct rational guards
are

    sin beta>beta-beta^3/6>.0734,
    cos beta>=1-beta^2/2=.997298875,
    R^2-(.0734)^2=.0074312084 <(.08621)^2=.0074321641.

Therefore

    rho>.997298875-.08621=.911088875>.911.

Each outer assigned arc is strictly longer than .147, so each
outer norm must exceed .911. This contradicts (U) and proves (1).
No global Q, middle variance, other marked fan or BC-lens credit
has been charged by this local zero-tau lemma.

## 7. Exact arbitrary-anchor interface and the remaining gap

For arbitrary actual tau, use the actual unit anchors, not the
unit-tip replacements. Define

    A=(q0+q1)/2,
    L=(q1-q0)/(2R),
    H=[(q0-c0)+(q1-c2)]/(2R),
    J=[(q1-c2)-(q0-c0)]/(2R).

The exact unit conditions and outer mean-square identity are

    A dot L=0, |A|^2+R^2|L|^2=1,
    (|c0|^2+|c2|^2)/2
       =1+R^2(|H|^2+|J|^2)-2R A dot H-2R^2 L dot J. (A)

Both endpoint disks cover each actual anchor, so
`|H|^2+|J|^2=(|q0-c0|^2+|q1-c2|^2)/(2R^2)<=1`.
The zero-tau formulas are specializations of (A); they cannot be
substituted into (A) at positive tau without an error argument.

The accepted anchor bound `|q_i-tip_i|/R<=tau_i/x_i^2` is too coarse
by itself for the desired strong constant. A useful sharper interface
keeps the lens direction. In base coordinates, with n pointing
toward the apex, write the actual lens displacement as

    q-tip=R[s e+v n], v=tau/x>=0,
    s^2+2x|s|+v^2<=2yv, y=sqrt(1-x^2).

For any signed linear direction `a e+b n`, this gives

    (a e+b n) dot(q-tip)/R
       <=[b+|a|(2y-v)/(sqrt(x^2+2yv-v^2)+x)]v
       <=[b+|a|y/x]tau/x.                    (Ldir)

Keeping b's sign can be substantially sharper than the Euclidean
norm bound. The parent is developing a separate exact lemma for
this interface. A whole-paid-domain norm bound from (A), with the
angular/skew budget and (Ldir) sharing one D, is still required.

In particular this document does not prove a general actual signed-A
strong-arc fee when tau>0, does not discard anchor slack as a
separately paid variable, and does not transfer a zero-tau result
to a global ordinary q2 profile. Its reusable outputs are the exact
general norm identity, the asymmetry/payment relaxation, and the
complete zero-tau conditional exclusion above.

## 8. Inputs, execution and audit boundary

The accepted sources are pinned as follows:

    actual_anchor_1_40_PROOF.md
       58dc829a24518ecacc44e3414c89973eb92b94573b72a631fc78e41d023173c1,
    actual_anchor_1_40_certificate.json
       c92f9879d82e80ca134eb93379d1e842189f836c5e5603d954361c95b7f286da.

Only their per-face deficit mechanism and peak enclosures are used.
Their three-face compatibility conclusion is not invoked as a q2
statement. No imported mathematical source or checker was executed.
The independent symbolic helper also performed no mathematical
program execution, old checker call, scan or file mutation.

This is a new analytic proof awaiting independent static review.
Any later nonzero-tau addendum or new computation must use a new
source/receipt and retain this frozen zero-tau version, rather than
silently widening its scope.
