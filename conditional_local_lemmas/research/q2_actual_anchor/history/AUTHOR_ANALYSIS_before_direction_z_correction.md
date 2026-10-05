# Two consecutive actual signed-anchor faces: positive conditional loss and its obstruction

Fix

    R=5661/50000=.11322,
    lambda=7/50,
    t=1071883/20000000.

Write `C=C_lambda(t)` for the already accepted complete independent
one-face support. Its acute branch is

    G(theta)=sin(theta)(1+2cos(theta))-4lambda sin(theta)^2+t theta.

The unique acute maximizer is `theta*`, with
`.676<cos(theta*)<.677`. Its accepted certificate is independent of R.
The original scalar or three-face checkers are not rerun here.

For a consistently oriented consecutive two-face fan
`(p,c0,c1),(p,c1,c2)`, define

    x_i=|c_(i+1)-c_i|/(2R),
    A_i=-det(c_(i+1)-c_i,q_i-c_i)/2,
    S_i=area_i/R^2+A_i/R^2-4lambda*x_i^2+t gamma_i.

Here `x_i` is the **base** half-sine. The actual unit anchors satisfy
`|q_i|=1`, `|q_i-c_i|<=R`, `|q_i-c_(i+1)|<=R`; either signed side is
allowed. Every triangle has circumradius at most R.

The assigned-arc parameter is different:

    xi_j=sin(ell_j/2)/R,

where `ell_j` is the length in radians of a connected unit-circle arc
entirely covered by the ball centered at `c_j`. A bound `xi>.6` must
refer to this parameter, not to the base `x`. No complete-cover or
topological conclusion is asserted in this note.

## Conclusions and scope

1. With no assigned-arc footprint constraint, the two actual signed-A
   terms can attain `S0+S1=2C`, even with all four centers in the unit
   disk and both base half-sines greater than `.6`. Thus no strictly
   positive universal two-face fee follows from that domain.
2. If `|p|<=1` and at least one outer center `c0` or `c2` covers a
   connected arc with `xi>=3/5`, there is an explicit conditional fee:

       S0+S1 < 2C-1/10000000.                 (P)

   This weak rational constant is proved below with generous exact
   bounds. It is much smaller than the requested global-row threshold.
3. That stronger arc condition still does not justify the requested
   fee strictly exceeding

       delta=178523206195158271774904671901
             /8011730250000000000000000000000
             approximately .02228272803.

   A fixed algebraic mirrored-skew configuration has both outer
   `xi>.6`, actual consecutive full covered arcs, and actual signed-A
   loss `<.014<delta`. This is a local obstruction, not a full
   100-circle covering or a validation of the high-M profile. The
   additional outer-Q/lens budget, degree-seven empty fan and full
   profile incidence remain separate obligations. In particular this
   note does not assume an unproved outer-Q or longer-arc bootstrap.

## 1. Exact zero-loss pair in the broad actual-anchor domain

Let `theta=theta*`, `y=cos(theta)`, `x=sin(theta)`,
`k=cos(theta/2)`, `l=sin(theta/2)`, `Q=1+2y`. Put

    H=Q l, T=Q k, u=sqrt(1-R^2 H^2),
    p=(u-R T,0),
    c0=p+2R k (y,-x),
    c1=p+2R k (1,0),
    c2=p+2R k (y,x),
    q0=(u,-R H), q1=(u,R H).

Both triangles are isosceles, their actual apex angles are theta,
their circumradii are R and their base half-sines are x. The two q's
are the actual outward lens tips. Their norms are exactly one by the
definition of u. Consequently skew loss and actual-anchor loss both
vanish, and each actual score equals G(theta)=C.

The outer center norms obey

    rho_outer^2=1-2R k u+R^2(2+y-2y^2),

approximately `.9035968^2`. They are too small for assigned `xi>.6`.
For any center of norm rho covering an arc of parameter xi, one must
have

    sqrt(1-R^2 xi^2)-R sqrt(1-xi^2)
          <=rho<=
    sqrt(1-R^2 xi^2)+R sqrt(1-xi^2).           (F)

At xi=.6 the lower root is about `.90711395`. This explains exactly
why the perfect zero-loss pair meets base `x>.6` while failing the
outer assigned-arc condition. The three-face unit condition does not
have a two-face analogue forcing loss: two points can be translated
onto the unit circle.

## 2. An explicit positive conditional fee at assigned xi>=.6

This proof uses Sections 1--3 of the accepted actual-anchor method;
those arguments apply to each face and sum over any number of faces.
It uses neither the three-anchor fit nor any three-face checker run.

Assume (P) fails and let

    d_i=C-S_i>=0, D=d0+d1<=10^(-7),
    theta_i=asin x_i, e_i=theta_i-theta*,
    z_i=gamma_i-theta_i,
    sigma_i=2x_i^2 sin(skew_i)^2/sin gamma_i,
    tau_i=x_i sqrt(1-x_i^2)-A_i/R^2.

The accepted full-domain scalar support, endpoint localization and
convex apex cut yield

    .7<theta_i<.96,
    2 sum e_i^2+sum z_i+sum sigma_i+sum tau_i<=D.     (B)

The actual gammas are acute. The same already accepted interval
constants are sufficient:

    sin theta in (.64,.82), cos theta in (.573,.765),
    Q_i=1+2cos theta_i<2.54,
    1-r_i<=(6/5)z_i,
    |U_i-Q_i|<=(11/4)z_i,
    skew-tip error <=.79 sqrt(sigma_i),
    actual-anchor error <=(625/256)tau_i,

where `r_i` is the normalized circumradius. Their original valid
budget was D<=1/40; here D<=10^(-7) is smaller. All zero-area limiting
faces are excluded before localization because their limiting score
is at most `1/2+t*pi<.7<C-10^(-7)`.

Rotate the shared actual ray `p c1` onto the positive real axis. Let
`w_i=(q_i-p)/R`. Its ideal values are

    w0*=(T,-H), w1*=(T,H),
    T=Q cos(theta*/2) in (2.15,2.16),
    H=Q sin(theta*/2) in (.94,.95).

The derivative norm of `Q(theta) exp(+-i theta/2)` is at most
`2sin theta+Q(theta)/2<3`. Comparing its actual tip to the acute tip,
then to its actual anchor, gives

    |w_i-w_i*|
       <=3|e_i|+.79 sqrt(sigma_i)+4.02(z_i+tau_i)
       <=sqrt(9/2+.79^2)sqrt(D)+4.02D
       <.00075=eta.                              (E)

The Cauchy step uses the one budget (B), not separate worst-case
angle and skew payments. The final rational square guard is
`(9/2+.79^2)*10^(-7)=.00000051241<.00072^2` and
`.00072+4.02*10^(-7)<.00075`.

Let `v=(w0+w1)/2`, `b=(w1-w0)/2`, `a=(q0+q1)/2`. Then

    |v-(T,0)|<eta, |b-(0,H)|<eta,
    a=p+R v, a dot b=0, |a|^2+R^2|b|^2=1.

Write `u=sqrt(1-R^2 H^2)`. We have u>.994 and |a|>.994.
The normalized perpendicular n to b with positive x component satisfies

    |n-(1,0)|<=2eta/(H-eta)<(107/50)eta.

Also

    ||a|-u|
      <=R^2(2H+eta)eta/(|a|+u)<eta/50.

Thus `|a-(+-u,0)|<(9/4)eta`. The negative sign is incompatible
with `|p|<=1`: it would give

    p_x<-.994+(9/4)eta-.11(2.15-eta)<-1.

For the positive sign, with `p*=(u-RT,0)`,

    |p-p*|<(9/4+.12)eta<(12/5)eta=.0018.       (T)

For an outer endpoint, its actual ray angle differs from its ideal
angle by at most `|e_i|+z_i`, and its normalized radial length is
`2r_i cos(gamma_i/2+-skew_i)`. Since

    |e_i|<1/4000,
    sqrt(sigma_i)<1/3000,
    |skew_i|<(9/5)sqrt(sigma_i)<3/5000,

the radial-and-direction error is at most

    3|e_i|+(17/5)z_i+2|skew_i|
       <3/4000+(17/5)10^(-7)+6/5000<1/500.

The skew bound follows from `sin^2(skew)<=625 sigma/512`,
`|skew|<pi/2` and `|skew|<=(pi/2)|sin(skew)|`; using
`pi/2<8/5` and `sqrt(625/512)<111/100` is sufficient.
Together with (T), each outer center differs from the corresponding
ideal outer center by less than

    .0018+.12/500=.00204<.0021.

The ideal outer norm is strictly less than `.904`. One deliberately
loose rational proof uses

    k>.915, u>.994,
    2+y-2y^2<1.764,
    .1132<R<.1133,

and the displayed norm formula from Section 1:

    rho_ideal^2
       <1-2(.1132)(.915)(.994)+(.1133)^2(1.764)
       =.81673120996<.904^2.

Hence the actual outer norms are `<.9061`.

On the other hand, covering an arc with xi>=.6 forces by (F)

    rho>=sqrt(1-.36R^2)-.8R
         >.9976-.8(.1133)=.90696>.9069.

The lower root is increasing in xi on the feasible interval; this also
follows directly from the cosine-law footprint inequality. The final
strict norm contradiction proves (P), including the endpoint D=10^(-7).
It requires an assigned covered arc at only one outer center and no
nearest-owner or reflected-center assumption.

This fee is included solely to make the positive conditional statement
quantitative. It is not proposed as sufficient for the global row.

## 3. Fixed algebraic skew obstruction to the requested threshold

Choose one fixed point, without scanning or optimizing:

    y=cos gamma=17/25,
    x=sin gamma=4sqrt(21)/25,
    k=cos(gamma/2)=sqrt(21)/5,
    l=sin(gamma/2)=2/5,
    cos d=899/901, sin d=60/901.

This is the rational half-angle choice `tan(d/2)=1/30`. Put

    Q=1+2y=59/25, W=2y-1=9/25,
    T=Q cos(d) k+W sin(d) l,
    H=Q cos(d) l-W sin(d) k,
    u=sqrt(1-R^2 H^2).

Equivalently, all required coordinates are exact radicals:

    a=(899sqrt(21)-120)/4505,
    b=(899sqrt(21)+120)/4505,
    T=(53041sqrt(21)+1080)/112625,
    H=(106082-540sqrt(21))/112625,
    v=(8091sqrt(21)+7080)/112625,

    p=(u-RT,0),
    c0=(u-Ra,-2Rbx),
    c1=(u-Rv,0),
    c2=(u-Ra,2Rbx),
    q0=(u,-RH), q1=(u,RH).

The left triangle has rays at angles `-gamma,0`, with radial side
lengths `2R cos(gamma/2-d)=2Rb` and
`2R cos(gamma/2+d)=2Ra`. The right triangle is its mirror. Thus
both circumradii are exactly R and both positive actual apex angles
are gamma. The two skew signs are opposite, with magnitude d.

To check actual tip orientation, in the left bisector frame the base
midpoint is

    R[(1+y)cos(d)-i(1-y)sin(d)].

The outward base normal is `exp(id)`. Adding its tip height Ry gives

    R[Q cos(d)+iW sin(d)].

After multiplying by `exp(-i gamma/2)` this is `R(T-iH)`, precisely
`q0-p`. The reflected expression gives q1. Therefore both q's are
actual covered outward lens tips, not formal replacement anchors;
their anchor slack tau is zero. Unit norms are exact from the
definition of u.

The actual signed-A joint score is exactly

    S0+S1=2G(gamma)-4x sin(d)^2,

so its loss from 2C is

    D_pair=2(C-G(gamma))+4x(60/901)^2.           (L)

One new exact Fraction-only helper guard, retained in the adjacent
receipts, used

    4.582575<sqrt(21)<4.582576,
    .99456<u<.994562,
    .9075<|c0|=|c2|<.9077,
    gamma>.823,
    C<1.473501.

It bounded (L) above by

    106697877127691/8118010000000000<.014<delta.

The geometry values `.7491269` for |p| and `.9501702` for |c1| are
supporting fixed-point evaluations. Coarse exact bounds suffice for
their unit-disk norm caps: `2.16<T<2.18`, `.38<v<.4`, and the
displayed u enclosure imply `0<p_x<.751`, `.949<c1_x<.952`.

Both outer footprint parameters exceed .6. Indeed `.9977^2` is
greater than `1-.36R^2`, so the lower root in (F) is less than
`.9977-.8R=.907124<.9075`; the upper root is greater than
`.99+.8R>1`, while rho_outer<.9077. The maximal connected footprint
half-sine therefore exceeds .6. The more precise supporting values
are outer `xi approximately .6061274` and full footprint angle
`approximately .13735945`.

These decimals are fixed-point evaluations, not optimization results.
The exact upper loss and outer norm guards suffice to obstruct delta.

## 4. Actual consecutive assigned arcs and nearest-owner scope

Let rho=|c0| and define the other unit-circle endpoint of the outer
left ball by radial reflection:

    w0=2[(q0 dot c0)/rho^2]c0-q0.

Define w2 as its mirror across the real axis. Reflection preserves
unit norm and distance to c0, so w0 is a genuine second endpoint of
`B(c0,R)` on the unit circle. The sign `det(c0,q0)>0` makes q0 the
counterclockwise endpoint of that outer footprint. Since c1 lies on
the positive real axis, q0 and q1 are its clockwise and
counterclockwise footprint endpoints. Reflection gives the right
outer endpoint order. The actual assigned arcs are consequently

    c0: counterclockwise arc [w0,q0],
    c1: counterclockwise arc [q0,q1],
    c2: counterclockwise arc [q1,w2].

Each is exactly the full ball/unit-circle footprint. Their interiors
are disjoint. Over the interior of any one of these arcs its center
is at distance <R, while the other two centers are at distance >R.
Thus these are actual connected canonical nearest arcs among the
three centers, with true equal-distance switches q0 and q1.
At the switches, the third center is farther than R. The center p has
norm below .751, so it cannot cover any unit-circle point at radius R.

This argument covers nearest interference among these four centers;
it does not silently prove that a completed degree-seven fan or an
unlisted whole-cover set has no interference. Additional boundary
cells can be placed on disjoint consecutive full-footprint intervals,
but their compatibility with the specified interior degree, empty
circumcircles, four-cell window, outer-Q and exactly 100 circles has
not been certified here. The parent's separately checked new
theta*-skew point is the primary chosen q2 counterexample artifact.

## 5. Execution ledger and independent audit boundary

This author's q2 work consists of analytic derivation and document
writing. I did not execute a q2 author guard, rerun an old verifier,
or call an optimizer. A targeted helper did execute new fixed-point
arithmetic before the parent's subsequent no-further-execution
instruction: four successful mathematical Python 3 calls and one
failed `python` invocation because that binary was absent. Those
calls must not be described as mathematical execution zero. Their
sources, stdout/stderr, order and correction are preserved under
`helper_receipts/`; no further mathematical execution is authorized
or planned here.

Those new evaluations included one exploratory comparison list of seven
prescribed angles, followed by a coarse fixed-angle bound, a fixed
skew-candidate numerical evaluation and its exact radical/score guard.
It would be misleading to call new exploration or scan zero. They did
not run old checkers, import old mathematical modules, resume a historical
search, or call an optimizer. They are independent of the already frozen
three-face package and its counts, which remain unchanged.

The radius-scope agent owns the separately authorized single new
radical guard of the chosen theta*-skew family. A fresh reviewer owns
its full arc order and nearest-interference audit. This note's fixed
algebraic point is an independent backup/analytic cross-check, not a
substitute for that chosen artifact or a claim of complete global
profile verification.
