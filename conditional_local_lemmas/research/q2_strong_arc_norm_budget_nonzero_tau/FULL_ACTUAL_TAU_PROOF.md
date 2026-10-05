# General actual signed-A two-face exclusion for full outer arcs above .147

This is a new author analytic proof, awaiting independent static review.
It is separate from the accepted zero-tau proof, the frozen working
interface, and the mean-envelope addendum. It does not use a webpage
model's provisional coefficient. No mathematical program, imported
mathematical module, checker, optimization or parameter scan was run.

Fix `R=5661/50000`, `lambda=7/50`, `t=1071883/20000000`, and the
accepted complete scalar support C. Let `(p,c0,c1)` and `(p,c1,c2)`
be consecutive consistently oriented positive-area triangles of
circumradius at most R. Let their actual covered common unit anchors
be q0 and q1. Use the same actual signed-A scores S_i, base half-sines
x_i, outward lens tips tip_i, and actual slack

    tau_i=x_i sqrt(1-x_i^2)-A_i/R^2 >=0.

Assume the actual middle center has
`|c1|<=sqrt(1-R^2)`, and each actual outer center covers a connected
assigned unit-circle arc of full length strictly greater than .147.
There is no additional base-x, symmetry, apex-position, global Q,
degree, profile, or middle-variance premise. Then

    D=2C-S0-S1 >223/10000=.0223.                 (T)

The parent's joint deficit includes a nonnegative middle half-arc
variance, so (T) also exceeds the parent's required joint payment.
This is a local conditional geometric theorem; its global profile
integration is external to this proof.

## 1. Single deficit budget and actual unit-anchor branch

Suppose `D<=.0223`, and put `W=tau0+tau1` and `B=D-W`.
The accepted per-face deficit mechanism gives

    2 sum(theta_i-theta*)^2+sum z_i+sum sigma_i+W<=D,
    theta_i=asin x_i, gamma_i=theta_i+z_i,
    sigma_i=2x_i^2 sin^2 d_i/sin gamma_i.

Its paid localization and limiting-face closure are those already
pinned in the accepted zero-tau proof; the primitive budget is
dimensionless and applies at this fixed R. Define

    s=(theta0+theta1)/2, Delta=theta0-theta1,
    eta=(d0-d1)/2, Z=(z0+z1)/2,
    alpha=3s/2+Z/2-eta, h=s/2+eta-Z/2.

Removing paid z only simultaneously with its skew coefficient, as
in that proof, gives

    4(s-theta*)^2+Delta^2+2 sum x_i sin^2 d_i+W<=D. (B0)

The resulting ranges, including the elementary eta guard, are

    x_i>.64, .75<s<.904, |Delta|<.15,
    |eta|<.095, 0<=Z<=.01115,
    1<alpha<1.46<pi/2, sin alpha>.8, cos alpha>.1,
    .274<h<.55.

For example `2 sum x_i sin^2 d_i>=2.56 sin^2 eta` forces
`|eta|<.095`, since `sin(.095)>.094` and
`2.56(.094)^2>.0223`. The positive base-angle definition has
`|d_i|<pi/2`, hence the eta branch lies inside the interval on which
this guard applies. All other angle ranges follow by substitution;
the lower cosine polynomial at .55 also gives `cos h>.85`.

Normalize actual anchor errors by

    epsilon_i=(q_i-tip_i)/R,
    U=(epsilon0+epsilon1)/2, V=(epsilon1-epsilon0)/2.

The accepted exact anchor bound `|epsilon_i|<=tau_i/x_i^2`
implies `|U|,|V|<=625W/512<.028`.
In the mean phase phi of the unit-tip directions from c1, write

    b0=exp(i(phi-alpha)), b1=exp(i(phi+alpha)),
    a0=exp(i(phi+h+Delta)), a1=exp(i(phi-h+Delta)),
    B_tip=(b0+b1)/2=cos alpha n0,
    L0=(b1-b0)/2=i sin alpha n0, n0=exp(i phi),
    H0=(a0+a1)/2=cos h exp(i(phi+Delta)).

Set `L=L0+V`, `n=-iL/|L|`, `u=sqrt(1-R^2|L|^2)`.
Here `|L|>.772`, and
`|n-n0|<=2|V|/(|L0|-|V|)<.073`, so `n dot n0>.99`.
Consequently `(B_tip+U) dot n>.1(.99)-.028=.071`.
Actual middle coverage gives

    |B_tip+U|^2+|L|^2
      =(|b0+epsilon0|^2+|b1+epsilon1|^2)/2<=1,

so `u>=sqrt(1-R^2)`. Unit-anchor equations force their mean A
to equal `+-u n`, and `c1=A-R(B_tip+U)`.
The negative branch would have
`c1 dot n<-u<=-sqrt(1-R^2)`, violating the actual middle norm cap.
The actual branch is therefore positive. This is a statement about
the same actual c1 and fan, with no independent radial reflection.

## 2. Exact rigid fan path from unit tips to the actual anchors

For `0<=v<=1`, put

    L_v=L0+vV, n_v=-iL_v/|L_v|,
    u_v=sqrt(1-R^2|L_v|^2), A_v=u_v n_v,
    c1(v)=A_v-R(B_tip+vU),
    f_i=(c_outer_i-c1)/R=b_i-a_i,
    c_outer_i(v)=c1(v)+R f_i,
    q0(v)=A_v-RL_v, q1(v)=A_v+RL_v.

Translate p by the same vector as c1. All fan shape parameters,
areas, base lengths and circumradii remain unchanged. The q's have
unit norm and their normalized endpoint offsets are exactly

    q_i(v)-c1(v)=R(b_i+v epsilon_i),
    q_i(v)-c_outer_i(v)=R(a_i+v epsilon_i).

The covered-lens intersections are convex, so these offsets stay
covered. In particular `|L_v|<=1` for every v, since the two middle
offsets have squared-norm average
`|B_tip+vU|^2+|L_v|^2<=1`.
The actual signed-A slacks are exactly `v tau_i` and
`D(v)=D-W+vW`. Thus the unit-tip shape at v zero has deficit B.
No preservation of outer arcs, all global norm caps, or middle
variance along this path is assumed. The endpoint v one is the
actual fan by Section 1.

Write

    H_v=H0+vU,
    J0=(a1-a0)/2, F=L0-J0,
    C_v=A_v-RH_v.

The physical outer half-difference is the constant RF, while their
mean is C_v. Hence their exact mean squared norm is

    M(v)=|C_v|^2+R^2|F|^2.                      (M)

The independently versioned `MEAN_BUDGET_ENVELOPE.md`, together
with its `ENDPOINT_CLARIFICATION.md`, proves the unit-tip shape bound

    M(0)<=N(B),
    N(B)=1+R^2-2R cos(.523+sqrt[(47/100)B]),
    N'(B)>.30 for 0<B<=.0223.                  (N)

Its weighted common-skew minimization pays Delta in the same B;
it is a mean-norm envelope, so no minimum-norm inequality is being
substituted here. The source and clarification are required inputs
to this proof, and require their own static review.

## 3. Exact directional supports and their unit-tip coefficients

Let `e_i` be the directed base tangent and `n_i` the base normal
toward p. In normalized lens coordinates,

    epsilon_i=s_i e_i+v_i n_i, v_i=tau_i/x_i,
    s_i^2+2x_i|s_i|+v_i^2<=2y_i v_i,
    y_i=cos theta_i.

Thus `|s_i|<=cot(theta_i)v_i`. For any vector g,

    g dot epsilon_i
       <=[g dot n_i+cot(theta_i)|g dot e_i|]tau_i/x_i.

The exact phase convention of Section 1 gives

    n0-cot(theta0)e0=-i b0/x0,
    n0+cot(theta0)e0=+i a0/x0,
    n1+cot(theta1)e1=+i b1/x1,
    n1-cot(theta1)e1=-i a1/x1.

In these four identities n0,n1 denote the two base normals, not
the separate mean-phase vector called n0 in Section 1. Define
the unambiguous support functionals

    ell0(g)=max((-i b0) dot g,(+i a0) dot g)/x0^2,
    ell1(g)=max((+i b1) dot g,(-i a1) dot g)/x1^2.

Then `g dot epsilon_i<=tau_i ell_i(g)`.
The four vectors inside the maxima are unit vectors, so

    |ell_i(g)-ell_i(g')|<=|g-g'|/x_i^2.         (Lip)

In the instantaneous frame of L_v, differentiation of (M) gives

    M'(v)=g0(v) dot epsilon0+g1(v) dot epsilon1,
    g0=-RC_v-(DA_v)^T C_v,
    g1=-RC_v+(DA_v)^T C_v.                     (G)

At v zero, let `r=sin alpha`, `u=sqrt(1-R^2r^2)`.
In the orthonormal mean-phase frame the gradients have components

    g0(0)=R(-P0,Q0), g1(0)=R(-P1,-Q1),
    P0=u-R cos h cos Delta+(u/r)cos h sin Delta,
    P1=u-R cos h cos Delta-(u/r)cos h sin Delta,
    Q0=R[r+cos h sin Delta-(Rr/u)cos h cos Delta],
    Q1=R[r-cos h sin Delta-(Rr/u)cos h cos Delta].

Each P is greater than
`.993-.114-.15/.8=.6915`, and each Q is less than
`.114(1+.15)<.132`.
The base tangent angles involve
`w0=alpha-theta0`, `w1=alpha-theta1`, with `.205<w_i<.628`.
Thus `P_i tan w_i>(.6915)(.205)>.132>Q_i`.
This proves `g0(0) dot e0<0`, `g1(0) dot e1>0`, so the b
branches, rather than the a branches, realize the unit-tip maxima.

Exact substitution simplifies those two maxima to

    ell0(g0(0))=R(u-R cos alpha)/x0^2
      [sin alpha(1-R cos h cos Delta/u)+cos h sin Delta],
    ell1(g1(0))=R(u-R cos alpha)/x1^2
      [sin alpha(1-R cos h cos Delta/u)-cos h sin Delta]. (C0)

Here `u-R cos alpha<.989`.
Also `cos h>.85`, `cos Delta>.988`, and `u<=1` imply
`1-R cos h cos Delta/u<.905`.
To bound the first bracket divided by x0 squared, put
`a=|Delta|<=.15`.
If Delta is nonnegative, its numerator is less than `.905+a`,
and

    2x0^2>1-cos(1.5+a)>.929+.996a>.905+a.

If Delta is nonpositive, its numerator is less than
`.905-.8466a`, since `cos h>.85` and `sin a>=.996a`, and

    2x0^2>1-cos(1.5-a)>.929-a>.905-.8466a.

The last positive gap is at least
`.024-.1534(.15)=.00099`.
These guards use `cos(1.5)<.071` from the eighth-order upper
cosine polynomial. On `[1.5,1.65]`, `sin z>.996`, since
`|z-pi/2|<.08` and `cos(.08)>1-(.08)^2/2>.996`.
At a zero endpoint the product comparison is non-strict; the
displayed numerator/denominator gap remains strict.
For the second coefficient, replace Delta by minus Delta and x0
by x1 in the same argument. Consequently

    ell_i(g_i(0))<2R(.989)=.22394916<.224.      (C1)

## 4. A derivative perturbation bound retains the factor R

For any L along the path define

    r=|L|, u=sqrt(1-R^2r^2), A(L)=u(-iL/r).

On the entire path `.772<r<=1`, `.993<u<=1`.
The eigenvalues in the radial and tangent directions, up to the
fixed orthogonal rotation, are `u'=-R^2r/u` and `u/r`.
Therefore

    ||DA||<1.3.

For completeness, writing `v=L/r`, the derivative of normalization
has norm at most `1/r`, and its second derivative norm is at most
`3/r^2`. Product differentiation of `A=u(-iv)` thus gives

    ||D^2A||<=|u''|+3|u'|/r+3u/r^2,
    u''=-R^2/u-R^4r^2/u^3.

Using `R^2<.013`, `u>.993`, `r>.772` gives
`|u''|<.014`, `3|u'|/r<.04`, `3u/r^2<5.04`.
In particular `||D^2A||<5.2`. These bounds hold on the whole
line segment joining L0 to every L_v, so

    |A_v-A_0|<=1.3 v|V|,
    ||DA_v-DA_0||<=5.2 v|V|.

Differentiate `|A|^2+R^2|L|^2=1` to obtain the exact cancellation

    (DA)^T A=-R^2L.

Substituting it into (G) is essential: it yields

    g0=R[-A+RH+RL+(DA)^T H],
    g1=R[-A+RH-RL-(DA)^T H].                  (Rfact)

Since `|H0|<=1` and `H_v-H0=vU`, the preceding bounds give

    |g_i(v)-g_i(0)|/R
      <=(1.3+R+5.2)v|V|+(1.3+R)v|U|
      <=6.614|V|+1.414|U|
      <=8.028(625/512)W.

After division by `x_i^2>.4096`, a fully explicit loose guard is

    |g_i(v)-g_i(0)|/x_i^2
      <(.114)(8.03)(1.221)(.0223)/.4096 <.061.

Indeed the numerator is `.024925330386`, while
`(.061)(.4096)=.0249856`.
At W zero the gradient difference is zero, so this strict fixed
positive bound still holds. This calculation is a paper rational
guard; no numerical evaluator was run.

By (Lip) and (C1), every support coefficient on the actual path
satisfies

    ell_i(g_i(v))<.224+.061=.285<.30.          (C2)

No sign premise on `g_i(v) dot e_i` at positive v is used: the
maximum includes both support branches, and (Lip) transports it.
Thus (G) and the lens supports prove the uniform bound

    M'(v)<=.285W  for every 0<=v<=1.           (Pay)

When W is positive its coefficient comparison is strict; when W
is zero the derivative and the payment are both exactly zero.

## 5. Closing the complete actual-anchor budget

Integrating (Pay), then using (N), gives

    M(1)<=M(0)+.285W
          <=N(D-W)+.285W<=N(D)<=N(.0223).

The middle inequality follows from `N'(B)>.30>.285` on positive
B and continuity at B zero. This uses W only once, as part of the
same D, and uses no BC-lens or middle-variance credit.
The lower cosine polynomial gives `cos(.626)>.81`, while
`.523+sqrt[(.47)(.0223)]<.626`; hence

    M(1)<1+R^2-2R(.81)=.8294023684.            (Upper)

A connected actual full assigned arc greater than .147 requires
the covering center norm to exceed `.911088875`, as proved by the
fixed rational guards in the accepted zero-tau proof:

    sin(.0735)>.0734,
    cos(.0735)>=.997298875,
    R^2-(.0734)^2<(.08621)^2.

In particular both actual outer norms would exceed .911, making
their mean squared norm greater than `.911^2=.829921`.
This contradicts (Upper) and proves (T), for arbitrary actual tau.

## 6. Exact scope, inputs, and execution record

Required sources are the accepted per-face input and peak enclosure
pinned by the zero-tau proof; that proof and its three preserved
clarifiers; `MEAN_BUDGET_ENVELOPE.md` and its separate zero-H
clarifier; and the path identities restated fully above. No accepted
three-face compatibility theorem is treated as a two-face theorem.
The general two-face proof has no independent-J relaxation and
charges no global Q, additional arc variance, or other marked fan.

This new author proof requires full independent static review,
including the phase convention, actual branch, both support branches,
mean-envelope payment, Hessian bounds, and all rational guards.
Earlier zero-tau review does not by itself accept these new steps.

Author mathematical program execution: 0. Mathematical source
imports, old checker replays, optimizer calls, parameter scans and
search: 0. Filesystem source writing and SHA256 hashing are the only
tool operations for this proof. Earlier exploratory calls in the
separate q2 directories retain their separate honest ledgers and
are not reclassified by this proof's zero-execution statement.
