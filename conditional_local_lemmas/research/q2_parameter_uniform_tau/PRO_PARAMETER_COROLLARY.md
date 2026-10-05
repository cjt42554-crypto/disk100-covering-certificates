# Separate continuous-parameter corollary of the Pro .0225 theorem

Status: new analytic author proof, awaiting independent static review.
The frozen local .0223 parameter proof is unchanged. This corollary
uses the separately accepted Pro .0225 norm-domain theorem, and
rederives the per-face .025 budget needed for parameter transport.
It does not use the local .0223 geometric conclusion.

No new mathematical program, source import, checker, numerical replay,
old 48-guard batch, optimization or parameter scan was executed.
The accepted Pro source's earlier root replay remains its ONE
separate receipted execution; it is not repeated or counted here.

## 1. Exact central theorem, independent domain and changed objective

Fix R=5661/50000. Let (p,c0,c1),(p,c1,c2) be consecutive positively
oriented nondegenerate actual triangles of circumradius at most R.
For each base, the corresponding actual unit anchor is covered by
both endpoint disks. Define actual signed A, normalized x, actual
apex gamma, and the matched score/support exactly as follows:

    A_i=-det(c_(i+1)-c_i,q_i-c_i)/2,
    x_i=|c_(i+1)-c_i|/(2R), theta_i=asin x_i,
    S_i(lambda,t)=area_i/R^2+A_i/R^2-4lambda*x_i^2+t gamma_i,
    C(lambda,t)=max_{0<=theta<=pi/2}{G(lambda,t;theta),Bo(lambda,t;theta)},
    G=sin theta+sin(2theta)-4lambda sin^2 theta+t theta,
    Bo=sin theta-4lambda sin^2 theta+t(pi-theta),
    D(lambda,t)=2C(lambda,t)-sum_i S_i(lambda,t).

These two support branches include the complete allowed apex domain;
their convex endpoint reduction is unchanged when lambda,t vary.
Actual A may be negative and actual tau is arbitrary. No symmetry,
equal-distance switch, closest owner, assigned arc, global Q, degree,
variance, or radial-inner positioning is assumed.

The geometric hypotheses are ONLY the original Pro norm interface:

    |c1|^2<=1-R^2, |c0|^2>=83/100, |c2|^2>=83/100.   (Norm)

These are the SAME actual three centers in the actual fan defining
the scores. They are not imposed after moving or reflecting centers.

At lambda0=7/50,t0=1071883/20000000, the separately accepted
Pro theorem gives

    D0=D(lambda0,t0)>dP=9/400=.0225.                 (Pro)

Its source was fully read, and the complete independent static audit
plus one already-completed, byte-matching isolated root replay are
pinned dependencies. That theorem covers (Norm), which is a different
domain from the local actual-outer-full-arc>.147 theorem. This source
retains (Norm) exactly; it does not infer it from an arc statement.

Put signed increments u=lambda-lambda0, v=t-t0 and a=|u|,b=|v|.
For the full CLOSED parameter rectangle

    a<=3/10000=.0003, b<=1/5000=.0002,              (PRect)

the new conclusion is

    D(lambda,t)>dP-(3/5)a-(69/400)b
                >=44571/2000000=.0222855.           (PT)

This is a local conditional theorem at exactly R22. It provides no
global covering or radius exclusion, and no automatic extension of
the Pro norm-domain theorem to another radius.

## 2. Rederive the central per-face budget up to .025

This is an explicit new extension of the scalar-budget implication
used in the Pro proof: it is valid for D0<=1/40, rather than only
its geometry-contradiction range D0<=9/400. None of the later Pro
two-face geometry or slab guards is invoked to establish it.

The complete central scalar support has a unique acute peak
theta* in (.828,.829). Write x*=sin theta*. The acquired Pro source's
Section 1 gives, with its acute branch called F there,

    C0>=G0(5/6)>7/5, Bo0(theta)<6/5,
    G0(5/6)-G0(.7)>.035,
    G0(5/6)-G0(.96)>.035,
    -G0''>4 on [.7,.96].                           (Scalar)

Its single-peak proof supplies the full acute monotonicity outside
the peak, not global concavity. The common frozen central scalar
peak enclosure additionally gives .676<cos theta*<.677; this is
a scalar premise identifying the same exact support maximizer,
not a use of any local .0223 geometric theorem.

Suppose D0<=1/40=.025. Every central face deficit is nonnegative
and at most .025, so S_i>=C0-.025>1.375>6/5. Convex apex endpoint
interpolation, with its obtuse branch below 6/5, therefore implies
G0(theta_i)>=S_i. Hence C0-G0(theta_i)<=C0-S_i<=.025.
The .035 endpoint gaps and acute monotonicity imply

    .7<theta_i<.96.                                (Locate)

For z=gamma-theta>=0 define exactly, as in the Pro source,

    L_t0(theta,z)=sin^2 theta[cot(theta/2)-cot((theta+z)/2)]-t0 z.

On the entire allowed apex interval,

    (L_t0)_z>=1-cos theta-t0>1-.77-.06>.17>0.

Here cos(.7)<.77 follows from its fourth-degree upper cosine
polynomial, and t0<.06. On 0<=z<=.1, the explicit elementary
ratio bound in the acquired proof gives

    sin theta/sin((theta+.1)/2)>=3871/2500,
    (L_t0)_z>=.5(3871/2500)^2-t0>1.               (Slope)

The ratio bound uses cos(theta/2)>=553/625,
cot(theta/2)<=20/7, sin(.05)<=.05, cos(.05)<=1;
it is unchanged by the deficit threshold. The final strict rational
slope is a frozen exact source guard, not a newly executed guard.

If z>=.1, integrating (Slope) to .1 and then the positive whole-apex
slope gives L_t0>.1, incompatible with a face deficit<=.025.
Thus z_i<.1 and L_t0(theta_i,z_i)>=z_i. The score identity is

    C0-S_i=C0-G0(theta_i)+L_t0(theta_i,z_i)+sigma_i+tau_i.

The actual sigma_i,tau_i are nonnegative, including actual A_i<0.
Curvature in (Scalar), integrated from the same interior theta*,
therefore yields the simultaneous budget

    2 sum_i(theta_i-theta*)^2+sum_i z_i
         +sum_i sigma_i+sum_i tau_i<=D0<=1/40.     (PBudget)

Every step above depends only on the complete scalar support, actual
positive-face formula, circumradius bound and endpoint lens. It uses
no symmetry and no independent charge of anchor tau, Q, BI or variance.
Its normalized formulas are radius-independent; the strong geometry
theorem (Pro), however, is currently invoked only at its accepted R22.

## 3. Support-plane comparison and bootstrap

Evaluation at the central peak gives

    C(lambda,t)>=C0-4u*x*^2+v theta*.

Exact score subtraction in the original unchanged actual fan yields

    D(lambda,t)>=D0+4u sum_i(x_i^2-x*^2)
                         -v sum_i(gamma_i-theta*).  (PCompare)

The common central cosine enclosure gives
.541671<x*^2<.543024, hence |x_i^2-x*^2|<.544 over the entire
base domain. Also 0<gamma_i<pi and .828<theta*<.829, so
pi<22/7 implies |gamma_i-theta*|<2.315. Thus before localization,

    D(lambda,t)>=D0-4.352a-4.63b.                 (PGross)

For (PRect) its maximum gross loss is

    4.352(.0003)+4.63(.0002)=.0022316<.0025.

If D0>=.025, (PGross) gives

    D(lambda,t)>=.025-.0022316=.0227684>dP,

so (PT) already holds. Only in the remaining branch D0<.025
do we invoke the newly rederived (PBudget).

## 4. One-budget localized comparison and strict endpoints

Put e_i=theta_i-theta*. The derivative bound
|(sin^2 theta)'|<=1 and gamma_i=theta_i+z_i give

    sum_i |x_i^2-x*^2|<=sum_i|e_i|,
    sum_i |gamma_i-theta*|<=sum_i|e_i|+sum_i z_i.

Cauchy against (PBudget) supplies

    sum_i |e_i|<=sqrt(D0), sum_i z_i<=D0.

Consequently

    D(lambda,t)>=(1-b)D0-(4a+b)sqrt(D0)=fP(D0).    (PTransport)

This bounds a score perturbation with the same D0 and leaves actual
sigma and tau in that original single budget. It does not spend
the two-face geometric loss twice or create an independent tau loss.

For d>=dP and d<.025, sqrt d>=.15, q=4a+b<=.0014,
and b<=.0002. In particular

    fP'(d)=1-b-q/(2sqrt d)
              >=1-.0002-.0014/.3>0.

By the accepted strict (Pro), D0>dP. Thus fP(D0)>fP(dP), and
sqrt(dP)=3/20 EXACTLY. It follows that

    fP(dP)=dP-.15(4a+b)-b dP
               =dP-.6a-.1725b.

This proves (PT) including every closed parameter endpoint. At u=v=0
the transport loss vanishes and all strictness comes from (Pro).
Zero tau, skew, z or parameter increments do not require strict
zero-product inequalities. No middle-halfarc variance is paid.

## 5. Explicit useful parameter rectangles and frozen-target comparison

The wider lambda rectangle (PRect) is exactly

    1397/10000<=lambda<=1403/10000,
    1067883/20000000<=t<=1075883/20000000,

with uniform fee

    D(lambda,t)>44571/2000000=.0222855.

The balanced rectangle a,b<=1/5000 is exactly

    699/5000<=lambda<=701/5000,
    1067883/20000000<=t<=1075883/20000000,

with the stronger uniform fee

    D(lambda,t)>44691/2000000=.0223455.

The earlier fixed parent target

    delta*=178523206195158271774904671901
                  /8011730250000000000000000000000

is strictly below .022283, by the same exact paper cross-product

    22283(8011730250000000000000000000000)
       -1000000(178523206195158271774904671901)
      =2178965591728225095328099000000>0.

Hence the wider lambda fee exceeds that frozen target by MORE THAN
.0000025, and the balanced fee by MORE THAN .0000625.
This is a local mathematical comparison only. At another lambda,t,
the actual complete global deficit requirement and conjugate budget
must be recomputed and matched; no changed global target is certified
by keeping delta* fixed here.

## 6. Supplementary exact conjugate control and future composition

Using the common central scalar endpoint gaps and peak curvature,
the same completed-square derivation as in the separate local source
gives, for the entire (PRect), a unique new acute support peak and

    |theta_new-theta*|<=q/4<=7/20000=.00035,
    0<=V:=C(lambda,t)-[C0-4u*x*^2+v theta*]
          <=q^2/8<=49/200000000=.000000245.

For completeness the required branch/localization test is
P=4a+(22/7)b<.00183, so .038-2P>.03434 and .9-2P>.89634,
using the frozen common scalar support certificate. The new acute
curvature satisfies -G(lambda,t)''>4-8a>=3.9976. At its new peak,
|G0'|<=q and the central curvature is >4, giving the displayed
movement. Finally G0(theta)<=C0-2(theta-theta*)^2 and the parameter
increment away from theta* has magnitude at most q|theta-theta*|;
maximizing -2e^2+q|e| gives exactly q^2/8. These supplementary
scalar statements use no local .0223 geometric fee.

The score-transport implication itself is radius-independent:
whenever a separately proved central actual-fan theorem supplies
D0>dP at a specified R with a specified actual geometric domain,
and its full scalar score/support plus (PBudget) are the ones here,
Sections 3-4 yield the same parameter fee in that unchanged domain.
Composing with any future radius theorem requires its full radius
and geometric hypotheses to be established first. In particular this
source does NOT declare the present outer squared norm threshold .83
valid up to R=.116, or acquire the middle cap by independent reflection.

The larger central margin .0225-delta* permits wider useful parameters
than the local .0223 transport. Neither rectangle is a physical
limitation on geometry or exact optimality. Actual same-fan admission,
its norm cap, outer norm threshold, global profile incidence and
once-only allocations remain separate obligations.

Author mathematical executions/imports/checker calls/replays/scans/
optimizers: 0. Original Pro bytes, old local proof, accepted reviews
and the earlier one-run exact guard receipt remain unchanged.
