# Continuous parameter transport of the general actual-tau q2 fee

Status: new analytic author proof, awaiting independent static review.
This source preserves all accepted predecessor bytes. No mathematical
program, mathematical import, checker, numerical replay, optimizer,
scan or old 48-guard batch was executed.

## 1. Exact geometric domain and matched support

The radius is exactly

    R=5661/50000=.11322.

Let (p,c0,c1) and (p,c1,c2) be two consecutive consistently oriented
positive-area triangles of circumradius at most R. For each directed
base c_i c_(i+1), the actual anchor q_i is a unit point covered by
both endpoint disks. Define

    x_i=|c_(i+1)-c_i|/(2R), theta_i=asin x_i,
    gamma_i=the actual positive apex angle at p,
    A_i=-det(c_(i+1)-c_i,q_i-c_i)/2,
    S_i(lambda,t)=area_i/R^2+A_i/R^2-4lambda*x_i^2+t*gamma_i.

The signed A_i may be negative. Actual lens slack
tau_i=x_i sqrt(1-x_i^2)-A_i/R^2 is nonnegative and arbitrary.
The SAME actual middle center satisfies

    |c1|<=sqrt(1-R^2),

and each ORIGINAL actual outer center covers a connected assigned
unit-circle arc of full length strictly greater than 147/1000.
There is no independent movement or reflection of c1, no additional
angle or base-x restriction, no symmetry, global Q, profile, degree,
or middle-variance premise.

For every lambda,t below, C(lambda,t) is the complete independent
full-apex support for that SAME objective:

    C(lambda,t)=sup F(lambda,t;x,gamma),
    F=x^2 cot(gamma/2)+x sqrt(1-x^2)-4lambda*x^2+t gamma,
    0<=x<=1, asin x<=gamma<=pi-asin x,

with limiting endpoints understood by continuous branch closure.
Equivalently, by the convex apex endpoint reduction,

    C(lambda,t)=max_{0<=theta<=pi/2} {G(lambda,t;theta),Bo(lambda,t;theta)},
    G=sin theta+sin(2theta)-4lambda sin^2 theta+t theta,
    Bo=sin theta-4lambda sin^2 theta+t(pi-theta).

The convexity F_gamma_gamma>0 is independent of lambda,t, so this
definition needs no new numerical conjugate evaluation. Every actual
S_i is at most the matched C(lambda,t). The theorem does not replace
C(lambda,t) by the old support value or a numerically approximated
support.

Write

    lambda0=7/50, t0=1071883/20000000,
    a=|lambda-lambda0|, b=|t-t0|,
    D(lambda,t)=2C(lambda,t)-S_0(lambda,t)-S_1(lambda,t),
    d*=223/10000=.0223,
    L(a,b)=d*-(3/5)a-(1723/10000)b.

For the complete CLOSED parameter rectangle

    a<=1/2500=.0004, b<=1/5000=.0002,                (Rect)

the reusable local conditional conclusion is

    D(lambda,t)>L(a,b)>=1101277/50000000=.02202554.   (T)

All hypotheses concern the original actual fan, and tau is included
once in the signed-A score. This is not a global radius exclusion.
No radius other than the specified R is asserted in this source.

## 2. Exact central inputs used without regeneration

The accepted full actual-tau theorem at lambda0,t0 supplies

    D0=2C0-S_0(lambda0,t0)-S_1(lambda0,t0)>d*.        (Central)

This holds for arbitrary actual tau under exactly Section 1's geometry.
Its proof and complete independent audit are pinned as dependencies.
Only this two-face theorem is used as a geometric exclusion; no
three-face compatibility conclusion is transferred to q2.

The separately accepted per-face scalar deficit mechanism states:
if D0<=1/40, each central face deficit is nonnegative and

    .7<theta_i<.96, gamma_i=theta_i+z_i, z_i>=0,
    2 sum_i e_i^2+sum_i z_i+sum_i sigma_i+sum_i tau_i<=D0,
    e_i=theta_i-theta*,                              (Budget)

where theta* is the unique central acute support peak. Its frozen
certificate gives

    .828<theta*<.829, .676<cos theta*<.677,
    -G0''>4 on [.7,.96].

It also gives C0-G0(.7)>1/40+.013 and
C0-G0(.96)>1/40+.013, both greater than .038, and an obtuse
support gap C0-Co>.9. Central G0 is increasing before its unique
peak and decreasing after it over its full acute branch. These are
pre-existing exact premises, not newly executed guards.

## 3. A centered support comparison before localization

Put signed increments u=lambda-lambda0 and v=t-t0, and
x*=sin theta*. Evaluating the new full support at the old peak gives

    C(lambda,t)>=C0-4u*x*^2+v*theta*.

Since the actual fan and anchors have not changed, exact score
subtraction therefore yields

    D(lambda,t)>=D0+4u sum_i(x_i^2-x*^2)
                         -v sum_i(gamma_i-theta*).   (Compare)

The accepted cosine bracket implies

    .541671<x*^2<.543024,
    |x_i^2-x*^2|<.544  for every 0<=x_i<=1.

Every positive triangle has 0<gamma_i<pi. The elementary bound
pi<22/7 and theta* in (.828,.829) imply

    |gamma_i-theta*|<2.315.

Consequently the uniformly valid global comparison is

    D(lambda,t)>=D0-4.352a-4.63b.                  (Gross)

The loss is at most

    (4.352)(.0004)+(4.63)(.0002)=.0026668

on (Rect). If D0>=1/40, (Gross) gives

    D(lambda,t)>=.025-.0026668=.0223332>d*>=L(a,b).

Thus (T) already holds for that entire branch. This step supplies
the bootstrap BEFORE the paid central localization is used. In the
only remaining branch, D0<1/40, (Budget) is available.

## 4. Localized perturbation uses one central deficit

The function sin^2 theta has derivative sin(2theta), of absolute
value at most one. Therefore

    sum_i |x_i^2-x*^2|<=sum_i |e_i|,
    sum_i |gamma_i-theta*|<=sum_i |e_i|+sum_i z_i.

From the SAME (Budget), Cauchy gives

    sum_i |e_i|<=sqrt(2 sum_i e_i^2)<=sqrt(D0),
    sum_i z_i<=D0.

Applying these two estimates in (Compare) gives the correct
simultaneous comparison

    D(lambda,t)>=(1-b)D0-(4a+b)sqrt(D0)=f(D0).    (Paid)

It introduces no extra geometric score credit. In particular tau
and sigma remain in the original single D0; they are not recharged,
reset, or converted to a different loss. The estimates on e and z
are used to bound a score perturbation, not to claim two independent
copies of (Budget).

On d<=D0<1/40 with d>=d*, sqrt(d)>.1, and on (Rect)
4a+b<=.0018, b<=.0002. Thus

    f'(d)=1-b-(4a+b)/(2sqrt d)
          >1-.0002-.0018/.2=.9908>0.

The accepted (Central) has strict D0>d*. Hence f(D0)>f(d*).
Since d*<(.15)^2,

    f(d*)=(1-b)d*-(4a+b)sqrt(d*)
           >=d*-.15(4a+b)-b*d*
           =d*-.6a-.1723b=L(a,b).

This proves the strict (T), including parameter-rectangle boundaries.
If a=b=0, the last comparison is equality and the needed strictness
comes solely from the accepted Central theorem. Zero tau, zero skew,
zero angle perturbation, and zero individual parameter increments
do not require a strict zero-product comparison.

## 5. Two explicit rectangles and one exact target comparison

The broad rectangle is exactly

    349/2500<=lambda<=351/2500,
    1067883/20000000<=t<=1075883/20000000,

and has uniform fee

    1101277/50000000=.02202554>.0220.

The narrower rectangle is exactly

    6999/50000<=lambda<=7001/50000,
    1071483/20000000<=t<=1072283/20000000.

Here a,b<=1/50000 and (T) gives

    D(lambda,t)>11142277/500000000=.022284554.     (Narrow)

For comparison only, the previously specified fixed parent target was

    delta*=178523206195158271774904671901
                 /8011730250000000000000000000000.

It is strictly below .022283: writing its numerator N and denominator
H, the exact paper cross-product difference is

    22283 H-1000000 N
       =2178965591728225095328099000000>0.

Thus the narrow fee exceeds delta* by more than

    .022284554-.022283=.000001554.

This comparison does not certify that a changed global lambda,t
has the same parent target. Its complete global signed budget,
profile incidence, and coefficient matching must be recomputed.
Any separately allocated nonnegative middle variance can increase
the established local deficit, but is not used in this proof.

## 6. Peak movement and support residual for the entire broad rectangle

These scalar statements are supplementary outputs, not assumptions
needed for Sections 3-4. They give explicit conjugate stability and
identify why the geometric argument was not rerun.

Let P=4a+(22/7)b. Uniformly over both full support branches,
|G(lambda,t)-G0|<=P and |Bo(lambda,t)-Bo0|<=P. On (Rect),

    P<=.0016+(22/7)(.0002)<.00223,
    .038-2P>.03354>0, .9-2P>.89554>0.

The old monotonicity and endpoint gaps imply every new acute
support maximizer lies in (.7,.96); the obtuse branch cannot win.
On that window,

    -G(lambda,t)''>4-8a>=3.9968>0.

The new support thus has a unique acute interior peak theta_new.
At it, |G0'(theta_new)|=|4u sin(2theta_new)-v|<=q, q=4a+b.
The old -G0''>4 and G0'(theta*)=0 imply

    |theta_new-theta*|<=q/4<=9/20000=.00045,

with equality interpreted nonstrictly at q=0.

The old curvature also gives, on the entire window,

    G0(theta)<=C0-2(theta-theta*)^2.

The parameter increment relative to its old-peak value is at most
q|theta-theta*|. Completing the square therefore bounds the exact
support residual

    V=C(lambda,t)-[C0-4u*x*^2+v*theta*],
    0<=V<=q^2/8<=81/200000000=.000000405.        (Residual)

The actual deficit identity may consequently be written exactly as

    D(lambda,t)=D0+4u sum(x_i^2-x*^2)
                        -v sum(gamma_i-theta*)+2V.

Only V>=0 was needed for the fee proof. Its upper bound is useful
for a later matched global dual budget; it does not evaluate or
replace the exact C(lambda,t). No new peak polynomial, 48 guards,
or pointwise numerical optimization has been executed.

## 7. Actual limits of this transport lemma

The broad rectangle's worst fee loss is .00027446. It is mainly
linear score transport, not loss of scalar concavity: new curvature
remains above 3.9968 and peak movement is at most .00045. The accepted
central geometry itself has a path/envelope coefficient gap
.30-.285=.015 and a final mean-norm contradiction gap exceeding
.829921-.8294023684=.0005186316. This transport proof does not
spend those geometric margins to raise the central fee.

Preserving the old fixed delta* with this particular method is
restricted by the small d*-delta* margin, about 17.27e-6, and the
explicit payment .6a+.1723b. The narrow rectangle resolves that
restriction. This does not prove that either the geometry or the
exact optimum is physically limited to that rectangle.

The geometric prerequisites remain decisive: admitted actual common
unit anchors, consecutive positive faces of circumradius at most R,
the cap on the SAME actual c1, and both ORIGINAL outer full arcs>.147.
The parameters change only the objective, so none of those geometric
conditions is weakened or inferred by this lemma. Radius changes,
alternative arc/norm thresholds, global incidence, once-only marked
face allocation and exact optimality all require separate proofs.

Author mathematical program executions, source imports, old checker
replays, scans and optimizer calls: 0. This paper derivation awaits
an independent full static audit before being described as accepted.
