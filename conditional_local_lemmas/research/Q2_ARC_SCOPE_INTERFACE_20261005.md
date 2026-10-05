# Two-face actual-anchor / assigned-arc interface

The parent supplies R=5661/50000,lambda=7/50,t0=1071883/20000000.
After its already-paid q2,d7 type slack, the required extra normalized
fee must be strictly larger than

    178523206195158271774904671901 /
    8011730250000000000000000000000.

This is a conditional integration target. The parent reports that the
R=.11322 profile reduction to a unique degree7 q2/q3 apex still depends
on the unproved Q<=.35754 radial premise. At that radius the actual
assigned full-arc lower bound presently established is a>.0597.
Neither xi>.6 nor the radial premise is treated as established here.

Use distinct symbols throughout:

    x_i=|c_(i+1)-c_i|/(2R)              (triangle half-base),
    xi_j=sin(a_j/2)/R                  (assigned-cell half-chord).

The three-face certificate uses x_i. The parent's proposed stronger
arc domain uses xi_j. Equal-score two-face examples can have x_i>.6
while their outer xi_j are about .55, so these conditions cannot be
substituted for one another.

## A precise necessary norm condition for a covered assigned arc

Let a=2beta be a single covered arc with endpoints
w_plus/minus=cos(beta)u +/- sin(beta)v. For the present small-arc domain
0<=beta<pi/2. Let c have rho=|c|, and suppose both endpoints are within
R of c. Their mean squared distance gives

    1+rho^2-2cos(beta)(c dot u)<=R^2.

Since c dot u<=rho,

    (rho-cos(beta))^2<=R^2-sin(beta)^2.

Therefore the assigned arc's xi=sin(beta)/R satisfies 0<=xi<=1 and

    rho>=rho_min(xi),
    rho_min(xi)=sqrt(1-R^2 xi^2)-R sqrt(1-xi^2).

For 0<R<1 this function is strictly increasing on 0<xi<1: its
derivative is R xi[1/sqrt(1-xi^2)-R/sqrt(1-R^2xi^2)]>0, because
1-R^2xi^2>R^2(1-xi^2). This proof only uses actual endpoint coverage.
It does not force radial placement, equal arcs, convexity, or nearest
ownership, and incurs no lens or score payment.

For xi>=3/5 the exact necessary bound is

    rho>=sqrt(1-(9/25)R^2)-(4/5)R.

The elementary norm bound is necessary, but actual arc endpoints and
their consistent cyclic placement must still be verified in any
counterexample or sufficient construction. It cannot itself establish
a two-face score fee, much less a whole-cover exclusion.

## Current research status

The original zero-skew two-face saturation is agreed by both parent
and local proofs, so it was not replayed here. New work instead uses
mirror skews to increase the outer covered arcs while paying a small
two-face score deficit. A fixed candidate uses sin(epsilon)=1/12 and
cos(epsilon)=sqrt(143)/12 at the accepted acute peak theta0. Its
prospective exact deficit is sin(theta0)/36. Acceptance requires the
new explicit algebraic arc calculation and independent static audit.

No claim that this candidate extends to a 100-circle cover or a
degree7 complete fan is made. Four-cell window, radial Q, full-cover
topology, other supports, and disjoint payments remain separate. In
particular, a local arc counterexample would refute a claimed local
fee on that stated domain, without refuting the existing global bound.
