# Independent static review of the relaxed-middle radius corollary

Reviewer: `/root/q2_uniform_review`. Verdict: **STATIC PASS** for
`COROLLARY_RELAXED_MIDDLE.md`, SHA256
`768a6e71756ad155331b2196b68394548f2eadb4affeb1352e210ef93d766e28`.
This is a new scope review independent of the author. It does not
automatically inherit a relaxed cap from the original single-point
audit or from a separate author's geometry note. Mathematical
executions, mathematical imports, checker calls, numerical replays,
scans and optimizers: **0**. Accepted predecessor bytes are unchanged.

The accepted conclusion is D>9/400 on the entire CLOSED radius band
[5661/50000,29/250], with fixed lambda=7/50 and t=1071883/20000000,
same actual consecutive positive faces of circumradius at most R,
actual covered unit anchors and arbitrary signed A/tau. The SAME
actual middle squared norm may now be at most 1+3/2000, while both
actual outer squared norms remain at least .83. In particular the
ordinary premise |c1|<=1 suffices. Alternatively, each SAME actual
outer disk covering a connected FULL unit-circle arc strictly longer
than .1561 suffices for the outer norm premises on this radius band.
No global covering exclusion, incidence rule or parameter rectangle
is established by this corollary.

## New branch selection for the changed cap

I read the frozen radius main and the full frozen geometry analysis.
I independently rederived the portions of that analysis used in this
composition. Low-D support/localization and endpoint-lens coverage
do not require any norm cap on c1. They give alpha<=11661/8000<1.46,
sin(alpha)>6/7, |g|,|d_b|<=27/1040 and |b|<=1. The last bound
comes from averaging the two covered anchor-minus-middle offsets:
|k|^2+|b|^2<=1. It therefore remains available before selecting a
translation sign, and along the relative lens path.

The exact sixth-order lower cosine guard at 73/50 was independently
checked by paper integer arithmetic:

    P6(73/50)-11/100=783848711/11250000000000>0.

Since the entire alpha interval is in (0,pi/2), cosine decreases
there and cos(alpha)>.11. With n=(b_y,-b_x)/|b|, the positive
numerator and |b|<=1 give

    n_x>6/7-27/1040=6051/7280,
    k dot n>(11/100)(6051/7280)-27/1040
             =47661/728000.

For R<=.116, u^2>=1-R^2>.993^2. Multiplication gives exactly

    2u(k dot n)>47327373/364000000>13/100.

The last comparison has the positive integer numerator margin
7373 at denominator 364000000. It is a genuine strict margin,
not a rounded numerical equality. The exact negative branch of the
two-unit-anchor translation therefore satisfies

    |c1|^2=u^2+2Ru(k dot n)+R^2|k|^2
          >1-R^2+.13R=1+R(.13-R).

Both factors in R(.13-R) are positive on the full radius band. The
coarse lower product R0(.13-.116)=39627/25000000=.00158508
exceeds .0015. Thus that branch contradicts the NEW cap
|c1|^2<=1.0015 and the actual branch must be positive. The argument
includes T=0; its strictness comes from the cosine and branch margins.
It concerns the same original scored c1 and performs no radial move,
reflection or replacement of the fan.

I checked the actual premise use in the radius main: its middle cap
is used solely to select the positive actual branch. The interpolated
positive formula does not need a cap along the path. Every subsequent
signed linear support, slab-transfer, Hessian/remainder, changing
baseline and mean-norm contradiction uses only the original angle/lens
domain, the established branch, and the two outer norm floors. The
new sign proof consequently replaces the old cap use completely.
No radius-dependent part of the main proof is weakened by this change.

## Optional full-arc premise for arbitrary actual centers

For an arbitrary center of norm rho covering a connected minor arc
with half-angle a, averaging endpoint squared distances yields
1+rho^2-2rho cos(a)<=R^2. No radial alignment is assumed. Its
necessary lower root is

    rho>=cos(a)-sqrt(R^2-sin^2(a)).

Both positive radial roots are retained. The lower root strictly
increases with a in the feasible interior and decreases with R;
continuity covers a=asin(R). An infeasible fixed half-angle simply
excludes the asserted longer covered arc. For u=sqrt(.83)<cos(a0),
the lower root at a0 exceeds u if
R^2<1.83-2u cos(a0). The square of cos(a0)-u is positive,
so this equivalence has no invalid squaring branch.

At a0=.07805 the source's fourth-order cosine bound gives
cos(a0)<.996956. Its norm root guard gives u<.91104336, and
cos(a0)>=1-a0^2/2>.99>u verifies the required branch. I checked
the rational product:

    1.83-2(.91104336)(.996956)=.01345971197568
                                          >.013456=(.116)^2.

Thus every R in the radius band and every covered connected FULL
arc longer than .1561 force the same center's squared norm above
.83. This proof concerns a full connected arc, not a halfarc, an
assigned-base sine, or a sum of disconnected arc components. Its
strict length input is stronger than the parent's old .147 premise;
the corollary does not assert that the global chain supplies it.

## Interface bottleneck and composition scope

For a fixed full length .147 the exact sharp one-disk interface is
Rcrit=sqrt(1.83-2sqrt(.83)cos(.0735)). The paper root and cosine
brackets give .113281<Rcrit<.113282. I checked their polynomial
and rational directions. Past Rcrit a radial center covering an arc
slightly longer than .147 can have squared norm below .83. This is
an obstruction to this arc-to-norm INTERFACE, not a global covering
or scored-fan counterexample. The sufficient lengths .1471, .1478
and .1561 correspond to separately checked upper cosine bounds and
cover radii .1133, .1135 and .116, respectively; none is asserted sharp.

The original middle cap is therefore removable within this new
low-D domain, while the outer floor is a material premise. The fee,
radius range, signed actual-anchor semantics and scalar support
remain those already independently accepted in PRO_RADIUS_REVIEW.
Changing lambda/t, D, outer norm floor or anchor incidence requires
its own proof. Actual global fan admission, full-cover arc bounds,
profile exclusions and once-only allocation remain external.

The corollary and its separate pin record, frozen main, required
geometry source and independent main review are checked by metadata
SHA256/byte count in PRO_RELAXED_CAP_SOURCE_PINS.json. No unresolved
mathematical defect was found in the frozen composition.
