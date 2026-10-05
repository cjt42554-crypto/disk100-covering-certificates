# Independent static review of actual-arc thresholds and middle-branch relaxation

Verdict: PASS for the complete frozen geometry interface in ANALYSIS.md,
SHA256 97283475b6176cd0d426759d73d931690a155dec69a86843d5f5a630ff6ebbdf.
The reviewing root is independent of the geometry author. The full new paper,
its source references and the accepted Pro actual-anchor geometry were read.
No mathematical program, checker, imported mathematical source, numerical
scan or old guard replay was used in this review.

The arbitrary-center endpoint average gives the necessary quadratic
rho^2-2rho cos(a)+1-R^2<=0. Both positive roots are retained. A nonempty
covered unit-circle set is a minor arc because its cosine threshold is
positive; AM-GM also bounds that threshold below sqrt(1-R^2), proving the
stated half-angle upper bound asin(R). Equality at either radial root proves
sharpness of the necessary lower norm interface, without imposing radial
placement on the original center. A connected arc crossing the angle cut
can be lifted inside this minor arc. A total of disconnected components
cannot replace the connected length premise.

The derivative in a is positive because cos(a)>sqrt(R^2-sin^2(a)); the
derivative in R is negative. Squaring the lower-root comparison is valid
only with sqrt(v)<cos(a). The paper states this condition and separately
handles infeasible arcs. The inverse full-length expression requires
1-sqrt(v)<R<sqrt(1-v), also stated explicitly. Beyond Rcrit, choosing a
slightly longer feasible radial arc produces the stated one-disk interface
obstruction; it is not a whole-cover or fan counterexample.

I checked the square brackets for sqrt(.83), the alternating P4/P6 cosine
bracket at .0735, their product directions, and the .113281/.113282 square
comparisons. Thus fixed full arc>.147 implies outer squared norm>.83
through R<=.11328; it cannot alone imply that floor at .1133. The three
sufficient arc lengths .1471/.1478/.1561 for radii .1133/.1135/.116 follow
from upper cosine and upper sqrt(.83) bounds, with strictly positive
rational margins. The separate .156 insufficient comparison has the
opposite bounding directions correct. None of these convenient rational
thresholds is claimed as the exact sharp critical length.

The new weaker norm floor .91075 for full arcs>.147 and R<=.1135 follows
from the cubic sine lower bound, cos(a)>=1-a^2/2, and the strictly positive
square guard .007490707763765625-.0074902851=.000000422663765625.
This is a new radius interval inequality, not a replay of the old R22 cut.
Its squared floor .8294655625 is weaker than .83 and cannot substitute
into the Pro fee without another mean-norm argument.

For the Pro low-deficit domain, B<=1 is a lens consequence independent of
any middle cap. The P6(73/50)-.11 numerator is exactly 783848711 over
11250000000000; alternating remainder and monotonic cosine give cos(alpha)
>.11 throughout the stated domain. The exact lower n_x>6051/7280 and
k dot n>47661/728000 follow without an absolute-value normal relaxation.
For R<=.116 the root u>.993 is justified by the radius square. The exact
product 2(.993)(47661/728000)=47327373/364000000 exceeds .13.
The negative unit-anchor branch therefore has squared middle norm
>1+R(.13-R), exceeding 1+.00158508 throughout [.11322,.116]. Hence
the SAME actual middle squared cap <=1+3/2000, including norm<=1,
selects the positive branch. T=0 retains strictness from the angle guard;
no strict zero-product step is needed. No original fan vertex is moved.

At fixed R22, replacing only the original Pro middle cap by this new cap
is justified because that hypothesis selects the branch and does not
enter its later direction, Hessian, or scalar closure. For other radii
this review accepts only the branch and arc interfaces. The full radius
fee requires its separate radius proof and independent review. This review
does not accept new lambda/t localization, any global fan transport,
profile matching, full-cover exclusion or exact optimum.

All reviewed predecessor bytes and histories remain unchanged. New
source byte pins accompany this verdict. Mathematical execution/import/
checker/replay/scan count for this independent review is zero.
