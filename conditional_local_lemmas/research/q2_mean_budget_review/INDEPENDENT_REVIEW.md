# Independent static review of the paid ideal mean-norm envelope

Verdict: **PASS for the shape-budget envelope in Sections 1–3**, with
the separately preserved ENDPOINT_CLARIFICATION. This review accepts
the following statements at `R=5661/50000`, under the accepted actual
per-face high-score localization, `D<=223/10000`, `W=tau0+tau1`,
`B=D-W>=0`, and at the positive common unit-tip branch:

    4E^2+Delta^2+(1666/625)sin^2 eta<=B,
    chi_eff<=.523+sqrt[(47/100)B],
    M(0)<=N(B)=1+R^2-2R cos(.523+sqrt[(47/100)B]),
    N'(B)>.30 for 0<B<=.0223.

The envelope is a mean-square estimate for the ideal unit-tip endpoint.
It does not by itself estimate the actual positive-tau endpoint. This
review does NOT accept the still-open actual-path bound `M'(t)<=.30W`,
its proposed support-branch signs, a general positive-tau strong fee,
or a global covering/radius exclusion. The conditional integration
stated in the source is correct if that additional uniform actual-path
bound is separately proved.

I read the complete MEAN_BUDGET_ENVELOPE, ENVELOPE_PINS, separate
endpoint clarification, original working interface, ledger and original
working manifest. The frozen zero-tau proof and its independently
accepted full-source audit are used as pinned dependencies, including
its separately recorded branch, peak-constant and endpoint clarifications.
The accepted scalar support/peak constants are pre-existing premises;
they were not regenerated. No mathematical source execution, import,
checker, numerical replay, scan or optimizer was performed. Reading
the working interface does not broaden this envelope verdict to all
of its path/derivative claims.

One supplemental reviewer independently checked the full envelope and
ENVELOPE_PINS, derived the same weighted skew and single-Delta payment,
and confirmed the N derivative guards. It also identified the same
zero square-root issue and the need to state the trivial B=0 maximum
case separately. Its scope and supplied summary are preserved in
SUPPLEMENTAL_STATIC_CHECK; it performed no execution or file mutation.

## One high-score payment and the weighted common skew

The accepted actual deficit budget includes W exactly once. Removing
paid apex slack with its simultaneous skew estimate, rather than
discarding z while retaining its skew credit, gives

    4E^2+Delta^2+2 sum x_i sin(d_i)^2<=D-W=B.

This is the same budget as the source's equation (1). The actual
positive-face localization retains `.75<s<.904`, `|Delta|<.15`
and the acute small-skew branch from the zero-tau review. Since
each x_i>.64 and
`sin^2(m+eta)+sin^2(m-eta)>=2sin^2 eta`, the budget implies
`2.56sin^2 eta<=D`. The explicit lower sine polynomial at .095
exceeds .094 and `2.56*(.094)^2=.02262016>.0223`, so |eta|<.095.
The inverse-sine step uses the already localized acute eta branch.

Put X=x0+x1>0 and Y=x0-x1. Expanding the two squared sines gives
exactly

    X-[X cos(2eta)cos(2m)-Y sin(2eta)sin(2m)].

The sinusoidal expression in m is bounded by its Euclidean amplitude
for EVERY m. This is a relaxation, not an assumed symmetric or
realizable minimizing fan. Writing that amplitude as
`X sqrt(1-u)`, with
`u=4x0x1 sin^2(2eta)/X^2 in [0,1]`, proves equation (2) using
`sqrt(1-u)<=1-u/2`. Thus no absolute-value common skew or separate
optimizing budget has been introduced.

The trigonometric product identities for
`x0=sin(s+Delta/2)`, `x1=sin(s-Delta/2)` give exactly

    8x0x1/X=4 sin s cos(delta)[1-cot^2 s tan^2 delta].

The stated uniform guards hold. For example
`cos(.075)>=1-(.075)^2/2=.9971875>.997`, and division of
`sin(.075)<.075` by that lower cosine gives `tan(.075)<.076`.
The previously accepted seventh-order lower sine at .75 exceeds .68;
the fourth-order upper cosine there is `.73193359375<.74`.
Hence `cot s<.74/.68<1.1`. The exact guarded weight factor is

    .997*[1-(1.1)^2(.076)^2]=.99003200688>.99.

Also `cos^2 eta=1-sin^2 eta>1-(.095)^2=.990975>.99`.
Their product exceeds `.9801>.98`. Multiplication by `4sin s>2.72`
therefore proves the nonstrict uniform coefficient
`4*.68*.98=2.6656=1666/625`, including eta=0. The remaining
angle terms are retained in the SAME B.

## Mean asymmetry is absorbed by its existing Delta payment

Removing Z increases chi, while replacing negative eta by zero
increases chi0 and lowers the payment. All derivative estimates
remain valid on the larger rectangle with `eta+=max(eta,0)<.095`:
`alpha0=3s/2-eta+>1.03` and `<1.356<pi/2`. Thus chi0_s<.6,
0<chi0_eta<1 over the entire connecting segment from the accepted
peak point. The lower sine polynomial gives
`eta+<=1.002sin eta+`, including its zero endpoint. The peak starting
value is strictly below .523, with both required cosine brackets
already explicitly supplied in the zero-tau clarifications.

Cauchy against the shape/skew payment remaining after Delta^2
uses the exact coefficient

    .6^2/4+1.002^2/(1666/625)=310977/666400<47/100.

The source's rational comparison is correct because
`.38*2.6656=1.012928>1.004004`. I flagged that the literal strict
square-root coefficient comparison at `H=B-Delta^2=0` is false.
The author retained the frozen main and added ENDPOINT_CLARIFICATION:
the uniformly valid chain has a strict FIRST step from the strictly
smaller peak starting value, followed by the nonstrict coefficient
comparison. This proves the desired strict scalar bound even at H=0;
no later argument requires strictness from the zero product.

At the positive unit-tip branch, the exact previously derived outer
mean-square identity is

    M(0)=1+R^2-2R cos chi cos Delta.

The actual phase chi is positive and acute throughout the paid domain;
chi<=chi0(s,eta+) and cos Delta>0 therefore permit this upper
norm relaxation. The new chi0 is between .375 and .662. Put
v=Delta^2 and t0=3v/2<.034. For v>0, the source's proof of

    cos chi0 cos Delta>=cos(chi0+3Delta^2/2)

uses valid positive cosines and `sin t0>.999t0`,
`tan chi0>chi0>.375`. Its coefficient
`.375*.999*1.5=.5619375>.5` pays the elementary cosine loss
`1-cos Delta<=Delta^2/2`. At v=0 the displayed inequality is
equality, as the source explicitly says. This transformation retains
the benefit of the mean-norm asymmetry while paying it from the same
Delta^2 already present in B; it is not an additional budget.

For `0<=v<=B`, rationalizing the square roots gives, when v>0,

    sqrt(kB)-sqrt(k(B-v))
       =kv/[sqrt(kB)+sqrt(k(B-v))]
       >=kv/[2sqrt(kB)]>1.5v,
    k=47/100.

Indeed `sqrt(kB)<.103<k/3`; at B=v=0 the claim is read by
continuity and equality. Therefore `sqrt(k(B-v))+1.5v` is
maximized at v=0. Combining this with the preceding mean identity
proves chi_eff<=.523+sqrt(kB) and the stated N(B). No separate
optimization of skew and angle asymmetry is needed.

## Uniform derivative of the envelope and its endpoint

For B>0, direct differentiation yields

    N'(B)=Rk sin(.523+y)/y, y=sqrt(kB)<.103.

The ratio decreases with y on this full interval because its
derivative numerator is `y cos(.523+y)-sin(.523+y)<0`:
the first term is at most .103, while the second exceeds
`sin(.523)>.49` by its cubic lower guard. Hence the ratio is
bounded below by `sin(.626)/.103`. The cubic lower sine at .626
is `.585114270666...>.585`, and

    R*.47*.585=.031129839>.30*.103=.0309.

This proves the strict derivative payment on the entire positive
budget interval, not just a sampled value. N is continuous at B=0;
its infinite right derivative there causes no missing endpoint.
Integrating the lower slope from a positive lower endpoint and taking
the continuous limit handles B=0 whenever a later actual-path bound
is available.

At D0, `kD0=.010481<(.103)^2=.010609`, so its angle is below
.626. The sixth-order lower cosine polynomial at .626 exceeds .81.
A separate coarse rational hand check uses z=(.626)^2=.391876,
`z^2/24>.00639`, `z^3/720<.000089`; its polynomial is then greater
than `.804062+.00639-.000089=.810363>.81`. Consequently

    N(.0223)<.8294023684.

The independently accepted original assigned-arc norm lower bound
`.911088875` has square larger than `.911^2=.829921`, so the
conditional final norm gap is correct.

## Unaccepted actual-path step and preserved source scope

If the fixed-fan actual path is proved to start at this positive
unit-tip branch, end at the original actual fan, and obey
`M'(t)<=.30W` uniformly with the same paid W, then the source's
conditional chain follows by integration:

    M(1)<=N(D-W)+.30W<=N(D)<=N(D0).

This last implication is elementary once those additional premises
are established. The present verdict supplies ONLY N and its budget
payment. It does not supply the actual path derivative coefficient,
whole-path support-branch signs, or a transport to a globally valid
actual fan. Section 4's proposed directional cancellations are useful
pending identities/support choices, not an accepted sufficient bridge.
The separately read WORKING_INTERFACE remains an independently
unfinished interface at this stage, despite being hash-pinned.

The original main, original working interface and manifests remain
unchanged. The endpoint correction is preserved separately. All
new review records and source pins are in q2_mean_budget_review;
the earlier tau-zero and small-tau frozen audits were not altered.
