# Exact new guards and fixed-radius scope

At `D0=19/500`, `L=7/10`, `H=24/25`, `zcut=1/10`, the old analytic
mechanism is reused with newly checked rational guards. The author ran
the new verifier exactly once, with no failed attempts; the receipt is
`attempts/author_001/receipt.json`. The original theorem/checker was not
run. The accepted lambda-14 scalar checker was not run.

The new guards use the accepted full-apex scalar certificate at exactly
`lambda=7/50`, `t=1071883/20000000`. Only the pinned rational
interval/Taylor definitions are imported. This proof needs concavity
on `[0,pi/3]` and the accepted decreasing tail, never concavity on the
entire acute endpoint interval.

The two localization gaps `C-G(L)-D0` and `C-G(H)-D0` have strict positive
lower enclosures, respectively about `.0007265081` and `.0010151067`.
Full-apex convexity, the cut ratios, and curvature `>4` give one budget

    2 sum e_i^2 + sum z_i + sum sigma_i + sum tau_i <= D <= 19/500.

The paid apex bound is `gamma_i<499/500=.998<pi/2`; this is obtained
before the acute `r,V` estimates. The exact `V` bound remains valid:

    Vlo=-686839/1250000=-.5494712,
    Vhi=-27/250=-.108,
    Vlo+11/20=661/1250000>0.

The new skew coefficient is `43/50=.86`, because

    skew_sq=3081055/4194304,
    (43/50)^2-skew_sq=13157649/2621440000>0.

The whole-box Jacobian bounds remain

    JB2=76960399529/20000000000,
    Jb2=77583868411/20000000000.

After the simultaneous angular/skew payment, their new coefficients are

    alphaB2=121336399529/40000000000,
    alphab2=92375868411/40000000000.

The exact shared linear-spending margins are

    alphaB2-4 LB^2 D0=23081167529/40000000000>0,
    alphab2-4 Lb^2 D0=49840036411/40000000000>0.

Thus both error functions have their maxima at `W=0`. The strict squared
error guards are

    (17/50)^2-alphaB2 D0=6608408949/20000000000000>0,
    (3/10)^2-alphab2 D0=44858500191/20000000000000>0.

They give `|B-B*|<.34`, `|b-b*|<.30`. The projection guard is

    1-(30/143)^2-(97/100)^2=3085359/204490000>0.

Consequently the exact unit-anchor necessary condition is contradicted
with strictly positive margin

    lhs=993/1250=.7944,
    rhs=60886704089/100000000000,
    lhs-rhs=18553295911/100000000000>9/50.

All degenerate limits are excluded by the universal limiting face score
`<=1/2+t*pi<7/10<C-19/500`, before localized divisions. Negative signed
anchors are included in `tau`, rather than being replaced by absolute
anchor area.

This is only a fixed-radius local actual signed-A joint-score theorem.
It does not cap `area/R^2+J_R(L)/R^2-lambda(L/R)^2+t gamma`, does not
spend actual BC-lens or BI credits, and does not supply a whole-cover
topological reduction. Any global replacement must use the same signed
shoelace terms once and disjoint eligible marked-face/anchor incidences.
A uniform extension to smaller radii is a separate theorem.

A targeted static helper read the original proof and source review and
checked the changed budget mechanism without code execution. This is
supporting analysis, not the full independent review of the frozen new
proof, checker and generated certificate. The parent will own that review.
