# Independent full-source static review

**Verdict: PASS for the local actual signed-anchor inequality
`T<3C-1/40` at the exact new quadratic parameters.** I read the
complete frozen main proof and checker after independently checking
the proposed analytic mechanism and fixed guards in `BLIND_REVIEW.md`.
This source-review stage ran no checker, import of an old mathematical
module, optimizer, candidate replay or new scalar guard calculation.
The author sources and author manifest were not modified. The root's
independent fresh checker replay remains a separate verification step.

## Frozen versions

| Source | SHA-256 |
| --- | --- |
| `PROOF.md` | `58dc829a24518ecacc44e3414c89973eb92b94573b72a631fc78e41d023173c1` |
| `verify_actual_quadratic.py` | `71811ec96c80d452028f9b434215a8899798e555056125b3adaab178d41bc86f` |
| `certificate.json` | `c92f9879d82e80ca134eb93379d1e842189f836c5e5603d954361c95b7f286da` |
| `replay.log` | `ee02905f57913c1dfa688d9f3faf0a0767869883ffb33b4d33e47c7029fed162` |

The independent blind review SHA-256 is
`b55f48fe9ef1f0b9bb8082856df72495767e6884c955fc5fca6f47fe38bfad02`.
The adjacent independent `SOURCE_REVIEW_SHA256.json` binds this review,
the blind review, all five author-manifest source pins, the author
manifest snapshot and the four accepted scalar dependency pins.

## 1. Exact objective and complete initial domain

The theorem uses exactly `R=11321/100000`, `lambda=7/50` and
`t=1071883/20000000`, with negative quadratic base potential
`-4lambda*x^2`. The actual signed anchor, actual triangle area and
actual apex potential are retained once. The constant `C` comes from
the already accepted full support for these same coefficients;
neither the old linear-objective fee nor its old conjugate is used.

Section 1 derives the area identity with
`sigma=2*x^2*sin(d)^2/sin(gamma)` and the signed lens constraint
directly. The two disk inequalities give `tau>=0` for either anchor
sign. Hence the exact one-face score is `F-sigma-tau` and is at most
the complete independent support `C` before any localization.
The sine-rule apex interval and convex endpoint reduction include
initially obtuse actual triangles.

The main source correctly fixes the new concavity issue identified
in the blind review: it uses concavity on `[0,pi/3]` plus the decreasing
tail, rather than false concavity on all `[0,pi/2]`. Its checker pins
the accepted conjugate proof, checker and certificate, checks the
exact coefficient fractions and accepted status, and reads that
certificate without executing the earlier conjugate verifier.

## 2. Full gamma cut and compatible deficit

For a hypothesized high-score triple, every individual deficit is
nonnegative and at most `D0`. The two new endpoint guards, the unique
acute peak and its monotone branches locate every theta in `(.7,.96)`.
Sine concavity and cosine monotonicity establish curvature greater
than 4 throughout this window. Its quadratic payment coefficient is
therefore 2, including the required factor one-half from integration.

The cut ratio has no interior minimum because its derivative
numerator has strictly negative derivative. Both rational endpoint
enclosures exceed `81/50`; this is stronger than the independent
blind proof's sufficient `3/2`. The cut stays within the full feasible
apex interval. Its value and the complete obtuse endpoint are both
strictly below the high-score threshold. Convexity then excludes
every gamma beyond the cut, rather than just a sampled obtuse or
nearby-apex subdomain.

After the cut, integration gives the same budget

    2 sum e_i^2 + sum z_i + sum sigma_i + sum tau_i <= D <= 1/40.

It further implies `gamma<.985<pi/2` before the circumradius and
cosine monotonicity estimates are used. No stronger slope or acute
assumption is used prematurely. There is no repeated angle, skew
or anchor payment.

## 3. Exact tips and signed error

The main source supplies the exact base midpoint and outward-normal
calculation for its complex tip coordinates. The common bisector
relations retain the actual shared rays and their lifted turns.
Changing the skew convention only reverses a signed imaginary term.

The outer-tip error `tau/x^2` is valid for negative signed anchor
height as well as positive height. The lower sine bound keeps its
denominator away from zero. The `U` derivative bound is global on
the localized rectangle. The `V` decomposition uses the established
acute gamma bound and the same `z<=D0` payment. Replacing the blind
review's `C0/S0` by the source's slightly larger bound `6/5` still
gives `Vlo>-.55`. The source's `cs=625/512`, `ct=625/256` and `.79`
skew bound have the correct normalization. Orthogonal real and
imaginary skew errors are combined before bounding their norm.

## 4. Full-box Jacobian and segment integration

The derivative columns stated in Section 4 agree with direct
differentiation of the three-angle acute-tip maps. The endpoint
imaginary component can change sign in this box, and the source and
checker explicitly retain the maximum of both possible magnitudes.
The middle real component has a positive whole-box lower bound;
the remaining middle components are bounded by differences or sums
of two positive-product intervals.

The squared Frobenius bounds match the independent exact arithmetic:

    JB2=76960399529/20000000000,
    Jb2=77583868411/20000000000.

These bound operator norms throughout `[L,H]^3`. The accepted peak
triple and actual triple both lie in that convex box, so the whole
connecting segment is controlled. This is a global localized map
bound, not a Hessian estimate at one ideal configuration.

## 5. Shared error spending

The actual combinations cancel the common apex exactly. Aggregating
the actual bisector errors gives the stated safe linear coefficients
`LB=4.02`, `Lb=2.645`, each dominating its anchor-error coefficient.
The squared norms of skew weights are `3/2` and `1/2`; Cauchy with
the existing angle coefficient 2 gives precisely

    alphaB2=114406399529/40000000000,
    alphab2=90065868411/40000000000.

Splitting the single budget into quadratic/skew spending and
`W=sum z+sum tau` proves both error functions in (15). The strict
guards `4Llinear^2D0<alpha^2` make them decreasing in `W`, so the
maximum is at `W=0`. For `D=0` the conclusion is immediate without
differentiation. The independently calculated squared fit margins
are respectively

    2233600471/1600000000000,
    2094131589/1600000000000.

They prove the strict `.27` and `.24` bounds simultaneously from
the same deficit. Those two consequences are not separately added
payments. No BC lens deficit occurs anywhere in this local estimate.

## 6. Unit-circle condition and degenerate closures

Section 6 now proves the three-unit-point necessary condition
self-containedly. The endpoint mean is perpendicular to their
half-chord; its norm relation and the displayed scalar square identity
give the exact factor `R` after normalization. The localized error
bound ensures a nonzero actual half-chord, so its perpendicular is
defined. Repeated anchors elsewhere in the original domain are not
excluded by assumption.

The ideal size guards, `.98` axis projection, lower left side `.9496`
and upper right side use the correct sign and scale. Their exact
independent rational difference is
`38496637031/100000000000>.38`. The strict contradiction therefore
also excludes the closed hypothesized endpoint `D=D0`.

The theorem itself assumes positive-area faces. Its limiting-degenerate
statement is also consistent. In addition to the explicitly discussed
`x=0` and `theta=pi/2` endpoints, a zero-area face with another
collapsed side cannot have high score: its normalized signed anchor
is at most `x*sqrt(1-x^2)<=1/2`, its negative quadratic term is
nonpositive, and its angle term is at most `t*pi`. Thus its score
is at most `1/2+t*pi<.67<C-D0`. This supplies a direct explanation
for any other zero-area closure without dividing by a vanishing
side or sine. Nondegenerate high-score faces are then subject to all
the positive denominators justified by localization.

I also independently read the author's subsequently separate
`DEGENERATE_CLOSURE_ADDENDUM.md`, SHA-256
`06ac07cdf3c1e55b32b8e06e31214390c0269ccfeee11831a4bae334bee53146`.
Its all-zero-area-sequence argument is correct: pass to convergent
apex angles, use the universal signed-anchor cap, drop the nonpositive
quadratic term, and obtain
`limsup score<=1/2+t*pi<46790713/70000000<7/10`.
The accepted `C>147/100` gives `C-1/40>289/200`, excluding any
such individual face from a high-score triple. Nonzero-area limiting
faces keep well-defined parameters and obey the original theorem.
This addendum is pinned separately; the original five author sources
and their manifest were not changed by either this review or the addendum.

## 7. Checker and global scope

Static inspection finds no optimizer, candidate search, old verifier
call, geometry enumeration or unproved numerical maximization in the
new checker. Its fixed scalar guards correspond to the analytic
whole-domain arguments above. Imported primitive bytes and accepted
scalar dependencies are checked before and after its verification.
This review did not execute its `verify()` function.

The local theorem covers arbitrary common covered unit anchors of
either signed side; no reflection, center norm, nearest-owner or
assigned-arc condition is used. It supports a direct replacement of
these same marked-area and actual-anchor terms in a valid signed
global area identity. Global incidence, disjoint charged face/edge
sets, matching coefficients, other profile rows and the once-only
surrounding budget remain separate obligations. The proof correctly
claims neither a complete global radius exclusion nor an exact optimum.

I found no mathematical blocker in the frozen local proof or checker.
