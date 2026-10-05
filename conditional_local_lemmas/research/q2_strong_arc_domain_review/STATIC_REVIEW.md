# Independent static review of the stronger two-face arc-domain analysis

Verdict: **PASS_RESTRICTED_MIRROR_FAMILY_AND_EXACT_GENERAL_INTERFACE**.

At the exact `R=5661/50000`, `lambda=7/50`,
`t=1071883/20000000`, and the stated exact target

```
delta=178523206195158271774904671901
     /8011730250000000000000000000000,
```

the entire exact-peak, mirrored-skew, circumradius-R, actual-outward-tip
family is incompatible with `D_joint<=delta`, the reflected norm cap,
and an outer actual assigned full arc greater than `147/1000`.
The general-fan identities (7)--(13), including directed correction
(8a) and its half-endpoint error weights, are also valid on their
explicitly stated geometric and actual-endpoint domain.

This verdict **does not accept a fee for every general two-face fan**.
The constrained general average-norm bound needed by the proposed
closure method remains unproved. No mirror extremality, Q bound,
global assignment/incidence, boundary-support allocation or whole-disk
exclusion is established. The accepted scalar peak and saved
dimensionless localization guards remain explicit premises; this
review does not regenerate their numerical enclosures.

## Sources and recorded correction

I read the complete final `../q2_strong_arc_domain/ANALYSIS.md` and
`EXECUTION_SCOPE.md`, the added directed-anchor paragraphs, source pins,
revision record and related-progress note. I also read the accepted
scalar proof/certificate and actual-quadratic support proof/certificate
identified in the source-pin map, the previous q2 symbolic analysis and
its independent review, and the necessary norm-cut note and review.
Only the scalar/one-face normalized facts are used as arithmetic
premises. The fixed `1/12` witness is not extrapolated by numerical
continuity or by transfer of its fixed-point guards.

The initial b>=0 display had `H<Qs`. I identified that as false at
zero skew, where `b=0`, `a=1`, and `H=Qs`. The author preserved the
prior source and revised that sign to `H<=Qs`. The final strict
rational bound `Qs<.948662<.949` remains valid, so no norm conclusion
or budget changed. The revision and prior SHA are retained in
`REVISION_RECORD.json`; this reviewer changed no author file.

The final author bytes reviewed here have these hashes, with paths
relative to `../q2_strong_arc_domain/`:

| File | SHA-256 |
| --- | --- |
| `ANALYSIS.md` | `9274e72efe3bee6c0aa4957e9e1bcc924bb217dda8558624a6a949f36846cac4` |
| `EXECUTION_SCOPE.md` | `5750407131f7821acb69b284131680dd2c15312560d9a8e6b83fa5ddc48b7eae` |
| `SOURCE_PINS.json` | `b8d082f0c9ea09bdf50c06af96754522ecba3b488c38f3786659715ab84d2bb4` |
| `REVISION_RECORD.json` | `c558acc376a4290d783afd8fb7952dc7d99a3ba9d6015682f8506890218e6be1` |
| `RELATED_PROGRESS.md` | `cfcecc10ed0391b4c97c2ef1db8da128221c33121b6e4b1996e51c5972b89a41` |
| `history/ANALYSIS_before_zero_skew_correction.md` | `26c456d69574a69a7a3c3633204148865fdee384c47c22a7b684abe087be5b03` |

The exact premise sources and the independent norm-cut review are
separately bound in this review directory's `SOURCE_PINS.json`.

## Actual covered arc gives the necessary norm cut

A radius-R ball with R<1 cannot cover an antipodal pair, so any
connected covered arc has full angle below pi. Endpoint separation
then gives `0<beta<=asin(R)<pi/2`. Adding the two endpoint disk
inequalities yields

```
rho^2+1-2 cos(beta)(c dot e)<=R^2.
```

Using `c dot e<=rho` produces a quadratic inequality in rho and its
lower-root bound `rho>=f(beta)`. The stated derivative is positive
on the feasible interior because `cos(beta)^2` exceeds the radicand
by `1-R^2`; continuity preserves strict order up to the feasible
endpoint. The comparison point z=.0735 is feasible since `sin z<z<R`.

The sine Taylor and cosine inequalities used are valid at this positive
z. Their exact finite-rational comparisons are correct:

```
z-z^3/6=.0734338224375>.0734,
cos z>=.997298875,
R^2=.0128187684,
R^2-(.0734)^2=.0074312084<.0074321641=(.08621)^2.
```

Thus `rho>f(z)>.911088875` for an actual full arc greater than .147.
This derivation requires neither radial assignment nor a complete
footprint. It uses an actual covered arc, never the triangle half-base
as an arc proxy. It establishes a necessary condition only, and does
not establish the parent full-arc hypothesis itself.

## Entire exact-peak mirrored family

For every admissible skew, `a^2+b^2=1`, `x=2sc`, and the defined
Q,V satisfy `V^2 c^2+Q^2 s^2=1`. These identities prove both displayed
distance-R anchor equations algebraically. Expanding also gives
`T+M=2k` and `T-k=2y(ca+sb)`, so the two actual ray pairs are the
stated mirrored triangles with a common middle side. Their bases have
length `2Rx`, their apex angles are theta*, and their circumradii
are R. The outward-tip formula, or its signed determinant, selects
the actual outer tips and zero signed-anchor slack.

The area identity is `area/R^2=x(1+y)-2x b^2`. Hence each score is
exactly `C-2x b^2`, and `D_actual=4x b^2` for either skew sign and at
zero skew. This is a symbolic identity, not a numerical property of
the former fixed `1/12` point. The nonnegative joint variance implies
`D_actual<=D_joint`. The integer comparison `1000A<23B` for the
stated exact target therefore gives `|b|<.09` and `a>.995` without
paying the variance twice.

Squaring the accepted cosine bracket and applying half-angle identities
gives all displayed rational x,c,s,Q,V bounds. For b>=0, the final
corrected bounds prove `k>.874155`, `l<.486`, `0<H<=Qs<.949` and
`u>.994`. In particular the displayed positive H lower estimate is
`.90691224>.9`. For b<=0 they prove `k>.910425`, `0<l<.403`,
`0<H<.978` and `u>.993`; the displayed lower l estimate `.31556`
is positive. Both sign cases include b=0 with the corrected weak
intermediate sign.

The exact norm identity (4) follows by expanding
`c0=q0-R(k,l)` and using `|q0|=1`, `k^2+l^2=1`.
The finite-rational products in (5) and (6) are correctly displayed:

```
1+R^2-2R(.994)(.874155)+2R^2(.949)(.486)
 =.8278871630468752,
1+R^2-2R(.993)(.910425)+2R^2(.978)(.403)
 =.8182098527881312.
```

Their gaps below `.910^2=.8281` are respectively
`.0002128369531248` and `.0098901472118688`. All quantities replaced
inside the negative and positive products have the required signs,
so the strict norm upper bounds follow for both epsilon signs,
including zero skew and budget equality. Reflection gives the same
norm for c2. Either outer strong arc therefore contradicts the
necessary norm cut on the positive translation branch.

On the negative translation branch the expansion changes only the
sign of the `2Ruk` term, giving a squared norm above one because
u,k,H,l are all positive in the target budget. This contradicts the
parent reflected cap `sqrt(1-R^2)<1`. The two unit anchors have fixed
vertical half-chord and nonzero H, so their midpoint must be horizontal
and have coordinate +u or -u. Up to a global rotation, these exhaust
the translations of this family. No branch is silently discarded.

## General-fan interface, including actual-anchor direction

Equation (7) follows from the sine-rule ray lengths in the accepted
skew convention. Both faces use the same positive middle ray length;
the equality is a real shared-side constraint. The lens-tip expressions
(8) use bisectors at `-gamma0/2` and `+gamma1/2` and the accepted
outward normal signs. They are valid without symmetry.

In normalized base coordinates, the exact lens inequality is
`v_tan^2+h^2+2x|v_tan|<=cos(theta)^2`. Substitution
`h=cos(theta)-Delta`, `Delta=tau/x` gives exactly (8a), including
the inward sign `eta=v_tan*T-Delta*N`. The tangent T=iN is the
directed base tangent for the stated positive fan convention.
The inequality automatically restricts `0<=Delta<=2cos(theta)`.
It covers negative-height anchors; eta is not a freely directed error
disk. The half-endpoint bounds for both mean and half-chord follow
from the triangle inequality and `x>.64`, with coefficient
`(1/2)(1/.64^2)=625/512`. E changes by minus the mean error and F
does not change. These bounds do not create a second available loss.

The unit equations for `a+Rb` and `a-Rb` give (9). For nonzero b,
the midpoint is `a=chi*u*n` with the stated u and a unit normal n.
Substituting `c0=a+R(E-F)`, `c2=a+R(E+F)` and averaging cancels the
F cross terms and proves (10) exactly. Both strong outer arcs imply
the necessary strict average lower bound (11). The displayed ideal
values of m,b,E,F are correct at zero skew, and (10) then reduces
to the known symmetric norm formula.

Equation (12) requires genuine distinct successive endpoints of a
connected short middle assigned arc. Its half-angle lies below pi/2,
so the principal arcsine follows from chord length `2R|b|`.
The document explicitly prohibits applying it to an arbitrary unit
anchor pair. Under that assignment, the actual variance is a function
of this same chord, rather than an independent free angle.

The one-face dimensionless localization and paid-budget argument applies
to each face because its individual nonnegative deficit is at most
`D_actual<=delta<.023<1/40`. This imports the scalar/one-face facts,
not the old fixed-radius three-face exclusion. The same argument
excludes zero-area or obtuse high-score faces before division; the
universal zero-area bound `1/2+t*pi<.7<C-delta` also excludes collapse
of any side. Summing the two paid face inequalities and adding the
one actual middle variance yields (13), with the exact current Knew.
It assumes the parent's once-only variance allocation. Actual anchor
tau, BC lens loss and variance remain distinct terms; none can be
spent again from an independent free budget.

## Execution boundary and unresolved claims

My mathematical execution count is **zero**. I performed no mathematical
program execution, Fraction calculation, project-source import,
checker/guard call, trigonometric evaluation, replay, optimizer or scan.
I used file reads, SHA-256 inspection, handwritten algebra and rational
comparisons, and writes confined to this separate review directory.
The author also records zero mathematical executions for this new
analysis. Previous witness or radius arithmetic executions remain
evidence belonging to those separate frozen artifacts.

`RELATED_PROGRESS.md` references a separate zero-tau norm proof awaiting
its own review. I read that reference's scope only; I did not audit or
accept the separate zero-tau theorem as part of this verdict. Its claimed
progress is not a premise needed by the accepted restricted-family
argument or the general identities.

After the recorded zero-skew display correction, no mathematical blocker
remains in this restricted-family result or exact interface. The positive-
tau constrained general norm estimate, general fee, global Q and full
cover reduction remain unproved outside this accepted scope.
