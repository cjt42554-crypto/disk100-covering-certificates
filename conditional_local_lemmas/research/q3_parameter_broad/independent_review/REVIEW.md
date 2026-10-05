# Independent static review of the broad actual-A rectangle

**PASS for the stated local theorem and restricted-curvature conjugate
cap.** On the actual positive-area, consistently oriented consecutive
common-apex three-face domain, with covered unit anchors of either
signed side and circumradii at most the same positive R, the proof gives

```
63/500 <= lambda <= 77/500,
791883/20000000 <= t <= 1351883/20000000,
0 < R <= 794400/5378209,
S(lambda,t) < 3C(lambda,t)-41/2000.
```

Both coefficient endpoints and the upper radius endpoint are included.
The normalized local fee is `.0205`. The independent support C has the
complete acute/obtuse endpoint definition from the referenced proof.
This review supplies neither an independent-J fee nor a whole-cover
radius exclusion, profile verification or exact optimum.

I read the complete final `PROOF.md` and `verify_fraction_guards.py`,
the final certificate/log/receipt, and the preserved failed receipt and
log. The complete final-to-failed proof and source diffs show only the
peak-displacement rational correction. I previously read the full
reference parameter proof, original and stronger actual-A proofs and
certificates, full-apex scalar proof/checker/certificate, accepted
stronger analytic/fresh audit, and the final `.038` radius addendum and
static audit during the immediately preceding independent tasks. For
this follow-up I checked their retained acceptance records and source
hashes; their frozen bytes and the earlier reviews remain unchanged.

No mathematical checker, Fraction calculation, mathematical module,
optimizer or geometry search was executed by this reviewer. File reads,
hashing and report writes are the only programmatic operations.
`SOURCE_PINS.json` binds the reviewed versions and acceptance records.

## Complete independent support and new peak location

The reference sine-rule identities hold for every admitted positive
triangle and every positive R. They give `theta=asin(x)`, the full
interval `theta<=gamma<=pi-theta`, and the exact actual score
`F(lambda,t;theta,gamma)-sigma-tau` with `sigma,tau>=0`. The signed
anchor identity and tau remain valid for negative anchor height.
Convexity in gamma is independent of lambda and t, so the full-apex
maximum is one of G and Bo before any acute or skew localization.
Consequently every individual actual score is at most C and all its
new individual deficits are nonnegative.

The broad lambda interval is positive. On `[0,pi/4]`, each term in
`-G_new''` is nonnegative and their sum is strictly positive. On
`[pi/4,pi/3]`, the displayed bound is
`7/10+4*(4/5)-4*(77/500)=821/250>0`. On the entire tail
`[pi/3,pi/2]`, the quadratic derivative term is nonpositive and
`G_new'<=-1/2+tmax<0`. This proves complete acute-branch unimodality
without incorrectly asserting concavity on `[0,pi/2]`.

The old peak lies strictly in `(.8,.85)`. Integrating the old curvature
bound on `[L,H]` gives `G0'(L)>.4` and `G0'(H)<-.44`. The new derivative
perturbation has magnitude at most `4|dl|+|dt|<=7/100`, leaving strict
endpoint signs `.33` and `-.37`. Thus the unique new global acute peak
lies strictly inside `(L,H)` for every parameter pair in the closed
rectangle.

Evaluation at the feasible old peak gives

```
C_new > 147/100-4*(7/500)-(7/500)*(24/25)
      = 17507/12500 > 7/5.
```

Here the old cosine/sine/peak facts are used only at the reference
coefficients. No old fixed-coefficient scalar C value is substituted
for the new C. The entire obtuse branch is at most `1+tmax*pi`, strictly
less than its rational `22/7` bound. Its gap below `7/5` exceeds the
claimed fee, so it cannot compete with the acute maximum or the later
high-score threshold. The proof includes the full initial apex interval,
including obtuse faces.

## The tight endpoint comparison is still strict

The perturbation function `-4dl*sin(theta)^2+dt*theta` is
`7/100`-Lipschitz, because `sin(theta)^2` is 1-Lipschitz. Comparing
the new endpoint value to the feasible old peak therefore loses at
most `(7/100)*(H-L)=.0182`. The saved reference endpoint gaps are
strictly larger than `.0387` and `.039`, respectively. Hence the
smaller new endpoint gap is strictly larger than

```
387/10000-(7/100)*(13/50)=41/2000.
```

The rational coarse-bound equality in this line is intentional and
valid. Its strictness comes from the already accepted strict old gap,
not from asserting that the coarse arithmetic difference itself is
larger than the fee. Thus an individual score at exactly `C-fee` is
also excluded outside the open theta window. The saved actual lower
enclosure leaves further positive room beyond that coarse strict bound.

If `D_new<=41/2000`, each face score is at least `C-41/2000`. The
global acute monotonicity, both endpoint gaps and the full obtuse gap
force every theta into `(L,H)`. The ratio guard has no hidden interior
minimum: its derivative numerator has strictly negative derivative,
so any interior extremum is a maximum. Both accepted endpoint ratios
exceed `81/50`. The new gamma-loss slope stays greater than
`(1/2)*(81/50)^2-tmax>1` up to the `.1` cut. That cut is feasible
because `2H+.1<pi`.

Its endpoint score loses at least `.1`. Full gamma convexity on the
remaining interval bounds it by the cut value and the obtuse endpoint,
both below `C-fee`. Thus every actual slack satisfies `0<=z_i<.1`.
No actual apex restriction was assumed before this argument.

## One old budget and centered sensitivity

On this newly established actual box, the reference curvature and
gamma-loss slope integrate to

```
2 sum(theta_i-theta0)^2 + sum z_i + sum sigma_i + sum tau_i
    <= D_old.
```

This does not assume that D_old is small. The broad argument uses only
this scalar budget to estimate sensitivity; it does not reuse old
small-deficit V or fit bounds outside their domains.

The separately accepted actual-A theorem and radius extension give
`D_old>19/500` for the same admitted triple at every radius in the
declared interval. The new geometry is in the expanded domain of that
accepted theorem; no center norm, assigned arc, reflection, ownership
or apex-location hypothesis is silently added.

The exact centered support comparison has the same signs as before:

```
D_new >= D_old
         +4dl*(sum x_i^2-3x0^2)
         -dt*(sum gamma_i-3theta0).
```

Cauchy gives `sum |theta_i-theta0|<=sqrt(3D_old/2)`, and
`sum gamma_i-3theta0=sum(theta_i-theta0)+sum z_i`. Applying the two
parameter bounds `|dl|,|dt|<=7/500` and the one budget gives exactly

```
D_new >= f(D_old),
f(D)=(493/500)D-(7/100)sqrt(3D/2).
```

For `D2>D1>=19/500`, rationalizing the square-root difference shows
that the secant slope is strictly greater than

```
493/500-(21/200)/(46/100)=4357/5750>0,
```

using `(23/100)^2<57/1000`. Thus f is strictly increasing on the entire
required ray. The accepted strict old deficit and the exact square
comparison `(191/800)^2-57/1000=1/640000>0` give

```
f(D_old)>f(19/500)
 > (493/500)(19/500)-(7/100)(191/800)
 = 41511/2000000
 > 41000/2000000=41/2000.
```

This contradicts the closed high-new-deficit threshold. The sensitivity
margin is `511/2000000`; the theta-localization comparison is the tighter
coarse guard. No new physical fee is charged for any budget estimate.

## Restricted-curvature cap and corrected displacement

Both peaks have already been placed in `(L,H)`. Changing lambda
perturbs the second derivative by at most `8*(7/500)`, while changing
t contributes no second derivative. Therefore

```
-G_new''>4-8*(7/500)=486/125=kappa
```

on the full segment between the peaks. At theta0 the new derivative is
`v=dt-4dl*sin(2theta0)`. The integral Taylor bound and completing the
square give the exact coefficient `1/(2*kappa)=125/972` in the displayed
upper cap. The full support is the acute peak already proved above,
so this bounds C itself, not only a restricted branch approximation.

Integrating the derivative gives displacement strictly less than
`|v|/kappa` when v is nonzero; for v zero, uniqueness makes the peaks
coincide. Since `|v|<=7/100`, the nonstrict uniform bound is

```
(7/100)/(486/125)=35/1944.
```

The failed author's assertion of equality to `7/388` was incorrect.
That fraction is slightly larger than the exact quotient, so the old
displayed upper bound was loose rather than an omitted geometric case;
the precise rational equality in the guard nevertheless failed and
was correctly repaired. The current proof and guard use the correct
exact constant. The cap can be outward evaluated from the retained
reference peak/cosine/C0 enclosures without isolating a new peak.

## Accepted dependencies, closure and computation record

The `.038` independent analytic/fresh audit genuinely accepts the
fixed-radius actual signed-A theorem. Its one fresh new guard call
matches the author JSON and serialized certificate, with historical
checks unchanged. The separate final uniform-radius proof and
independent static acceptance justify the same normalized premises
on `0<R<=794400/5378209`; this endpoint is closed by strict actual-fit
inequalities. Fixed-radius certificate metadata has not been relabeled
as a range certificate. The zero-skew product qualification is explicit
in the accepted addendum/clarification and is compatible with strict
final fit errors. Previously accepted scalar root/value enclosures
remain an explicit arithmetic premise, not a new computation here.

At each fixed positive R, a zero-area limiting face has score at most
`1/2+tmax*pi<1`, while `C_new>7/5` and the fee is less than `.4`.
Such a face cannot belong to a high-score triple, because the other
two face scores are at most C. This covers all collapsed-side limits
before localized divisions. Positive-area limits retain their defined
parameters and satisfy the proof. R is used consistently in all
coverage, circumradius, area and base normalizations; `R=0` is excluded.
Physical area tending to zero with a vanishing radius is not by itself
the fixed-positive-radius zero-area closure.

The final guard imports only Fraction and JSON and rejects disabled
assertions. Its arithmetic comparisons agree with every new constant
above. It explicitly declares that geometry and accepted dependencies
are external; its success is not treated as a substitute for their
proofs and reviews. Its complete final certificate agrees with the
source's returned fields.

The preserved attempt-001 log/receipt record one failed new author
Fraction call at the incorrect exact-displacement assertion. The final
receipt records two total new author calls, one failure and one success.
The source/proof diff is limited to replacing `7/388` by `35/1944`.
The receipt also separately preserves one patch-format failure before
any write or mathematical execution; it is not an extra Fraction call.
There were zero old checker calls/imports and zero optimizer or geometry
search calls reported. This reviewer adds zero mathematical executions.

No mathematical blocker was found in the final local theorem. A
whole-cover use still requires the same dual coefficients, actual
eligible ordinary fan incidence, disjoint charged face/anchor sets
and every remaining support/profile/reduction. Assigned-cell lenses,
terminal BI credits and variance are not charged by this broad result.
The previously frozen four trees and their ZIP were not modified;
this new directory is a separate later result.
