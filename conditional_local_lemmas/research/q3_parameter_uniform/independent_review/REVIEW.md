# Independent static review of the parameter rectangle theorem

**PASS for the stated local theorem.** For
`27/200 <= lambda <= 29/200`,
`1071883/20000000-1/200 <= t <= 1071883/20000000+1/200`, and
`0 < R <= 231/2000`, the proof establishes

```
sum_i [area_i/R^2 + signedA_i/R^2 -4lambda*x_i^2+t*gamma_i]
    < 3C(lambda,t)-16027/800000,
```

on its actual common-apex, consecutive, consistently oriented,
three-positive-face domain with covered unit anchors of either signed
side. The fee exceeds `1/50` by `27/800000`. This is an independent
static mathematical verdict. It supplies neither a complete covering
radius exclusion nor the exact optimum.

I read the complete new proof; the frozen actual quadratic proof,
checker, certificate and all-zero-area closure addendum; the full-apex
conjugate proof, checker and certificate; and the uniform-radius
addendum and independent static review. I subsequently read the complete
new Fraction-only guard, its dependency manifest, certificate, author
log and author receipt. No verifier, imported mathematical module,
optimizer, search, or numerical scalar guard program was executed by
this reviewer. `SOURCE_PINS.json` binds the reviewed versions.

## Full single-face support, including signs and all apex angles

For every admitted positive triangle the sine rule supplies
`theta=asin(x)` and `theta<=gamma<=pi-theta`. The exact area identity
has nonnegative skew loss `sigma`. The two endpoint disk inequalities
give `tau=x*cos(theta)-signedA/R^2>=0`, including negative anchors.
Consequently the actual score is exactly `F-sigma-tau` with the stated
`F`. The quadratic and tilt changes introduce no geometric dependence
into `sigma,tau`.

Convexity of `F` in the complete actual apex interval puts its maximum
at `gamma=theta` or `gamma=pi-theta`; the endpoint values are precisely
the displayed `G` and `Bo`. The compact endpoint support therefore
gives every actual score at most `C(lambda,t)`. The support is exact as
an independent single-face support: the isosceles area and outer lens
tip are feasible, and the whole single-face geometry may be translated
so that the tip is a unit point. No constraint on the actual shared
apex is introduced by this independent maximization. The three-face
compatibility is treated separately by the accepted actual-anchor fee.

The positive `theta=pi/2` endpoint has a singleton feasible gamma
interval and is included. The endpoint functions are continuous at
`theta=0`. No obtuse or skew branch is discarded by an initial acute
assumption, and the signed anchor is never replaced by its absolute
area or by a standalone independent lens allowance.

## Uniform new support location and stronger-fee localization

The old peak is strictly inside `(.8,.85)`, hence inside `(L,H)` with
`L=.7,H=.96`. The old curvature bound `-G0''>4` on this box and the
new derivative perturbation bound

```
|-4(lambda-lambda0)*sin(2theta)+(t-t0)| <= 1/40
```

give the stated strict new derivative signs at `L,H`. In particular,
the right endpoint magnitude from the coarse old peak bracket is
`4*(H-.85)-.025=.415`; the proof's strict sign statement is valid.

For every lambda in the rectangle, `-G''` is positive on `[0,pi/4]`.
On `[pi/4,pi/3]` its explicit lower estimate is
`7/10+4*(4/5)-4*(29/200)=83/25>0`. On `[pi/3,pi/2]`, the positive
quadratic coefficient and `sin(2theta)>=0` permit
`G'<=1/2-1+t<0`. Thus the entire acute branch has one maximum inside
`(L,H)` and is monotone on either side. Full-interval concavity, which
would be false at these coefficients, is not used.

The old peak is a feasible comparison point for the new support.
The displayed reference lower bound is greater than `7/5`, whereas
the full obtuse branch is bounded by

```
Bo <= 1+tmax*pi < 1+tmax*(22/7)
   = 82890713/70000000.
```

The entire obtuse gap is therefore greater than
`15109287/70000000`, which exceeds `16027/800000` by a wide margin.
This proves acute dominance without having to locate the new obtuse
maximum.

The perturbation of an endpoint difference is at most
`(1/40)*(H-L)=13/2000`. The copied old endpoint margins consequently
give the explicit lower gaps

```
C-G(L) > 1/40 + 137/10000 - 13/2000 = 161/5000,
C-G(H) > 1/40 + 14/1000  - 13/2000 = 13/400.
```

Both are strictly larger than `16027/800000`. Thus the stronger-fee
contradiction, not only the initial `1/50` contradiction, forces
every `theta_i` into `(L,H)`.

The endpoint ratio guard has no interior-minimum gap: differentiation
of its derivative numerator gives the stated strictly negative
quantity. Both endpoint ratios exceed `81/50`. Up to the cut
`gamma=theta+1/10`, the denominator sine is increasing, so the new
loss slope is strictly greater than

```
(1/2)*(81/50)^2-tmax = 25072117/20000000 > 1.
```

The cut is feasible since `2H+1/10<3<pi`. Its loss exceeds `1/10`,
and convexity on the remainder of the full apex interval bounds that
remainder by its cut and obtuse endpoints. Both endpoint gaps exceed
the stronger fee. It follows that all actual `z_i=gamma_i-theta_i`
lie in `[0,1/10)`.

## Old budget after the new localization

The proof uses only the old curvature and old gamma-loss slope on
the newly established actual box. Integrating the old curvature from
the old peak and the old slope from each acute endpoint supplies

```
2 sum e_i^2 + sum z_i + sum sigma_i + sum tau_i <= D_old,
e_i=theta_i-theta0,
D_old=3C0-S(lambda0,t0).
```

This implication does not assume `D_old<=1/40`. The normalized
old budget remains valid on this box even when the old total deficit
is larger. In particular, the new proof does not attempt to reuse the
old fit bounds, such as `|V|<.55`, outside their old small-deficit
domain. It needs only the budget inequality for sensitivity.

The separate accepted radius extension supplies `D_old>1/40` for
the whole admitted radius interval. Its exact upper endpoint exceeds
`231/2000`. The same actual triple and the same positive radius are
used at both parameter pairs; no relabeling or invariance of a fixed
configuration under radius changes is assumed.

## Centered sensitivity and the strict final fee

The new score identity is exact. Evaluating the new independent
support at the old acute peak gives the lower support comparison
used in the proof. Subtracting score from three support values has
the correct signs:

```
D_new >= D_old
         +4*dl*(sum x_i^2-3*x0^2)
         -dt*(sum gamma_i-3*theta0).
```

The function `sin(theta)^2` is 1-Lipschitz, and Cauchy with the
budget gives `sum |e_i|<=sqrt(3D_old/2)`. Since
`sum gamma_i-3theta0=sum e_i+sum z_i` and `sum z_i<=D_old`, the
two centered sensitivity estimates are valid for arbitrary
`D_old>1/40`. Taking the two independent parameter deviations at
their rectangle extremes yields

```
D_new >= f(D_old),
f(D) = (199/200)*D-(1/40)*sqrt(3D/2).
```

No differentiability of the new maximum or equality in the old-peak
comparison is needed. The support comparison is one-sided and
suffices even if a maximizing branch changed.

For `D2>D1>=1/40`, rationalizing the square-root difference shows
that the secant slope is strictly greater than

```
199/200 - (3/80)/(2*(19/100)) = 1703/1900 > 0,
```

because `(19/100)^2<3/80`. Thus `f` is strictly increasing on the
entire required ray. The upper square-root guard is exact:

```
3873^2-15000000=129>0,
sqrt(3/80)<3873/20000.
```

Therefore

```
f(1/40) > 199/8000-3873/800000
        = 16027/800000.
```

The accepted strict `D_old>1/40` makes the final bound strict. Under
the stronger hypothetical threshold all preceding localization
arguments still apply, so the contradiction proves the displayed
strict stronger fee, including both parameter endpoints.

## Degeneracies, radius and allocation

At each fixed positive radius a zero-area limiting face has normalized
area tending to zero. Its anchor score is at most `1/2`, its quadratic
term is nonpositive, and its tilt is at most `tmax*pi`. Hence its
limsup score is at most

```
1/2+tmax*pi < 47890713/70000000 < 1.
```

Since `C>7/5`, this is well below `C-16027/800000`. The other two
scores are at most `C`, so no high-score triple can contain such a
face. This treats also collapse of a side different from the base.
Positive-area limits have well-defined parameters and obey the same
proof. `R=0` is excluded, and the fixed-positive-radius closure must
not be applied merely because physical area tends to zero while a
varying radius tends to zero.

The radius addendum's only physical radius dependence is the
three-unit-point condition. Its strict fit implies the contradiction
for all `0<R<=949600/4987489`, including the closed upper endpoint.
The parameter theorem uses only the subset `0<R<=231/2000`.

All sensitivity losses are estimates from the same old deficit, not
additional physical charges. A whole-cover use still requires the
same lambda and t in all supports, a correct signed area identity,
eligible fan incidence and disjoint charged ordinary face/anchor
sets. Assigned-cell lenses, variance, ear allowances and reserved
BI credits each retain their separate single allocations. Nothing
here converts this theorem to an independent-J fee.

## Static guard consistency and remaining scope

The new guard imports only `Fraction`, paths, hashing and JSON. It
reads pinned accepted certificates and does not import or execute
any frozen mathematical source. Its peak brackets, endpoint margins,
old curvature and ratio enclosures, new concavity and tail bounds,
acute and obtuse comparisons, cut slope, degeneracy bound,
square-root comparisons and sensitivity fee all match the analytic
steps above. The complete emitted certificate is consistent with
these formulas. The author preserved the first successful source,
log, receipt and certificate under `attempts/author_001`. After this
static review suggested two rerun-hygiene improvements, the final
guard was changed to reject disabled assertions explicitly and to
check the radius status, coefficients, fee, exact endpoints and
open/closed flags. I read the full original guard and the complete
final diff. Those changes are correct and change no theorem or
numeric guard value.

The second author receipt reports one necessary successful new
Fraction call on the revised final source and two total successful
new author calls. Its log/certificate is unchanged from the first
successful result. The exact computational accounting for this
parameter theorem is therefore two new author Fraction guard calls,
zero reviewer guard calls, zero failed calls reported, zero frozen
checker calls or frozen mathematical imports, and zero optimizer or
geometry-search calls. This review supplies no additional runtime
receipt or numerical replay. Source hashing and static byte comparisons
are the only programmatic operations performed by this reviewer.

No mathematical blocker was found. Global profile rows, matching
surrounding supports and full-cover optimality remain outside this
review's scope.
