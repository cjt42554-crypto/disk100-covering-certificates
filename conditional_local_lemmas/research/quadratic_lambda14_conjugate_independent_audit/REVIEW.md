# Independent audit of the fixed lambda=7/50 full-apex conjugate

Analytic verdict: PASS for the corrected frozen source. The initially
supplied acute second-derivative formula had a factor-two error, found before
any independent replay. The owner preserved it as superseded history and
corrected the formula and rational guard. The accepting runtime outcome and
actual counts are recorded in the separate audit certificate and manifest.

Corrected proof SHA256:
`aaaad3b1629ba98fd44907735b75ae111fecaa45c8d77f078f4e28acf4713ea3`.
Corrected checker SHA256:
`e95db0a6ead5ace03884ab54743b6238f8fe4cd242b1563617b789fce9cb57f4`.
The saved numerical certificate remains unchanged at
`f325d42bfe651599881b7dc652567729a1586389ee385c38060840bc41696bc2`.

## 1. Actual circumradius and fixed-base endpoint reduction

For any nondegenerate triangle with base `L=2Rx` and actual apex angle
gamma, its sine rule gives `rho=L/(2sin gamma)`. Thus `rho<=R` forces
`sin gamma>=x`, or `arcsin x<=gamma<=pi-arcsin x`, even for an obtuse
apex or an actual circumradius strictly smaller than R. At that fixed
base and angle the other two angles have fixed sum. Their sine product
is maximized by equality, which gives normalized triangle area at most
`x^2*cot(gamma/2)`. This is the height/area cap for the actual triangle,
not an assumption that it is isosceles.

The exact one-face score in the source adds the cap anchor
`J_R(L)/R^2=x*sqrt(1-x^2)`, subtracts `4lambda*x^2`, and adds `t*gamma`.
The second derivative with respect to gamma is
`(x^2/2)*csc(gamma/2)^2*cot(gamma/2)>0` for positive x. Its maximum
therefore lies at either permitted endpoint, for every real lambda and
t. The endpoints themselves are attained by isosceles triangles with
circumradius R, so this computes the exact unrestricted single-face
cap score, rather than just a safe upper bound. The x=0 and x=1 cases
follow in the closure; x=0 includes the obtuse `t*pi` endpoint.

With `x=sin a`, `a in [0,pi/2]`, the endpoint branches are exactly

```
f(a)=sin a+sin(2a)-4lambda*sin(a)^2+t*a,
h(a)=sin a-4lambda*sin(a)^2+t*(pi-a).
```

This definition uses the full cap anchor J_R. It does not establish any
actual shared unit-circle anchor, three-face compatibility fee, annular
projection condition, or signed outer-area budget.

## 2. The whole discarded tail and the retained interval

For `lambda=7/50` and `t=1071883/20000000`, the entire acute tail
`a>=pi/3` has `f'(a)<=t-1/2<0`. For the obtuse derivative
`h'=cos a*(1-8lambda*sin a)-t`, a nonpositive bracket proves negativity
immediately. If its bracket is positive, the tail bounds `cos a<=1/2`
and `sin a>=sqrt(3)/2` give
`h'<=1/2-2lambda*sqrt(3)-t<3/125-t<0`. Hence the complete positive-tilt
domain can be restricted to `[0,pi/3]` without assuming global
concavity of h.

Direct differentiation gives

```
f''=-sin a-4sin(2a)-8lambda*cos(2a),
h''=-sin a-8lambda*cos(2a).
```

The initially supplied source had `-8sin(2a)` in f''; the child audit
identified this factor-two error before any fresh replay. The correct
coefficient is -4. On `[0,pi/4]` both second derivatives are strictly
negative. On `[pi/4,pi/3]`,

```
h''<=-sqrt(2)/2+4lambda<0,
f''<=-2sqrt(3)+4lambda<0.
```

These use the exact guards `(14/25)^2<1/2` and `sqrt(3)>17/10`.
Thus the intended strict-concavity conclusion remains true with the
corrected formula and the corrected guard `-2*(17/10)+4lambda<0`.
The strict endpoint derivative signs give exactly one retained
interior maximum per positive-tilt branch. No claim from the old
module's declared `lambda<=1/8` concavity domain is extended here.

## 3. Root and value certification over the full domain

The new checker bisects the rational derivative signs on
`[0,PI_lower/3]`; its strict negative right derivative puts the unique
root below that rational endpoint. Strict concavity on the complete
retained interval and the proved discarded-tail decrease make this
root the global branch maximum. Thus the slight difference between
the rational right endpoint and pi/3 leaves no unverified gap.

The pinned primitive's sinq/cosq are rational Taylor sums with absolute
Lagrange-remainder bounds. These remain valid for the doubled derivative
arguments up to 2pi/3; the checker does not use the primitive's narrower
monotone sin()/cos() wrappers. Its interval trig helper evaluates the
midpoint Taylor enclosure and enlarges by the argument radius, using
the global 1-Lipschitz property of sine and cosine. This encloses the
score of every angle in the final root bracket. A lower enclosure at
the bracket midpoint bounds the actual maximum below from the correct
side. The upper interval evaluation bounds that root value above.

Branch dominance is established by strict disjoint outward value
enclosures. The bisection tolerance alone is not an acceptance argument;
the sign bracket, concavity, tail proof and outward value enclosure are
all needed together. Any unresolved derivative sign before the fixed
tolerance raises an error.

## 4. Zero tilt and exact obtuse maximum

The acute tail still decreases at t=0, and the corrected acute
strict-concavity proof applies unchanged. The entire zero-tilt obtuse
branch is `s-(14/25)*s^2` for `s=sin a in [0,1]`. Its vertex is
`s=25/28`, which belongs to this interval, and its exact value is
`25/56`. This argument needs neither the positive-tilt obtuse tail
test nor any unproved global angle concavity. Strict outward comparison
with the acute maximum then gives the zero-tilt conjugate.

## 5. Evidence and scope

The final audit performs exactly one newly authorized
`verify_conjugate.verify()` call after the corrected source is frozen,
and compares its full return dictionary with the saved certificate.
All target files, their manifest, and the sole arithmetic primitive
are checked byte-for-byte before and after. Importing that primitive
does not call its historical verifier. The child audit is analytic only.

No old batch, optimizer, local incidence, covering construction, joint
epsilon, arc bootstrap, or global lower radius is established here.
The theorem supplies only the unrestricted full-apex cap-anchor
single-face conjugate at the fixed lambda and positive tilt, plus its
zero-tilt reference.
