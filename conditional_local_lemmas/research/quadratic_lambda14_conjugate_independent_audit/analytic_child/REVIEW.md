# Independent static review of the corrected lambda-14 conjugate

**Verdict: PASS for the corrected fixed-parameter full-apex single-face
support theorem and its scalar enclosure logic.** This is an analytic
and source review. Child verifier calls, batches, optimizers, parameter
searches, and source mutations are all **zero**. Actual execution and
certificate matching are owned separately by the parent.

## Reviewed source identities and repair

| Source | SHA-256 |
| --- | --- |
| `research/quadratic_lambda14_conjugate/PROOF.md` | `aaaad3b1629ba98fd44907735b75ae111fecaa45c8d77f078f4e28acf4713ea3` |
| `research/quadratic_lambda14_conjugate/verify_conjugate.py` | `e95db0a6ead5ace03884ab54743b6238f8fe4cd242b1563617b789fce9cb57f4` |
| `research/coupled_edges/verify_support.py` | `90fd6e3426d8c113bacedd55f9a13f5e47981e9a48d21bfef48183ec4602efec` |

The initial proof pin was
`266dc36cf80d739dd2f99bcf1703fb77589a88644228f3adde64128b780913e4`
and initial checker pin was
`a5573868e6af596d48eefdd3fc0a15fc9f53285e4a30de42b9614622f8cae960`.
That proof was not accepted: differentiating `sin(a)+sin(2a)` gives the
second derivative `-sin(a)-4sin(2a)`, rather than the displayed
`-sin(a)-8sin(2a)`.

The owner repaired the coefficient, the consequent curvature bound, and
the checker's scalar guard. Static diffs against
`history_factor2_correction/` confirm precisely these three replacements:

```
f'' = -sin(a)-4sin(2a)-8lambda cos(2a),
f'' <= -2sqrt(3)+4lambda on [pi/4,pi/3],
guard: -2*(17/10)+4lambda < 0.
```

The score, first derivatives, root isolation, and output expressions
were unchanged. The initial snapshot is explicitly labeled
`SUPERSEDED_ANALYTIC_GUARD_FACTOR2`. This verdict concerns the corrected
source identities above, not the initial erroneous formula.

## 1. Every actual triangle is included

For a nondegenerate triangle, let its actual circumradius be rho<=R,
the chosen opposite side be `L=2Rx`, and the opposite apex angle be gamma.
The sine rule gives

```
L=2rho sin(gamma),  sin(gamma)>=x,
arcsin(x)<=gamma<=pi-arcsin(x).
```

In particular `0<x<=1`. If the other angles are beta and eta, their sum
is `pi-gamma`, and

```
area = (L^2/2) sin(beta)sin(eta)/sin(gamma),
sin(beta)sin(eta) <= cos(gamma/2)^2.
```

Therefore `area/R^2 <= x^2 cot(gamma/2)`. Equivalently, the perpendicular
height above the base is at most `Rx cot(gamma/2)`. Equality holds when
beta=eta. This is valid for arbitrary actual triangle skew and for
acute or obtuse gamma; it does not substitute a radius-R triangle for
an actual triangle without proving a majorant.

Since `J_R(L)/R^2=x sqrt(1-x^2)`, the normalized area-plus-J score with
quadratic penalty and angle tilt is bounded by exactly the displayed
function

```
x^2 cot(gamma/2)+x sqrt(1-x^2)-4lambda x^2+t gamma.
```

At fixed positive x, its gamma-second derivative is
`(x^2/2)csc(gamma/2)^2 cot(gamma/2)>0`. The feasible gamma interval is
closed, so its maximum is at one of its two endpoints, for every real
lambda and t. At either endpoint the isosceles equality triangle has
rho=R, and thus is feasible. This proves exact support, not merely an
upper bound obtained from incompatible shapes.

Write `a=arcsin(x)`. At gamma=a the area-plus-J term is
`sin(a)(1+2cos(a))`; at gamma=pi-a it is `sin(a)`. These give precisely
f and h in the source. At x=1 the endpoints coincide at gamma=pi/2 and
the two branch values agree. At x=0 the closure gives branch values
0 and `t*pi`, covering the limiting linear gamma support
`max(0,t*pi)`. The normalized area and J vanish in this limit. These
degenerate closure cases do not omit any positive-area triangle.

## 2. The fixed positive-tilt tails are both covered

Here `lambda=7/50` and `T=1071883/20000000`, with
`3/125<T<1/2`. For `a in [pi/3,pi/2]`, one has
`cos(a)<=1/2`, `cos(2a)<=-1/2`, and `sin(2a)>=0`. Thus

```
f'(a)=cos(a)+2cos(2a)-4lambda sin(2a)+T <= T-1/2<0.
```

For `h'(a)=cos(a)(1-8lambda sin(a))-T`, a nonpositive bracket gives
strict negativity immediately. A positive bracket may be bounded using
`cos(a)<=1/2` and `sin(a)>=sqrt(3)/2`, giving

```
h'(a)<=1/2-2lambda sqrt(3)-T < 3/125-T<0.
```

The last comparison follows from `sqrt(3)>17/10`. This sign split is
valid over the entire tail and does not assume global concavity of h.
Both positive-tilt maxima therefore lie in `[0,pi/3]`.

## 3. A new strict-concavity proof at lambda=7/50

The corrected second derivatives are

```
f''=-sin(a)-4sin(2a)-8lambda cos(2a),
h''=-sin(a)-8lambda cos(2a).
```

On `[0,pi/4]` both are strictly negative, including a=0, where their
value is `-8lambda`. On `[pi/4,pi/3]`,
`sin(a)>=sqrt(2)/2`, `sin(2a)>=sqrt(3)/2`, and `cos(2a)>=-1/2`. Hence

```
h'' <= -sqrt(2)/2+4lambda <0,
f'' <= -2sqrt(3)+4lambda <0.
```

The first strict inequality is certified by `(14/25)^2<1/2`; the
second follows already from `sqrt(3)>17/10`. These are valid new bounds
for lambda=.14, outside the old module's declared lambda<=1/8 range.
The tilt is linear and changes neither second derivative.

The left derivatives are `3+T` and `1-T`, both positive. The tail
bounds give negative derivatives at pi/3. Strict concavity therefore
provides exactly one interior critical point for each positive-tilt
branch. The same acute argument applies at zero tilt.

## 4. The zero-tilt obtuse branch is global

The positive-tilt obtuse tail proof is not reused at t=0. Instead,
`s=sin(a)` ranges over the entire `[0,1]` as a ranges over `[0,pi/2]`,
and

```
h_0=s-4lambda s^2.
```

Its vertex is `s=1/(8lambda)=25/28`, strictly inside `[0,1]`, and its
global value is `1/(16lambda)=25/56`. This includes any zero-tilt obtuse
maximum outside the retained positive-tilt interval.

Branch dominance also has a simple analytic check at these fixed
parameters. Since T>0 and sqrt(3)>17/10,

```
f_T(pi/6) >= f_0(pi/6) = 9/25+sqrt(3)/2 >121/100.
h_T(a) <=25/56+T*pi <25/56+33/175 <121/100,
```

using `T<3/50` and `pi<22/7`. Thus the acute branch dominates the whole
obtuse branch at both the prescribed tilt and zero tilt. This reinforces
the source's stricter interval comparisons without a numerical run.

## 5. Static enclosure logic

The checker uses the reviewed, source-pinned rational interval class,
Machin pi enclosure, and sine/cosine Taylor enclosures. The helper
`trig` encloses an interval argument through a midpoint enclosure plus
the interval radius, using the global 1-Lipschitz property of sin and
cos. Natural interval operations then enclose the branch scores.

The derivative formulas match the analytic first derivatives. Bisection
starts at 0 and `PI.a/3`, with strict outward derivative signs checked.
Because `PI.a/3 < pi/3`, this upper endpoint lies inside the proven
concavity interval. Its negative derivative and concavity cover the
small remaining segment to pi/3, while the tail proof covers the rest
to pi/2. There is no discarded gap in the angle domain.

Strict concavity makes every resolved midpoint sign select the side
containing the unique root. An unresolved sign causes failure, not an
assumed decision. The midpoint score is a valid lower bound on the
maximum, and the score interval over the final root bracket is a valid
upper bound. The width and acute/obtuse dominance comparisons are
explicit checks. The exact zero-tilt obtuse value replaces any attempted
zero-tilt obtuse root search.

Importing the primitive constructs its arithmetic constants and pi
enclosure; its old `verify()` is protected by the main guard and is not
called. The new source uses no old interval catalogue or optimizer. This
child only read these definitions and did not import or execute them.

## Scope and remaining interfaces

No analytic blocker remains for the corrected fixed-lambda single-face
support and its zero-tilt reference. The theorem does not prove a joint
fan deficit, outer-Q condition, assigned-arc bootstrap, eligible q3
incidence, same-score epsilon, or global radius exclusion. Parent runtime
verification and certificate matching remain separate evidence; no such
execution is attributed to this child.
