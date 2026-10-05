# Independent static review of the stronger assigned-arc norm cut

**PASS for the necessary norm restriction.** At
`R=5661/50000=.11322`, a center covering an actual connected
unit-circle arc of full length `a>.147` must have

```
|c|>.911088875>.911.
```

This is a local necessary geometric condition. It is not a two-face
score fee, a proof that any covering has `a>.147`, a verification of
the reported gamma=1.85 support, or a global radius exclusion.

I read `Q2_STRONG_ARC_NORM_CUT_20261005.md` completely and checked its
argument against the independently reviewed arc-endpoint norm lemma.
No checker, Fraction calculation, trigonometric evaluator, mathematical
module, optimizer, search or historical replay was executed by this
reviewer. Only file reads, source hashing and this new review directory
were used; no root source or previously frozen artifact was changed.

## Arc domain and monotonic comparison

Since R is less than one, a covered connected unit-circle arc cannot
contain a pair of antipodal points: their distance is two and exceeds
the diameter `2R` of the covering ball. Hence its full length is less
than pi. With `beta=a/2`, this gives `0<beta<pi/2`.
Endpoint coverage and the triangle inequality give
`2sin(beta)<=2R`, so `beta<=asin(R)`.

The accepted necessary norm lemma therefore applies in exactly its
small-arc domain:

```
rho>=f(beta),
f(beta)=cos(beta)-sqrt(R^2-sin(beta)^2).
```

Its strict order in xi is equivalent to strict order in beta on this
feasible interval, because `xi=sin(beta)/R` is strictly increasing
there. The finite continuous endpoint at `beta=asin(R)` is also
covered; no derivative at a zero square-root denominator is needed.
Thus `beta>z=147/2000` gives `rho>=f(beta)>f(z)`.

The fixed comparison point is feasible: `sin(z)<z=.0735<R`, so its
square-root argument is strictly positive. The proof does not assume
the center is radial, inside the unit disk, or the nearest owner of
the whole arc; actual endpoint coverage suffices.

## Taylor bounds and exact rational comparisons

For the positive small z, alternating sine Taylor terms give the strict
bound `sin(z)>z-z^3/6`. The omitted term signs and decreasing magnitudes
are valid at this `z<1`. Direct finite-rational simplification gives

```
z-z^3/6=.0734338224375>.0734.
```

The elementary global inequality `cos(z)>=1-z^2/2` gives the stated
exact lower bound `.997298875`. All decimal constants here are exact
finite rational numbers rather than floating-point enclosures.

The four square comparisons in the source are correct:

```
R^2=.0128187684,
(.0734)^2=.00538756,
R^2-(.0734)^2=.0074312084,
(.08621)^2=.0074321641.
```

The last square exceeds the preceding radicand bound by
`.0000009557>0`. Since `sin(z)>.0734>0`, the actual radicand is
strictly smaller still. Positive square-root monotonicity therefore
gives

```
sqrt(R^2-sin(z)^2)<.08621.
```

Combining the nonstrict cosine bound with this strict square-root
bound yields

```
f(z)>.997298875-.08621=.911088875.
```

The margin over `.911` is exactly `.000088875>0`. Hence the claimed
strict norm cut follows, including cases where the actual covered
arc reaches its feasible maximum length.

## Scope and dependence

The parent supplies `a>.147` as a new conditional-domain hypothesis.
This review does not establish that hypothesis, the gamma=1.85
boundary support, or the reported still-unproved radial Q premise.
No numerical scalar check of those external statements can be inferred
from this norm derivation.

Any previously studied configuration with an outer norm below this
cut cannot belong to the new covered-arc domain. That implication
follows from the necessary norm lemma. It does not independently
audit the old configurations or extend their local obstruction to
the new domain. Their acceptance records, actual arc conditions and
score deficits remain separate evidence.

The distinction between triangle half-base x and assigned-arc
half-chord xi remains in force. The norm cut charges no triangle
score, actual anchor, assigned-cell lens, BI or variance credit.
A future joint proof must verify its actual arcs and fan incidence,
match its dual coefficients, and retain the parent's single
allocation of middle-half-arc variance and all other credits.

No mathematical blocker was found in this necessary norm cut. This
review belongs to a new directory and is not part of a previously
frozen q2 or four-tree ZIP.
