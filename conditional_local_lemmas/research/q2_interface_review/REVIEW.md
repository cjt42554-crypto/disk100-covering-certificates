# Independent static review of the q2 assigned-arc interface

**PASS for the necessary norm lemma and the stated scope distinctions.**
The interface does not establish a two-face fee or a whole-cover bound.
The reported mirror-skew candidate remains pending: its arc algebra,
prospective deficit and admissibility were not accepted or audited here.

I read `Q2_ARC_SCOPE_INTERFACE_20261005.md` completely. This review uses
analytic geometry only. I executed no checker, Fraction calculation,
mathematical module, optimizer or geometry search, and changed no root
interface or previously frozen source. Source hashing and these separate
report writes are the only programmatic operations.

## Endpoint coverage gives the claimed necessary lower root

Take u as the unit arc-bisector direction and v as its perpendicular
unit direction, as implied by the unit-circle endpoint parametrization.
For the stated small-arc domain `a=2beta`, `0<=beta<pi/2`, the two
unit endpoints are

```
w_plus = cos(beta)*u+sin(beta)*v,
w_minus= cos(beta)*u-sin(beta)*v.
```

Both squared distances being at most `R^2` implies their average is
at most `R^2`. Expanding that average cancels the opposite v terms
and gives exactly

```
1+rho^2-2*cos(beta)*(c dot u) <= R^2,
rho=|c|.
```

Because u is a unit vector, `c dot u<=rho`. Since `cos(beta)>0`,
replacing `c dot u` by rho makes the left side smaller, so the valid
consequence is

```
1+rho^2-2*cos(beta)*rho <= R^2,
(rho-cos(beta))^2 <= R^2-sin(beta)^2.
```

The squared left side is nonnegative; hence the right side is
nonnegative. For positive R this proves
`0<=xi=sin(beta)/R<=1`. Taking the lower root gives

```
rho >= cos(beta)-sqrt(R^2-sin(beta)^2)
     = sqrt(1-R^2*xi^2)-R*sqrt(1-xi^2).
```

The square-root choice for cosine uses the small-arc domain. No
radial placement of c, nearest ownership or equal arc assumption
was used. For `0<R<1` the lower root is positive: the difference
between the squares of its two nonnegative terms is `1-R^2>0`.

The edge cases are also included. At `xi=0` the expression is
`1-R`; at `xi=1` it is `sqrt(1-R^2)`. The original squared inequality
and its lower root remain valid when the discriminant is zero.
The bound is necessary and nonstrict; it does not assign a score
payment or by itself certify cyclic placement of several assigned arcs.

## Strict monotonicity and the `.6` threshold

On `0<xi<1`, both derivative denominators are positive, and direct
differentiation yields

```
rho_min'(xi)=R*xi*[1/sqrt(1-xi^2)
                   -R/sqrt(1-R^2*xi^2)].
```

The square comparison in the interface is exact:

```
(1-R^2*xi^2)-R^2*(1-xi^2)=1-R^2>0.
```

Therefore `sqrt(1-R^2*xi^2)>R*sqrt(1-xi^2)` and the bracket is
strictly positive. This proves the stated strictly increasing
function on the open interval. Its continuous endpoint values also
extend strict order to any two distinct points of `[0,1]`; no
derivative at the singular endpoint xi=1 is needed.

Thus `xi>=3/5` gives the correct necessary bound

```
rho >= rho_min(3/5)
     = sqrt(1-(9/25)*R^2)-(4/5)*R.
```

The interface correctly keeps this inequality nonstrict. A strict
xi threshold would yield a strict norm comparison with this value,
but xi equal to `3/5` is not excluded by the lemma.

## Base x, assigned-arc xi and conditional q2 scope

The two normalized quantities have different geometric inputs:

```
x_i = |c_(i+1)-c_i|/(2R),
xi_j= sin(a_j/2)/R.
```

The first comes from a triangle base joining two centers; the second
comes from two unit-circle endpoints assigned to a cell. The same
normalizing radius does not make these lengths equal or imply an
ordering. In particular the three-face certificate's x condition
cannot be silently replaced by an assigned-arc xi condition, or
conversely. The interface correctly states this distinction.

For the small-arc range the exact extra domain condition `xi>3/5`
is equivalent to

```
a > 2*asin(3R/5).
```

The reported `a>.0597` lower bound at `R=5661/50000` does not imply
that stronger condition. As a comparison of the scalar assumptions,
`a=.06` is above `.0597` but has
`sin(a/2)<a/2=.03<3R/5`, hence xi below `3/5`. This comparison
does not purport to supply an entire-cover or q2 counterexample.

The exact rational extra-fee requirement is correctly presented as
a conditional integration target supplied by the parent, not as an
accepted conclusion of this norm proof. The root interface explicitly
retains the unproved radial `Q<=.35754` premise and does not claim
that `xi>.6` has been established. This audit does not independently
verify the numerical q2,d7 profile row or its required fee.

The mirror-skew candidate with `sin(epsilon)=1/12` and its stated
prospective deficit remains pending. No conclusion about its actual
arc lengths, admissibility, local fee obstruction, four-cell window,
complete degree-seven fan or extension to a 100-circle cover is made
by this audit. Unknown guard results have not been treated as evidence.

The norm lemma consumes no assigned-cell lens, triangle score,
terminal BI, variance or other credit. Any later local obstruction
or two-face support must use its exact domain and actual assigned
endpoints. Whole-cover incidence, radial control, all other supports,
consistent cyclic placement and disjoint charges remain separate
obligations. No mathematical blocker was found in the stated
necessary norm lemma.
