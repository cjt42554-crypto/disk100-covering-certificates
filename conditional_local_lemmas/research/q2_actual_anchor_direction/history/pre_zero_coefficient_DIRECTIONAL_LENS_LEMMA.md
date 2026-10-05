# Exact directional support for actual signed-anchor slack

This lemma applies to one actual base and its actual anchor in the
intersection of the two endpoint radius-R disks. It is independent
of a unit-anchor equation, another face, assigned-arc ownership,
apex/circumradius, or any near-peak localization. It does not prove
a two-face fee or a complete-cover exclusion.

Let R>0, 0<x<1, y=sqrt(1-x^2)>0. Choose an oriented orthonormal basis
(e,n) with det(e,n)=1 and actual base endpoints

    c0=m-R*x*e, c1=m+R*x*e.

Write the actual anchor q and the outward lens tip q* as

    q=m+R*t*e+R*(-y+v)*n,
    q*=m-R*y*n.

The signed external area is

    A=-det(c1-c0,q-c0)/2=R^2*x*(y-v).

Thus the exact nonnegative signed-anchor slack tau=x*y-A/R^2 equals
x*v. In particular v=tau/x. The endpoint disk conditions are exactly

    (t+x)^2+(-y+v)^2<=1,
    (t-x)^2+(-y+v)^2<=1.

Combining them gives the equivalent lens constraint

    t^2+2*x*|t|+v^2<=2*y*v.                     (1)

It follows that 0<=v<=2*y and

    |t|<=T_x(v)=sqrt(x^2+2*y*v-v^2)-x.          (2)

Conversely, each v in [0,2*y] and each t with |t|<=T_x(v) obey both
endpoint disk constraints. Hence (2) is an exact fixed-slack bound,
not a relaxation that assumes the anchor lies on an outer boundary.
At v=0 or v=2*y, T_x(v)=0; signed A can be negative when v>y.

For any real coefficients a,b, the exact largest directional anchor
error at that fixed slack is

    max [a*e+b*n] dot (q-q*)/R
      =|a|*T_x(v)+b*v.                         (3)

The equality is attained with t having the sign of a and magnitude
T_x(v). Formula (3) retains the sign of the normal coefficient b;
a separate norm bound need not retain it.

By rationalizing (2), for v>0,

    T_x(v)=v*(2*y-v)/(sqrt(x^2+2*y*v-v^2)+x).

The same identity holds at v=0 by continuity. Consequently

    [a*e+b*n] dot (q-q*)/R
    <= [b+|a|*(2*y-v)/(sqrt(x^2+2*y*v-v^2)+x)]*tau/x
    <= (b+|a|*y/x)*tau/x.                      (4)

The last inequality uses 0<=v<=2*y, so its nonnegative tangential
ratio is at most y/x. No absolute value of b was introduced. When
b<0, actual motion toward the apex can decrease the chosen functional.
This potentially improves a later two-unit-anchor norm estimate, but
such improvement must be established for that estimate's actual
coefficient direction rather than inferred from (4) alone.

As a corollary,

    |t|<=y*v/x,
    |q-q*|/R<=v/x=tau/x^2,                     (5)

recovering the accepted nondirectional anchor-error bound. Another
exact implication of (1) is

    |q-q*|^2/R^2+2*x*|t|<=2*y*tau/x.           (6)

Equations (3)--(6) spend the same tau. A later argument must not also
count tau as an assigned-cell lens loss in the global Q ledger.
Those two slacks have different geometric definitions.

Degenerate endpoints are simple: at x=1, y=0 the intersection of
the two endpoint disks is the single midpoint, so q=q* and tau=0.
At x=0 the signed area vanishes and tau cannot control tangential or
normal anchor motion; formulas dividing by x are not asserted there.
Near-support localization in an application may bound x away from 0,
but that is a separate premise.

All equations above follow directly by algebra from the two actual
endpoint inequalities. No new mathematical program, guard, optimizer,
parameter scan, frozen source import, or historical replay was run to
write this lemma. File reads, writes and hashing are metadata work.
The unresolved two-face step is to apply these signed directional
errors to the exact two-unit-anchor translation and actual outer norms
with one deficit budget, including the real middle HALF-arc variance.
