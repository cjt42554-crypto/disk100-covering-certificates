# Replacing three marked allowances by an actual-anchor quadratic support

This is a symbolic substitution rule at the exact same coefficients as the
new local theorem. It is not a new unrestricted radius exclusion. The
ordinary reduction, the applicable type/angle window, the same-radius
boundary support, all remaining types, and all nongeneric limits remain
obligations of the whole-cover proof.

## 1. Start before the independent anchor majorant

For a retained ordinary boundary, write its actual counterclockwise centers
as c_i, its actual covered unit switches as q_i, and the assigned short arc
from q_(i-1) to q_i as 2a_i, with sum a_i=pi. Use the signed terms

    A_i = -det(c_(i+1)-c_i, q_i-c_i)/2,
    L_i = a_i + [det(q_i,c_i)+det(c_i,q_(i-1))]/2.

The signed center polygon area P satisfies the exact cancellation identity

    pi = P + sum A_i + sum L_i.                         (1)

If the center triangulation is retained and embedded, P is its actual sum
of triangle areas. Formula (1) follows by adding shoelace determinants;
it holds with either anchor sign. No independent lens or anchor maximum
has yet been substituted. On an ordinary boundary edge of length l_i=2Rx_i,
coverage of q_i gives

    A_i/R^2 <= J_i = x_i sqrt(1-x_i^2),
    tau_i = J_i-A_i/R^2 >= 0.                           (2)

The new theorem includes tau_i in its one paid deficit. It does not charge
the assigned-cell lens displacement in L_i.

## 2. The quadratic base and lens term stay together

Fix lambda and the actual marked-face angle coefficients t_i. Each selected
ordinary marked face has an interior apex, its own boundary base, and is
allocated once by the old proof. For that face set

    T_i = area_i/R^2 + A_i/R^2 - lambda*(l_i/R)^2 + t_i*gamma_i.

Peeling selected faces and their selected anchors from (1) is an exact
identity. Adding back their lambda*l_i^2 and subtracting their
R^2*t_i*gamma_i leaves precisely the old remaining area budget. In the
usual complete ordinary partition, the complete base/lens term is

    Q_actual = sum L_i + lambda*sum l_i^2.              (3)

The reflected radial pair support estimates (3), including the effects of
the actual radial/tangential displacements in L_i and in l_i. Its standard
form is

    Q_actual <= B*TARGET + 2*pi*MU
                - GAMMA*sum(a_i-pi/B)^2.              (4)

It uses the actual centers, assigned arcs, lens terms and base lengths.
The anchor majorant J_i and its loss tau_i are absent from this functional.
Using (4) once together with a theorem about T_i therefore does not spend
tau_i or a lens loss twice. Any lens or variance payment already used to
prove (4) remains inside (4); it cannot be added as another fee. If a
different whole-cover implementation incorporates tau_i into its boundary
support, this substitution needs a new audit of that implementation.

This defined Q_actual is not an identification of the still unspecified
"outer-Q <= .35748752" from the latest parent proposal. No numerical
boundary value is inferred here.

## 3. The exact local replacement

At R=11321/100000, lambda=7/50, and t=1071883/20000000 the accepted full-apex
conjugate C_lambda(t) bounds each old independent score

    area_i/R^2 + J_i - 4lambda*x_i^2 + t*gamma_i.

It also bounds T_i by (2). The new actual shared-anchor theorem states,
for three consecutive common-apex ordinary faces and their actual covered
unit anchors,

    T_0+T_1+T_2 < 3*C_lambda(t)-1/40.                  (5)

Both sides use the same lambda,t, triangle areas, base lengths and actual
angle sum. Thus a whole-cover proof that previously allocated 3C to these
three faces may replace that allowance by 3C-1/40, decreasing its area
upper bound by R^2/40. The original angle-star identities, all unmarked
corner terms, ear handling, q0 terminal BI credits and boundary support
stay in their existing locations and are each used once.

Only ordinary one-base faces are eligible for (5). When ears occur elsewhere,
their two anchors and base penalties remain in their existing separate ear
cap. An ear allowance is not an eligible three-face allowance. All unselected
marked or ear anchors can still use A_i<=R^2 J_i as in the old proof.

Several selected fans provide additive payments when their marked-face
and boundary-anchor sets are disjoint. A marked face has one interior
apex, so distinct eligible apexes give disjoint such sets. Sharing a
terminal boundary center causes no duplicate face or anchor charge. The
new theorem imposes no assigned-minimum-arc or reflected norm condition;
those may nevertheless be needed by (4) and the rest of the global proof.

The strict rational comparison

    1/40 > 19931551/1000000000 = .019931551

exceeds the parent's stated local target for the actual-anchor replacement.
This does not prove the proposed target for the independent J score.
The explicit 37-site diagram disproves that standalone target in the
listed local domain. Passing from that score to (5) keeps its existing
tau_i instead of replacing it by zero.

## 4. Exact remaining interfaces

Before using (5) in a radius exclusion, establish that each chosen fan is
actually consecutive, has circumradii at most the same R, has its common
covered unit anchors on the correct boundary edges, and has exactly the
same common angle coefficient t in the global dual. Recheck any one-face
or type support outside the old declared lambda<=1/8 interval at lambda
7/50. The frozen new conjugate handles these three marked faces only;
it does not automatically extend every other component of an old proof.

A same-radius boundary/domain bootstrap, the allocation of all eligible
fans, all other type/topology ranges and closure are still required.
The .1132 archive audit supplied with this interface concerns the structure
of (3)-(4), and is not a fresh replay of its full global theorem or a
transfer of its numerical constants to R=.11321. No optimal polynomial,
matching cover, or unrestricted global optimality follows from this note.
