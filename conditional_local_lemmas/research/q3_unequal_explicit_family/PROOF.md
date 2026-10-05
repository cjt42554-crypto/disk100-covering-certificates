# An algebraic obstruction to the proposed standalone three-face fee

At R=11321/100000, lambda=7/50 and t=1071883/20000000, the listed local
reflection, actual-switch, assigned-arc, empty-circumcircle and q3,d7
conditions do not imply a same-score loss >.019931551. This note gives
one explicit algebraic diagram with loss approximately .0183173655881.
It has 37 sites, 33 boundary cells whose entire assigned unit-circle
arcs are covered, and an interior degree-seven q3 apex. Its entire unit
circle is covered. It is NOT a 100-site cover of the whole closed unit
disk, a counterexample to a global bound, or a certificate of outer-Q.

The distinction is essential: a local support argument using only the
conditions verified here cannot reach the proposed fee. The global
profile, remaining faces, or a compatible allocation of existing lens
loss may still exclude the diagram's occurrence in a 100-site cover.

## 1. Exact algebraic centers and the seven-face star

Put A=359/500, B=1691/2500, sA=sqrt(1-A^2), sB=sqrt(1-B^2),
a=acos(A), b=acos(B), eta=1/10^7 and

    rho=sqrt(1-R^2(2-B^2))-R(1+B)-eta,
    p=(rho,0).

All square roots mean their positive real branches. The four consecutive
boundary sites are

    c0=(rho+R(2A^2-1+B), -R(2A*sA+sB)),
    c1=(rho+R(1+B), -R*sB),
    c2=reflection of c1 in the x axis,
    c3=reflection of c0 in the x axis.

These are a symmetric but unequal three-face fan. Its successive apex
angles are a,b,a. The central spokes have normalized length 2cos(b/2),
and the end spokes length 2cos(a-b/2). The three circumcenters are

    O0=p+R(B,-sB), O1=p+R(1,0), O2=p+R(B,sB).

Every incident site is exactly at distance R from its corresponding O.
The outer-face asymmetric payment is exactly

    sigma=sin(a)[1-cos(a-b)].

Thus, with G(z)=sin(z)(1+2cos(z))-.56sin(z)^2+t*z, the actual whole
three-face score is

    S=2G(a)+G(b)-2sin(a)[1-cos(a-b)].

This is not the restricted equal-isosceles case a=b. It attains equality
in the newly derived shared-spoke Qmean lower bound: its differences
are (d,0,-d) and outer weights are equal. Qmean cannot be spent a second
time to exclude it.

To complete the apex degree to seven, let k be the unique root in
(.85,.87) of

    4k^3-3k+(2A^2-1)=0,

and put s=sqrt(1-k^2). Equivalently k=cos(u), u=(pi-2a)/3. Define three
additional interior neighbors

    d0=p+R(2k(1-2k^2),4k^2*s),
    d1=p+R(-2k,0), d2=reflection of d0.

The cyclic star is c0,c1,c2,c3,d0,d1,d2. Its remaining circumcenters are

    O3=p+R(2A^2-1,2A*sA),
    O4=p+R(-k,s), O5=p+R(-k,-s),
    O6=reflection of O3.

Their incident distances are exactly R. To check the nontrivial identity,
the cubic gives 2A^2-1=3k-4k^3; positivity gives
2A*sA=s(4k^2-1), since both sides have the same squared value. Then
d0-O3=R(-k,s) and d0-O4=R(2A^2-1,2A*sA), and the other identities
follow by reflection. The seven positive apex angles sum to 2pi:
the known four neighbor half-angle parameters sum to 2a, and the
three new half-angle parameters sum to pi-2a. All faces are positively
oriented; their rays sweep once around p.

The checker isolates k using exact rational signs and monotonicity of
the cubic on [.85,.87]. It checks every nonincident site lies strictly
outside all seven exact circumcircles, after adding the boundary sites
below. The minimum squared exclusion margin is >.0227. Hence the seven
faces are strict Delaunay faces, and form the complete apex star.

## 2. A 33-cell covered boundary and genuine interior apex

Let n=sqrt(1-R^2). Add 29 sites

    n*(cos(j*pi/16),sin(j*pi/16)), j=2,3,...,30.

Their coordinates are algebraic by repeated half-angle construction;
their norm is exactly n. Combined with c0,c1,c2,c3 these are 33 distinct
boundary sites. Every c_i has strictly positive slack in both
(1-R)^2<|c_i|^2<1-R^2. The four interior sites p,d0,d1,d2 satisfy
|site|<1-R, so none can be closest on a circle covered by radius-R disks.

Each added ring site's circle footprint has half-angle beta=asin(R).
Consecutive ring spacings pi/16 are shorter than 2beta, since the exact
checker proves sin(pi/32)<R. The right footprint endpoint of the ring
site at angle -pi/8 is the algebraic unit point

    w=(n*cos(pi/8)+R*sin(pi/8), R*cos(pi/8)-n*sin(pi/8)).

It is strictly covered by c0; its reflected point is covered by c3.
Both transitions between c0,c1,c2,c3 are strictly covered too. Each
short arc between consecutive covered endpoints is contained in its
site's disk: throughout the relevant angular window (-.3,.3), site
dot unit-point is positive and concave, so its minimum is at an endpoint.
Together the 29 ring footprints and these four arcs cover the entire
unit circle, with strict overlap at the joins. The checker records the
join slack >.0000172 in squared distance.

All added ring sites are strictly outside all seven star circumcircles.
The covering of the unit circle also ensures p and d0,d1,d2 have cells
disjoint from it, since their minimal distance to the circle exceeds R.
Thus p is a genuine interior apex, and exactly four of its seven
neighbors are boundary sites. The three consecutive faces with two
boundary vertices are the marked faces; the other four are unmarked.

## 3. Actual nearest-owner switches and strictly long assigned arcs

For v=c1-c0, ell2=|v|^2 and delta=(|c1|^2-|c0|^2)/2, define

    q0=(delta*v-sqrt(ell2-delta^2)*Jv)/ell2,
    q1=(1,0), q2=reflection of q0.

Here J is counterclockwise rotation by pi/2. The branch has
v dot Jq0>0, selecting the ordinary CCW owner transition. It satisfies
|q0|=1 and v dot q0=delta identically. The other algebraic intersection
lies in the upper-left quadrant. The three q_i are covered by their
two incident boundary disks with strict distance slack. They are true
nearest-owner transitions in the complete 37-site diagram.

For a strictly longer terminal owned segment use the rational unit rotation

    h=461/10000, sin(alpha)=2h/(1+h^2),
    cos(alpha)=(1-h^2)/(1+h^2), Qleft=Rot(-alpha)q0.

Its alpha is strictly between .092 and .1. The checker proves Qleft is
covered by c0 with squared slack >.0000339, and the reflected point
by c3. All ring sites are at distance >R throughout the protected
interval from Qleft to its reflection: the nearest ring directions are
the two first sites at +/-pi/8, their distances are minimized at the
corresponding endpoint, and both endpoint squared excesses exceed
.0000349. Interior sites also have circle distance >R. Thus within this
interval the nearest owners are among c0,c1,c2,c3 alone.

For the negative outer arc, c0 versus c1 changes sign exactly at q0;
the difference derivative is positive throughout (-.3,0), since both
coordinates of c1-c0 are positive. The c0 versus c2 derivative is likewise
positive and its value is already strictly favorable at q0. The c0/c3
equal-norm comparison has the sign of the unit point's y coordinate.
Hence c0 owns the entire segment Qleft to q0. For the negative middle
arc q0 to q1, c1 beats c0 by the same sign change, beats c2 by the
equal-norm vertical comparison, and beats c3 by the outward interval
inequality in nearest_guard13. Reflection proves the positive arcs.
No arbitrary common-covered point is substituted for a true switch.

The middle arc lengths are strictly between .092 and .2. The two
terminal cells each contain an owned segment of length alpha>.092.
All assigned full arcs are covered because the entire circle is covered
and a nearest disk is no farther than any covering disk.

For completeness, the reflected norm cap guarantees connected assigned
arcs here. If distinct sites x,y have norm <=n, their bisector chord has
distance |delta|/|y-x| <= |(x+y)/2| < n from the origin. Its two unit
intersection points have shorter angular separation >2asin(R).
A single radius-R circle footprint of any site has length <=2asin(R).
Thus every pair comparison crosses at most once in a site's connected
footprint, and its favorable subset is an interval, a prefix or suffix.
Since circle coverage confines a site's nearest-cell arc to its footprint,
intersection over comparisons gives a single interval. The local four
footprints lie in (-.3,.3). Other ring cells have at least their outward
half-cell interval pi/32>.092, unaffected by those local footprints.
Consequently ALL 33 boundary assigned arcs have length >.092.

## 4. Strict whole-score obstruction and the separate existing lens loss

The accepted full-apex conjugate at these exact lambda,t is read from
its immutable certificate; its checker is not re-executed. Direct outward
calculation of the actual three triangle areas and their exact base
terms gives

    .01 < 3C_lambda(t)-S < .019 < .019931551.

More precisely the enclosed loss is about .018317365588123673, leaving
margin about .001614185411876326 below the requested target. The centers
and true switches are algebraic; all comparisons are outward rational
interval comparisons, with analytic exact equalities as specified above.
This is a counterexample to the standalone local fee in the listed
domain, even with q3,d7 and 33 covered boundary arcs. It is not a full
100-disk covering counterexample, and does not satisfy or refute an
undefined global outer-Q condition.

The checker also computes, without spending, the EXISTING two middle
assigned-cell lens losses. For one middle cell, let its half-angle be
delta, s=sin(delta), m=cos(delta), h=m-c1 dot midpoint_direction and
hmax=sqrt(R^2-s^2). Its old H-area slack is P=s(hmax-h). The two cells
have (P1+P2)/R^2 approximately .90240502584. Therefore a joint score
S-(P1+P2)/R^2 is a materially different and potentially useful target.
The P_i must be reserved from their existing boundary H-area budget
exactly once, with a global incidence allocation. They cannot be appended
as an independent new fee to a budget which already consumed them.

## Reproducibility and scope

Only one explicit algebraic candidate was built. A preceding one-point
floating diagnostic used a=.77,b=.828; it was not an optimizer or scan.
Root checker history contains four actual attempts: two plumbing errors
(Fraction wrapper and final str/Path identity guard), then a successful
exact run, then the final successful run with additional footprint guards.
Both failure sources/logs and the first successful certificate are retained.
No parameter campaign, graph catalogue, old checker or C verifier was run.
The final certificate, exact source identity and independent review record
are material evidence. This note cannot complete the original optimum.
