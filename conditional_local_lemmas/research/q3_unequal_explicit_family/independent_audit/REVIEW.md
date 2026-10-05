# Independent audit of the explicit unequal q3,d7 local obstruction

PASS for the single exact algebraic diagram and the stated local obstruction
at R=11321/100000, lambda=7/50, t=1071883/20000000. It satisfies the
reflected norm caps, actual canonical nearest-owner switches, connected
covered assigned arcs longer than .092 for all 33 boundary cells, and a
complete interior q3,d7 seven-face star with strict empty circumcircles of
radius R. Its entire unit circle is covered. Its retained whole three-face
score has loss strictly below the proposed .019931551 payment.

This does not assert a cover of the whole closed unit disk, a 100-site
cover, the required global profile, all other face circumradius bounds,
outer-Q or a counterexample to a global radius bound. Those additional
conditions could still exclude this diagram from the global problem.

## Source identities and exact construction

The reviewed frozen identities match the requested final source:

| File | SHA256 |
|---|---|
| PROOF.md | 6d168f09df53368e088080b51fc59849c704291a24919eae6bc09134a6890525 |
| verify_explicit.py | c64f86cbc2fe2657ee00df18534d2901188b30356ea7268a5f9f7f7ca9688efb |
| certificate.json | e5db4d60907d510bcc606edc89243b579036dc8d63c2004fceb3703580c1db75 |

A=359/500 and B=1691/2500 are rational cosine values. All positive
square-root branches and the cubic root are specified algebraically.
The ring coordinates arise from positive half-angle square roots for
pi/8 and pi/16 followed by exact rotation addition, so they are algebraic
sites rather than floating trigonometric samples. Interval expressions in
the checker enclose these exact real sites; correlation loss only enlarges
their enclosures.

With a=acos(A), b=acos(B), the four p-relative fan vectors factor as

```
c0-p = 2R*cos(a-b/2)*exp(-i*(a+b/2)),
c1-p = 2R*cos(b/2)*exp(-i*b/2),
c2-p = 2R*cos(b/2)*exp(+i*b/2),
c3-p = 2R*cos(a-b/2)*exp(+i*(a+b/2)).
```

Their actual apex angles are a,b,a, their normalized bases are
sin(a),sin(b),sin(a), and they are nondegenerate. The circumcenter
formulas O0,O1,O2 give exact incident radius R. Direct determinant
calculation independently gives outer normalized area
sin(a)[cos(a)+cos(a-b)] and middle area sin(b)(1+cos(b)).
Their whole retained objective is consequently exactly

```
S = 2G(a)+G(b)-2*sin(a)*(1-cos(a-b)),
G(z)=sin(z)*(1+2cos(z))-(14/25)*sin(z)^2+t*z.
```

This is the requested area+J-quadratic-base+angle objective, rather than
an actual signed-anchor objective or an equal-isosceles approximation.

## Cubic completion and the complete actual degree-seven star

The cubic f(k)=4k^3-3k+(2A^2-1) has opposite exact rational signs at
.85 and .87, and f'>0 throughout that interval. Its bisection therefore
isolates a unique actual root. For u=(pi-2a)/3, cos(u) lies in this
interval and satisfies the cubic, identifying k=cos(u).

Write s=sqrt(1-k^2), x2=2A^2-1 and y2=2A*sqrt(1-A^2). The cubic gives
x2=3k-4k^3. The positive quantities y2 and s(4k^2-1) have the same
squared value, so y2=s(4k^2-1). Thus d0-O3=R(-k,s) and
d0-O4=R(x2,y2); d1-O4=R(-k,-s), with reflected identities for d2.
Together with the displayed O0 through O6, these check every incident
distance exactly R. The checker does not numerically test a false
interval equality: these are analytic exact identities for the algebraic
sites enclosed by its intervals.

The seven consecutive ray gaps around p are

```
a, b, a, a+u-b/2, 2u, 2u, a+u-b/2.
```

All are positive and below pi, and their sum is 4a+6u=2pi. This supplies
the full once-around argument behind the proof's compressed half-angle
sentence; positive triangle orientations alone would not suffice.

For each of the seven circles, the checker tests all 34 nonincident
sites, including every added ring site. All 238 squared-distance
exclusions are strict. The minimum certified squared exclusion margin
is approximately .022766513499028283. Therefore these seven faces are
strict Delaunay faces, and their once-around union exhausts p's star,
giving full degree seven rather than only seven selected triangles.

The four interior sites p,d0,d1,d2 have norms strictly below 1-R.
Their distances to every unit-circle point exceed R. Once the circle
is covered, none can own a circle point. In particular p's convex
Voronoi cell cannot extend outside the unit disk without meeting its
circle, so its complete star is retained in the clipped diagram.
Its four other neighbors are actual boundary owners. The three
successive switched edges c0c1,c1c2,c2c3 mark precisely the first three
faces. The other four have at least one additional interior vertex.
This is actual q=3,d=7 incidence in this diagram.

## Whole unit-circle coverage and all 33 boundary sites

Every ring center has exact norm n=sqrt(1-R^2). Its radius-R circle
footprint has half-angle beta=asin(R). The exact guard
sin(pi/32)<R implies pi/16<2beta, so all consecutive ring footprints
j=2,...,30 overlap and cover their entire long circle chain.

The point w listed in the proof is the right footprint endpoint of the
ring center at -pi/8: it is exactly the rotation of that direction by
beta, because cos(beta)=n and sin(beta)=R. c0 strictly covers w;
c3 covers its reflected partner. Their squared overlap slack exceeds
.0000172837. The three canonical fan transitions are each covered by
both adjacent local disks.

For each local endpoint-to-endpoint joining arc, the checker gives
positive site dot unit-point throughout (-.3,.3). Its second angular
derivative is the negative of that dot product. Hence the dot product
is concave there, with minimum at an endpoint, so endpoint coverage
covers the whole intervening short arc. These four covered arcs bridge
the exposed ring endpoints, proving coverage of the entire unit circle,
including all joins. This argument does not claim coverage of its disk.

The footprint-containment guards put all four local footprints strictly
inside (-.3,.3). A negative local site's polar angle is above -.3 and
below zero; exclusion at both +/- .3 therefore encloses its connected
footprint, and reflection supplies the positive sites. This containment
is material for the ring cells' minimum assigned arcs.

For ring indices j=2,...,29, use the outward owned interval
[j*pi/16,j*pi/16+pi/32]; for j=30 use
[30*pi/16-pi/32,30*pi/16], with angles unwrapped as needed. Each lies
outside (-.3,.3), is contained in its ring disk because pi/32<beta,
and is nearest among the equal-norm ring sites. The missing ring sites
at indices 31,0,1 do not compete. Local disks lie outside their own
footprints and interior sites are farther than R, so none can beat
that ring owner. Every one of the 29 ring cells thus contains an owned
interval of length pi/32>.092.

## Genuine switches, terminal ownership and connected full arcs

The first bisector equation is v dot q=delta. Its displayed algebraic
q0 has unit norm and v dot Jq0=sqrt(ell2-delta^2)>0, selecting the
ordinary CCW transition. The other root lies in the distant upper-left
quadrant, outside all four local footprints. q1=(1,0) and q2 are the
reflection-related actual transitions.

On the negative protected interval, both coordinates of c1-c0 and
c2-c0 are positive. Their comparison functions increase with angle.
The c0/c1 comparison changes sign exactly at q0, while nearest_guard02
remains favorable to c0 for the whole incoming segment. The equal-norm
reflection comparisons have the sign of the point's y coordinate.
On the negative middle segment, c1 beats c3 because its horizontal
difference is negative, its vertical difference positive, and the
outward nearest_guard13 is negative for every angle in [-.2,0].
Reflection proves the positive segments.

The rational half-tangent rotation defines a unit rotation with
.092<alpha<.1. q0's angle has magnitude between .092 and .2.
The terminal owned segment from Qleft to q0 therefore has length
alpha>.092 and stays inside (-.3,0). Its coverage slack exceeds
.0000339319. Over the entire protected interval, the nearest ring
directions are +/-pi/8 and their distances are minimized at the
corresponding endpoint. The certified squared excess over R^2 is
approximately .0000349204 at each endpoint, so no ring site can
compete anywhere in that interval. Interior sites cannot compete either.
Thus the terminal and middle segments are actual owned segments,
rather than arbitrarily chosen common-covered points.

For completeness, the connected-arc lemma is correct for all 37 sites.
For distinct sites x,y with norms at most n, their bisector's distance
from the origin is

```
abs(delta)/|y-x| <= |(x+y)/2| < n.
```

Strictness follows from strict convexity of the norm ball and the sites'
distinctness. Its two circle intersections have shorter separation
2acos(distance)>2asin(R). Any one radius-R disk's footprint has length
at most 2asin(R): minimizing its cosine threshold
(1+rho^2-R^2)/(2rho) gives sqrt(1-R^2). A pair comparison therefore
has at most one crossing in a site's connected footprint, and its
favorable portion is an interval, prefix, suffix, all or empty.

Full circle coverage confines every nearest-owned arc to its owner's
footprint. Intersecting these favorable intervals proves the complete
assigned arc is connected. Nearest ownership and coverage then ensure
the entire assigned arc is covered by its owner. The four local cells
contain the owned segments above, and all 29 ring cells contain the
half-cell intervals. Thus all 33 distinct boundary cells have connected,
covered full assigned arcs of length strictly greater than .092.

## Whole-score obstruction and existing payments

The checker reads only the accepted frozen complete conjugate at the
exact lambda and positive tilt; it does not rerun the conjugate checker.
Direct interval evaluation of the actual three areas and their base
terms gives

```
3C-S = .018317365588123673... < .019 < .019931551,
target-(3C-S) = .0016141854118763264... > 0.
```

Thus the listed local conditions, even with the actual full q3,d7 star
and 33 covered boundary arcs, do not imply the proposed standalone
whole-S loss. This conclusion is consistent with the earlier tiny-fee
theorem; that theorem is not re-executed or contradicted.

The shared-spoke skew equality is also exact. The normalized lengths
are (V,U,U,V), with signed half-differences (d,0,-d) and adjacent mean
differences (d,-d). Equal outer weights make the completed common square
vanish. Consequently Q_mean equals the already included skew sum
2*d^2*cot(a/2)=2*sin(a)*(1-cos(a-b)). It cannot be charged again.

For the two middle owned arcs, the midpoint and half-angle formulas
give exactly P=sin(delta)*(hmax-h), the existing old H-area slack.
The normalized sum is approximately .9024050258434848 and positive.
The checker records it without subtracting it from S. A future joint
target using this slack must reserve its existing boundary payment once
and supply compatible global incidence; it is not an additional free
local fee after that loss has already been consumed.

## Execution history and final scope

The parent made four earlier attempts on this one algebraic candidate:
two retained plumbing failures, one successful corrected run, and the
final successful run adding footprint guards. The preserved sources,
logs, status records and first successful certificate remain unchanged.

This independent audit made exactly one new verify_explicit.verify()
call in one fresh isolated Python process, with assertions enabled,
bytecode writing disabled and BLAS/parallel-library threads set to one.
The entire returned JSON matched the saved certificate. The fresh
certificate has the identical SHA256
e5db4d60907d510bcc606edc89243b579036dc8d63c2004fceb3703580c1db75.
All 29 current, history and dependency files matched their before/after
hashes. The source's __main__ was not executed and no frozen file was
modified. The resulting total is five checker attempts: two plumbing
failures and three successes, including this audit.

The two analytic subreviews were read-only and made zero verifier calls.
Conjugate, primitive and latest tiny-fee verify calls: zero. Optimizer,
parameter campaign, graph catalogue, old batch and additional candidate
calls: zero. The earlier one-point floating diagnostic is retained
history and was not repeated by this audit.

No analytic source error or execution blocker was found. Acceptance is
one exact local obstruction within the enumerated conditions. Whole-disk
coverage, a 100-site global diagram, the remaining face circumradius
conditions, the global profile and outer-Q remain uncertified. This is
not a global counterexample or an exact optimum result.
