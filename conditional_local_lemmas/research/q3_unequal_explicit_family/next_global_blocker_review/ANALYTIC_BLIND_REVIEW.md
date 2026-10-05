# Independent analytic blind review of the 37-site obstruction

Verdict: **PASS for the stated local geometric domain and the standalone
fee obstruction.** This is a source-specific analytic review. No candidate
verifier, conjugate verifier, optimizer, parameter scan, graph catalogue,
old bound batch, or cover checker was executed. All those call counts are
zero. Numeric replay and complete certificate matching belong to the
separately assigned portable reviewer, not to this review.

Reviewed frozen source identities:

| Source | SHA256 |
| --- | --- |
| `PROOF.md` | `6d168f09df53368e088080b51fc59849c704291a24919eae6bc09134a6890525` |
| `verify_explicit.py` | `c64f86cbc2fe2657ee00df18534d2901188b30356ea7268a5f9f7f7ca9688efb` |
| `certificate.json` | `e5db4d60907d510bcc606edc89243b579036dc8d63c2004fceb3703580c1db75` |
| `SHA256.json` | `c3a299cd60ca949372d8f708ae916cfa7047cce1da55f07af0c4c9f69511e00f` |

The saved certificate was read as a pinned owner execution artifact;
it was not regenerated or accepted by a fresh run in this review.

## 1. Algebraic definitions and the unequal score are consistent

The positive radical definitions give A=cos(a), B=cos(b) and the stated
four sites. In complex notation relative to p,

```
c0-p=R(exp(-2ia)+exp(-ib)),
c1-p=R(1+exp(-ib)),
c2-p=R(1+exp(ib)),
c3-p=R(exp(2ia)+exp(ib)).
```

Their successive actual apex angles are a,b,a. End-spoke and middle-spoke
lengths are respectively `2Rcos(a-b/2)` and `2Rcos(b/2)`. Each outer
base length is `2Rsin(a)`, and the middle one is `2Rsin(b)`. The outer
normalized actual area is `sin(a)[cos(a-b)+cos(a)]`, so its shortfall
from the acute isosceles area cap is exactly
`sigma=sin(a)[1-cos(a-b)]`. Hence the original retained-J score is
exactly `2G(a)+G(b)-2sigma`, with no extra endpoint or lens payment.

The shared-spoke signed skews are `(d,0,-d)` and the two outer weights
are equal. The exact completing-square Qmean minimum therefore has
middle skew zero and is attained. This does not yield a second payable
skew loss.

The current checker computes actual triangle cross-product areas and
these exact base terms, reads the immutable accepted C certificate,
and checks the strict loss interval. Its score expression matches the
analytic one; it does not rerun C or maximize over fan parameters.

## 2. All seven exact faces fill the actual apex star

The isolated cubic root k belongs to (.85,.87), where the polynomial is
strictly increasing. With `u=acos(k) in (0,pi/3)`, its cubic identity
gives `3u=pi-2a`; this follows from cosine injectivity on [0,pi]. The
positive branch also gives `sin(2a)=sqrt(1-k^2)(4k^2-1)`. Substitution
then proves every stated incident site is exactly distance R from its
stated circumcenter, including d0/d1/d2. No interval overlap is being
used to assert an exact equality.

The seven neighbor ray angles can be written as

```
-(a+b/2), -b/2, b/2, a+b/2,
pi-2u, pi, pi+2u.
```

Their positive consecutive gaps are

```
a,b,a,a-b/2+u,2u,2u,a-b/2+u.
```

They sum to 2pi. Positivity follows from b<2a, which is immediate from
`B>2A^2-1` and decreasing cosine, and from u>0. Each gap is <pi.
Thus the listed star sweeps once around p, rather than selecting seven
faces from a possibly larger star.

There is an additional self-contained complete-cell proof. In cyclic
order the seven circumcenters O_j all lie on the circle of radius R
about p. Every neighbor satisfies `neighbor=O_previous+O_next-p`.
The p/neighbor perpendicular bisector is therefore the chord joining
those two O's. These seven p-side halfplanes intersect in the convex
seven-gon `conv{O_j}`; p is strictly inside it because every consecutive
circumcenter-direction gap is <pi. The checked strict exclusion of
each other site at every O means its p-versus-site affine comparison
halfplane contains the whole polygon strictly. It cannot clip any
edge or add a facet. Thus the actual Voronoi cell of p is exactly this
seven-gon, and its actual Delaunay degree is exactly seven.

The source's empty-circle loop includes every one of the 37 sites
except each face's three incident sites, for all seven circles. Strict
outside guards and positive orientations are sufficient. No assumption
about the remaining triangulation or its face circumradii is imported.

## 3. The whole unit circle really is covered

The 29 new sites are exactly `sqrt(1-R^2)` times the unit directions
`j*pi/16`, j=2,...,30. The positive half-angle construction and repeated
rotation give these exact algebraic directions, and their indices are
all distinct. The four local sites have strict norm below that ring;
the four interior sites have norm below 1-R. This proves all 37 sites
are distinct and inside the unit disk in the stated classes.

Each ring footprint has half-width `beta_R=asin R`. The guard
`sin(pi/32)<R` proves strict overlap of consecutive ring footprints.
The two remaining ring/fan joins are exactly the displayed unit points
at angles `+/-[pi/8-beta_R]`, and the guarded local covering distances
are strict.

These joins lie on the correct side of the fan switches, not merely
somewhere in the protected interval. Indeed `beta_R>pi/32` and
`beta_R<.15` (because `sin(.15)>.15-.15^3/6>R`) give

```
-.3 < beta_R-pi/8 < -.2 < arg(q0) < 0.
```

Here use pi<22/7, pi>3 and the source's middle-arc guard. The reflected
ordering holds on the right. The exact endpoint guards and positive-dot
window then cover the four short arcs between the joins and q0,q1,q2:
site dot unit-point is concave throughout this window, so its minimum
is at an endpoint. Together with the ring chain this covers the entire
unit circle. No coverage of its interior is inferred.

Because every circle point is covered by some boundary disk and each
of p,d0,d1,d2 has minimum circle distance >R, none of those four can be
nearest on the circle. Their convex Voronoi cells cannot extend outside
the unit disk without crossing its circle, so these are genuine I sites.

## 4. The three switches are true nearest-owner switches

The formula for q0 solves equal distance and selects the positive CCW
distance-difference derivative. The other unit root is guarded in the
far upper-left quadrant. Reflection and the equal-norm central pair
give q1=(1,0) and q2, with the correct positive transition signs.

The rational rotation has its actual angle alpha in (.092,.1), by the
strict sine/cosine comparisons and its known positive principal branch.
The protected interval from Qleft to its reflection lies in (-.3,.3).
The closest ring directions in that interval are +/-pi/8; their minimal
distances occur at its corresponding outer endpoints. The strict two
exclusion guards therefore rule out every ring site throughout that
interval, not only at the switches. Interior sites are excluded as above.

Among the four local sites, c0 versus c1 has positive derivative over
the negative window and crosses at q0. The c0/c2 comparison likewise
increases and is still strictly favorable to c0 at q0; it remains
favorable on the entire negative outer arc. Equal-norm c0/c3 and c1/c2
comparisons have the required sign from the unit point's y coordinate.
On the negative middle arc, c1 versus c3 satisfies
`q dot v13-delta13<=v13_x cos(.2)-delta13<0`, since v13_x<0,
v13_y>0, cos(theta)>=cos(.2), and sin(theta)<=0. This is exactly the
source's guard13. Reflection supplies both positive arcs.

Thus c0 owns Qleft to q0, c1 owns q0 to q1, c2 owns q1 to q2 and
c3 owns q2 to the reflected terminal point. All are covered, and each
has a strict owned segment longer than .092. No freely chosen covered
point has been substituted for an actual switch.

## 5. Connected assigned arcs and ALL 33 minimum arcs

The connected-arc lemma is sound and applies to every pair of these
distinct sites. For sites x,y of norms <=n=sqrt(1-R^2), let
`v=y-x`, `delta=(|y|^2-|x|^2)/2`. Then

```
|delta|/|v| <= |(x+y)/2| < n.
```

Strictness follows from
`|(x+y)/2|^2=(|x|^2+|y|^2)/2-|x-y|^2/4<V`.
The two bisector unit roots therefore have shorter angular separation
`2acos(|delta|/|v|)>2asin R`. A radius-R disk's connected circle
footprint has length at most `2asin R`, by the inequality
`r+(1-R^2)/r>=2sqrt(1-R^2)` for its center norm r>0. Thus a comparison
has at most one root inside any one site's footprint; its favorable
subset is an interval, prefix, suffix, all or empty. Whole-circle
coverage confines each site's nearest-owned set to that footprint.
Intersecting all comparison intervals gives one assigned interval.
Ties at endpoints do not create positive-length extra components.

Every ring site retains an outward half-cell of length pi/32, directed
away from the local gap near zero. The ring sites have equal norm, so
their ordinary angular Voronoi comparisons leave this half-cell intact.
It lies outside (-.3,.3), where none of the four local sites can cover
or beat the covered ring disk. Interior sites cannot beat a covered
circle point either. Middle ring sites can retain their full cells;
the two end ring sites retain the outward half. Since `pi/32>.092`,
all 29 ring assigned arcs have the required strict minimum length.
Combined with Section 4, this proves the claim for ALL 33 BC sites.

## 6. Honest scope and consequence for the proposed local fee

The three BC/BC marked pairs are actual boundary-owner neighbors, and
their Voronoi bisector segments connect the corresponding exact star
circumcenters to the actual unit switches: squared-distance comparisons
are affine and hold at both ends. The apex's other neighbors are the
three I sites. Therefore p really has the stated q3,d7 local star.

The source-pinned owner certificate encloses the original score loss
near .0183173655881 and below .019; this is below .019931551. Its exact
score formula, guards, normalizations and dependencies were reviewed
above. The portable reviewer is responsible for actual fresh replay.
Consequently the proposed standalone fee cannot follow from just these
listed local conditions. Restricting to an equal-isosceles fan, charging
Qmean again, or using arbitrary switch points cannot repair that claim.

The large recorded middle-cell P loss belongs to an existing boundary
lens budget. It is not an additional fee here. A combined proof can use
it only with one compatible incidence allocation and without spending
it twice. A precisely defined global condition, remaining global faces
or the global profile may also rule out this local diagram's occurrence.

This is NOT a 100-site cover of the closed disk, an unrestricted-radius
counterexample, an outer-Q certificate or refutation, or an optimum.
The owner's history of two plumbing failures and two successful calls
is retained and reported honestly; none of those calls is attributed
to this reviewer. No analytic blocker remains within the stated scope.
