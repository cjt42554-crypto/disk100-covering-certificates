# Static independent analysis of the two-face unit-anchor compatibility family

The proposed fee `>.02228272803` cannot follow from the local two-face
common-apex/unit-anchor hypotheses plus normalized centers and three
assigned boundary arcs having `xi=sin(fullarc/2)/R>.6` alone. The
mirrored skew family below has all those local properties and loss
`x/36<37/1800=.020555...`. Its external caps are genuine full
ball/unit-circle footprints, are consecutive with the middle cap,
and have unique nearest owners in their interiors among the four
specified sites. This is a local counterexample to that proposed
strength, not to every possible smaller positive fee, nor to a
theorem with additional degree-7/whole-cover assumptions.

This analysis is entirely static. No checker, optimizer, search,
mathematical module or new guard program was executed. The scalar
peak is the already accepted peak for the exact coefficients below;
its saved enclosure is an explicit accepted dependency. The final
chosen family is the radius_scope author's `sin epsilon=1/12`
point. The other author's rational-angle skew point obeys the same
geometric endpoint argument; I do not substitute its approximate
scalar loss or footprint values for a separate exact certificate.

## Parameters and distinction between base x and assigned-arc xi

Fix

    R=5661/50000=.11322,
    lambda=7/50,
    t=1071883/20000000.

Let `theta=theta*` be the accepted unique acute maximizer of

    G(theta)=sin theta+sin(2theta)-4 lambda sin(theta)^2+t theta,

so `G(theta*)=C_lambda(t)`. Write

    y=cos theta, x=sin theta,
    c=cos(theta/2), s=sin(theta/2),
    Q=1+2y, V=1-2y.

Here `x` is the triangle base parameter `|C_i-C_(i+1)|/(2R)`.
It is different from `xi=sin(assigned fullarc/2)/R`.
The accepted bracket is `.676<y<.677`; in particular

    .73<x<.74,
    .9154<c<.9157,
    .4018<s<.4025.

The half-angle bounds follow directly by squaring their rational
endpoints using `c^2=(1+y)/2`, `s^2=(1-y)/2`. They are not the
result of a newly executed mathematical calculation.

Choose the final skew point

    B=sin epsilon=1/12,
    A=cos epsilon=sqrt(143)/12.

Then `A^2+B^2=1`, `A>.9965`, `0<epsilon<theta/2`.
Define

    j=c A+s B=cos(theta/2-epsilon),
    k=c A-s B=cos(theta/2+epsilon),
    l=s A+c B=sin(theta/2+epsilon),
    T=Q c A-V s B,
    H=Q s A+V c B,
    M=V c A-Q s B=2k-T,
    u=sqrt(1-R^2 H^2).

The exact algebraic identities

    k^2+l^2=1,
    M^2+H^2=1,
    T-2j y=k,
    2j x-H=l

will be used below. The second follows from
`V^2 c^2+Q^2 s^2=1`; the A/B cross terms cancel.

## Actual sites and tips

Place the common apex and three neighboring centers at

    p=(u-RT,0),
    C0=p+2Rj(y,-x),
    C1=p+(2Rk,0)=(u+RM,0),
    C2=p+2Rj(y,x).

The actual unit anchors are

    q0=(u,-RH),
    q1=(u,RH).

The two triangles `(p,C0,C1)` and `(p,C1,C2)` have positive areas,
successive apex angles `theta`, and skew parameters `-epsilon`
and `+epsilon` respectively. The sine-rule side lengths give
circumradius exactly R for both faces. The shared ray `pC1` is
exactly the same side of both triangles, so no independent-length
compatibility has been suppressed.

The outer lens tips of those faces relative to their actual apex
are `R(T,-H)` and `R(T,H)`. This follows either from the accepted
complex actual-tip formula at bisectors `-theta/2,+theta/2`, or
directly from the midpoint and outward-normal construction. Thus
q0/q1 are their actual outward lens tips, not free surrogate points.
In particular

    |q0|=|q1|=1,
    C0=q0-R(k,l),
    q0-C1=R(-M,-H),

with the mirrored identities for C2 and q1. The identities above
give

    |q0-C0|=|q0-C1|=|q1-C1|=|q1-C2|=R.

Hence every required anchor is genuinely common to the two endpoint
disks. It is on the outward side of its actual base and has zero
actual-anchor slack `tau`. The normalized actual signed anchors
are `xy` on both bases.

## Useful elementary bounds and normalized-center condition

The accepted y bracket and the A/B choice imply the following loose
but sufficient rational bounds:

    2<T<2.2,
    .9<H<1,
    .87<k<1,
    .476<l<.5,
    -.42<M<0,
    .99<u<1,
    .113<R<.114.

For example,
`l>.4018*.9965+.9154/12>.476` and
`H>2.352*.4018*.9965-.354*.9157/12>.9`.
These are plain rational product inequalities after the displayed
radical endpoint bounds. They are stated transparently so a future
standalone new guard can check them without any old checker.

They immediately give

    0<p_x=u-RT<1-R,
    0<(C1)_x=u+RM<1.

For C0, put `rho=|C0|`. Since `C0=q0-R(k,l)`,

    rho^2=1+R^2-2R(uk-RHl).

The parenthesized quantity exceeds
`.99*.87-.114*.5=.8043`, so `rho<1`. In fact

    rho^2<1+.114^2-2*.113*.8043=.8312242<.92^2,

and hence `rho<.92`. Also `rho>1-R`: equality in the reverse
triangle inequality would require parallel q0 and `(k,l)`, whereas
their determinant is `ul+RHk>0`. The same statements hold for C2.
All four actual sites p,C0,C1,C2 therefore lie strictly in the unit
disk. The common apex cannot cover any unit-circle point, because
for every such z,

    |z-p|>=1-|p|>R.

## Actual full-cap endpoints and orientation

Let

    beta=asin(RH),
    q0 angle=-beta,
    q1 angle=+beta.

Both angles lie in `(−pi/2,pi/2)` since u>0 and H>0.

The middle center C1 lies on the positive real axis. Since q0,q1
are its two unit-circle intersections, its entire closed footprint
is exactly the short CCW cap

    I1=[q0,q1]=[-beta,+beta].

This can also be checked without any direction convention:
`|z-C1|<=R` is equivalent to `cos(arg z)>=u` because C1>0
and equality holds at q0/q1. Thus q0 is its lower endpoint and
q1 its upper endpoint.

For the lower outer center,

    det(C0,q0)=R(lu+RkH)>0,
    q0 dot C0=1-R(uk-RHl)>0.

Therefore q0 is counterclockwise of the radial direction of C0,
with an offset strictly between 0 and pi/2. Because
`1-R<rho<1`, the disk centered at C0 meets the unit circle in one
short nontrivial cap. Consequently q0 is its upper cap endpoint,
not its lower endpoint and not an interior point.

Its other genuine unit-circle endpoint is the exact reflection
of q0 across the radial C0 axis:

    w0=2[(q0 dot C0)/rho^2] C0-q0.

Reflection preserves both the unit norm and distance from C0.
It reverses the determinant sign. If `eta` is the cap half-angle,
the angles are

    arg(q0)=arg(C0)+eta=-beta,
    arg(w0)=arg(C0)-eta=-beta-2eta.

The full C0 footprint is exactly

    I0=[w0,q0].

Let w1 be the reflection of w0 across the real axis. Mirroring
gives the full C2 footprint

    I2=[q1,w1].

No wrapping or long-arc ambiguity occurs. Indeed
`sin eta=det(C0,q0)/rho<=R`, while eta is acute. Both beta and
eta are less than `asin R`; since R<.114, all displayed endpoints
are within the right half-plane. Thus the consecutive CCW order is

    w0, q0, q1, w1,

and the three entire caps have disjoint interiors, share just the
specified adjacent endpoints, and form one continuous closed arc.

## Complete assigned arcs and nearest-site interference

The full footprints I0,I1,I2 are valid complete assigned arcs among
the four specified sites. In the interior of I0 the distance to C0
is strictly less than R; the point lies outside I1 and I2, so its
distances to C1 and C2 are strictly greater than R. Its distance to
p is also greater than R by the preceding radial bound. Thus C0
is the unique nearest site there. The same argument holds in I2,
and C1 is the unique nearest site in the interior of I1.

At q0 the only nearest tie among these sites is C0/C1, both at
distance R; at q1 it is C1/C2. The opposite outer site has its
cap entirely on the other side of the positive middle cap and
cannot tie or intrude there. A canonical tie convention may assign
those individual endpoint points either way without changing any
arc length. At w0/w1 only the corresponding outer site is at
distance R among the four sites.

These facts resolve the possible local Voronoi objection. They do
not say that unspecified additional sites in a complete cover
would preserve these arcs. The present example only asserts the
local four-site nearest/assigned configuration. A global embedding
must separately control all additional sites.

## Assigned-arc xi values and the target loss

The actual middle full cap has

    xi_middle=sin beta/R=H>.9>.6.

For each outer full cap, the determinant expression gives

    xi_outer=sin eta/R=(lu+RkH)/rho.

Using the loose bounds established above,

    lu+RkH>.476*.99+.113*.87*.9=.559719,
    rho<.92.

Since `.559719>.608*.92`, this proves statically

    xi_outer>.608>.6.

The radius_scope author intends a sharper standalone enclosure
around `.62`; it is not needed for the local counterexample here.
All three actual complete assigned arcs satisfy the proposed
strict `.6` condition. The triangle bases independently have
`x=sin theta*>.73>.6`.

The actual area loss on each face is

    sigma=2x^2 sin(epsilon)^2/sin theta=2x B^2.

The apex potential and quadratic base potential are unchanged,
the actual tip anchor has tau=0, and theta is the accepted scalar
peak. The exact pair score is therefore

    2C-4x B^2=2C-x/36.

Its total deficit satisfies

    D=x/36<.74/36=37/1800<.02228272803.

Thus a universal fee above `.02228272803` under only the stated
local arc/base/normalization assumptions would exclude this
explicit admissible pair and is false. This conclusion retains
the accepted scalar-peak dependency rather than pretending its
saved interval was independently regenerated.

## Actual local Delaunay fan check

The circumscribed-circle centers of the lower and upper faces are
respectively

    O0=p+R(k,-l),
    O1=p+R(k,l).

Each has radius R and contains its three stated vertices. The
other outer site is strictly outside that circle; for example

    |C2-O0|^2-R^2=8R^2 j l x>0.

This follows from
`j=cos(theta/2-epsilon)=ky+lx`. The mirrored identity holds
for C0 and O1. Thus the two triangles are also locally compatible
with the empty-circumcircle condition for these four sites, and
their shared common-apex fan is not merely an algebraic collection
of unrelated triangles.

The common apex has no degree-7 certificate, and these four disks
do not cover the whole unit disk. The example does not establish
a permissible global profile or any whole-cover topology. Such
additional hypotheses could still force a stronger conditional
two-face loss.

## Exact two-anchor curvature interface for a future lemma

For arbitrary actual anchors with relative vectors
`v_i=(q_i-p)/R`, define

    K=(v0+v1)/2,
    L=(v1-v0)/2.

The two unit equations give the exact conditions

    (p+RK) dot L=0,
    |p+RK|^2+R^2|L|^2=1.

For L nonzero, these determine

    p=+/-sqrt(1-R^2|L|^2) n-RK,

where n is a unit normal to L. They do not furnish the three-point
sagitta condition used in the accepted three-face proof.
The shared middle center, with `W=(C1-p)/R`, additionally obeys

    |W-K|^2+|L|^2+2|(W-K) dot L|<=1,

equivalent to both anchor disk inequalities. Any other center
with relative vector Wj has exact physical norm

    |Cj|^2=1-R^2|L|^2
          +/-2R sqrt(1-R^2|L|^2) n dot(Wj-K)
          +R^2|Wj-K|^2.

This is a concrete interface for testing a future extra radial,
degree or incidence hypothesis. The mirrored skew example shows
that this interface plus `xi>.6` alone remains too permissive for
the parent's target fee. The unreviewed external Q bound is not
used anywhere in this analysis.

## Corollary: zero local assigned-cap lens slack

All three assigned arcs above are exactly their disks' complete
short caps. Let a denote one cap's HALF angle, e its unit radial
midpoint direction, and rho its center norm. Its center is `rho e`,
and both cap endpoints have distance R from that center. Therefore

    (rho-cos a)^2=R^2-sin(a)^2.

To select the inner radial branch one must establish
`rho<cos a`. Center normalization `rho<1` alone is insufficient
in general: an outer branch close to the maximum possible cap
width can also have center norm below 1.

For the middle cap in this family,

    rho_middle=u+RM<u=cos beta,

because M<0 and its half-angle is beta. For each outer cap,
`rho_outer<.92`, whereas `sin eta<=R<.114` and eta is acute,
so `cos eta>=sqrt(1-R^2)>.99`. Hence all three actual centers
are strictly on the inner branch:

    rho=cos a-sqrt(R^2-sin(a)^2).

The signed area of the curvilinear wedge bounded by the two
center-to-anchor segments and the CCW unit arc is exactly

    L_actual=a-rho sin a.

One may derive it by adding the unit-sector arc integral a and
the two straight-segment shoelace contributions `-rho sin a`.
Its independent inner-radial support at the same assigned
half-angle is

    F_R(a)=a-sin a cos a+sin a sqrt(R^2-sin(a)^2).

Their exact difference is

    F_R(a)-L_actual
       =[rho-(cos a-sqrt(R^2-sin(a)^2))] sin a=0.

Thus all three assigned-cell lens slacks of this particular
complete-cap configuration are exactly zero. The two marked
edge anchors independently satisfy `A_i=J_R(L_i)` and `tau_i=0`,
since each is the outward tip of its base's radius-R intersection
lens. There is no unspent positive local assigned-cap lens loss
in this example that could rescue the rejected standalone
actual-A target fee under the same local domain.

This corollary does not consume or analyze any global arc
variance, window, radial-Q or profile constraint.

## Scope relative to the parent's new joint objective

The parent subsequently proposed the distinct joint deficit

    D_joint=D_actual+K(a_mid-pi/33)^2,
    K=(33/20)(33/32)/R^2,

where a_mid is the middle assigned arc's HALF angle. In this
family `a_mid=beta=asin(RH)`. The counterexample analysis above
concerns D_actual alone. It does not bound D_joint below the
parent's target and does not refute that joint claim. Zero
assigned-cap lens slack likewise says nothing about the separate
arc-variance payment. Additional global arc/window/Q/degree
hypotheses remain independent obligations. Any different fixed
diagnostic point for that new joint objective requires its own
static analysis and fresh guard evidence; no such execution was
performed here.
