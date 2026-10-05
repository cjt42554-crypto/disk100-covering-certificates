# Independent static review of the continuous-radius Pro-domain lemma

Reviewer: `/root/q2_uniform_review`. Verdict: **STATIC PASS** for
`research/q2_radius_uniform_tau/PROOF.md`, SHA256
`692a4526d88dd3e856aa96b60961b43c1726d9ba029803199b90e8d7b0c9da2c`,
with its separately frozen analytic slab-transfer note. This verdict
is independent of both the main and slab authors. Mathematical
executions, mathematical imports, old/new checker runs, parameter
scans and optimizers: **0**. All operations were source reading,
new review writing and metadata hashes; frozen inputs were unchanged.

The accepted theorem fixes lambda=7/50 and t=1071883/20000000 and
covers the complete CLOSED radius interval [5661/50000,29/250].
The two faces must be the same consecutive positive oriented actual
triangles, with circumradii at most the current R and actual unit
anchors covered by the corresponding endpoint disks. On the SAME
actual centers it assumes |c1|^2<=1-R^2 and both outer squared norms
at least 83/100. Under these hypotheses the matched signed-A deficit
satisfies D>9/400, including arbitrary actual tau and negative A.
No variance, Q, BI, independent J, or additional copy of tau is paid.
The theorem is conditional local geometry, not a global radius bound.

## Radius-independent support and actual branch

I read the complete accepted Pro proof and its independent review,
then independently checked every radius-dependent argument in this
new source. The normalized full-apex support, unique peak, paid
localization, skew convention and single deficit budget contain no R.
They can therefore be retained as the pinned scalar premises under
D<=9/400, without transferring a pointwise geometric estimate.

The phase convention alpha=3s/2+Z/2+eta, h=s/2-Z/2-eta matches the
source's 2k cos(gamma/2 +/- d) side-length convention. The exact
normal/tangent lens error and its signed support retain negative A.
The common anchor equations give c1=+/-u n-Rk; this eliminates
translation exactly rather than choosing radial-inner centers.

Every segment v e_i remains in each convex covered endpoint lens.
Its two middle-offset vectors have squared-norm average
|k(v)|^2+|b(v)|^2<=1, so |b(v)|<=1 on the whole segment. The lower
bound 6/7-27/1040>83/100 and current-R bound u^2>=1-R^2>.99^2
hold uniformly. Also n_x=(B0+d_y)/|b|>=B0-|d_b|>.83, since
|b|<=1 and the numerator is positive. Thus k dot n>.083-27/1040>0.
On the negative unit-anchor branch the actual norm satisfies
|c1|^2=u^2+2Ru(k dot n)+R^2|k|^2>u^2>=1-R^2,
contradicting the SAME actual cap. The positive branch is forced.
No actual arc or cap is presumed preserved along the comparison path.

## Signed support and transfer of all eight complete slabs

The four exact linear gradient components remain algebraic identities
at the variable R. The upper x-component and opposite y-component
bounds yield the sufficient tangent-sign margin

    R[163/1000-(27/20)R]>0

on R<=.116. The sign proof uses sin A_i>1/5 and the x-component
negative bound, so it retains the smaller directional lens support
instead of substituting a gradient norm.

I checked the continuous slab transfer separately from the point
certificate. At fixed a in the full angular range put
rho=u-R cos a, u=sqrt(1-R^2 sin^2 a). The endpoint enclosures
remain valid for R<=.12 because

    rho_a=R sin a(1-R cos a/u)>0,
    (rho sin a)_a=cos a(1-2R^2 sin^2 a)/u-R cos(2a)>0.

The first strict sign uses u^2-R^2 cos^2 a=1-R^2; the second uses
a in (pi/4,pi/2) and 1-2R^2 sin^2 a>0. Hence each full s slab,
with its entire omega interval, has the stated A-, Q-, P-, P+ bounds
at every radius. Dropping 1/u from the positive Q term is safe.
The complete Delta range is enclosed by the min cosine and max
absolute-value expression, and reversing Delta covers the other face.

Direct differentiation gives H=(Rrho)'=u-2R cos a-R^2 sin^2 a/u
between zero and one, and (R^2rho)'=R[2u-3R cos a-R^2 sin^2 a/u]>0
throughout [.11322,.12]. Therefore A- loses at most one per unit R,
Q- is nondecreasing, and each P endpoint is 1-Lipschitz. Absolute
values and maxima preserve that last bound. The map
q -> min(q,q cos(.15)) is nondecreasing, including when q crosses
zero. These facts prove, for each ENTIRE slab,

    L_j(R)>=L_j(R0)-(23/20)(R-R0).

The original outward rational certificate bounds L_j(R0) from below,
so its accepted margin >1/4000 is a valid frozen starting premise.
No outward-rounding recipe is implicitly transported to another R.

Crucially, the new coefficient is added to the EXACT residual:
G_i(k)=G_i(.205)+(k-.205)x_i^2. It is not added to a redefined
conservative slab function. With k=.213 and x_i^2>13/30, all eight
slabs have strict residual

    1/4000-(23/20)(139/50000)+(1/125)(13/30)
      =1559/3000000>0.

Thus kappa_i<.213 holds on the full closed radius band. The sharper
interface .205+(69/26)(R-R0) follows from dividing the same residual
loss by x_i^2>13/30. The slab child additionally gives a valid
9/8-Lipschitz version: H<1-R0/5 makes its loss coefficient smaller
than 9/8. These are analytic bounds on continuous functions, not
new sample evaluations of the slabs.

## Radius-dependent Hessian and complete tau remainder

For f(b)=u Jcw b/|b| the radial/tangent derivative singular values
are R^2|b|/u and u/|b|. The path domain proves ||Df||<=100/83
uniformly. Writing w=a e+v t, direct second differentiation gives

    D^2f[w,w]=-Jcw[(R^2 a^2/u^3+v^2/(u|b|^2))e
                                      +2av/(u|b|^2)t].

The elementary square inequality
sqrt(v^4+4a^2v^2)<=(2/sqrt(3))(a^2+v^2) gives the displayed
directional Hessian bound. On the current band R^2<.014,
u^3>.97 and u|b|^2>17/25. Therefore

    R^2/u^3+2/(sqrt(3)u|b|^2)
       <14/970+175/102<7/4.

Expanding the exact mean-norm identity after its linear part leaves
-R^2|d_b|^2+R^2|g|^2+2R[r_f dot E*-(f(b)-f(b*)) dot g].
The favorable first term may be dropped; the other terms all appear
in the source's remainder. Taylor's one-half factor and |E*|<=1
produce the stated coefficient. Its entire-radius upper bound is

    [.014+.116*(104/25)]*(225/169)
      =111726/169000<2/3.

Consequently M<=M*+.213T+(2/3)T^2, with T<=9/400, and the single
tau cost is at most .228T<=.23T. This cost is nonstrict at T=0;
the final norm gap supplies strictness without a zero-product error.

## New baseline and the decreasing norm envelope

Discarding nonnegative Z increases chi, since its exact derivative
is -1/2+R cos alpha/(2u)<0. On the entire connecting rectangle
alpha>1, u>.99, R<.12, so q=R cos alpha/u<13/198 and
0<g_s<.6, -1<g_eta<0. The original .522 baseline cannot simply
be kept at another R. The new source instead integrates
g_R(theta*,0)_R<100/99 and uses

    b(R)=.522+(100/99)(R-R0).

Cauchy with the same paid shape/skew budget gives
chi<=b(R)+.71 sqrt(E). At the largest E and R its angle is
124999/198000<.632. It is positive and acute; the quadratic cosine
lower bound at .632 exceeds .8. Thus the cosine/sine monotonicity
used later is justified throughout the band.

The Delta correction costs at most R Delta^2. Since R<.23,
Delta^2 and T are merged from the SAME remaining budget D-E.
Phi_R(E) therefore uses .23(d-E), not two separate worst-case
copies of d. The frozen scalar derivative gap >.04 for a .22 cost
loses exactly .01 and remains >.03. Continuity handles E=0.
The maximum is E=9/400, giving

    F(R)=1+R^2-2R cos[.6285+(100/99)(R-R0)].

Its derivative is strictly below
.232-1.6+.168=-1.2: the last term comparison is exactly
229100<259875. Hence F decreases on the COMPLETE radius interval.
Its frozen left endpoint is below .83-.0003. Both actual outer
squared norms >=.83 force their mean >=.83, giving the contradiction.
There is no numerical maximization over R in this argument.

## Review scope and optional guard clearance

The new proof and slab note are independently sound. Metadata pins
for the author manifest, source record and all thirteen accepted
dependencies are verified separately. The proposed finite rational
GUARD_PLAN is consistent with the proof and is cleared as an OPTIONAL
new isolated uniform-domain guard. This review does not execute it,
require replay of the old 48 guards, or count an old replay as new.

The current verdict retains the middle cap <=sqrt(1-R^2). It does
not yet accept a replacement by <=1; that needs its new branch proof.
Nor does fixed full arc>.147 imply the required outer norm floor
throughout R<=.116. The source correctly states the arc-to-norm root
decreases with R and must be handled separately. No parameter
rectangle, changed geometric domain, global fan admission, global
budget or exact optimum is inferred by this review.
