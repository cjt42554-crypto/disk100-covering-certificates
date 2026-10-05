# A continuous-radius actual signed-anchor q2 lemma, with a unit-disk middle cap

This is a NEW author paper proof, awaiting independent review. All old
single-radius sources, reviews and packages remain frozen. The theorem
changes both the radius domain and the middle-center hypothesis; the
new branch argument below explicitly proves the changed hypothesis.
It does not assert a global radius exclusion or an exact optimum.

Fix lambda=7/50 and t=1071883/20000000. Let C be the accepted complete
full-apex scalar support for these exact parameters. For every

    5661/50000 <= R <= 227/2000,       i.e. .11322 <= R <= .1135,

let (p,c0,c1), (p,c1,c2) be consecutive consistently oriented positive
faces, of circumradii at most R. Their ACTUAL common unit anchors q_i
are covered by both endpoint disks. Define

    x_i=|c_(i+1)-c_i|/(2R),
    A_i=-det(c_(i+1)-c_i,q_i-c_i)/2,
    S_i=area_i/R^2+A_i/R^2-4lambda*x_i^2+t*gamma_i,
    tau_i=x_i sqrt(1-x_i^2)-A_i/R^2 >=0,
    D=2C-S0-S1.

Assume the SAME actual middle center satisfies |c1|<=1. Assume each
original actual outer center covers a connected assigned unit-circle
arc of full length strictly greater than .147. Then

    D > 223/10000 = .0223.                         (T)

In particular the theorem applies when the stronger actual cap
|c1|<=sqrt(1-R^2) holds. Neither cap is imposed by moving a center
within the proof. The theorem includes negative actual signed A,
arbitrary actual anchor slack, arbitrary compatible apex slack and
skew. It requires no symmetry, fixed circumradius ratio, apex-position
bound, global Q, degree profile or middle half-arc variance.

All radius-dependent inequalities below are proved on the WHOLE
closed radius interval. The inherited scalar support and paid
localization are dimensionless fixed-parameter statements; their
old certificate's displayed radius is not a radius premise for those
particular scalar statements. Their three-face compatibility conclusion
is not used here. The proof uses no Pro skew convention or Pro .0225
theorem.

## 1. Radius-free high-score budget, and its acute branch

Suppose D<=D0=.0223. The exact score before localization is

    S_i=F(x_i,gamma_i)-sigma_i-tau_i,
    F=x^2 cot(gamma/2)+x sqrt(1-x^2)-.56x^2+t gamma,
    sigma_i=2x_i^2 sin^2(d_i)/sin(gamma_i),
    theta_i=asin x_i, gamma_i=theta_i+z_i, z_i>=0.

The actual sine-rule domain is theta_i<=gamma_i<=pi-theta_i and
|d_i|<(pi-gamma_i)/2<pi/2. The accepted full-apex support gives
S_i<=C before any acute restriction. Thus every individual deficit
is nonnegative and at most D0<1/40. The accepted endpoint cut,
curvature -G''>4 on .7<theta<.96, and paid gamma slope >1 imply

    .7<theta_i<.96,
    2 sum(theta_i-theta*)^2+sum z_i+sum sigma_i+W<=D,
    W=tau0+tau1, B=D-W>=0,
    .828<theta*<.829, .676<cos(theta*)<.677.     (1)

Here G(theta)=sin(theta)+sin(2theta)-.56sin(theta)^2+t theta.
Every formula and constant in this step is independent of R after
normalizing the actual bases and areas by R. The complete support
and localization premises are the pinned accepted input identified
in SOURCE_PINS.json. No regenerated numerical premise is needed.
Limiting zero-area scores are at most 1/2+t*pi<.7<C-D0, so they
cannot enter this high-score contradiction before the divisions.

Put

    s=(theta0+theta1)/2, E=s-theta*, Delta=theta0-theta1,
    Z=(z0+z1)/2, eta=(d0-d1)/2, m=(d0+d1)/2.

Since x_i>.64, the budget first gives
sin^2 d_i< (625/512)D0 <.028. The positive base angles select the
acute inverse-sine branch: sin(.18)>.17 and .17^2>.028 give
|d_i|<.18. Therefore |eta|<.18 and cos(2eta)>0 BEFORE using

    sin^2(m+eta)+sin^2(m-eta)>=2sin^2 eta.

As in the frozen proof, sin(theta+z)<=sin(theta)+z yields
sigma_i>=2(x_i-z_i)sin^2 d_i. Removing only the positive leftover
z_i[1-2sin^2 d_i] terms from (1) gives the SINGLE budget

    4E^2+Delta^2+2 sum x_i sin^2 d_i+W<=D.     (2)

It implies 2.56 sin^2 eta<=D0, hence |eta|<.095, because
sin(.095)>.094 and 2.56(.094)^2=.02262016>.0223. Also

    .75<s<.904, |Delta|<.15, 0<=Z<=.01115,
    alpha=3s/2+Z/2-eta in (1.03,1.456575),
    h=s/2+eta-Z/2 in (.274425,.547).

Consequently 1<alpha<1.46<pi/2, sin(alpha)>.8,
cos(alpha)>.1 and cos(h)>.85. These trigonometric guards are
radius-free lower Taylor bounds. The acute skew ordering and every
zero-product endpoint agree with the pinned clarifications.

## 2. A NEW actual-branch proof under only |c1|<=1

In the original fan's mean phase phi, let b_i be the unit direction
from c1 to outward tip_i, and a_i from the corresponding outer
center to tip_i. The SAME actual skew convention used throughout
this local proof gives

    b0=exp(i(phi-alpha)), b1=exp(i(phi+alpha)),
    a0=exp(i(phi+h+Delta)), a1=exp(i(phi-h+Delta)),
    n0=exp(i phi),
    B_tip=(b0+b1)/2=cos(alpha)n0,
    L0=(b1-b0)/2=i sin(alpha)n0,
    H0=(a0+a1)/2=cos(h)exp(i(phi+Delta)).

These are identities of the actual shared-ray triangles, rather
than independently chosen formal anchors. For actual errors set

    epsilon_i=(q_i-tip_i)/R,
    U=(epsilon0+epsilon1)/2, V=(epsilon1-epsilon0)/2,
    L=L0+V, n=-iL/|L|, u=sqrt(1-R^2|L|^2).

The exact covered-lens geometry gives
|epsilon_i|<=tau_i/x_i^2, whence |U|,|V|<.028. Since |L0|>.8,
|L|>.772, and

    |n-n0|<=2|V|/(|L0|-|V|)<.073,
    n dot n0=1-|n-n0|^2/2>.99,
    k=(B_tip+U) dot n >.1*.99-.028=.071.       (3)

Both ACTUAL anchors are covered by the middle disk, so

    |B_tip+U|^2+|L|^2
       =(|b0+epsilon0|^2+|b1+epsilon1|^2)/2<=1.

Thus u>=sqrt(1-R^2)>.993 uniformly: R<=.1135<.114,
R^2<.013 and .993^2=.986049<1-.013=.987.
The two actual unit anchors have mean A=+-u n and c1=A-R(B_tip+U).
If the sign were negative, then the SAME ACTUAL c1 would satisfy

    |c1|^2=u^2+2Ru k+R^2|B_tip+U|^2
       >1-R^2+2R(.993)(.071)
       =1+R(.141006-R)>1,                    (4)

since R<=.1135<.141006. This contradicts |c1|<=1. The actual
branch is therefore positive, A=u n. This explicitly re-proves the
branch after widening the old sqrt(1-R^2) cap. It uses no radial
reflection, separate movement or presumed fan transport.

## 3. Exact fixed-fan unit-anchor path on the complete radius interval

For 0<=v<=1 define

    L_v=L0+vV, n_v=-iL_v/|L_v|,
    u_v=sqrt(1-R^2|L_v|^2), A_v=u_v n_v,
    c1(v)=A_v-R(B_tip+vU),
    f_i=(c_outer_i-c1)/R=b_i-a_i,
    c_outer_i(v)=c1(v)+Rf_i,
    q0(v)=A_v-RL_v, q1(v)=A_v+RL_v.

Translate p and all fan vertices rigidly with c1(v). This preserves
shape, areas, bases, positive orientation and circumradii. The q's
are unit because A_v is perpendicular to L_v and
u_v^2+R^2|L_v|^2=1. The endpoint offsets are exactly

    (q_i(v)-c1(v))/R=b_i+v epsilon_i,
    (q_i(v)-c_outer_i(v))/R=a_i+v epsilon_i.

They remain in the covered-lens intersection by convexity. Their
middle squared-norm average gives |L_v|<=1 for EVERY v. Therefore

    .772<|L_v|<=1, .993<u_v<=1               (5)

on the whole path, for every admitted R. The v=1 endpoint is the
original actual fan by (4); v=0 is its positive unit-tip comparison.
The actual signed-height slack is v tau_i, so D(v)=B+vW exactly.
The original outer arcs are needed only at v=1. No cap or assigned
arc is assumed transported along the path.

Let J0=(a1-a0)/2, F=L0-J0, H_v=H0+vU, C_v=A_v-RH_v. The physical
outer half-difference is the constant RF. Thus the EXACT physical
mean squared norm is

    M(v)=(|c0(v)|^2+|c2(v)|^2)/2
         =|C_v|^2+R^2|F|^2.                 (6)

The chord error cancels; no separate mean and chord losses are paid.

## 4. Re-proving the ideal mean envelope uniformly in R

For delta=Delta/2, X=x0+x1, Y=x0-x1, minimization over the actual
common skew, used as a relaxation only, yields

    2 sum x_i sin^2 d_i
      >=X-sqrt(X^2 cos^2(2eta)+Y^2 sin^2(2eta))
      >=(8x0x1/X)sin^2 eta cos^2 eta,
    8x0x1/X=4 sin s cos(delta)[1-cot^2 s tan^2(delta)].

The radius-free guards cos(delta)>.997, |tan(delta)|<.076,
cot(s)<1.1, cos^2 eta>.99, sin(s)>.68 give

    4E^2+Delta^2+(1666/625)sin^2 eta<=B.       (7)

Indeed .997[1-1.1^2*.076^2]=.99003200688>.99, and the two factors
>.99 have product>.9801>.98. This is one budget, not another fee.

For eta+=max(eta,0), put

    chi0_R(s,eta+)=s/2+eta++asin(R sin(3s/2-eta+)).

Actual Z removal increases the norm-driving chi, while replacing
negative eta by zero increases chi and lowers its payment. Since
alpha0=3s/2-eta+ lies in (1.03,1.356), its derivatives on the whole
connecting rectangle satisfy

    0<chi0_eta<1,
    chi0_s<1/2+3(.114)(.55)/(2(.993))<.6.

The radius-dependent peak constant MUST be re-proved. The accepted
upper and lower cosine bounds give

    sin(3theta*/2)=(1+2cos(theta*))sin(theta*/2)
       <2.354*.403=.948662.

Now uniformly at R<=227/2000,

    R*.948662 <=.1135*.948662=.107673137
       <.107790048=.108-(.108)^3/6<=sin(.108).

Hence chi0_R(theta*,0)<.4145+.108=.5225<.523. This replaces the
old SINGLE-radius .1075 guard; that old guard is not inherited.
Also eta+<=1.002 sin eta+ on [0,.095]. Cauchy using (7) gives

    chi0_R(s,eta+)
       <.523+sqrt[(310977/666400)(B-Delta^2)]
       <=.523+sqrt[(47/100)(B-Delta^2)].       (8)

The second comparison is NONSTRICT at B-Delta^2=0; the first is
strict from the peak margin at all endpoints.

At the positive tip endpoint,
M(0)=1+R^2-2R cos(chi)cos(Delta). Throughout the paid rectangle,
.375<chi0_R<.662: the upper bound follows also from
asin(R)<.114, because sin(.114)>=.113753076>.1135.
For Delta!=0, sin(3Delta^2/2)>.999(3Delta^2/2) and
tan(chi0_R)>.375 imply

    cos(chi0_R)cos(Delta)
       >=cos(chi0_R+3Delta^2/2).              (9)

The coefficient .375*.999*1.5=.5619375 exceeds .5, which pays
1-cos(Delta)<=Delta^2/2. At Delta=0, (9) is equality.

For 0<=v<=B<=D0, sqrt(.47(B-v))+1.5v is maximized at v=0.
For B>0 this follows by rationalizing the square-root loss, whose
coefficient is at least .47/(2sqrt(.47B))>1.5, since
sqrt(.47B)<.103<.47/3. When B=0, only v=0 is admissible.
Equations (8)-(9) therefore prove, for every admitted R,

    M(0)<=N_R(B),
    N_R(B)=1+R^2-2R cos(.523+sqrt[(47/100)B]). (10)

For B>0, with y=sqrt(.47B)<.103,

    N'_R(B)=R(.47)sin(.523+y)/y >.30.         (11)

The ratio sin(.523+y)/y decreases on (0,.103], and its endpoint
is sin(.626)/.103>.585/.103 by the cubic sine lower bound.
The minimum radius gives
R(.47)(.585)>=.031129839>.0309=.30*.103. Thus (11) is uniform
on the radius interval and on the positive budget interval. N_R
is continuous at B=0; no finite derivative at zero is claimed.

## 5. Both exact lens support branches, with a new uniform tip bound

Use the consistently oriented base tangent e_i and normal n_i toward
the apex. In each actual base frame,

    epsilon_i=s_i e_i+v_i n_i, v_i=tau_i/x_i,
    s_i^2+2x_i|s_i|+v_i^2<=2y_i v_i, y_i=cos theta_i.

Hence |s_i|<=cot(theta_i)v_i. The exact phase identities give

    n0-cot(theta0)e0=-i b0/x0,
    n0+cot(theta0)e0=+i a0/x0,
    n1+cot(theta1)e1=+i b1/x1,
    n1-cot(theta1)e1=-i a1/x1.

Here n0,n1 in these four identities are BASE normals, distinct
from the mean-phase n0 of Section 2. Define without ambiguity

    ell0(g)=max((-i b0) dot g,(+i a0) dot g)/x0^2,
    ell1(g)=max((+i b1) dot g,(-i a1) dot g)/x1^2.

For arbitrary signed g and actual tau,

    g dot epsilon_i<=tau_i ell_i(g),
    |ell_i(g)-ell_i(g')|<=|g-g'|/x_i^2.       (12)

Differentiation of (6) gives

    M'(v)=g0(v) dot epsilon0+g1(v) dot epsilon1,
    g0=-RC_v-(DA_v)^T C_v,
    g1=-RC_v+(DA_v)^T C_v.                   (13)

At v=0, in the mean-phase orthonormal frame, put r=sin(alpha),
u=sqrt(1-R^2r^2). The EXACT gradients are

    g0(0)=R(-P0,Q0), g1(0)=R(-P1,-Q1),
    P0=u-Rcos(h)cos(Delta)+(u/r)cos(h)sin(Delta),
    P1=u-Rcos(h)cos(Delta)-(u/r)cos(h)sin(Delta),
    Q0=R[r+cos(h)sin(Delta)-(Rr/u)cos(h)cos(Delta)],
    Q1=R[r-cos(h)sin(Delta)-(Rr/u)cos(h)cos(Delta)].

Uniformly, P_i>.993-.114-.15/.8=.6915 and Q_i<.114(1+.15)<.132.
For w_i=alpha-theta_i, .205<w_i<.628, so
P_i tan(w_i)>.6915*.205=.1417575>.132. This selects the b branch
at the TIP endpoint, without claiming any positive-v tangent sign.
Exact substitution gives

    ell0(g0(0))=R(u-Rcos alpha)/x0^2
       [sin alpha(1-Rcos h cos Delta/u)+cos h sin Delta],
    ell1(g1(0))=R(u-Rcos alpha)/x1^2
       [sin alpha(1-Rcos h cos Delta/u)-cos h sin Delta]. (14)

The factor u-Rcos(alpha) is positive and strictly below .989,
since R>=.11322 and cos(alpha)>.1. Also
R cos(h)cos(Delta)/u>.11322*.85*.988=.095082156>.095,
so the common bracket factor is positive and below .905.

Put a=|Delta|<=.15. For Delta=a>=0, the first numerator is below
.905+a, whereas

    2x0^2>1-cos(1.5+a)>.929+.996a>.905+a.

For Delta=-a<=0, its numerator is below .905-.8466a and

    2x0^2>1-cos(1.5-a)>.929-a>.905-.8466a.

The smallest latter gap is .024-.1534*.15=.00099>0. These are
radius-free guards. They use cos(1.5)<.071, sin(z)>.996 on
[1.5,1.65], and sin(a)>=.996a. The sine interval follows from
concavity and seventh-order LOWER Taylor polynomials at both
endpoints, as fully recorded in the pinned SINE_INTERVAL_GUARD.
At a=0 the zero-product step is nonstrict; the starting gap remains
strict. Reflecting Delta's sign proves the identical x1 bound.

Thus each bracket in (14), divided by x_i^2, is below 2. On the
NEW whole radius interval the tip bound is

    ell_i(g_i(0))<2R(.989)<=2(.1135)(.989)
       =.224503<.225.                       (15)

The old .224 coefficient is NOT claimed above R=.11322.

## 6. Full-path operator estimates, re-proved on the radius interval

For L on the full path, write r=|L|, u=sqrt(1-R^2r^2) and
A(L)=u(-iL/r). By (5), .772<r<=1 and .993<u<=1. The radial
and tangent singular values of DA, up to an orthogonal rotation,
are |u'|=R^2r/u and u/r. Hence ||DA||<1.3: the latter is at most
1/.772<1.3, while the former is below .013/.993<.014.

The normalization w=L/r has ||Dw||<=1/r and

    D^2w[h,k]=-[ (w dot h)Pk+(w dot k)Ph
                         +w(Ph dot Pk)]/r^2,
    P=I-ww^T,

so ||D^2w||<=3/r^2. The scalar radial Hessian has operator norm
at most |u''|+|u'|/r. Product differentiation therefore gives

    ||D^2A||<=|u''|+3|u'|/r+3u/r^2,
    u''=-R^2/u-R^4r^2/u^3.

Every guard below is uniform using R^2<.013, rather than an old
fixed-radius Hessian receipt:

    |u''|<.013/.993+.013^2/.993^3<.014,
    3|u'|/r<3(.013)/.993<.04,
    3u/r^2<3/(.772)^2<5.04.

Their sum is below 5.094<5.2. These are operator bounds along
the ENTIRE segment from L0 to every L_v. Consequently

    |A_v-A_0|<=1.3v|V|,
    ||DA_v-DA_0||<=5.2v|V|.

Differentiating |A|^2+R^2|L|^2=1 gives (DA)^T A=-R^2L exactly.
This retains the essential factor R in both gradients:

    g0=R[-A+RH+RL+(DA)^T H],
    g1=R[-A+RH-RL-(DA)^T H].                 (16)

Since |H0|<=1 and H_v-H0=vU, splitting the product difference
using H0 and DA_v gives

    |g_i(v)-g_i(0)|/R
      <=(1.3+R+5.2)v|V|+(1.3+R)v|U|
      <6.614|V|+1.414|U|
      <=8.028(625/512)W.

Using R<.114, W<=D0 and x_i^2>.4096 proves the following NEW
uniform interval estimate with the same loose rational guard:

    |g_i(v)-g_i(0)|/x_i^2
      <(.114)(8.03)(1.221)(.0223)/.4096<.061.

The exact numerator .024925330386 is below
.061*.4096=.0249856. No separately paid U or V budget was introduced;
their coefficients were summed against the SAME W. When W=0,
the gradient difference is zero and the fixed positive .061 bound
still holds strictly. By (12) and the NEW tip bound (15),

    ell_i(g_i(v))<.225+.061=.286<.30,
    M'(v)<=.286W                             (17)

on the whole path and radius interval. No positive-v tangent sign
is presumed; both support branches are transported by (12).
When W=0, both errors and M' vanish. Its payment is equality then.

## 7. One-budget integration and a NEW full-arc lower bound

Integrating (17), then using (10)-(11) and continuity at B=0, gives

    M(1)<=M(0)+.286W
       <=N_R(D-W)+.286W<=N_R(D)<=N_R(D0).     (18)

The angular/skew/asymmetry budget is B; the actual anchor budget
is W; D=B+W. Every variable is spent once. No BC-lens, global Q,
middle variance or other fan fee appears in (18).

Since .47D0=.010481<.103^2, its effective angle is below .626.
The sixth-order cosine LOWER Taylor polynomial gives cos(.626)>.81.
Therefore

    N_R(D0)<G(R)=1+R^2-1.62R.

On the entire interval G'(R)=2R-1.62<0. Its worst value is thus
at Rmin=.11322, and

    M(1)<G(.11322)=.8294023684.               (19)

The old radius-specific outer norm floor .911088875 CANNOT be
inherited at R=.1135. Here is a new direct geometric proof.
For beta=.147/2=.0735, any center of norm rho covering a connected
unit-circle arc of full length at least .147 must satisfy

    rho>=cos(beta)-sqrt(R^2-sin^2(beta)).     (20)

To justify (20) without assuming a radially centered assigned arc,
consider its two endpoints, of arguments z+-beta after choosing
a length-.147 subarc. Their larger squared distance from a center
of argument psi is at least 1+rho^2-2rho cos(beta): their minimum
cosine is at most cos(beta). For example their average cosine is
cos(z-psi)cos(beta)<=cos(beta). Endpoint coverage then gives
rho^2-2rho cos(beta)+1-R^2<=0, whose lower root is (20).

The NEW explicit Taylor guards are

    sin(beta)>beta-beta^3/6=.0734338224375>.07343,
    cos(beta)>=1-beta^2/2=.997298875,
    R^2-(.07343)^2<=.01288225-.0053919649
       =.0074902851
       <(.086548875)^2=.007490707763765625.

Thus (20) gives rho>.997298875-.086548875=.91075 uniformly at
every admitted R. Each original actual outer arc satisfies this
bound, so

    M(1)>(.91075)^2=.8294655625.              (21)

The positive margin between (21) and (19) is

    .8294655625-.8294023684=.0000631941>0.

This contradicts D<=D0 and proves (T) for the entire continuous
radius interval, with arbitrary actual tau and only |c1|<=1.

## 8. Scope, quantitative bottlenecks and proof record

The dominant limitation in this deliberately loose interval proof
is the combination of the worst-case norm upper bound at Rmin with
the worst-case arc lower bound at Rmax. The .0000631941 squared-norm
gap is far smaller than the derivative-payment reserve .30-.286=.014.
It is not evidence of a sharp mathematical radius obstruction.

Other individually identifiable interval guards are:

* The peak guard needs R*.948662<.107790048. Its numerator margin
  at Rmax is .000116911. A different peak constant requires its own
  proof, rather than changing .523 by inspection.
* The bracket guard <.905 needs Rmin*.85*.988>.095; its strict
  margin at Rmin is .000082156. This currently limits radius
  extension downward more than the full-tau payment does.
* The uniform Hessian used R<.114 and R^2<.013. The tip coefficient
  was changed to .225, while its entire-path perturbation remains
  <.061. None of these steps is a sampled radius check.
* The envelope slope >.30 follows whenever the adopted guards hold
  and R>206/1833, because R*.47*.585>.30*.103. This exact sufficient
  threshold is below Rmin and is not a new global lower bound.
* The NEW negative-branch exclusion under |c1|<=1 only additionally
  needs R<.141006 once the adopted angle/error guards hold. Its
  proof is (3)-(4); the middle cap remains a hypothesis on the
  actual fan rather than an asserted transport theorem.

Relaxing connected outer arcs, the actual endpoint coverage, the
fixed-parameter complete support, or the actual middle unit-disk
cap changes the theorem and requires new proof. General finite
covering incidence and loss allocation remain external. Adding a
separately allocated nonnegative middle variance can only increase
the proven D, but no particular joint target at another R is
certified here. No global bound is inferred from this lemma alone.

Author mathematical program executions, mathematical imports,
old checker replays, optimizers, parameter scans and searches: 0.
Only source reading, new document writing and SHA256 metadata
operations were used. In particular the old 48 Pro guards and the
accepted single-radius local checker were not run. All inherited
scalar certificates remain exact pinned premises, and all new
continuous estimates above are paper inequalities awaiting a new
independent full-source audit.
