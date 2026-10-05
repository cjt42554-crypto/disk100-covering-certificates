# Actual-tau mean-norm envelope and a two-outer-arc sum interface

New root-authored analytic lemma, awaiting independent static review. It
does not modify any accepted source. No mathematical program, imported
source, old guard replay, scan or optimizer was run. All finite decimals
below are exact rationals. This lemma removes a conservative individual
outer-norm premise by exposing the quantitative mean envelope already
proved before that premise is used.

Fix lambda=7/50, t=1071883/20000000, R0=5661/50000, and

    R0<=R<=29/250=.116,   d*=9/400=.0225.

Let the two actual consistently positive consecutive triangles have
circumradii<=R and their own actual unit anchors covered by their base
endpoints. Use their actual signed A, arbitrary nonnegative actual tau,
the matched complete support C, and D=2C-S0-S1. Impose on the SAME actual
middle center only

    |c1|^2<=1+3/2000.

No outer norm or assigned-arc premise is needed for the first theorem.
The radius theorem, its independent review, and the independently
reviewed relaxed-branch interface are the frozen inputs in SOURCE_PINS.
They supply the dimensionless budget, exact positive actual branch,
signed directional supports and complete remainder on this entire R
interval. Their individual-outer-floor conclusion is not assumed here.

Define for 0<=d<=d*

    b_R=261/500+(100/99)(R-R0),
    beta_R(d)=b_R+(71/100)sqrt(d),
    H_R(d)=1+R^2-2R cos(beta_R(d)).                    (1)

The new quantitative theorem is

    D<=d*  ==>  M=(|c0|^2+|c2|^2)/2 < H_R(D).         (2)

In particular it needs neither individual outer floor .83 nor equality
of the two outer norms. At d=d* the inherited uniform envelope implies
H_R(d*)<=H_R0(d*)<.8297, but (2) retains its actual R and D dependence.

## 1. Reclose the actual budget using D, rather than the cutoff d*

Assume D<=d*. Complete scalar support gives D>=0 before localization.
The frozen radius proof establishes, without any outer norm hypothesis,

    4(s-theta*)^2+c_eta eta^2+Delta^2+2Z+T<=D,
    c_eta=223532401/90000000, T=tau0+tau1,
    M<=1+R^2-2R cos(chi)cos(Delta)+(23/100)T,
    0<chi < b_R+(71/100)sqrt(E),
    E=4(s-theta*)^2+c_eta eta^2.                       (3)

The strict chi bound follows from the strictly smaller peak starting
value; it persists at E=0. The coefficient comparison and all zero
products themselves are nonstrict. All involved angles lie between 0
and .632<pi/2, so cosine monotonicity has the required sign.

The relaxed middle cap is used solely to select the actual positive
unit-anchor translation. Its independent proof excludes the negative
branch by |c1|^2>1+.00158508>1+.0015. This selection uses D<=d*,
the original actual lenses and the radius-free angle box; it does not
use either outer norm. The cap need not hold on every interpolated fan.

Using 1-cos(Delta)<=Delta^2/2, then R<23/100 and the ONE remaining
budget Delta^2+T<=D-E, (3) gives

    M < Phi_R,D(E)
      =1+R^2-2R cos[b_R+(71/100)sqrt(E)]
                         +(23/100)(D-E), 0<=E<=D.     (4)

Strictness in (4) comes from chi, not multiplication by T or Delta.
For 0<E<=D<=d*,

    Phi'_R,D(E)
      =R(71/100)sin[b_R+(71/100)sqrt(E)]/sqrt(E)-23/100
      >=R0(71/100)sin(261/500)/(3/20)-23/100
      >3/100.                                        (5)

The final rational/trigonometric input is the frozen radius proof's
already accepted lower derivative guard. Nothing in (5) depends on
replacing D by d*. Continuity covers E=0, including D=0, where only
E=0 is admissible. Thus Phi_R,D is maximized at E=D, proving (2).
Actual tau, angular and skew terms have been paid once. No middle
HALF-arc variance, Q, BI, separate J or assigned-cell lens is charged.

## 2. A sharp connected-arc norm interface, averaged over two centers

For a covered connected full arc of length L, put a=L/2. At R<1,
0<=a<=asin(R). The independently reviewed arbitrary-center endpoint
quadratic yields

    |c|>=rho_R(a)=cos(a)-sqrt(R^2-sin^2(a)).            (6)

The statement includes the lower and upper radial branches; no radial
alignment of the given center is assumed. A zero length uses the
continuous endpoint rho_R(0)=1-R. Multiple disconnected components of
one center cannot be combined into L.

The squared lower-root function f_R(a)=rho_R(a)^2 is strictly convex
on [0,asin(R)], with a continuous endpoint extension. For w=sqrt(R^2-
sin^2(a)) in the interior,

    rho'_R=sin(a)[cos(a)/w-1]>=0,
    rho''_R=cos(a)[cos(a)/w-1]
                  +sin^2(a)(1-R^2)/w^3>0,
    f''_R=2(rho'_R)^2+2rho_R rho''_R>0.               (7)

Here rho_R>0 and cos(a)>w follow from R<1. The a=0 derivative
extends regularly; the upper endpoint is treated by continuity and
interior approximation. Thus if the TWO DIFFERENT actual outer disks
cover connected full arcs L0 and L2, then

    M>=(f_R(L0/2)+f_R(L2/2))/2
      >= f_R((L0+L2)/4).                             (8)

Summing these two connected lengths is legitimate. It does not sum
disconnected assigned components of a single center, and it does not
assert that either actual arc survives an auxiliary translation.

## 3. Explicit adaptive fee from an outer-norm or outer-arc mean

For any proposed 0<=d<=d*, let beta=beta_R(d), U=sqrt(H_R(d)), and

    a_R(d)=atan[R sin(beta)/(1-R cos(beta))].           (9)

The denominator is positive, 0<beta<.632, and R<cos(beta), so this
atan selects the positive acute branch. It is not an unspecified
principal-angle identity. Direct algebra gives

    cos(a_R)=(1-R cos(beta))/U,
    sin(a_R)=R sin(beta)/U,
    U^2+1-2U cos(a_R)=R^2,
    U-cos(a_R)=R(R-cos(beta))/U<0.                   (10)

Thus U is exactly the LOWER norm root at a_R. Also
U^2-sin^2(beta)=(cos(beta)-R)^2>0, hence sin(a_R)<R and
0<a_R<asin(R); no infeasible arc threshold has been introduced.

Equations (2), (8) and (10) prove both reusable implications:

    actual M>=H_R(d)                    ==> D>d;
    actual L0+L2>=4a_R(d)               ==> D>d.        (11)

This permits unequal outer norms and unequal arc lengths. At the closed
threshold in (11), strictness still follows from (2). A simpler strictly
larger arc input also works without relying on that equality endpoint.

More detailed unequal arc information can retain the first line of (8)
and compare its exact mean with H_R(d). Alternatively, for a prescribed
mean squared norm V in the admissible range, solve the monotone function
H_R(d) by

    d=[(acos((1+R^2-V)/(2R))-b_R)/(71/100)]^2,          (12)

provided the acos argument is in [-1,1], its chosen angle is at least
b_R, and 0<=d<=d*. Outside these conditions (12) is not asserted.
H_R is strictly increasing on [0,d*]; its continuous zero endpoint and
positive derivative on d>0 justify the inverse. Truncating to d* is
permitted only when V>=H_R(d*). There is no new global target budget here.

## 4. Two new paper-certified arc-sum corollaries

The critical a_R(d*) increases with R. To see this, its tangent
t_R=R sin(beta)/(1-R cos(beta)) has positive partial derivative in R,
and partial derivative in beta equal to
R(cos(beta)-R)/(1-R cos(beta))^2>0. Since beta_R increases with R,
the whole derivative is positive. Consequently interval arc thresholds
are certified at their RIGHT radius endpoint, in contrast to the mean
norm's decreasing endpoint. This is a monotonic proof, not a sample.

First, for all R0<=R<=.116, the exact upper beta is
124999/198000<.6314. The sine upper alternating polynomial gives

    sin(.6314)<.5903.

An explicit paper certificate is x^3=.251717687144>.2517 and
x^5<(.25172)(.399)=.10043628<.10044, so
P5(x)<.6314-.2517/6+.10044/120=.590287<.5903.
Also beta>=.6285 and the cosine upper polynomial at .6285 gives
cos(beta)<.809 (its fourth power is <.1561). Therefore

    t_R < (.116)(.5903)/(1-(.116)(.809))
        =.0684748/.906156 <.07559.

The last strict rational comparison follows from
.906156*.07559-.0684748=.00002153204>0. Since atan(t)<t,

    4a_R(d*)<4(.07559)=.30236.

We can improve this sufficient bound without evaluating atan: tan(.07545)
>.07545+(.07545)^3/3=.075593171467875>.07559. Hence

    4a_R(d*)<4(.07545)=.3018.                         (13)

Thus actual connected outer arc sum >=.3018 proves D>.0225 on the
ENTIRE R interval through .116, including unequal arcs. This replaces
the stronger old requirement that EACH arc be >.1561 (.3122 sum).
It is a sufficient rational threshold, not a sharp optimum.

Second, retain the original total .294 (e.g. both arcs>.147), but now
take R<=.11364. At its right endpoint

    beta<=.6285+(100/99)(.11364-.11322)<.62893.

For x=.62893, x^3=.248775113635957>.24874 and
x^5<(.24878)(.39556)=.0984074168<.0986, hence

    sin(beta)<P5(x)<.62893-.24874/6+.0986/120
                        =.588295<.5883.

Keep cos(beta)<.809. Put t0=.0735+(.0735)^3/3=.073632355125;
tan(.0735)>t0, proved by differentiating tan(x)-x-x^3/3 and using
tan(x)>x for x>0. The exact rational margin is

    t0-(.11364)(.5883+.809 t0)=.000008570228348355>0.   (14)

Thus t_R<t0<tan(.0735), giving a_R(d*)<.0735 and

    4a_R(d*)<.294 for every R0<=R<=.11364.              (15)

Therefore the SAME actual middle cap and only total connected outer
arc length >=.294 give the STRONG fee D>.0225 through R=.11364.
In particular both original arcs>.147 suffice there. No individual
outer floor .83 is required. The earlier critical radius bracket
(.113281,.113282) is a sharp obstruction to the particular .83 floor
interface, not a sharp obstruction to an actual-tau fee. The correlated
mean threshold H_R decreases with R and removes that unnecessary floor.

## 5. Matching parameter transport across the radius interval

The independently reviewed Pro parameter proof compares the SAME
actual fan against the central lambda0,t0 support. Its bootstrap,
curvature and apex-loss budget for central D0<=.025 are dimensionless;
all those scalar statements are independent of R. If geometry supplies
the central D0>.0225 by (11), (13) or (15), then exactly the same score
subtraction proves, with a=|lambda-lambda0|, b=|t-t0|,

    D(lambda,t)>.0225-.6a-.1725b,                     (16)

on the closed rectangle a<=.0003, b<=.0002, for EVERY admitted R.
Indeed D0>=.025 gives the gross bound >=.0227684>.0225; otherwise
the central single budget gives
Dnew>=(1-b)D0-(4a+b)sqrt(D0), increasing for D0>=.0225.
No new-parameter angular localization or branch is presumed: the
central geometry theorem is applied to the unchanged physical fan,
then only the objective is transported. The exact matched new support
C(lambda,t) is retained. Thus (16) includes fee>.0222855 on that
rectangle, or >.0223455 on a,b<=.0002. A changed global target still
must be recomputed independently; this lemma does not reuse its old
delta, profile, window, graph or allocation assumptions.

## 6. Scope and remaining mathematical obligations

The whole mean envelope is valid in the central low-deficit range
D<=.0225. A fee above that cutoff needs new localization/support/tau
guards; (12) cannot silently extend it. Larger radius ranges also need
new R-dependent payment and derivative proofs. Both connected outer
arcs must be covered by the same actual centers scored in this fan,
and actual unit anchors retain their original coverage semantics.
The middle cap remains an actual-fan premise. Global normalization,
incidence, once-only marked face allocation, complete disk cover and
exact optimum are external. The parent has independently published
the unconditional .11322 bound; none of the new conditional outputs
is itself a global lower-bound certificate.

Root mathematical execution/import/replay/scan count for this paper: 0.
Only writing, full-source reading, paper reasoning and metadata hashing
were used. New rational guards above may be checked in ONE fresh guard
after independent static clearance; the old 48 guards are not rerun.
