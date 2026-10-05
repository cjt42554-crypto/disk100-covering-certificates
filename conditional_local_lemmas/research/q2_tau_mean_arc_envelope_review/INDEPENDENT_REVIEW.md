# Independent static review of the actual-tau mean and two-arc envelope

Reviewer: `/root/q2_uniform_review`, independent of the root author.
Verdict: **STATIC PASS** for the frozen root paper PROOF.md, SHA256
`570696466a5e90d0d06e1738e41a4a94571ec7a8c45a66559fea27f2bd3cd663`.
Mathematical executions, mathematical imports, checker calls, numerical
replays, scans and optimizers: **0**. The author paper and its inputs
were not edited. Metadata hashes and byte-count checks are the only
computational operations in this review.

I read the whole new paper, its source pins, ledger and manifest, and
the cited frozen radius main, relaxed-middle corollary, geometry source
and reviews, and the Pro parameter source and review. The new strict
mean envelope, convexity, root selection, two continuous-radius arc
thresholds and radius-uniform objective transport were independently
derived below. They are not accepted merely because an earlier local
fee or its audit passed.

The central geometric domain is the full closed radius interval
[5661/50000,29/250], fixed lambda=7/50 and t=1071883/20000000,
two same actual consecutive positive nondegenerate faces with
circumradii at most the current R, and their actual covered unit base
anchors, permitting signed A of either sign and arbitrary actual tau.
The SAME actual middle squared norm is at most 1+3/2000. There is
NO outer norm or arc hypothesis for the mean-envelope conclusion.
Under D<=9/400 it establishes

    (|c0|^2+|c2|^2)/2 < H_R(D),
    H_R(d)=1+R^2-2R cos[b_R+.71 sqrt(d)],
    b_R=.522+(100/99)(R-.11322).

For actual connected outer full arcs L0,L2, it then proves the
closed-threshold implication L0+L2>=4a_R(d) => D>d, for 0<=d<=.0225.
In particular, a sum >=.3018 implies D>.0225 throughout R<=.116;
a sum >=.294 implies D>.0225 throughout R<=.11364. These are
conditional local conclusions, not global covering certificates.

## Removing the individual outer-floor hypothesis is justified

I traced every use of the old outer squared-norm floor in the frozen
radius proof. It is used only in the final contradiction with the
mean upper envelope. The normalized support/localization, one deficit
budget, signed directional support, whole-slab radius transfer,
complete tau remainder and changing angle baseline require no outer
norm premise. The relaxed-middle proof also selects the actual positive
translation branch using only the same middle cap, angle localization
and actual endpoint lenses. Its negative-branch lower squared norm
exceeds 1+.00158508 independently of either outer norm.

Consequently all pre-contradiction mean estimates genuinely remain
valid with both outer floors removed. This is not inference of a new
geometric domain merely from the old theorem's final fee statement.
The cited source formulas explicitly supply the intermediate estimates
needed for the new theorem.

## Actual D, strictness and the zero-deficit endpoint

Complete scalar support gives each actual face score at most C and
thus D>=0. The low-D budget is exactly

    E+Delta^2+2Z+T<=D,
    E=4(s-theta*)^2+c_eta eta^2.

The peak starting angle is strictly below b_R. Integrating the shape
derivative bounds and then applying nonstrict Cauchy therefore gives
chi<b_R+.71 sqrt(E), including E=0. This strictness is preserved
even if all increments vanish; it does not come from multiplying
a strict coefficient inequality by a zero budget.

All angles are positive and below .632<pi/2. From the exact ideal
mean and the single tau cost,

    M<=1+R^2-2R cos(chi)cos(Delta)+.23T
      <=1+R^2-2R cos(chi)+R Delta^2+.23T.

Here 2R cos(chi)(1-cos(Delta))<=R Delta^2 since cos(chi)<=1.
Replacing chi by its STRICT upper bound strictly raises the displayed
norm bound. Since R<.23 and Delta^2+T<=D-E from the SAME budget,
one obtains M<Phi_R,D(E) exactly as the author states. There is no
additional tau, angular, skew, Q, BI, J or variance debit.

For 0<E<=D<=.0225 the derivative lower bound is greater than .03
on the whole radius interval. Its frozen scalar guard was established
for the larger cutoff; decreasing the actual D changes only the
constant term in Phi and not this derivative. Continuity handles E=0.
The maximum is E=D, giving M<H_R(D). When D=0 only E=0 exists,
and the same strict peak-baseline argument supplies M<H_R(0).
Thus the actual-D conclusion does not require a positive-deficit
assumption or replacing D by the cutoff in its right-hand side.

H_R is continuous and strictly increasing on [0,.0225]. Its derivative
for d>0 is positive, and its right-hand square-root endpoint preserves
strict monotonicity. Therefore M>=H_R(d) contradicts D<=d, while
the branch D>.0225 already implies D>d. This validates the norm
implication for every proposed d in the specified closed interval.

## Convex squared radial root and two actual connected arcs

At R<1, a covered connected unit-circle arc is a minor arc, with
half-length in [0,asin(R)]. Averaging its endpoint distances gives
the necessary quadratic norm interval without radial alignment of
the actual center. Both positive roots are retained. For
rho=cos(a)-sqrt(R^2-sin^2(a)), w=sqrt(R^2-sin^2(a)) in the interior,
I independently differentiated and found

    rho'=sin(a)[cos(a)/w-1],
    rho''=cos(a)[cos(a)/w-1]+sin^2(a)(1-R^2)/w^3.

Because cos(a)>w>0 and rho>0, rho''>0 and
(rho^2)''=2(rho')^2+2rho rho''>0. At a=0 the derivative is regular
and rho''=1/R-1>0. At asin(R) the function is continuous; taking
interior approximations extends its convexity and Jensen inequality
to that endpoint. Thus the claimed strict convexity on the closed
interval is valid, although only nonstrict Jensen is needed.

For each ORIGINAL actual outer disk its own connected covered full
arc provides |c_i|^2>=f_R(L_i/2). Averaging and Jensen give
M>=f_R((L0+L2)/4). This combines two different actual disks and
does not merge disconnected components belonging to one center.
The two arcs need not be equal or disjoint: their sum here is a
scalar input to a norm inequality, not an assertion of disjoint area
coverage. They refer to the original physical fan; no actual arc
is presumed preserved along the auxiliary translated path.

## Explicit critical root and closed-threshold implication

Put beta=b_R+.71 sqrt(d), U=sqrt(H_R(d)), and
a=atan[R sin(beta)/(1-R cos(beta))]. The denominator is positive,
the ratio is positive and finite, and R<cos(beta) on the complete
domain. Hence atan selects its acute branch. Direct normalization
of the vector (1-R cos(beta),R sin(beta)) gives the author's two
sine/cosine expressions and the equation

    U^2+1-2U cos(a)=R^2.

Moreover U-cos(a)=R(R-cos(beta))/U<0, so U is specifically the
LOWER radial root. The identity
U^2-sin^2(beta)=(cos(beta)-R)^2>0 shows sin(a)<R; thus
0<a<asin(R). This confirms both the root's position and feasibility.
It also proves f_R(a)=H_R(d) with no unresolved sign after squaring.

If L0+L2>=4a, their averaged half-angle lies within the feasible
convex domain and is at least a. The nondecreasing f_R and Jensen
then give M>=H_R(d), contradicting D<=d by the STRICT upper mean
bound. Equality at the arc threshold is allowed. The explicit inverse
formula for d is valid under its stated acos-branch and range tests:
the nonnegative difference from b_R chooses the positive square-root
branch. Truncation at .0225 is justified only at or above H_R(.0225),
as the source states. No claim is made beyond the low-D localization.

## Exact uniform sum .3018 through R=.116

The tangent of the critical angle has partial derivative
sin(beta)/(1-R cos(beta))^2 in R, and partial derivative
R(cos(beta)-R)/(1-R cos(beta))^2 in beta. Both are positive. Since
beta_R(.0225) increases with R, a_R(.0225) increases on the whole
interval. Right-endpoint thresholds therefore control all smaller R.

I checked each new polynomial and rational comparison on paper.
At the right endpoint beta=124999/198000<.6314. On this acute
interval sine increases; the alternating upper P5 gives
sin(beta)<.5903, using

    (.6314)^3=.251717687144>.2517,
    (.6314)^5<.25172*.399=.10043628<.10044,
    .6314-.2517/6+.10044/120=.590287<.5903.

Beta>=.6285 and the upper P4 at .6285 give cos(beta)<.809:
its square is .39501225 and its fourth power is below .1561,
so P4<.808998041667<.809. All denominator comparisons therefore
retain positive signs. The critical tangent satisfies

    t_R<.0684748/.906156<.07559,

with the exact positive margin
.906156*.07559-.0684748=.00002153204. Also

    .07545+(.07545)^3/3=.075593171467875>.07559.

For x>0 below pi/2, differentiating tan(x)-x-x^3/3 gives
tan(x)^2-x^2>0 because tan(x)>x. Therefore tan(.07545) exceeds
the displayed cubic, so a_R(.0225)<.07545. Multiplication by four
gives the strict upper threshold .3018, uniform through .116.
Thus the CLOSED sufficient input L0+L2>=.3018 establishes D>.0225.

## Exact uniform sum .294 through R=.11364

At the second right endpoint beta<.62893. Its square is
.3955529449<.39556 and its cube is
.248775113635957, between .24874 and .24878. The upper sine
polynomial then gives

    x^5<.24878*.39556=.0984074168<.0986,
    sin(beta)<.62893-.24874/6+.0986/120=.588295<.5883.

The same acute cosine bound <.809 applies. The cubic tangent lower
bound at .0735 is t0=.073632355125. I independently checked the
author's long multiplication:

    t0-.11364(.5883+.809t0)=.000008570228348355>0.

Equivalently .11364*.5883/(1-.11364*.809)<t0, with positive
denominator. Hence the actual critical tangent is below t0 and
below tan(.0735); its critical angle is below .0735. Multiplying
by four proves 4a_R(.0225)<.294 throughout the closed radius band
through .11364. Thus total actual connected outer length >=.294,
including the sufficient example of both lengths>.147, gives D>.0225.
The old sharp obstruction to the single .83 floor does not obstruct
this different correlated mean-norm interface.

## Parameter transport preserves the physical fan uniformly in R

I independently checked the composition with the accepted Pro
parameter source. Its full-apex score/support, gross comparison,
central .025 bootstrap budget and curvature statements are all
dimensionless and hold at every admitted current R. The new central
mean/arc theorem supplies D0>.0225 on the same physical fan whenever
the applicable d=.0225 threshold is established. The parameter change
alters the objective only; the actual centers, original outer arcs,
anchors, circumradius and middle cap do not change.

The broad gross loss 4.352a+4.63b is at most .0022316, so D0>=.025
already gives Dnew>=.0227684. On the other branch the scalar budget
gives Dnew>=(1-b)D0-(4a+b)sqrt(D0), increasing for D0>=.0225.
Since sqrt(.0225)=.15 exactly, the resulting fee is
.0225-.6a-.1725b. All constants are independent of R. Thus the
full closed rectangle a<=.0003,b<=.0002 yields Dnew>.0222855
at every admitted R; the balanced a,b<=.0002 yields >.0223455.
The matched NEW complete support is retained. No unproved new-parameter
angle localization or translation branch is inferred, and no original
arc is transferred to a different fan.

## Scope, source pins and remaining action

All new mathematical implications in this frozen paper are sound.
This review clears the paper's new rational inequalities for optional
verification in one fresh guard; no guard implementation has yet been
reviewed or executed here and the old 48 guards remain unrerun.
Source and author-manifest identities are recorded in SOURCE_PINS.json.

The result remains conditional on actual-fan admission, actual anchor
coverage, original connected arcs when used, and the same middle cap.
It supplies no global incidence, marked-face allocation, complete-cover
classification or exact optimum. A changed global parameter target must
be recomputed independently. Raising d above .0225, changing anchor
semantics or enlarging the radius interval requires new localization
and payment proofs. No unresolved mathematical defect was found in
the frozen source.
