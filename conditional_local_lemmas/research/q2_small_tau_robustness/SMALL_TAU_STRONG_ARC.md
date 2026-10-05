# A small-positive-anchor-slack extension of the strong-arc norm cut

This is a new analytic addendum awaiting independent static review.
It invokes the separate zero-tau proof's norm estimate, not a global
covering conclusion. No mathematical program, numerical guard, source
import, optimizer or historical replay was executed to write it.

Fix R=5661/50000, lambda=7/50,t=1071883/20000000. Use two actual
consecutive positive faces (p,c0,c1),(p,c1,c2), circumradii<=R,
actual unit anchors q0,q1 in their endpoint R-disk intersections,
and the same actual middle center norm cap

    |c1|<=sqrt(1-R^2).

Both actual outer centers cover connected unit-circle arcs of full
length>.147. All score and anchor-slack definitions are those of
q2_strong_arc_norm_budget/TAU_ZERO_NORM_PROOF.md. Let

    D=2C-S0-S1, T=tau0+tau1.

Claim, conditional on that proof's independently accepted dimensionless
localization and zero-tau norm estimate:

    D<=223/10000 implies T>1/12500.              (1)

In particular the actual signed-A two-face deficit in this strong
arc/norm-cap domain is strictly greater than 1/12500=.00008, because
T<=D. This smaller general actual-anchor fee does not reach the parent's
required deficit .022282728..., and does not claim a global radius
exclusion or that the upper-layer norm cap can be imposed on an already
fixed fan without rebuilding its actual geometry.

## 1. Paid ranges and exact ideal tips

Suppose D<=D0=223/10000 and T<=T0=1/12500. Full scalar support gives
D>=0. Its accepted per-face localization gives .7<theta_i<.96,
x_i=sin theta_i>.64, z_i>=0 and the single budget

    2 sum(theta_i-theta*)^2+sum z_i+sum sigma_i+T<=D.

Thus the zero-tau proof's paid ranges remain valid: .75<s<.904,
|Delta|<.15, |eta|<.18, Z<=.01115. Its exact tip-direction construction,
which depends only on the actual triangles and not on the actual
anchor location, gives

    alpha=3s/2+Z/2-eta in (.945,1.55),
    b*=(tip1-tip0)/(2R)=i sin(alpha) n*,
    (tip0+tip1)/2-c1=R cos(alpha) n*,

where n* is a unit vector represented in complex notation. All original
vertices and apex remain actual. No independent midpoint or side is
substituted. In particular .8<|b*|<=1 and cos(alpha)>.02. The first
bound follows from sin(.945)>=.945-(.945)^3/6>.8. For the second,
alpha<1.55 and cosine is decreasing on this acute interval. Direct
alternating Taylor evaluation gives the purely rational guard

    cos(1.55)>=1-(1.55)^2/2+(1.55)^4/24-(1.55)^6/720
                 +(1.55)^8/40320-(1.55)^10/3628800>.02.

This uses no additional pi enclosure or unbundled arithmetic input.

## 2. Constructing the nearby unit-tip translate

Set epsilon_i=(q_i-tip_i)/R. The exact actual-anchor lens lemma gives

    |epsilon_i|<=tau_i/x_i^2<=(625/256)tau_i.

Define e=(epsilon0+epsilon1)/2,f=(epsilon1-epsilon0)/2 and

    epsilon=(625/512)T.

Both |e| and |f| are at most epsilon; they do not receive separate
budgets. Since T<=1/12500, epsilon<1/10000.
The actual normalized half chord is b=b*+f, so it is nonzero and
|b|<1.001. Choose n perpendicular to b with its orientation consistent
with n*; the normalization inequality gives

    |n-n*|<=2|f|/|b*|<(5/2)epsilon.

The actual unit-anchor mean has exactly the form

    a=(q0+q1)/2=k*u*n, k in {+1,-1},
    u=sqrt(1-R^2|b|^2).

Put u*=sqrt(1-R^2|b*|^2). The exact radius and |b|<1.001 give
u,u*>.993, since R^2<.013 and .013*(1.001)^2<1-(.993)^2.
Rationalizing the square roots gives

    |u-u*|
      <=R^2 (|b|+|b*|) |f|/(u+u*)
      <[.013*2.001/1.986]epsilon<.014epsilon.

Translate the entire fan, with no rotation of its relative vertices,
so that its ideal tips have unit-anchor mean a*=k*u* n*. The resulting
apex p' and centers c_j' differ from the originals by one common
translation delta_p. Actual mean a equals ideal tip mean plus R e;
hence

    |delta_p|<=|u*n-u* n*|+R|e|
      <(2.5+.014+.114)epsilon
      =2.628*(625/512)T<3.21T<=.0002568.         (2)

Both new ideal tips are exactly unit points, since their half chord
is R b* and their mean is perpendicular to it with norm u*. They
remain the outward endpoint-disk lens tips of the translated actual
bases. All relative triangle side lengths, orientations, apex angles,
areas and circumradii are unchanged.

The branch k=-1 is impossible: the translated middle center would
have norm

    |c1'|=u*+R cos(alpha)
      >=sqrt(1-R^2)+R*.02
      =sqrt(1-R^2)+.0022644.

Its distance from the actual c1 is <.0002568 by (2), contradicting
the assumed actual middle norm cap. Thus k=+1. For this branch,

    |c1'|=u*-R cos(alpha)>0,
    |c1'|^2-(1-R^2)=2R cos(alpha)(R cos(alpha)-u*)<0.

The ideal unit-tip translated fan therefore obeys the same middle
norm cap. This is an auxiliary translation for a proof comparison,
not a claim that its assigned arcs or its full-cover incidence are
preserved.

## 3. Score and norm comparison with one budget

The translated unit tips have tau_i'=0. All other score terms are
unchanged, and actual signed areas increase by R^2*tau_i. Therefore

    D'=2C-S0'-S1'=D-T, 0<=D'<=D<=D0.

Nonnegativity follows from the complete scalar support applied to
the translated faces. Section 5 of the zero-tau proof, before it invokes
any outer covered arcs, gives

    min_j in {0,2}|c_j'|^2<.8294023684<(.9108)^2.

That norm estimate only needs D'<=D0, the two unit tips, the paid
shape variables and the middle norm cap; it does not require that the
translated outer arcs satisfy the original cut. Consequently (2)
implies at least one original actual outer norm obeys

    |c_j|<.9108+.0002568=.9110568<.911088875.

But both original actual outer covered arcs>.147 independently force
|c_j|>.911088875. This contradiction proves (1). Taking D<=T0 would
also give T<=D<=T0 and D<=D0, so (1) proves D>T0 for the entire
stated actual-anchor domain, without a zero-tau assumption.

The result is deliberately weaker than the desired .022282728... fee.
It establishes a concrete nonzero-tau neighborhood and a small uniform
actual signed-A fee. It does not spend BC assigned-cell lens loss,
global Q, or any middle HALF-arc variance. A sharper general tau
payment remains necessary for the parent's larger target.
