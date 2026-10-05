# Blind analytic review of a new actual quadratic signed-anchor support

**Verdict: PASS for the proposed proof mechanism and constants, with the
global-concavity correction stated below.** The main new source was not
yet frozen when this review was prepared. This is an independent analytic
check of the supplied constants, not a full-source replay or certificate.
No old checker, imported interval module, optimizer or parameter search
was executed. A new standalone Fraction/Taylor calculation checked the
fixed scalar guards and rational arithmetic. Old sources were not changed.
A separate read-only child review independently confirmed the full-apex
convex cut and paid localization.

The proposed parameters and objective are

    R=11321/100000, lambda=7/50, t=1071883/20000000,
    T=sum_i[area_i/R^2+A_i/R^2-(14/25)*x_i^2+t*gamma_i],
    D=3C-T, D0=1/40.

Each `A_i` is the actual opposite-apex-side signed anchor, which may
be negative. Each anchor is a unit point covered by both incident
radius-R disks. There is one actual consistently oriented three-face
fan and one common apex. No min-arc, reflection, individual center
norm or canonical nearest-owner hypothesis is needed for this expanded
common-unit-point domain.

The complete one-face constant, unique acute peak and obtuse endpoint
bound are the accepted frozen `quadratic_lambda14_conjugate` inputs.
Their certificate SHA-256 is
`f325d42bfe651599881b7dc652567729a1586389ee385c38060840bc41696bc2`.
The exact area/skew and signed outer-tip formulas are geometric
identities from the already reviewed actual-anchor mechanism; their
use below does not transfer the old linear-objective fee.

## 1. Full gamma cut and one paid budget

Use `theta=asin x`, `z=gamma-theta>=0`, boundary-angle skew `d`,
and the two existing nonnegative shape and anchor deficits

    sigma=2*x^2*sin(d)^2/sin(gamma),
    tau=x*sqrt(1-x^2)-A/R^2.

The exact actual score is `F-sigma-tau`, where

    F=x^2*cot(gamma/2)+x*sqrt(1-x^2)-(14/25)*x^2+t*gamma,
    G(theta)=sin(theta)+sin(2theta)-(14/25)*sin(theta)^2+t*theta.

The complete support gives each `C-score_i>=0`. Thus `D<=D0`
implies each score at least `C-D0`. Convexity of `F` on the full
feasible interval `[theta,pi-theta]` and `C-Co>.91>D0` force
the acute endpoint `G(theta)>=C-D0`, including for initially
obtuse actual faces.

**Do not copy the old global strict-concavity sentence.** For this
new quadratic potential, `G''(pi/2)=.12>0`. The frozen conjugate
instead proves strict concavity on `[0,pi/3]` and strict decrease
on `[pi/3,pi/2]`. Those statements give the unique acute peak and
permit localization using two endpoint guards.

Let `L=.7`, `H=.96`. The peak is in `(.82,.84)`. On `[L,H]`,

    -G''=sin(theta)+4*sin(2theta)+(28/25)*cos(2theta)>4.

One fully rational endpoint mechanism is
`sin(.7)>16/25`, `sin(1.4),sin(1.92)>469/500` and
`cos(1.92)>-7/20`. Sine is concave on `[1.4,1.92]`, so its
minimum is at an endpoint; the other two terms use monotonicity.
The resulting lower expression is exactly
`16/25+4*(469/500)-(28/25)*(7/20)=4`, with strict guards.
These endpoint guards follow from alternating Taylor sums at
the stated rational arguments.

Since both endpoint distances from the peak exceed `.12`, Taylor's
theorem along the corresponding segment gives
`C-G(L),C-G(H)>2*(.12)^2=.0288>D0`. The frozen branch monotonicity
then excludes the entire outer theta domain and yields `L<theta<H`.
A direct independent Taylor enclosure gives still larger endpoint
margins after subtracting `D0`, exceeding `.01372` and `.01401`.

Set `zcut=.1` and
`r(theta)=sin(theta)/sin((theta+zcut)/2)`. The sign of its derivative
is that of

    N=cos(theta)*sin(u)-(1/2)*sin(theta)*cos(u),
    u=(theta+zcut)/2,
    N'=-(3/4)*sin(theta)*sin(u)<0.

Thus `r` has no interior minimum. At both endpoints, the simple
bounds `sin(v)>v-v^3/6` and `sin(u)<u` give `r>3/2`. Consequently
on the whole cut strip

    -F_gamma>(9/8)-t>1,

with last margin `1428117/20000000`. Also `2H+.1=2.02<pi`,
so the cut lies inside the feasible gamma interval. At the cut,
`F<G-.1<=C-.1<C-D0`. Between that cut and `pi-theta`, convexity
puts the maximum at the cut or the obtuse endpoint, both strictly
below `C-D0`. This excludes all remaining actual gamma values.
Integration and the curvature bound therefore give the single budget

    2*sum(theta_i-theta_*)^2+sum z_i+sum sigma_i+sum tau_i
                                                        <=D<=1/40.

This sequence handles the full gamma domain before localizing it;
it does not assume acute actual apices in the initial support.

## 2. Exact tips, signed anchors and skew

The fixed rational angle window guards are

    S0=.64<sin(theta)<.82=S1,
    C1=.573<cos(theta)<.765=C0,
    Qmin=2.14<Q_i=1+2*cos(theta_i)<2.54=Qmax.

The exact outer tip, in the actual bisector phase, is
`U*cos(d)+i*V*sin(d)`, with

    U=sin(theta)*cot(gamma/2)+cos(theta),
    V=sin(theta)*tan(gamma/2)-cos(theta).

Its radial change from the acute tip satisfies
`|U-Q_i|<=cot(L/2)*z<2.75*z`. On the localized acute gamma strip,
`1-q<=cot(theta)*z<=(C0/S0)*z`, where
`q=sin(theta)/sin(gamma)<=1`. The exact decomposition of `V`
then gives

    V>=1-2*C0-(C0/S0)*(1-C1)*D0>-.55,
    V<=1-2*C1+D0<.55.

Also `sin(d)^2<=cs*sigma` with `cs=625/512=1/(2*S0^2)`.
Because `|d|<pi/2`, `1-cos(d)<=sin(d)^2`. Hence the squared
complex skew displacement is at most

    cs*[.55^2+Qmax^2*cs*D0]*sigma < .79^2*sigma.

The strict rational margin in the coefficient is positive (greater
than `.01449`). This combines orthogonal real and imaginary tip
errors before spending the original `sigma` payment.

The signed lens inequality gives the distance from the actual anchor
to the unique outer tip at most `tau/x^2`, in units of R. Thus its
coefficient is `ct=625/256`. This holds even for negative signed
anchors: their nonnegative `tau` is larger. Taking an absolute
anchor area or selecting another tip according to its sign would
change the objective and invalidate this payment argument.

## 3. Full-box Jacobian and the changed imaginary sign

Fix the actual middle bisector phase. For independent acute angles
in `[L,H]^3`, set

    a=(theta0+theta1)/2, b=(theta1+theta2)/2,
    Bacute=(Q0*exp(-ia)+Q2*exp(ib))/2-Q1,
    bacute=(Q2*exp(ib)-Q0*exp(-ia))/2.

The endpoint derivative imaginary component lies between
`S0^2-Qmax*C0/4` and `S1^2-Qmin*C1/4`. Its lower possible
value is negative here. Accordingly the absolute bound must be

    u=max(S1^2-Qmin*C1/4,Qmax*C0/4-S0^2)=.365845.

It is invalid to reuse the old positivity assumption. With this
correction, define the proposed component bounds

    tJ=S1*C0+Qmax*S1/4=1.148,
    v=2*S1-Qmin*S0/2=.9552,
    w=(Qmax*C0-Qmin*C1)/4=.17922,
    vb=(Qmax*S1-Qmin*S0)/4=.1783,
    wb=Qmax*C0/2=.97155.

The middle real column of B remains positive because its lower
bound is `2*S0-Qmax*S1/2=.2386>0`; its upper magnitude is `v`.
The middle imaginary difference is bounded by `w`. The corresponding
b middle difference and sum are bounded by `vb` and `wb`. All
endpoint real magnitudes are bounded by `tJ`, and endpoint imaginary
magnitudes by `u`.

The resulting squared Frobenius bounds are exactly

    JB2=76960399529/20000000000,
    Jb2=77583868411/20000000000.

They bound the operator norms on the entire closed angle box.
The segment from the ideal peak triple to the actual triple stays
in that box, so integration supplies the stated global acute-map
errors, not merely a local linearization at the peak.

## 4. One shared budget for every correction

For the actual normalized combinations

    B=((w0+w2)/2-w1)/R, b=(w2-w0)/(2R),

the common apex cancels exactly. Their ideal values are
`B*=(-A,0)` and `b*=(0,Hhalf)`. The combined angular/skew bounds
by Cauchy have squared coefficients

    alphaB2=JB2/2+(3/2)*.79^2
            =114406399529/40000000000,
    alphab2=Jb2/2+(1/2)*.79^2
            =90065868411/40000000000.

The factors `3/2` and `1/2` are the squared norms of the respective
three tip-combination weights. Actual bisector displacements are
`(z0+z1)/2` and `(z1+z2)/2`. Safe linear coefficients are
`LB=2.75+Qmax/2=4.02` and `Lb=(2.75+Qmax)/2=2.645`.
They dominate the anchor coefficients `ct` and `ct/2`.

Let `E=2*sum e_i^2+sum sigma_i` and
`W=sum z_i+sum tau_i`. The same budget gives `E+W<=D`; hence

    |B-B*|<=sqrt(alphaB2)*sqrt(D-W)+LB*W,
    |b-b*|<=sqrt(alphab2)*sqrt(D-W)+Lb*W.

Both exact guards `4*L^2*D0<alpha^2` hold, so the tangent
inequality for the square root makes the maximum occur at `W=0`.
This is an analytic one-budget bound, not independent worst-case
spending on each kind of error. The strict fit margins are

    .27^2-alphaB2/40=2233600471/1600000000000>0,
    .24^2-alphab2/40=2094131589/1600000000000>0.

Thus both error bounds `.27` and `.24` hold simultaneously.

## 5. Ideal sizes and strict unit-circle contradiction

The frozen peak gives `.676<y=cos(theta_*)<.677`. The rational
guarded ideal sizes are

    .76<A=(1+2y)*(1-y)<.763,
    1.73<Hhalf=(1+2y)*sqrt(1-y^2)<1.74.

The first follows from monotonicity of `1+y-2y^2` on this interval;
the second can be checked by squaring positive lower and upper
product bounds. All four strict margins are positive.

For any three actual unit points, the exact necessary condition is

    2*|B_perp_to_b|<=R*(|B|^2+|b|^2).

For completeness, let their physical endpoint mean be `a` and
physical half-chord be `b0`. Then `a dot b0=0` and
`|a|^2+|b0|^2=1`. Along a unit normal to `b0`, writing the
coordinates of `a` and the middle unit point as `s,t` gives
`|s-t|<=1-s*t`, since
`(1-s*t)^2-(s-t)^2=(1-s^2)*(1-t^2)>=0`.
This proves the condition after scaling. The near-ideal half-chord
bound ensures `b` is nonzero, so its perpendicular is defined.

The horizontal fraction of the actual `b` is at most
`.24/(1.73-.24)=24/149`. Therefore the ideal horizontal axis
projects onto its perpendicular by more than `.98`.
The necessary left and right sides have the strict bounds

    lhs>2*(.76*.98-.27)=.9496,
    rhs<R*[1.033^2+1.98^2].

Their exact positive rational difference is
`38496637031/100000000000`. This contradicts the unit-circle
condition and proves the proposed local support `T<3C-1/40`,
provided the final main document states precisely these hypotheses
and proves/includes the listed guards.

## Scope and next audit

The mechanism works for the new quadratic objective in the expanded
actual-common-unit-point domain, without min-arc or reflection input.
It pays each marked area and signed anchor once and creates no new
copy of a BC lens loss. No old `.05` or other linear-objective fee
has been transferred. Applying the result to a retained-J score
would still require the separate same-objective anchor-slack bridge
and a single-use allocation of the existing middle-cell deficits.

A static full-source review should follow once the new main proof,
checker and certificate are frozen. This pre-freeze review does not
claim that replay or a global covering exclusion.
