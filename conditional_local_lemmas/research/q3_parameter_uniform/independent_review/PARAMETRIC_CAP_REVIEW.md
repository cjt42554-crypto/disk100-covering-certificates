# Independent static review of the moving-conjugate cap

**PASS.** I read `PARAMETRIC_CONJUGATE_CAP.md` completely and checked
its analytic argument against the accepted reference conjugate and
the already reviewed parameter-rectangle localization. No scalar
maximizer, guard or checker was executed. This review adds no
whole-cover support or global radius claim.

The rectangle theorem has already established acute dominance and
one global maximizer of G inside `(L,H)`. The accepted reference
peak is also inside that box. Thus the straight interval between
the two peaks remains in the domain of the old strict curvature
bound `-G0''>4`; no claim of concavity on the entire `[0,pi/2]`
is needed.

Changing lambda alters `G''` by `-8*dl*cos(2theta)`, of magnitude
at most `8/200=1/25`. The tilt change has zero second derivative.
Consequently `-G_new''>99/25=kappa` throughout the box. Since the
old derivative vanishes at theta0, the new derivative there is
exactly

```
v=dt-4*dl*sin(2theta0).
```

Taylor's integral formula in either direction between the peaks
therefore gives the stated upper quadratic expression. Maximizing
`v*h-(kappa/2)*h^2` gives `v^2/(2*kappa)=(25/198)*v^2`, with no
need to assume the unconstrained quadratic maximizer lies in the
box. Evaluation at theta0 gives the lower support comparison.
Thus both sides of the displayed cap hold for the exact full-apex
support, rather than only for an approximate acute maximum.

Integrating the strictly decreasing derivative between its zero
at the new peak and its value v at theta0 gives
`|theta_new-theta0|<|v|/kappa` when v is nonzero. If v is zero,
uniqueness puts the new peak at theta0, so the displayed nonstrict
bound also covers this case. The uniform constants are correct:

```
|v|<=4|dl|+|dt|<=1/40,
(25/99)*(1/40)=5/792,
(25/198)*(1/40)^2=1/12672.
```

The reference score is exactly
`C0-4*dl*x0^2+dt*theta0`. Replacing C by this upper cap in an
already accepted strict three-face fee keeps the result strict,
since `S<3C-E<=3*cap-E`. It creates no separate fee or geometric
charge. Use of the stronger E and wider radius domain still depends
on the separate stronger corollary acceptance.

The stated outward evaluation from the saved reference enclosures
is legitimate: `x0^2=1-y0^2`, and
`sin(2theta0)=2*y0*sqrt(1-y0^2)>0` because the accepted reference
cosine lies strictly between zero and one. Exact interval operations
must retain coefficient signs and squared-interval extrema; the
appendix asserts no unevaluated numerical certificate for a specific
whole-cover row. No mathematical blocker was found.
