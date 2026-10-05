# Full-apex conjugate at lambda 7/50 and the new positive tilt

This is a scalar marked-face support theorem, not a three-face compatibility
loss, a radius exclusion or an exact optimum. It supplies a new justification
outside the old implementation's declared lambda interval [0,1/8]. No old
implementation or certificate is modified.

Let the triangle circumradius be at most R, its opposite base be L=2Rx,
its actual apex angle be gamma, and J_R(L)=(L/2)sqrt(R^2-L^2/4).
For 0<x<=1 put a=arcsin(x). The circumradius condition implies
a<=gamma<=pi-a. At fixed x and gamma the height of the apex above the base
is at most Rx cot(gamma/2), with equality for the isosceles triangle.
Consequently its normalized score is at most

    x^2 cot(gamma/2)+x sqrt(1-x^2)-4 lambda x^2+t gamma.

Its second derivative with respect to gamma is

    (x^2/2)csc^2(gamma/2)cot(gamma/2)>0.

Therefore the maximum is at one of the two endpoints gamma=a or pi-a.
This endpoint reduction is valid for every real lambda and t. For x=0
the endpoint result follows by continuity, including t*pi. Degenerate
angles are included in this closure; no acuteness assumption is made.
The exact full-apex support is

    C_lambda(t)=max_{0<=a<=pi/2} {f(a),h(a)},
    f(a)=sin(a)(1+2cos(a))-4lambda sin(a)^2+t a,
    h(a)=sin(a)-4lambda sin(a)^2+t(pi-a).

Now set lambda=7/50 and t=1071883/20000000. Both maxima can be restricted
to [0,pi/3]. For a>=pi/3,

    f'(a)=cos(a)+2cos(2a)-4lambda sin(2a)+t <= t-1/2<0.

Also h'(a)=cos(a)(1-8lambda sin(a))-t. If the bracket is nonpositive,
this is negative. Otherwise, cos(a)<=1/2 and sin(a)>=sqrt(3)/2 give

    h'(a)<=1/2-2lambda sqrt(3)-t < 3/125-t<0,

where sqrt(3)>17/10 and 2lambda=7/25. These estimates cover the whole
tail, without an assumption that h is globally concave.

Both branches are strictly concave on [0,pi/3]. Their second derivatives
are

    f''=-sin(a)-8sin(2a)-8lambda cos(2a),
    h''=-sin(a)-8lambda cos(2a).

On [0,pi/4] they are negative. On [pi/4,pi/3], the second expression
is at most -sqrt(2)/2+4lambda<0, since (14/25)^2<1/2. The first is
at most -4sqrt(3)+4lambda<0. Endpoint derivative signs are strict,
so each restricted branch has one interior maximum. The checker
isolates both roots of the derivative with outward rational Taylor
enclosures and encloses their values; the acute maximum is larger.

The checker also encloses C_lambda(0). The same acute tail and
concavity estimates apply for t=0. The entire obtuse branch in that
case is the concave quadratic s-4lambda s^2 in s=sin(a), so its maximum
is exactly 1/(16lambda)=25/56. The acute maximum strictly dominates
this value, giving a safe zero-tilt reference without assuming the
old lambda<=1/8 implementation applies.

The source-pinned arithmetic primitive is the previously reviewed
`../coupled_edges/verify_support.py`; only its interval/Taylor definitions
are imported. No old verify function, interval catalogue, optimizer or
covering classification is run. This theorem does not establish the
conditional outer-Q bound, assigned-arc bootstrap, eligible q3 incidence
or the required same-score three-face epsilon > .019931551 at R=.11321.
