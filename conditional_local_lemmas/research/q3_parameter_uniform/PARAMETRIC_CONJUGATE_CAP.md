# An explicit quadratic cap for the moving full-apex conjugate

Use the parameter rectangle and exact full-apex support C(lambda,t) from
PROOF.md. Write dl=lambda-lambda0, dt=t-t0, where lambda0=7/50 and
t0=1071883/20000000. Let theta0 be the accepted acute maximizer at
the reference parameters, x0=sin(theta0), and C0=G0(theta0).

The following explicit cap and peak displacement hold throughout the
closed rectangle:

    C0-4dl*x0^2+dt*theta0 <= C(lambda,t)
      <= C0-4dl*x0^2+dt*theta0
           +(25/198)[dt-4dl*sin(2theta0)]^2,

    |theta(lambda,t)-theta0|
      <=(25/99)|dt-4dl*sin(2theta0)| <=5/792.

Here theta(lambda,t) is the unique global acute support maximizer proved
in Section 2 of PROOF.md. These statements do not assume global G
concavity on [0,pi/2] and do not execute any new scalar maximization.

Both reference and moving peaks lie in [L,H]. The accepted old curvature
bound gives -G0''>4 there, while changing lambda perturbs G'' by at most
8|dl|<=1/25; changing t leaves it unchanged. Thus, on the entire segment
between these peaks,

    -G(lambda,t;theta)'' > 99/25 = kappa.

At theta0 the new derivative is exactly v=dt-4dl*sin(2theta0), since
G0'(theta0)=0. Taylor's integral formula gives

    G_new(theta)<=G_new(theta0)+v(theta-theta0)
                                  -(kappa/2)(theta-theta0)^2.

The right side has maximum G_new(theta0)+v^2/(2kappa), proving the
upper cap; evaluation at theta0 proves the lower bound. Integrating the
strict derivative decrease between the two peaks proves the displacement
bound (with a strict bound when v is nonzero). The nonstrict statement
also covers v=0, when the peaks coincide. Finally |v|<=4|dl|+|dt|<=1/40
and (25/99)(1/40)=5/792.

This cap is useful with the actual-A joint theorem because it retains
the same reference-peak first-order terms. Once the stronger rectangle
fee E*=25473/800000 is accepted, for its stated radius domain one has

    S(lambda,t) < 3[C0-4dl*x0^2+dt*theta0
          +(25/198)(dt-4dl*sin(2theta0))^2] - E*.

The cap correction is at most 1/12672 per face in the rectangle. For a
particular rational lambda,t, a whole-cover verifier can outward-bound
the displayed expression using the already pinned reference peak,
cosine and C0 enclosures; x0^2=1-cos(theta0)^2 and
sin(2theta0)=2cos(theta0)sqrt(1-cos(theta0)^2)>0. It need not isolate
a new conjugate peak for each parameter pair. The exact support still
has the full-domain endpoint definition, and this does not certify any
entire-cover row, boundary support or global optimum by itself.
