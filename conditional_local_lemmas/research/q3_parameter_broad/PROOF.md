# Broader actual-A support: a closed coefficient rectangle

Use the actual common-apex positive-area three-face fan, covered unit
anchors of either signed side, score S(lambda,t) and exact full-apex
support C(lambda,t) defined in ../q3_parameter_uniform/PROOF.md. This is
a local geometric theorem; it is not an independent-J fee or a global
covering bound. The accepted stronger fee19/500 and its independently
reviewed radius extension are premises.

For the closed rectangle

    lambda in [63/500,77/500]=[.126,.154],
    t in [791883/20000000,1351883/20000000]=t0+[-.014,.014],
    0<R<=794400/5378209,

one has S(lambda,t)<3C(lambda,t)-41/2000. The fee is .0205; the radius
endpoint is included. No center norm, minimum arc, reflection, ownership
or common-apex location restriction is imposed.

Put eps=7/500, dl=lambda-7/50,dt=t-t0, L=.7,H=.96. The accepted old
acute peak theta0 is in(.8,.85); old curvature -G0''>4 on[L,H], and
saved old endpoint gaps are greater than .0387 and .039 respectively.
The new derivative perturbation is bounded by4|dl|+|dt|<=5eps=.07.
Thus G_new'(L)>.4-.07>0, G_new'(H)<-.44+.07<0.

The full acute branch is strictly concave on[0,pi/3]: on[0,pi/4] every
term of -G_new''=sin(theta)+4sin(2theta)+8lambda cos(2theta) is
nonnegative and their sum is positive; on[pi/4,pi/3] their sum exceeds
.7+4(.8)-4(.154)>0. On[pi/3,pi/2], G_new'<=-.5+tmax<0.
Consequently its unique global maximum lies in(L,H). The reference
evaluation gives C_new>1.47-4eps-eps*H=1.40056>7/5. The entire obtuse
branch satisfies Bo_new<=1+tmax*pi<1+tmax*22/7, more than41/2000 below
7/5. Thus C_new is the full-domain acute maximum, without asserting
concavity of G on the entire[0,pi/2].

Using sin^2 being1-Lipschitz and comparing at theta0, each new endpoint
gap is strictly larger than .0387-5eps(H-L)=.0205. Hence a hypothetical
total new deficit Dnew<=41/2000 localizes every theta_i to(L,H): each
individual deficit is nonnegative by the full convex-apex endpoint
support, and both branches outside that window are below its threshold.
On theta<=gamma<=theta+.1, the accepted ratio guard>81/50 gives
-F_new,gamma>(1/2)(81/50)^2-tmax>1. The cut is feasible since2H+.1<pi;
its loss is at least.1. Convexity on the rest of the apex interval and
the obtuse gap therefore force z_i=gamma_i-theta_i in[0,.1).

On this actual domain the old scalar support inequalities imply the
same single reference budget, without assuming Dold small:

    2 sum(theta_i-theta0)^2+sum z_i+sum sigma_i+sum tau_i <= Dold,
    Dold=3C0-S(7/50,t0)>19/500.

The last inequality is the separately accepted stronger actual-A fit at
the same radius, valid through the displayed endpoint. The centered
parameter identity and Cauchy now give

    Dnew >= f(Dold),
    f(D)=(493/500)D-(7/100)sqrt(3D/2).

For D>=19/500 this f is strictly increasing: sqrt(57/1000)>.23 makes
every secant slope exceed493/500-(21/200)/(46/100)>0. The strict square
bound sqrt(57/1000)<191/800 gives

    f(19/500)>(493/500)(19/500)-(7/100)(191/800)
             =41511/2000000>41000/2000000=41/2000.

This contradicts the hypothetical new deficit, including equality and
all coefficient/radius endpoints. A zero-area limiting face has score
at most1/2+tmax*pi<1, whereas C_new>7/5, so it cannot enter such a
high-score triple. All divisions follow localization. This uses neither
assigned-cell lens loss nor terminal BI/variance payments.

The same reference-peak Taylor proof additionally gives

    C_new<=C0-4dl sin(theta0)^2+dt theta0
                +(125/972)(dt-4dl sin(2theta0))^2,
    |theta_new-theta0|<=(125/486)|dt-4dl sin(2theta0)|<=35/1944,

because restricted curvature is -G_new''>4-8eps=486/125 on[L,H].
The formula can be outward evaluated from the old pinned peak enclosures;
no new peak was isolated. Whole-cover use still requires the matching
dual coefficients, actual ordinary fan incidence, disjoint charged sets
and every remaining global support/reduction. Scalar arithmetic alone
does not replace the accepted premise proofs or their independent audits.
