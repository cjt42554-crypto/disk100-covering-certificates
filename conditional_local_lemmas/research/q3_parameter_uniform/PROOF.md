# A uniform rectangle of actual-anchor three-face supports

This is a local theorem. It is not an independent-J fee, a whole-cover
radius exclusion, or a proof of the exact 100-disk optimum.

Put lambda0=7/50 and t0=1071883/20000000. For

    27/200 <= lambda <= 29/200,
    t0-1/200 <= t <= t0+1/200,
    0 < R <= 231/2000 = .1155,

consider three positive, consistently oriented, consecutive triangles
(p,c_i,c_(i+1)) with a common actual apex p and circumradii at most R.
Each base has an actual unit anchor q_i covered by both endpoint R-disks.
Set x_i=base_i/(2R), gamma_i=the positive apex angle, and
A_i=-det(c_(i+1)-c_i,q_i-c_i)/2. Signed A_i of either sign is allowed.
Define

    S(lambda,t)=sum[area_i/R^2 + A_i/R^2 -4lambda*x_i^2+t*gamma_i].

Let C(lambda,t) be the exact full-apex independent support, equivalently
the maximum on 0<=theta<=pi/2 of the two endpoint functions

    G(lambda,t;theta)=sin theta+sin(2theta)-4lambda sin(theta)^2+t theta,
    Bo(lambda,t;theta)=sin theta-4lambda sin(theta)^2+t(pi-theta).

Then

    S(lambda,t) < 3C(lambda,t)-16027/800000,

and in particular S(lambda,t)<3C(lambda,t)-1/50. The excess over 1/50
is 27/800000. No norm, minimum arc, nearest-owner, reflected-center,
canonical-switch, or restriction on p is imposed.

## 1. Accepted inputs and the same physical domain

The accepted fixed-parameter theorem at (lambda0,t0) gives

    D0_actual=3C0-S(lambda0,t0)>1/40.

Its new parameterized-radius addendum permits all
0<R<=949600/4987489, which contains the present interval.
It derives this from the old proof's strict unit-anchor fit and is not
a replay of any frozen checker. All numerical guards below are read as
exact enclosures from the pinned accepted fixed-parameter certificate.

For a nondegenerate triangle, write theta=asin x, z=gamma-theta and let
sigma>=0 be its normalized skew area loss. The two covered endpoints
give tau=x cos(theta)-A/R^2>=0, even for negative signed anchors. The
actual score is exactly

    F(lambda,t;theta,gamma)-sigma-tau,
    F=sin(theta)^2 cot(gamma/2)+sin theta cos theta
        -4lambda sin(theta)^2+t gamma,
    theta<=gamma<=pi-theta.

F is convex in gamma; thus its maximum is one of G or Bo. Endpoint
limits also give this support. Therefore every individual actual score
is at most C(lambda,t), and all three individual new deficits are
nonnegative. This equivalence requires no assertion of global G concavity.

## 2. Uniform location of the new acute support

Fix L=7/10, H=24/25. The accepted peak theta0 of G0 is in (.8,.85),
inside (L,H); C0>1.47. Its saved endpoint gaps give

    C0-G0(L)>1/40+.0137,
    C0-G0(H)>1/40+.0140.

On [L,H] the accepted bound is -G0''>4. The derivative perturbation
has magnitude at most 4|lambda-lambda0|+|t-t0|<=1/40. Integrating
from theta0 shows G'(L)>.4-.025>0 and G'(H)<-.44+.025<0.

For every lambda in the present rectangle, G is strictly concave on
[0,pi/3]. Indeed -G''=sin theta+4sin(2theta)+8lambda cos(2theta).
On [0,pi/4] it is positive. On [pi/4,pi/3] it exceeds
7/10+4(4/5)-4(29/200)>0, using sin theta>7/10 and
sin(2theta)>=sqrt(3)/2>4/5. On [pi/3,pi/2],

    G'<=1/2-1+t<0.

Thus G has one global maximum in (L,H), increases up to its peak and
decreases thereafter, and every theta outside (L,H) has
G(theta)<=max(G(L),G(H)).

At theta0, G(lambda,t;theta0)>=C0-1/50-H/200>7/5. For any theta,
Bo<=1+t*pi<1+(t0+1/200)(22/7)<7/5-1/50. Hence Bo cannot be the
full support maximum, and its entire branch is more than 1/50 below C.

The endpoint gaps of the new G are compared to its value at theta0,
not to a numerical approximation of its new maximizer. Since sin^2
is 1-Lipschitz and |theta0-L|,|theta0-H|<H-L=13/50,

    C-G(L or H)
      >=G(theta0)-G(L or H)
      >[old endpoint gap]-(1/40)(13/50)>1/50.

This proves a uniform full-domain localization without isolating a new
maximum for each parameter value.

## 3. A high new score localizes the actual apex, then uses the old budget

Suppose for contradiction D_new=3C(lambda,t)-S(lambda,t)<=1/50.
Each individual new score is at least C-1/50. Section 2 and convex
endpoint reduction force L<theta_i<H for every face.

For theta<=gamma<=theta+1/10, the accepted endpoint ratio guard gives

    sin theta / sin((theta+1/10)/2)>81/50,
    -F_gamma>(1/2)(81/50)^2-(t0+1/200)>1.

The cut is feasible since 2H+1/10<pi. At that cut F<=G-1/10<=C-1/10.
Convexity on the remaining interval and the uniform Bo gap exclude
all gamma>=theta+1/10. Thus each actual z_i is in [0,1/10).

On this actual theta/apex domain, the OLD fixed-parameter inequalities
apply directly: C0-G0(theta_i)>=2e_i^2, with e_i=theta_i-theta0,
and G0(theta_i)-F0(theta_i,gamma_i)>=z_i. It follows, without any
assumption that the total old deficit is small, that

    2 sum e_i^2+sum z_i+sum sigma_i+sum tau_i <= D_old,
    D_old=3C0-S(lambda0,t0).

The old budget is used here only to estimate parameter sensitivity.
There is no additional subtraction from the global area budget and no
use of assigned-cell lens losses or terminal BI credits.

## 4. Centered sensitivity cancels the independent-support first derivative

Write dl=lambda-lambda0, dt=t-t0, and x0=sin theta0. The exact identities
and the feasible old peak give

    S_new=S_old-4dl sum x_i^2+dt sum gamma_i,
    C_new>=C0-4dl*x0^2+dt*theta0.

Consequently

    D_new>=D_old+4dl(sum x_i^2-3x0^2)-dt(sum gamma_i-3theta0).

The same paid old budget and Cauchy imply

    |sum x_i^2-3x0^2| <= sqrt(3D_old/2),
    |sum gamma_i-3theta0| <= sqrt(3D_old/2)+D_old.

Therefore throughout the rectangle

    D_new >= f(D_old),
    f(D)=(199/200)D-(1/40)sqrt(3D/2).

The accepted actual-anchor theorem gives D_old>1/40. For D>=1/40,
f is strictly increasing: by rationalizing the square-root difference,
its secant slope is at least

    199/200-(3/80)/(2sqrt(3/80))
      >199/200-(3/80)/(38/100)>0.

Also sqrt(3/80)<3873/20000, by a strict integer-square comparison, so

    f(1/40)>199/8000-3873/800000=16027/800000>1/50.

This contradicts D_new<=1/50. The same calculation actually establishes
the stronger displayed fee: if D_new<=16027/800000, all localization
steps still hold after replacing 1/50 by 16027/800000 (the saved endpoint,
Bo and convex-cut gaps are larger), and the strict final inequality
again contradicts that assumption. Every guard for this substitution
is included in verify_parameter_uniform.py.

## 5. Closure, allocation and limits of this result

As theta tends to zero or pi/2, the endpoint score support remains
continuous. A zero-area limiting triangle has actual normalized
area+anchor score at most 1/2+t*pi<1; all nondegenerate scores are <=C.
Since C>7/5, a degenerate face cannot occur in a hypothetical triple
deficit <=16027/800000. Thus all divisions and acute-tip estimates are
preceded by uniform localization, and degeneracy cannot restore the
excluded high-score triple.

The full-apex support C is defined above exactly; this theorem does
not provide its algebraic expression or a numerical certificate for
an entire-cover row at each parameter value. Whole-cover use requires
matching lambda/t throughout the dual, circumradius and fan incidence,
and disjoint charged marked-face/anchor sets. The retained term
Q_actual=sum assigned-cell lens+lambda sum base^2, other supports,
variance and reserved BI credits must each be allocated once.

The frozen fixed-parameter .025 theorem, parameterized-radius addendum,
and their reviews are dependencies. The guard below checks only new
rational sensitivity/localization comparisons against those accepted
inputs. It does not rerun a frozen local/global checker or optimizer.
