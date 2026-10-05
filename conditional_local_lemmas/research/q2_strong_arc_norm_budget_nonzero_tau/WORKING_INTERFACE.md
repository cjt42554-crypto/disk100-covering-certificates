# Nonzero actual-anchor slack: exact translation path and directional derivative interface

This is a working analytic interface, not a proved strong-arc fee.
The frozen general zero-tau proof is
`../q2_strong_arc_norm_budget/TAU_ZERO_NORM_PROOF.md`, SHA
`c0077a0fa8ca16bf56c4c13a93571d048082c11ba60f1da4c4f7a6b3ffa549a1`.
It must not be widened silently to positive tau.

Assume the same R22, lambda, t, consecutive two-face domain, actual
unit anchors, circumradius cap and middle reflected norm cap. Suppose
the actual-score deficit `D=2C-S0-S1<=223/10000`.
The later joint deficit with the nonnegative middle variance also
implies this inequality whenever it is below the parent's target.
The existing per-face budget applies to the actual anchors:

    2 sum e_i^2+sum z_i+sum sigma_i+W<=D,
    W=tau0+tau1.

Here `tau_i` is actual signed-A slack; no BC-lens, middle variance,
outer-Q or additional fan credit is used in these formulas.

## 1. Exact positive unit branch at high score

Retain the zero-tau shape angles

    s=(theta0+theta1)/2,
    E=s-theta*, Delta=theta0-theta1,
    Z=(z0+z1)/2, eta=(d0-d1)/2,
    alpha=3s/2+Z/2-eta.

The same simultaneous skew estimate as in the frozen proof gives

    D>=4E^2+Delta^2+2 sum x_i sin(d_i)^2+W.

Because `x_i>.64` and
`sin^2 d0+sin^2 d1>=2sin^2 eta`,

    2.56sin^2 eta<=D<=.0223,
    |eta|<.095.

The last conclusion uses the explicit elementary bound
`sin(.095)>.094` and `2.56(.094)^2>.0223`.
Thus `.75<s<.904`, `Z<=.01115`, and

    1<alpha<1.456575<1.46<pi/2,
    cos alpha>.1, sin alpha>.8.

These trigonometric guards follow from lower alternating Taylor
polynomials at 1 and 1.46. They are paper inequalities here, not
new program evaluations.

Let `a_i` be the normalized unit vector from the outer center to its
outer lens tip, and `b_i` the normalized unit vector from c1 to that
same tip. These are the exact phase directions in the frozen proof.
Let `epsilon_i=(q_i-tip_i)/R`. Put

    B=(b0+b1)/2=cos(alpha)n0,
    L0=(b1-b0)/2=i sin(alpha)n0,
    U=(epsilon0+epsilon1)/2,
    V=(epsilon1-epsilon0)/2,
    L=L0+V,
    n=-iL/|L|,
    u=sqrt(1-R^2|L|^2).

The accepted anchor bound gives

    |U|,|V|<=625W/512<.028.

Hence `|L|>.8-.028=.772`, and the projection of n0 onto n is
greater than `.99`. For example the normalized-vector difference
bound `|n-n0|<=2|V|/(|L0|-|V|)<.073` suffices, since
`n dot n0=1-|n-n0|^2/2>.99`.
Therefore

    (B+U) dot n >(.1)(.99)-.028=.071>0.

Both actual anchors are covered by the common middle disk, so

    |B+U|^2+|L|^2
       =(|b0+epsilon0|^2+|b1+epsilon1|^2)/2<=1.

Thus `u>=sqrt(1-R^2)`. Unit-anchor equations give their mean
`A=(q0+q1)/2=+-u n`, while

    c1=A-R(B+U).

The negative branch would have
`c1 dot n=-u-R(B+U) dot n<-u<=-sqrt(1-R^2)` and violate the
middle norm cap. Therefore the actual branch is positive, `A=u n`.
This handles the branch before perturbative differentiation.

## 2. Exact unit-anchor translation path at fixed fan shape

For `0<=t<=1`, define

    L_t=L0+tV,
    n_t=-iL_t/|L_t|,
    u_t=sqrt(1-R^2|L_t|^2),
    A_t=u_t n_t,
    c1(t)=A_t-R(B+tU),
    o0=c0, o1=c2,
    f_i=(o_i-c1)/R=b_i-a_i, i=0,1,
    o_i(t)=c1(t)+R f_i, i=0,1,
    q0(t)=A_t-RL_t, q1(t)=A_t+RL_t.

All fan vertices, including p, are translated rigidly with c1(t).
Their actual angles, areas, circumradii and base lengths do not
change. Both constructed q's have unit norm at every t. Their
offsets from the fixed-shape endpoint disks are exactly

    (q_i(t)-c1(t))/R=b_i+t epsilon_i,
    (q_i(t)-o_i(t))/R=a_i+t epsilon_i.

The intersection of the two endpoint unit disks is convex, so those
offsets remain covered throughout the path. Each signed actual
anchor slack is exactly `t tau_i`, since the signed height is linear
in this displacement. The score deficit along the path is

    D(t)=D(0)+tW,
    D(0)=D-W.

The path does not preserve the given outer assigned arcs or every
global center norm; such preservation is unnecessary for a norm
derivative bound and must not be assumed. The endpoint t=1 is the
original actual configuration by the positive-branch argument.
The endpoint t=0 is a legitimate zero-anchor-slack fan with unit
tips and the same shape.

## 3. Mean norm square, exact derivative and signed lens support

Write

    H=(a0+a1)/2, J=(a1-a0)/2,
    F=L0-J,
    H_t=H+tU,
    C_t=A_t-RH_t.

Then C_t is the actual mean of the two outer centers at time t,
their physical half-difference is RF, and

    M(t)=(|c0(t)|^2+|c2(t)|^2)/2
         =|C_t|^2+R^2|F|^2.                    (M)

This cancels the half-difference anchor error exactly; it is better
than separately bounding the apex translation.
Let `r_t=|L_t|`, `ell_t=i n_t`. Decompose a variation dL as
`xi n_t+zeta ell_t`. The exact derivative of A is

    dA=-(u_t/r_t)xi ell_t
            -(R^2 r_t/u_t)zeta n_t.             (DA)

Consequently

    M'(t)=g_U dot U+g_V dot V,
    g_U=-2R C_t,
    g_V=2(DA)^T C_t.

For `C_t=c_n n_t+c_l ell_t`,

    g_V=-(2u_t/r_t)c_l n_t
             -(2R^2 r_t/u_t)c_n ell_t.

Set `g0=(g_U-g_V)/2`, `g1=(g_U+g_V)/2`. Then

    M'(t)=g0 dot epsilon0+g1 dot epsilon1,

with exact component forms

    g0=(-R c_n+(u_t/r_t)c_l)n_t
          +(-R c_l+(R^2r_t/u_t)c_n)ell_t,
    g1=(-R c_n-(u_t/r_t)c_l)n_t
          +(-R c_l-(R^2r_t/u_t)c_n)ell_t.       (G)

No independent mean/chord losses have been added here. The actual
endpoint lens directions retain their signs. For base i, let e_i
be the directed unit base tangent and n_i the normal toward the
apex. Write

    epsilon_i=s_i e_i+v_i n_i,
    v_i=tau_i/x_i>=0,
    s_i^2+2x_i|s_i|+v_i^2<=2y_i v_i,
    y_i=sqrt(1-x_i^2).

The exact directional lens support gives

    g_i dot epsilon_i
      <=[g_i dot n_i+|g_i dot e_i| y_i/x_i]tau_i/x_i. (L)

Thus a rigorous whole-domain bound on the bracketed coefficients
in (L), integrated along the path, would give the missing linear
tau payment for M. A Euclidean norm replacement discards a useful
signed radial cancellation and may be too loose.

The minimum outer norm square is exactly

    min(|c0(t)|^2,|c2(t)|^2)
       =M(t)-2R|C_t dot F|.                    (Min)

The asymmetry penalty in (Min) is relevant: bounding only M can
give an artificial benefit from rotating the chord. A successful
general proof should combine (L), (Min), the frozen angular/skew
payment, and W within one budget. No uniform sufficient coefficient
for that combination is proved in this working note.

## 4. Status and execution

The positive branch and the path/derivative identities are analytic
interfaces intended for independent checking. The full strong-arc
actual-A bound at positive tau is still open in this note. No new
mathematical program, imported mathematical module, old checker,
optimizer or scan was run. The parent has a separately frozen
directional lens lemma; this note does not alter or replace it.

Any future quantitative bound or certificate must retain a separate
source and full attempt ledger rather than treating (G)--(L) as an
already verified strong-fee calculation.
