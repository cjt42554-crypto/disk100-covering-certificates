# Continuous-radius actual-anchor two-face lemma

Author: `/root/q2_radius_uniform`. Status: paper proof, awaiting independent static review. Mathematical executions and mathematical imports for this extension: **zero**. No old checker was rerun, imported, or sampled. This extension uses the already accepted exact scalar enclosures from the acquired Pro proof and its one isolated root replay as frozen inputs; all changes of radius are justified below by inequalities on entire intervals.

## 1. Precise theorem

Fix lambda=7/50 and t=1071883/20000000. Let

- R0=5661/50000=.11322 and U=29/250=.116;
- R be any real number in the **closed continuous interval** [R0,U];
- (p,c0,c1) and (p,c1,c2) be the same two consecutive positively oriented nondegenerate actual triangles, each with circumradius at most R;
- for each actual base (ci,ci+1), qi be an actual unit anchor covered by both endpoint R-disks;
- the **same actual** centers satisfy |c1|²≤1−R² and |c0|²,|c2|²≥83/100.

Keep actual signed A_i=−det(ci+1−ci,qi−ci)/2, including negative A. Write x_i=|ci+1−ci|/(2R), gamma_i for the actual apex angle, and

S2=Σ[(|Ti|+A_i)/R²−4 lambda x_i²+t gamma_i].

C is the complete unrestricted apex/skew single-face support at these lambda,t, including both apex endpoint branches; its normalized expression is independent of R. Then

**D=2C−S2 > 9/400**

uniformly for every R in [R0,U]. No middle-halfarc variance, Q, BI, independent J, or another copy of tau is charged. This is a conditional two-face theorem and is not a global 100-disk certificate.

Neither equal endpoint distances, radial-inner position, equality of the two angles/circumradii, reflection symmetry, nearest-owner conditions, convexity of the boundary chain, nor actual A≥0 is assumed. Parameters lambda,t are fixed: this file proves no parameter perturbation result.

## 2. Frozen radius-independent input

The acquired source is `research/pro_sources_20261005/sources/two_anchor_tau_lemma/two_anchor_tau_lemma/PROOF.md`, SHA256 `8d2289a8816aefd05bea1e12783d5d78e925c31773a0d6222fee60fba3476c6e`. Its scalar certificate is `verification.json`, SHA256 `b83b9882560b3d3b6c98ef077e1e9e15482b0a671b43c9cf4a63c4faade8f30e`. `research/pro_tau_review/STATIC_REVIEW.md` independently audits its derivation, and records the root's already completed, single, matching local replay. Exact input pins are in this directory's SOURCE_PINS.json.

The normalized support/localization of source Sections 1–2 contains no R. Therefore for a hypothetical D≤d=9/400, at every R in the new interval it gives, without changing a guard,

- the unique support peak .828<theta*<.829;
- theta_i=asin x_i, z_i=gamma_i−theta_i, sigma_i≥0, tau_i≥0;
- 2Σ(theta_i−theta*)²+Σz_i+Σsigma_i+Σtau_i≤D;
- s=(theta0+theta1)/2∈(.753,.904), Delta=theta0−theta1 with |Delta|≤3/20;
- Z=(z0+z1)/2∈[0,9/800], T=tau0+tau1∈[0,9/400];
- .72<theta_i<.94, x_i²>13/30;
- d_i is the source skew convention (side lengths 2k cos(gamma/2+d_i), 2k cos(gamma/2−d_i)), eta=(d0−d1)/2, |eta|<12/125;
- Σsigma_i≥c_eta eta², c_eta=223532401/90000000.

Hence the only budget is

4(s−theta*)²+c_eta eta²+Delta²+2Z+T≤D.                 (2.1)

The source's phase convention is retained exactly:

alpha=3s/2+Z/2+eta, h=s/2−Z/2−eta,
A0=alpha−theta0, A1=alpha−theta1.

The R-independent angular box is

alpha∈[2067/2000,11661/8000], h∈[2199/8000,137/250],
A_i∈[411/2000,5029/8000].

In particular, all these angles are in the first quadrant, B0=sin alpha>6/7, cos alpha>1/10, sin A_i>81/400, and alpha>1. The source's Cauchy scalar inequality also remains unchanged:

9/100+1/c_eta < (71/100)².                           (2.2)

## 3. Exact geometry remains valid at every new R

Source Sections 3–6 derive identities, not radius-point fits. For every 0<R<1 they give outward-tip directions exp(−i alpha), exp(i alpha) after a common rotation and normals N0=exp(−iA0), N1=exp(iA1), T_i=iN_i. Actual relative error is

e_i=(qi−qi*)/R=−h_i N_i+ell_i T_i, h_i=tau_i/x_i,
(x_i+|ell_i|)²+(cos theta_i−h_i)²≤1.                  (3.1)

Consequently |ell_i|≤h_i cot theta_i and |e_i|≤tau_i/x_i². Put

g=(e0+e1)/2, d_b=(e1−e0)/2, b=(0,B0)+d_b, B=|b|,
n=(b_y,−b_x)/B, u=sqrt(1−R²B²).

The two unit anchors exactly eliminate common translation:

c1=±un−R[(cos alpha,0)+g].                          (3.2)

Along the straight error path v e_i, 0≤v≤1, the endpoint lenses remain valid by convexity. Both anchor-minus-middle vectors have length at most one, so B(v)≤1. With P=(15/13)T≤27/1040,

|g|,|d_b|≤P,
B(v)≥6/7−27/1040>83/100,
u(v)²≥1−U²> (99/100)².                             (3.3)

Thus all denominator domains used in the old derivation are valid uniformly. Also n_x>.83, and k=(cos alpha,0)+g satisfies

k·n>(1/10)(83/100)−27/1040>0.

On the negative branch, |c1|²=u²+2Ru(k·n)+R²|k|²>u²≥1−R², contradicting the stated same-actual-middle cap. Hence the positive branch in (3.2) is forced at each R. No cap is assumed along the interpolated fan.

With

E*=−cos h(cos Delta,sin Delta),
F=(−sin h sin Delta,B0+sin h cos Delta), E=E*−g,

the source identities are

M=(|c0|²+|c2|²)/2=u²+2Ru n·E+R²(|E|²+|F|²),
M*=1+R²−2R cos chi cos Delta,
chi=h+asin(R sin alpha).                            (3.4)

Every expression still refers to the same actual triangles and anchors.

## 4. Signed linear support and continuous transfer of all eight slabs

The gradients L_i and kappa_i are the exact source formulas (6.1), (6.3), with their occurrences of R now variable. The directional lens support requires L0·T0<0<L1·T1. The source's sufficient bound reduces to

[R(99/100)−R²−7R/40]/5−(23/20)R²
=R[163/1000−(27/20)R]>0,                            (4.1)

because R≤29/250=.116<163/1350. Thus the same signed support is valid on the whole new interval. It is not replaced by a norm estimate.

For k0=41/200, the exact residual corresponding to kappa0≤k0 is

G_k0=k0 sin²(theta0)−numerator(kappa0)
=A+Q cos Delta+P sin Delta,                         (4.2)

where A,Q,P are source (6.5); kappa1 uses Delta reversed. The source's eight **whole** s slabs cover [753/1000,113/125], retaining the full omega=Z/2+eta interval [−12/125,813/8000] and the full Delta interval [−3/20,3/20]. For any slab let a_−,a_+ and h_−,h_+ be its source angular endpoints; put

rho_R(a)=sqrt(1−R²sin²a)−R cos a.

The actual lower enclosure used in the original certificate is the mathematical function

A_−(R)=k0/2−R rho_R(a_+)sin a_+,
Q_−(R)=−k0 cos(2s_−)/2+R²rho_R(a_−)sin a_− cos h_+,
P_−(R)=k0 min{sin(2s_−),sin(2s_+)}/2−R rho_R(a_+)cos h_−,
P_+(R)=k0/2−R rho_R(a_−)cos h_+,
L(R)=A_−(R)+min{Q_−(R),Q_−(R)cos(3/20)}
     −max{|P_−(R)|,|P_+(R)|}sin(3/20).              (4.3)

The omission of 1/u from the positive term in Q_− is the original safe relaxation. These endpoint enclosures remain valid throughout the enlarged R interval: rho_alpha=R sin(alpha)[1−R cos(alpha)/u]>0 follows from u²−R²cos²(alpha)=1−R²>0, and (rho sin alpha)_alpha=−R cos(2alpha)+cos(alpha)[1−2R²sin²(alpha)]/u>0 follows from alpha∈(pi/4,pi/2) and R≤.116. The cosine and sine slab endpoint comparisons have no R dependence. Radius-dependent trigonometric endpoint sampling is absent.

Here is a uniform analytic bound on (4.3), valid even on the larger radius interval R0≤R≤3/25. Let c=cos a, beta=sin a and v=sqrt(1−R²beta²). Then

H=∂_R(R rho_R)=v−2Rc−R²beta²/v.

Since v>.99, R≤.12, c≤1,

H≥.99−.24−(.12)²/.99>0,
H≤v≤1.                                               (4.4)

Also

∂_R(R²rho_R)=R[2v−3Rc−R²beta²/v]>0                  (4.5)

by the lower estimate R[1.98−.36−(.12)²/.99]>0. Consequently A_− has derivative ≥−1, Q_− is nondecreasing, and both P endpoints are 1-Lipschitz. Absolute value, max and min preserve those Lipschitz/monotonicity statements; the multiplier cos(3/20) is positive. Therefore for **each of the eight complete slabs**, throughout R≥R0,

L(R)≥L(R0)−[1+sin(3/20)](R−R0)
    ≥L(R0)−(23/20)(R−R0).                           (4.6)

The second inequality is strict when R>R0; at R=R0 it is equality. The subsequent strict residual comes from the already strict frozen margin, so no strict coefficient is incorrectly multiplied by zero.

The original frozen certificate gives a rigorous lower enclosure L(R0)>1/4000 in every slab. Its original outward rounding is only used at R0; (4.4)–(4.6) are exact analytic inequalities for the functions, so no unproved rounding transfer is involved.

For new k=213/1000, increase the **exact pointwise residual**, rather than claiming a fresh conservative slab enclosure increases by this amount:

G_k=G_k0+(k−k0)sin²(theta_i)
   >1/4000−(23/20)(R−R0)+(1/125)(13/30).

As U−R0=139/50000, all eight slabs now satisfy the uniform strict bound

G_k > 1/4000−3197/1000000+13/3750
    =1559/3000000>0.                                 (4.7)

Hence kappa0,kappa1<213/1000 for every R∈[R0,U]. A sharper radius-dependent interface, if useful, is

kappa_i < 41/200+(69/26)(R−R0).                      (4.8)

No old slab was rerun and no new radius points were tested.

## 5. Uniform Hessian, complete remainder and one tau cost

For f(b)=sqrt(1−R²|b|²) Jcw b/|b|, the independently reviewed exact directional Hessian identity is unchanged. The source bounds give

||Df||≤100/83,
||D²f[b][w,w]||≤[R²/u³+2/(sqrt(3)u|b|²)]|w|².

Uniformly, R²≤U²<14/1000, u³>(99/100)³>97/100 and u|b|²>(99/100)(83/100)²>17/25. Using 2/sqrt(3)<7/6,

R²/u³+2/(sqrt(3)u|b|²)
<14/970+175/102<7/4.                                (5.1)

The radial first derivative is <(.014)/.99<100/83, and the tangential one is u/|b|<100/83. Thus both bounds hold along the whole path, not only at the endpoints.

The complete nonlinear remainder is the source expression

−R²|d_b|²+R²|g|²+2R[r_f·E*−(f(b)−f(b*))·g],

where |r_f|≤(7/8)|d_b|² and |E*|≤1. Dropping only its negative first term yields

remainder ≤[R²+(7/4)R+(200/83)R](15/13)² T².

Since 200/83<241/100 and R≤29/250,

coefficient <[14/1000+(29/250)(104/25)](225/169)
            =111726/169000 <2/3.                    (5.2)

The reusable radius-uniform actual-tau estimate is therefore

M≤M*+kappa0 tau0+kappa1 tau1+(2/3)T²,
M≤M*+(213/1000)T+(2/3)T²
 ≤M*+(57/250)T≤M*+(23/100)T,                        (5.3)

because T≤9/400 and 213/1000+(2/3)(9/400)=57/250=.228<.23. At T=0 these are non-strict bounds and do not assert a false strict product. The final contradiction supplies strictness independently.

## 6. Radius-dependent angular envelope

As in the exact source identity,

chi_Z=−1/2+R cos alpha/[2sqrt(1−R²sin²alpha)]<0,

because u0²−R²cos²alpha=1−R²>0. Thus discarding nonnegative Z only increases chi. On the connecting shape rectangle define

g_R(s,eta)=s/2−eta+asin(R sin(3s/2+eta)).

Here alpha>1, R≤.116<3/25 and u0>.99. The elementary cosine bound cos(1)≤1−1/2+1/24=13/24 gives

q=R cos alpha/u0 < (3/25)(13/24)(100/99)=13/198,
0<∂_s g_R=1/2+(3/2)q<79/132<3/5,
−1<∂_eta g_R=−1+q<0.                              (6.1)

These inequalities cover the entire rectangle from (theta*,0) to (s,eta), including both skew signs; no symmetry is used.

The frozen source baseline gives g_R0(theta*,0)<261/500. At fixed theta*,

∂_R g_R(theta*,0)=sin(3theta*/2)/sqrt(1−R²sin²(3theta*/2))<100/99.

Integration on the entire radius interval yields

g_R(theta*,0)<b(R):=261/500+(100/99)(R−R0).           (6.2)

This replaces the original constant .522 baseline, which cannot simply be inherited over the enlarged interval. For E=4(s−theta*)²+c_eta eta², (2.2), (6.1) and Cauchy give

chi≤b(R)+(3/5)|s−theta*|+|eta|
   ≤b(R)+(71/100)sqrt(E).                            (6.3)

Throughout D≤9/400 this angle is positive and at most

a(R)=1257/2000+(100/99)(R−R0)
 ≤124999/198000<79/125=.632<pi/2.                    (6.4)

Thus all cosine monotonicity steps below are justified.

## 7. Close the uniform fee using one budget

From (3.4), (5.3), and 1−cos Delta≤Delta²/2,

M≤1+R²−2R cos chi+R Delta²+(23/100)T.

R<23/100 allows Delta² and T to be combined with just the remaining budget in (2.1), giving

M≤Phi_R(E):=1+R²−2R cos[b(R)+(71/100)sqrt(E)]
             +(23/100)(9/400−E), 0≤E≤9/400.         (7.1)

No term is charged twice. For 0<E≤9/400,

Phi_R'(E)≥R0(71/100)sin(261/500)/(3/20)−23/100.

The frozen source exact scalar guard says the same expression with subtraction 11/50 is >1/25. Subtracting the extra 1/100 proves

Phi_R'(E)>3/100>0.                                  (7.2)

Continuity covers E=0. Therefore the maximum occurs at E=9/400 and

M≤F(R):=1+R²−2R cos a(R).                           (7.3)

This last envelope improves, rather than worsens, as R increases. Indeed (6.4) gives cos a(R)>4/5 and sin a(R)<79/125. Consequently

F'(R)=2R−2cos a(R)+2R sin a(R)(100/99)
 <29/125−8/5+(29/125)(79/125)(100/99)
 <29/125−8/5+21/125=−6/5.                           (7.4)

The last rational comparison is 229100<259875. Thus the whole continuous interval is bounded by its **left** endpoint. The frozen, already replayed exact norm guard gives

F(R)≤F(R0)=1+R0²−2R0 cos(1257/2000)<83/100−3/10000.

This contradicts the same actual two outer norms ≥83/100. The assumed D≤9/400 is impossible, proving the uniform theorem. Closed inequalities in the circumradius, endpoint coverage, middle cap and outer floors are allowed. Degenerate limits are excluded from the low-deficit domain by the source's radius-independent localization; its closure argument is unchanged.

## 8. Interfaces, limitations and validation request

The source's optional halfarc identity and signed lens enclosure remain available at the current R if the actual anchors really are the endpoints of the middle connected assigned arc. In particular a_mid=asin(R|b|). They are unnecessary for the fee. No variance has been consumed here.

**The outer floor is an explicit premise.** Full assigned arc length >.147 does not automatically imply |c|²≥.83 for every R≤.116. Its radius-dependent lower root is

cos(.0735)−sqrt(R²−sin²(.0735)),

which decreases in R. A separate continuous arc-to-norm lemma or a stronger arc threshold is needed to transport this theorem into an enlarged global-radius argument. Likewise the same actual middle cap must be established for the actual fan to which the score is assigned; norm-capping a different fan is not a proof of this premise.

The genuinely radius-sensitive constants in this derivation are the signed gradient support, the eight-slab linear coefficient, the Hessian/remainder, and baseline chi. Each has been re-proved on the whole interval above. The numerical support peak, localization and eight initial slab margins are unchanged normalized inputs, not newly executed tests. The ultimate strong norm contradiction becomes stronger with R under this envelope.

Suggested optional new exact-guard validation, **not run here**: check the finite rational inequalities (4.1), (4.4)–(4.7), (5.1)–(5.2), (6.1), (6.4), and (7.4), and check unchanged frozen certificate pins/margins by parsing metadata. This can be an independent new uniform-domain guard; it must not execute or import the old 48-guard checker. A local run is not needed to interpret the paper proof and is not claimed by this author.
