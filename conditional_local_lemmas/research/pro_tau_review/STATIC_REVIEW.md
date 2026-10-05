# Independent static audit: Pro two-actual-anchor tau theorem

Status: **STATIC PASS — cleared for one isolated local exact-guard replay.**
This audit does not constitute a local numerical replay. The acquired proof is a conditional two-face theorem, not a global 100-disk certificate. No mathematical source was executed or imported during this audit. The source's reported 48 remote PASS guards are retained as source evidence and are not represented as local executions.

Subsequent root replay: **PASS**, one mathematical guard execution in `research/pro_tau_local_replay_20261005`, exit 0. Its receipt records 48 guards, eight full slabs, JSON byte MATCH, source log exact text MATCH, and all six acquired source hashes unchanged. I inspected the receipt and compared these output bytes independently as metadata, without executing/importing the checker. This is the root's one execution, not a second audit execution. `REPLAY_METADATA_REVIEW.json` pins that evidence.

## Frozen source and theorem interface

Fully read the six files in `research/pro_sources_20261005/sources/two_anchor_tau_lemma/two_anchor_tau_lemma/`: `PROOF.md`, `RESULT_FOR_HANDOFF.md`, `README.md`, `verify.py`, `verification.json`, and `verification.log`, together with `TAU_SOURCE_ACQUISITION.json`. Exact pins are in `SOURCE_PINS.json`.

Parameters are R=5661/50000, lambda=7/50, t=1071883/20000000. C is the unrestricted apex/skew single-face support at those same parameters. Two consecutive positively oriented nondegenerate actual triangles share the actual edge pc1; their circumradii are at most R. Each actual unit anchor is covered by both endpoints of its own base. Actual signed A is retained; A may be negative, the two distances need not be equal, and no radial-inner, mirrored, equal-angle, nearest-owner, or convex-boundary-chain assumption is used.

The additional conditions are the **same actual** middle center |c1|²≤1−R² and both **actual** outer centers |c0|², |c2|²≥83/100. Under precisely this domain, the source proves D=2C−S2>9/400. The exact difference from the requested fee is

`1740724429841728225095328099 / 8011730250000000000000000000000 > 1/5000`.

The variance coefficient 1221/(640R²) is only added as a nonnegative corollary. Neither the main proof nor its guard checker pays middle-halfarc variance, Q, BI, or an independent J−A fee. Transport of the actual norm cap and the outer-arc-to-norm facts into a selected global fan remains an external premise to be checked by the parent proof.

## 1. Full support, localization, and single budget

The score identity (1.1), sigma≥0 from fixed-base/apex area maximization, and tau≥0 from the endpoint lens hold for actual signed A. Circumradius≤R supplies theta≤gamma≤pi−theta. The support expression is convex in gamma, so both endpoint branches F and G must be included. The source includes both; G<6/5 while F(5/6)>7/5, hence C=max F. The proof of F's unique peak uses separate angular ranges and does not assume global concavity. F'(0)>0, F''<0 through pi/3, and F'<0 afterwards supply the required uniqueness.

The .7/.96 localization uses the full two-branch support, not merely a low-apex subdomain. On this localized domain, L_z≥1−cos(theta)−t>0 for every permitted z. Up to z=.1 the explicit cotangent bounds give L_z>1; a z≥.1 would already exceed D≤9/400. The curvature lower bound on [.7,.96] is valid because sin(theta)≥sin(.7), sin(2theta)≥sin(1.92), and cos(2theta)≥cos(1.92). Thus the sole budget (1.3) is valid, with 2 times the squared theta displacements plus z, sigma, and tau paid exactly once.

The subsequent theta, s, Delta, Z, and T box follows from this budget. The sigma coefficient increases with theta and decreases with z on the stated acute domain. The skew range first uses |d|<pi/2 and sin²d≤9/500, then the Taylor lower bound for |sin d|/|d|. Minimizing d0²+d1² at fixed eta=(d0−d1)/2 gives the correctly recorded c_eta=223532401/90000000. No skew term is spent twice.

## 2. Phase convention and exact geometry

I rederived the source convention: with sides 2k cos(gamma/2+d) and 2k cos(gamma/2−d), the base is `2 i x exp(i(m−d))`. Its midpoint and outward tip are `exp(im)[k(1+cos gamma)cos d + i k(1−cos gamma)sin d]` and `exp(im)[U cos d + i V sin d]`. Hence the half-separation of the two tip-minus-middle directions is

`alpha=3s/2+Z/2+eta`, with `h=2s−alpha=s/2−Z/2−eta`.

The source signs are therefore correct for its stated d convention. Reversing both d signs changes the convention, not the geometry; mixing conventions would be an error. After one common rotation the normals are exp(−i(alpha−theta0)) and exp(i(alpha−theta1)). Direct addition of the two base vectors gives exactly the E* and F of (5.4). This derivation leaves the two faces asymmetric and does not require a mirrored optimizer or shared-edge-relaxation equality beyond the original actual geometry.

At tau=0, the two unit-anchor translations give c1=±u0−R cos(alpha). The positive branch yields the two individual outer norms in (3.2) and the mean M* in (3.3). The sign of chi_Z is exact: u0²−R² cos²(alpha)=1−R²>0. Discarding Z later increases chi but does not recreate a second budget.

## 3. Actual tau, lens, and positive branch

For each actual anchor, `e=−h N+ell T`, h=tau/x≥0, has the correct outward-normal sign. Endpoint coverage is exactly `(x+|ell|)²+(cos(theta)−h)²≤1`, including negative A. This implies the stated signed tangent range and |e|≤tau/x². The latter is used only for the quadratic remainder and coarse localization; first-order support retains the sign.

The exact two-unit-anchor translation is c1=±un−Rk with b=(a1−a0)/2. Along the linear error path, both anchor-minus-middle relative vectors remain in their own endpoint lenses by convexity. Their lengths remain at most one, so |b|≤1. Meanwhile |b|>6/7−27/1040>83/100 and u>.99. Since b_y>.83 and |b|≤1, n_x>.83. Thus k·n>.1·.83−27/1040>0. On the negative branch,

`|c1|²=u²+2Ru(k·n)+R²|k|² > u² ≥ 1−R²`,

contradicting the actual middle cap. This proves the required positive actual branch. The proof does not assume that the original middle cap persists at every interpolated fan; the homotopy is used to bound the smooth positive-branch formula, whose lens and denominator domains are already established.

## 4. Linear support and eight continuous slabs

The gradient (6.1) was independently rederived, as was the signed lens support producing (6.3). The tangent signs follow from the stated whole-box bounds on L_ix, L_0y, L_1y, and sin(A_i)>.2025; the single exact sign guard is sufficient after these analytic reductions. Algebraically, multiplying kappa0≤k by sin²(theta0), with theta0=s+Delta/2, gives the exact A+Q cos(Delta)+P sin(Delta) expression (6.5). Kappa1 is the same expression with Delta reversed.

All eight s slabs are closed intervals, obtained by equal subdivision of [753/1000,113/125], with adjoining endpoints identical. They cover the entire localized s range. In **each** slab, the whole omega range [−12/125,813/8000] is retained; Delta ranges over all [−3/20,3/20]. There are no numerical sample points substituted for a continuous domain.

For alpha in this first-quadrant box, rho(alpha) increases because rho'=R sin(alpha)(1−R cos(alpha)/u)>0. Also

`(rho sin alpha)'=−R cos(2alpha)+cos(alpha)(1−2R² sin²(alpha))/u>0`.

The code uses these monotone endpoints for rho and rho B. The Q lower bound safely omits division by u≤1 from a positive term, uses cos(h_max), and uses cos(2s_lower) for its negative term. P is enclosed with the minimum of sin(2s) at the slab's endpoints, upper bound one, and the independent rho/cos h bounds. The final min/max reduction in (6.6) handles either sign of Q and both signs of Delta. The certificate reports every continuous-slab lower bound strictly above 1/4000, with its minimum in slab 1. Reproducing those numerical intervals locally remains the next step.

## 5. Hessian and complete remainder

A separate static reviewer independently rederived this part. Remove the fixed clockwise rotation and write f(b)=u(r)e_r, r=|b|. For w=a e_r+v e_t,

`D²f[w,w]=−[R² a²/u³+v²/(u r²)] e_r−[2av/(u r²)] e_t`,

up to that fixed rotation. Here u''=−R²/u³ and u'/r−u/r²=−1/(u r²). Consequently

`|D²f[w,w]|≤R² a²/u³ + sqrt(v⁴+4a²v²)/(u r²)`

and `sqrt(v⁴+4a²v²)≤(2/sqrt(3))(a²+v²)`. This is (6.7); the code's upper estimate using u>.99, r>.83, and 2/sqrt(3)<7/6 is valid. Df has radial norm R²r/u and tangential norm u/r, so M1=100/83 is also valid throughout the path.

Let r_f=f(b)−f(b0)−Df(b0)d_b. Expanding the actual mean M completely, after subtracting its linear part, leaves

`−R²|d_b|²+2R r_f·E*−2R(f(b)−f(b0))·g+R²|g|²`.

Discard the first negative term, use |E*|≤1, |r_f|≤(7/8)|d_b|², and the first-derivative path bound. This yields exactly (6.8), rather than omitting a positive quadratic term. The reported scalar guard makes its coefficient <13/20. Together with kappa_i<41/200 and T≤9/400 it gives the non-strict uniform bound M≤M*+(11/50)T. At T=0 there is no unjustified strict multiplication: the final contradiction comes from the independent strict norm envelope.

## 6. One-budget closure and degeneration

The mean-shape budget is `4(s−theta*)²+c_eta eta²+Delta²+2Z+T≤D`. Dropping Z increases chi. The mean-value estimate for chi uses the full rectangular domain connecting (theta*,0) to (s,eta), which lies within the already certified alpha bounds. Cauchy gives chi≤.522+.71 sqrt(E) with E=4(s−theta*)²+c_eta eta². The Delta correction is at most R Delta²; since R<11/50, Delta² and T are merged under **one remaining budget** D−E. Thus the Phi bound (7.3) is valid.

Phi'(E)>1/25 for E>0 and continuity at E=0 cover the square-root endpoint. The right endpoint is E=9/400, giving angle 1257/2000. The source's claimed exact upper interval is below .83−.0003, so it contradicts the two actual outer squared norms ≥.83. Equality in the input cap, circumradius, coverage, or outer norms is allowed. T=0 is included. At x=0 the supported score tends to at most t pi, and x=1 belongs to F/G; neither can have D≤9/400. Within the contradiction domain, all the x, B, u, sine, and angular denominators are bounded strictly away from zero.

The optional HALFarc interface (Section 8) is relevant only when the two actual anchors really are the endpoints of the middle connected assigned arc. Such an R-covered unit-circle arc with R<1 has length below pi, which selects the minor-arc branch of asin(R|b|). This optional identification is not needed for the fee. Section 9's arc-to-norm estimate follows from endpoint averaged distances and monotonicity of the radial inner root; it requires no radial alignment and uses no area debit.

## 7. Checker safety, proof linkage, and remaining action

The entire checker was read. Standard-library imports are `__future__`, `fractions`, `math`, `pathlib`, and `json`. AST parsing found no exec/eval/process/network/old-checker symbols. No mathematical source import was performed. The scalar interval class uses correct floor/ceil outward rounding, including signed multiplication and reciprocal; its repeated interval powers may over-enclose but remain valid. Square roots use integer floor/ceil roots. Sine degree 79 and cosine degree 78 are rational Taylor polynomials with globally valid Lagrange remainders of orders 80 and 79 respectively. Every point argument is within the declared absolute bound three.

The recorded JSON has 48 guard entries, eight continuous-slab entries, and the matching source log. The program's geometry is in the proof, not justified merely by guard names: the phase formulas, actual branch, lens support, derivative/Hessian reductions, and single-budget closure were independently checked above.

`__main__` overwrites `verification.json` beside `verify.py`. Therefore the safe reproducibility step is to copy the frozen six-file bundle to a dedicated replay directory and run `python3 verify.py` there once; keep stdout, exit code, exact invocation count, and compare the new JSON with the frozen source JSON. No old checker or geometry search is called. This static audit clears that replay but does not pretend it already happened. No unresolved defect was found in the main conditional theorem.

The root has now completed that one isolated-copy replay, as recorded above. No further replay is required by this audit. Conditional theorem acceptance is supported by the independent static derivation plus that matching local exact guard; global integration and fan transport remain external.
