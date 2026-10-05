# Actual two-face fee: interface for a global covering proof

Author: `/root/q2_radius_uniform`. Status: a scope extraction of the new mean/arc paper, awaiting its independent review; **not a global covering certificate**. This note makes no mathematical execution/import/replay. Its mathematical source is `research/q2_tau_mean_arc_envelope/PROOF.md`, SHA256 `570696466a5e90d0d06e1738e41a4a94571ec7a8c45a66559fea27f2bd3cd663`. The independently audited continuous-radius and relaxed-cap inputs are separately pinned by that source. No frozen source is changed.

## Reusable continuous statements

Put R0=5661/50000, lambda0=7/50, t0=1071883/20000000. At each **actual** R in the continuous closed interval [R0,29/250], use two consecutive consistently positive nondegenerate actual triangles (p,c0,c1),(p,c1,c2), each of circumradius≤R. Its actual unit anchor qi must be covered by both endpoints of its own base. Keep signed actual A, including negative A, and the complete matched single-face support C. Define D0=2C(lambda0,t0)−S0−S1 using that exact fan and its anchors. The same actual middle center must satisfy |c1|²≤1+3/2000; |c1|≤1 is sufficient.

Under these common premises, the mean paper supplies the quantitative low-deficit implication

D0≤9/400 ⇒ M=(|c0|²+|c2|²)/2 < H_R(D0),
H_R(d)=1+R²−2R cos[b_R+(71/100)sqrt(d)],
b_R=261/500+(100/99)(R−R0), 0≤d≤9/400.

Thus actual M≥H_R(d) suffices to conclude D0>d. It permits unequal outer norms. If the two **different actual outer disks** cover connected full unit-circle arcs of lengths L0,L2, the convex lower-root norm function gives

M≥[rho_R(L0/2)²+rho_R(L2/2)²]/2
 ≥rho_R((L0+L2)/4)²,
rho_R(a)=cos a−sqrt(R²−sin²a).

Consequently the closed arc-sum condition

L0+L2≥4 atan[R sin beta/(1−R cos beta)],
beta=b_R+(71/100)sqrt(d),

also implies D0>d. The acute atan branch and the feasible connected-arc domain are proved in the source. Closed thresholds are legitimate because its mean upper bound is strict, including d=0. Multiple disconnected components of one disk must not be summed into one L.

At d=9/400 the source gives useful rational sufficient interfaces:

| Radius interval | Actual two-outer connected full-arc sum | Central fee |
|---|---|---|
| [R0,.11364] | L0+L2≥.294 | D0>.0225 |
| [R0,.116] | L0+L2≥.3018 | D0>.0225 |

In particular both original outer arcs>.147 satisfy the first line. The weaker requirement is their sum; neither individual arc must reach a common threshold.

Once central D0>.0225 is established in that physical fan, the same objective may be transported to lambda,t without altering its geometry. For a=|lambda−lambda0|≤.0003 and b=|t−t0|≤.0002, the source's radius-independent central bootstrap and score subtraction yield

D(lambda,t)>.0225−.6a−.1725b.

The support must be the exact matched C(lambda,t), including both apex branches. The full closed rectangle has uniform fee>.0222855; the balanced rectangle a,b≤.0002 has fee>.0223455. These are continuous R-and-parameter statements after the central geometric admission is supplied. They do not infer new-parameter angular localization or a different actual fan. A changed global required deficit, dual/conjugate budget, and profile need independent recomputation.

## What the global proof must supply

All objects used in a score and in its arc/norm comparison must belong to one admitted covering configuration. The global argument must establish:

1. **An actual fan.** The two faces share the actual pc1 edge and common apex p, have the asserted cyclic orientation and nondegeneracy, and have actual circumradii≤R. A relaxed independent pair of faces or a recomputed fan with old scores is inadmissible.
2. **Actual anchors and arcs.** Each qi is a unit point covered by both actual base endpoints; signed A is computed from it. The two connected full arcs used in the sum are covered by the same actual c0,c2 appearing in the scored fan. They may be unequal and need no radial alignment. A HALFarc, an x_i=base/(2R), or a union of disconnected assigned pieces is not a full arc length.
3. **Same actual middle cap.** The cap is a condition on c1 in that fan. Projecting centers into the closed unit disk can establish |c|≤1 for a normalized covering configuration, but fan selection, geometric adjacency and structural arc ownership must then be derived in that configuration. See the separate projection lemma.
4. **A disjoint global charge.** This local conclusion bounds the combined deficit of exactly its two scored faces. If the global proof sums such fees, each marked face and its actual A term must enter at most once, or an explicit weighted allocation must have total weight≤1 on every face. Overlapping stars/pairs cannot each debit the same face at full fee.
5. **One slack ledger.** Actual tau, angular loss and skew have already funded the local mean envelope once. This theorem consumes no middle-halfarc variance, Q, BI or independent J debit. Those quantities may be used elsewhere only according to the parent’s own disjoint global accounting; tau/J−A cannot be credited again merely under another name. Any optional variance coefficient must use the actual middle HALFarc and be counted once.
6. **Complete global admission and exclusion.** The full disk cover, boundary partition/profile, graph incidence, tail/topology cases, target bound and actual transfer from configuration to these marked faces are external. Local fee statements do not exclude a global radius by themselves.

These are requirements of mathematical composition; the note does not prescribe a particular global triangulation algorithm or claim that a reported highM branch supplies them already.

## The .83 floor is a conservative interface, not the fee bottleneck

The earlier condition |c0|²,|c2|²≥.83 allowed a simple uniform contradiction, but it discarded the correlation between R, the actual deficit D and the mean norm. Full arc>.147 ceases to force that individual floor above its one-disk critical radius in (.113281,.113282). That is a sharp limit of the **arc-to-.83 implication**, not a counterexample to the actual-tau fee.

The new source compares the actual arc-induced mean directly with H_R(d). It consequently retains the full original total arc .294 through R=.11364, avoiding the unnecessary individual floor. The useful global bottleneck is therefore the admissible actual fan and its connected arc sum, together with once-only marking and the remaining global deficit requirement. Raising the cutoff d above .0225, enlarging the radius range beyond .116, or changing anchor/arc semantics still requires fresh support, localization and payment proofs.
