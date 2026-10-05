# Radius-uniform fee with a relaxed actual-middle cap

Status: new composed paper corollary, **awaiting independent static review of both inputs and this composition**. Author mathematical executions/imports/replays: zero. The original main PROOF.md, SHA256 `692a4526d88dd3e856aa96b60961b43c1726d9ba029803199b90e8d7b0c9da2c`, remains byte-for-byte frozen.

The new geometry input is `research/q2_tau_geometry_thresholds/ANALYSIS.md`, SHA256 `97283475b6176cd0d426759d73d931690a155dec69a86843d5f5a630ff6ebbdf`. Its author status is paper proof awaiting independent review; this file does not turn an author result into accepted evidence by citation. Source and composition pins are separate in COROLLARY_PINS.json.

## 1. Precise changed premise

Retain every premise and definition of the main radius-uniform theorem: fixed lambda=7/50, t=1071883/20000000, R in the full closed interval [5661/50000,29/250], two consecutive positive nondegenerate **actual** triangles of circumradius≤R, and actual covered unit anchors on their respective bases, allowing signed negative A and unequal endpoint distances. Retain the same two actual outer norm floors |c0|²,|c2|²≥83/100.

Replace only the actual middle cap |c1|²≤1−R² by the weaker cap

**|c1|²≤1+3/2000.**

Then the same fee **2C−S2>9/400** follows on the entire radius interval. In particular, the convenient condition |c1|≤1 suffices. All centers here belong to the same scored actual fan; no norm transport, fan reselection, radial move or symmetry operation is performed.

## 2. Re-prove the sole use of the changed cap

Assume D≤9/400. The complete normalized support and localization in main Section 2 do not use any middle norm hypothesis. They give

alpha≤11661/8000<73/50, sin alpha>6/7,
P=(15/13)T≤27/1040, |g|,|d_b|≤P,
b=(0,sin alpha)+d_b, B=|b|≤1.

The latter bound is already implied by actual endpoint coverage: the two actual anchor-minus-middle vectors each have norm≤1; their half-difference b thus has norm≤1. It is also valid along the relative lens path, without a cap on its translated fan.

Use the new geometry source's explicit alternating Taylor guard

P6(73/50)−11/100=783848711/11250000000000>0.

All alpha lie in the first quadrant, and cosine decreases there, so cos alpha>11/100. With n=(b_y,−b_x)/B and k=(cos alpha,0)+g,

n_x>6/7−27/1040=6051/7280,
k·n>(11/100)(6051/7280)−27/1040
    =47661/728000.                                    (2.1)

Uniformly R≤29/250 and B≤1 imply

u=sqrt(1−R²B²)>.993,
2u(k·n)>(993/500)(47661/728000)
        =47327373/364000000>13/100.                    (2.2)

The exact two-unit-anchor translation is still c1=±un−Rk. On its negative branch,

|c1|²=u²+2Ru(k·n)+R²|k|²
      >1−R²+(13/100)R
      =1+R(13/100−R)
      ≥1+R0(13/100−29/250)
      =1+39627/25000000
      =1+.00158508>1+3/2000.                          (2.3)

This contradicts the **new** actual middle cap. Hence the actual positive branch is forced under that cap, including T=0. The claimed stronger cap is not inferred by reusing the old negative-branch inequality: (2.1)–(2.3) are the new argument.

After this positive sign is established, the cap is not used again in main PROOF.md. Sections 4–7 rely only on the original normalized angle and lens domain, B bounds, the actual positive-branch formula, and the outer floors. In particular their new radius-dependent signed linear support, eight whole-slab transfer, complete remainder, changing baseline chi and scalar norm contradiction are unaffected by the cap change. No middle cap is required along the interpolated positive-branch formula. Consequently the whole main proof closes under the stated changed premise.

## 3. Optional full-arc version on the entire radius interval

One may replace the explicit two outer norm floors by the following stronger, transparent sufficient geometric premise: each **same actual outer** R-disk covers a connected unit-circle arc of full angular length **strictly greater than 1561/10000=.1561**. This is a covered connected full arc, not a halfarc or a sum of disconnected assigned components.

For an arbitrary center of norm rho covering an arc with half-angle a, endpoint averaged squared distances give

rho≥cos a−sqrt(R²−sin²a).

No radial alignment is assumed. The lower root increases in a and decreases in R. The geometry source's paper bounds

sqrt(.83)<.91104336,
cos(.07805)<.996956,
1.83−2(.91104336)(.996956)=.01345971197568
                                  >(.116)²=.013456

prove that at a0=.07805 the lower root exceeds sqrt(.83), provided the arc is feasible. The needed branch condition cos(a0)>sqrt(.83) follows from cos(a0)≥1−a0²/2>.99>.91104336; thus squaring the positive difference cos(a0)−sqrt(.83) is legitimate. If a0 were infeasible the assumed longer covered arc could not exist. Thus every longer arc and every R≤.116 has rho²>.83. This yields both required outer floors for the exact same centers.

Combining this with the new actual-middle cap in Section 1 proves the same conditional fee on [R0,.116] with no explicit outer norm premise. The middle restriction can again be stated as |c1|≤1. This optional corollary uses the stronger full-arc premise; it does not assert that the parent's existing full-arc floor>.147 has already been strengthened.

## 4. Actual bottleneck and scope

The main score/tau theorem now extends to R=.116 at fixed lambda,t, and this composition removes the tighter middle cap. The remaining simple geometry interface is the two actual outer floors. If one uses only connected full arcs>.147 to obtain .83, the sharp one-disk radius threshold is

Rcrit=sqrt[1.83−2sqrt(.83)cos(.0735)]

and the geometry source proves .113281<Rcrit<.113282. Beyond it, the arc premise alone permits an inner radial center with squared norm below .83; that is a genuine failure of that interface, not a global-cover counterexample. Convenient strengthened full lengths .1471, .1478 and .1561 suffice respectively up to radii .1133, .1135 and .116. None is claimed sharp.

The parent still must establish an admissible actual fan, applicable actual arc floors or norm floors, and full-cover graph/profile exclusions. Neither this lemma, its relaxed cap, nor its continuous radius range proves the exact 100-disk optimum or a new global lower bound. Parameter perturbation, a different outer floor, a changed D threshold and new anchor semantics require separate proofs.
