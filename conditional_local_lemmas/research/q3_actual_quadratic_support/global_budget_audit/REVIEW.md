# Independent read-only audit of the quadratic actual-anchor budget interface

## Result and execution scope

**Analytic interface PASS, conditional on the local theorem and its same-coefficient global prerequisites.** The reflected radial boundary functional in the downloaded .1132 archive does not contain the marked-edge anchor slack J−A. A valid three-face theorem for the actual signed anchors can therefore replace the corresponding three independent marked conjugate allowances, retaining the compatible boundary support and its variance term once. This replacement starts from the signed area identity before A<=J is used.

This audit does not prove the new local 1/40 theorem, replay the .1132 unrestricted lower-bound chain, transfer its coefficients to .11321, or prove a new radius exclusion. Source imports, old verifier executions, new local guard executions, conjugate evaluations, optimizers, scans and graph catalogue executions: **all 0**. Only ZIP/text/JSON reads and provenance hashing were performed; no ZIP member was extracted or modified.

The current archive is attachments/[REDACTED_SESSION]/disk100_lower_bound_0p1132_20261004.zip, 2,013,535 bytes, SHA256 8194bbc46ee1e8a4c7df894e8cd55f7035609706b007b3e3b5307c8399f27806. Member references below are relative to its disk100_lower_bound_0p1132_20261004/ root. SOURCE_PINS.json records their exact hashes. The .11319 ZIP was only listed before the parent supplied this latest archive; none of its verifiers was executed.

The parent's ROOT_QUADRATIC_BUDGET_INTERFACE.md, final reviewed SHA256 9e451e083dfc5316988f4ee3d08e809ff84b5f1ae680555f5b62c8127cd7f4c4, also passes this symbolic review, including its ordinary-face/ear qualification at lines 88–91. Its assertions about the new local theorem are external premises, not results of this budget audit.

## 1. Signed identity before independent anchor maxima

Let c_i be the actual counterclockwise boundary centers, q_i the actual covered switch between c_i and c_(i+1), and 2a_i the assigned short arc of c_i from q_(i−1) to q_i, with sum a_i=pi. Define

~~~
A_i = −det(c_(i+1)−c_i,q_i−c_i)/2,
L_i = a_i+[det(q_i,c_i)+det(c_i,q_(i−1))]/2,
P   = sum det(c_i,c_(i+1))/2.
~~~

A_i is the signed boundary-edge anchor triangle. L_i is the signed center-to-two-switch triangle plus its circular segment. Direct expansion gives

~~~
2A_i = det(q_i,c_(i+1))+det(c_i,q_i)−det(c_i,c_(i+1)).
~~~

In P+sum A_i+sum L_i, the center-edge terms cancel, det(c_i,q_i) cancels det(q_i,c_i), and the remaining switch-center terms cancel by an index shift. Thus

~~~
pi = P+sum A_i+sum L_i.                              (I)
~~~

For an embedded retained center triangulation, P is the sum of its positively oriented triangle areas, even for a nonconvex boundary-center polygon. Exterior anchor and cell terms remain signed.

Coverage of q_i by both centers gives ell_i<=2R and

~~~
A_i <= J_R(ell_i),
J_R(ell)=(ell/2)*sqrt(R²−ell²/4),
tau_i=J_R(ell_i)/R²−A_i/R²>=0.                        (II)
~~~

The actual bisector switch is available in the application; the upper anchor bound itself only needs the two distance constraints and the consistent sign convention.

Traceable archive evidence: COUPLED_BOUNDARY_LOSS_PROOF.md:76–82 gives the separate signed cell outer area, :93–109 defines the center-face and edge-anchor majorants, and :167–182 allocates center faces, edge anchors and cell outer pieces separately. lower_bounds_extension_20261004/COMMON_APEX_FAN.md:7–13 proves actual boundary-face orientation; :78–109 performs signed patch cancellation. This review derives (I) directly.

## 2. Exact quadratic split and selected replacement

Use one common lambda for every boundary-edge penalty and define

~~~
Q_actual = sum_i L_i+lambda*sum_i ell_i².
T_i = area_i/R²+A_i/R²−lambda*ell_i²/R²+t_i*gamma_i.
S_i^J = area_i/R²+J_R(ell_i)/R²
        −lambda*ell_i²/R²+t_i*gamma_i.
T_i = S_i^J−tau_i.
~~~

With no ears, the exact ordinary partition is

~~~
pi = R²[sum_unmarked area_f/R²
         +sum_marked(T_i−t_i*gamma_i)]+Q_actual
   = R²[sum_unmarked area_f/R²
         +sum_marked(S_i^J−t_i*gamma_i)
         −sum_marked tau_i]+Q_actual.                (III)
~~~

The independent proof normally drops the nonnegative sum tau_i while applying the conjugates to S_i^J. A selected actual-anchor theorem restores this previously discarded term for those anchors. One must rebuild (III), rather than replace a larger score by a smaller score in an already truncated inequality without its signed justification.

For three selected ordinary common-apex faces with the identical lambda and common t,

~~~
T_0+T_1+T_2 <= 3C_lambda(t)−delta
~~~

replaces their independent 3C allowance and decreases the physical area upper bound by R² delta. Delta already contains the anchor slack paid by the new theorem; tau cannot then be deducted again. Unselected anchors may still use A<=J. Thus not every face needs a new actual-anchor theorem.

## 3. The archived Q and variance have not used J−A

Write s_i=sin a_i, h_i=sqrt(R²−s_i²), rho_i=cos a_i−h_i, and c_i=rho_i u_i+Delta_i u_i+tang_i v_i. The actual cell term is

~~~
L_i=H_R(2a_i)−s_i Delta_i,
H_R(2a_i)=a_i−s_i rho_i.
~~~

RADIAL_PENALIZED_STAR_0p11309_PROOF.md:25–33 gives the endpoint-lens inequality and signed L_i; :37–61 expands the full actual edge-square cycle and pays its displacement remainder from the cell outer loss; :96–110 states the affine boundary total and cancellation with the negative marked-base penalties.

The current radial_stability_backend.py:48–54 builds only

~~~
[a−s_a rho_a+b−s_b rho_b]/2+lambda*ell_rad(a,b)²
−MU*(a+b)+GAMMA*[(a−pi/B)²+(b−pi/B)²]/2.
~~~

Its site() at :31–36 supplies s, cos a, h and rho; CENTER is pi/B at :29. There is no marked face, apex, actual edge-anchor height, J_R, A_i or tau_i in this functional.

radial_corrected_stability_backend.py:58–60 adds half of each pair endpoint's displacement correction. The reflected override at reflected_radial_stability_backend.py:14–21 is the positive residual allowance

~~~
P_ref(a)=(sqrt(1−R²)−rho(a))*max(F_K(a),0).
~~~

RADIAL_NODE_SOFT_CORRECTION.md:14–27 proves the one combined cyclic functional and demands one correction per node. REFLECTED_CENTER_RADIAL_CORRECTION.md:34–58 strengthens the displacement cap and expressly disallows a second use of that outer-loss allocation.

Summing the pair support gives

~~~
Q_actual <= B*TARGET+2pi*MU−GAMMA*V,
V=sum_i(a_i−pi/B)².                                  (IV)
~~~

Each half-arc appears in two pairs, so the affine terms total 2MU*sum a_i=2pi*MU, the variance halves total GAMMA V, and correction halves count P_ref once per node. REFLECTED_RADIAL_B33_0p1132_COMPONENTS.md:20–30 states (IV). verify_deep_face_B33_0p1132.py:109–111 uses it verbatim.

Therefore neither Q=B*TARGET+2pi*MU nor −GAMMA V has spent the marked-edge tau=J−A in this implementation. They can be retained with a compatible actual-anchor replacement. This statement is source-specific: another boundary implementation which incorporates anchor slack needs its own allocation audit.

## 4. Three losses must not be conflated

The geometric assigned-cell lens loss is

~~~
P_lens,i=H_R(2a_i)−L_i=s_i Delta_i>=0.
~~~

The old proof already combines −sum P_lens with the increase in actual edge squares while proving (IV). Keeping its Q−GAMMA V therefore does not authorize another subtraction of P_lens. A retained lens payment needs a new joint boundary support or an explicit allocation between valid inequalities.

The backend P_ref is a nonnegative upper correction, included with a plus sign; it is not an available negative fee. Both P_lens and P_ref concern a cell-to-two-switch piece. Tau concerns a boundary-edge switch triangle. These are different terms in (I).

A theorem for the three T_i includes three center faces and three edge anchors, with no full cell L_i. Sharing a terminal boundary center is consequently harmless for this allocation. A different theorem which absorbs the entire four-cell fan patch would need explicit outer-piece weights before it could coexist with (IV).

## 5. Centered angles, BI credits and the one variance budget

The q0 face support is for area_f/R²+sum l_v(alpha_fv−pi/3). At an interior vertex with degree d, marked incidence p and label l_v, its unmarked centered-angle sum is

~~~
l_v[(6−d+p)pi/3−sum_marked gamma].
~~~

With boundary label l_B=−b, the unmarked boundary angle sum is Gamma+(e−2)pi and there are U boundary q0 corners. Moving these potential terms gives the marked tilt t_v=l_v+b and the remaining term −l_v(6−d+p)pi/3. These are exactly the signs at verify_deep_face_B33_0p1132.py:79–93; the final boundary-centered term is −b[(2−e)pi+U pi/3] at :108. Identical selected t_v preserves all complete-star cancellations. A locally convenient different tilt cannot replace it.

The q0 checks reserve terminal BI credits at :69–74, and the vertex score pays the matching two terminal incidences at :91. Those credits come from q0 face supports, not unspent marked-anchor losses. RADIAL_FACE_B33_0p11316_PROOF.md:65–86 states the incidence and centered score; DEEP_FACE_B33_0p1132_PROOF.md:60–66 preserves it. An actual-anchor replacement leaves these reservations in place once.

Variance remains governed by the applicable old branch. Current no-ear low-M rows use the one per-p3 charge at verifier :80,92; ear-low rows instead use one ear charge at :110–111. The initial high-M support's GAMMA=1/4 is reserved only for proving the improved arc domain. It is not added to the different final high-M support, whose GAMMA=0. See DEEP_FACE_B33_0p1132_PROOF.md:21–44,75–80 and the component note :54–75. The anchor improvement creates no second variance budget.

## 6. Disjoint objects and ear scope

Each selected ordinary face has one boundary edge and one interior apex. Each boundary edge owns one switch-anchor triangle. Distinct ordinary faces therefore have distinct anchors, and distinct eligible apex fans have disjoint face and anchor sets, even when they share a terminal center. Overlapping three-face subfans at one longer run cannot receive multiple fees without a disjoint selection.

Eligibility requires three actual consecutive one-base faces with a common apex, actual covered switches on the correct edges, positive orientation, and circumradii <= the same R. Three abstract marked incidences do not imply this geometry. research_appendix/SHARED_THREE_FACE_INTERFACE.md:9–31,56–70 records the distinction; :44–54 identifies the q0 BI payment and disjoint marked-face ownership.

An ear has two boundary edges and two anchors. It is treated once with both anchors and both penalties and may retain its original separate ear cap. The local replacement applies only to selected ordinary one-base faces; an ear cap is not a 3C allowance. Every boundary edge still appears once in Q_actual and its corresponding negative penalty cancels once. See RADIAL_FACE_B33_0p11316_PROOF.md:57–63,114–127 and the current total at verify_deep_face_B33_0p1132.py:102–111.

## 7. Current coefficients and remaining obligations

This archive supplies structural evidence at R=283/2500=.1132. Its final low-M support has lambda=.056,GAMMA=1.5; final high-M has lambda=.11,GAMMA=0 (the corresponding entry points :7–10). Its high-degree final high-M tilt is

~~~
t=(−36702176+35686535)/100000000
 =−1015641/100000000.
~~~

These are not the new R=.11321,lambda=.14,t=1071883/20000000. The archived marked-conjugate module declares and asserts lambda<=1/8 at exact_quadratic_marked_conjugate.py:1–4,31–34. The separately accepted fixed lambda14 conjugate validates its fixed marked support, not every old type, different tilt, ear cap or boundary certificate at .14.

Using a conditional delta=1/40 theorem in a new radius proof still requires:

1. The complete local theorem and its independent acceptance.
2. A same-radius, same-penalty support for the defined Q_actual and its actual domain. No unspecified numeric outer-Q is identified here.
3. A global type/angle row with exactly the same selected common t, and valid supports for every other marked/unmarked/ear type.
4. Proof that the required eligible fans actually occur in the actual diagram.
5. One q0 BI reservation and one variance allocation, every remaining count/topology branch, and endpoint closure.

Subject to these conditions, the parent's final symbolic interface is correct and the R²/40 decrease is valid. It is a result for the actual-anchor score, not for the independent J score. The previously reviewed 37-site local diagram blocks the stated standalone independent-score target in its local domain; it does not invalidate a distinct bound retaining actual J−A.

## Provenance

SOURCE_PINS.json records archive/member hashes, the preliminary and final root interface versions, and zero execution counts. Root-owned and coupled-owned source files were not edited. Current .1132 bytes are only structural evidence; there was no fresh full-chain replay or claimed new global endpoint.
