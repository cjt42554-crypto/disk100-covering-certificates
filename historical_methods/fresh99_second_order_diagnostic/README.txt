PUBLICATION SCOPE: CURATED HISTORICAL RESEARCH RECORD
This is a lean, privacy-edited copy. Original replay/rebuild instructions below describe the complete historical archive, not a ready-to-run public reproduction. Large run snapshots, inventories, caches and some dependency trees are intentionally omitted. No complete replay of these finite runs is claimed here; no program was executed for publication. Read ../README.md and this route's PUBLICATION_SCOPE.json. Original files were retained unchanged.

Three fresh99-binding cores: physical second-order numerical diagnosis
====================================================================

Result
------
All three new99-binding cores have positive numerical t=R^2 Lagrangian curvature
on their core critical tangent after removing only simultaneous rotation:

    source/free    active contacts   core dimension   min eigenvalue
       0/24             144               54        .0015733843036035775
       1/73             148               50        .0015684482813187349
       2/80             143               55        .0015605909277838254

All normalized contact multipliers are numerically positive; minima are
.0005023875193511698/.0001625769741998483/.0004422949411362035. There are no
zero-stress contacts in these three diagnoses. Adding the genuine free-center
coordinates gives exactly2 zero modes, rather than positive curvature. Before
the rotation quotient there are3 numerical zero modes: rotation plus the2 free
coordinates. No trustworthy negative direction was found, so no negative-step
probe or new descent proposal was generated.

This supports a bad stable99-core surrounded by a free-center feasible lens.
Merely moving the free center within that lens need not reduce the core radius.
It is useful evidence against treating nominal source minima as easy saddles.
It is not a strict local-minimum proof: there is no interval root or Hessian
enclosure. In particular the full100 configuration cannot be called strict
because it retains two physical flat directions.

Source binding and numerical root correction
--------------------------------------------
Read-only pinned sources are fresh_rattler_geometry_bridge/inputs/basin00.npy,
basin01.npy, basin02.npy, with respective file SHA:

    e053e7012512613df84d2329e549f6a9722d263eb7154116a6bd5771b79017ce
    2119c57582b6142b7650285b8329cbcf7f132b752cc79b18150a3096f9f13db5
    e482be7b4d0f6eb8d94316eee003c371d1017bea70efa9e3ba74548e7e591367

Those sources and the entire frozen80-file geometry bridge are never written.
Raw active counts are stable for squared-gap thresholds1e-13 through1e-9.
Raw stationarity residuals are ~4.92e-9/3.31e-9/1.58e-9; their artificial
smallest contact singular values ~5.12e-8/3.45e-8/1.64e-8 make a direct
threshold-based exact-rank claim unsafe. One small Newton step solves the
numerical core KKT system, with a fixed affine rotation slice based on the raw
core rotation vector. Maximum99-center correction is only
1.11948094e-7/5.65055981e-8/5.16159777e-8. Corrected stationarity norms are
1.80e-15/9.16e-16/3.23e-16. The free coordinates remain bitwise unchanged.

The correction is a local equality/stationarity root solve, not a full100
optimizer or a branch search. It removes no physical translation, dilation or
additional core degree of freedom. Excluding the free coordinates from this
root calculation is legitimate because every active branch is independent of
them; the final Hessian restores those columns as exact zero rows/columns.
No artificial anchor or penalty is used to give the rattler positive curvature.

Corrected actual Delaunay/boundary topology equals the source topology. Whole-D
numeric oracle radius differs from model sqrt(t) by at most1.67e-14. Model
radius changes only at ~1e-14 scale; these are not new upper-bound proposals.
Original and corrected coordinates and every Newton checkpoint are retained.

Correct physical curvature
--------------------------
Use smooth actual internal3-center and boundary2-center branches. At contact k,
g_k(C) is its squared radius, with objective multipliers lambda_k>=0 and sum1.
Internal hole p solves the two equal-distance equations; boundary hole q solves
one equal-distance equation and q.q=1. Their analytic IFT derivative Q=Dx/Dc
is computed by a2x2 solve, with nonzero determinants checked numerically.

Let w_ki be the circumcenter barycentric weights, or boundary equal-distance
weights from the tangent formula. Set alpha_ki=lambda_k*w_ki. For a boundary
contact, beta_k=-lambda_k*(1-(sum_i w_ki*c_ki).q_k). Then the pulled-back physical
Lagrangian Hessian has quadratic form

    H[v,v] = 2*sum_ki alpha_ki*||Dx_k[v]-v_ki||^2
             +2*sum_boundary beta_k*||Dq_k[v]||^2.

The boundary beta multipliers are negative and their curvature is retained.
The formula equals Hessian(sum_k lambda_k*g_k(C)) because the equal-distance
and unit IFT definitions plus hole stationarity cancel the second derivatives
of auxiliary hole coordinates. Details and conditions are in cone_review.txt.
The formulas use t multipliers throughout; radius-objective scaling is not mixed.

An independent central difference of the analytic branch gradients at steps
1e-5,1e-6,1e-7 checks the entire200x200 physical Hessian. Maximum projected
operator discrepancies over all three sources/steps are at most2.08e-10, over
seven million times smaller than the smallest positive core eigenvalue. Hessian
asymmetries, every spectrum and every step are saved, rather than hidden by
averaging alone. This is strong numerical consistency, not a certified error
bound on unknown exact roots.

Critical cone and what zero multipliers would change
----------------------------------------------------
At an exact KKT root, for epigraph constraints g_k(C)-t<=0 and a rotation gauge,
critical directions satisfy dt=0, gauge'=0, Dg_k[v]=0 on lambda_k>0 and
Dg_k[v]<=0 on lambda_k=0. Thus all-active equalities suffice here only if the
numerical positive multipliers can be validated at an actual root. Zero-stress
constraints cannot silently be imposed as equalities: positive all-active
curvature alone would not test the larger critical cone. cone_review.txt gives
this distinction, negative-direction sufficiency under regularity, multiplier
nonuniqueness counterexamples, and second-order-restoration requirements.

The numerical core epigraph+gauge Jacobian is full row rank, giving a consistent
unique normalized multiplier diagnosis. The active contact matrix has rank
m-1 after correction. All active nonsupport squared-clearance gaps are positive
(.007854/.006295/.006226 minima); inactive squared-radius gaps are
5.54e-6/3.22e-7/1.11e-6. Hole IFT determinants and barycentric weights are
positive away from degeneracy. These are eligibility checks, not interval claims.

Files and reproduction
----------------------
RESULT.json contains all contacts, weights, geometry eligibility, Newton steps,
analytic spectra and three-step errors. matrices/source00..02.npz stores full
analytic and difference Hessians, core bases, free flat bases, and six low core
modes per source. AUDIT.json independently checks rank, rotation and free zero
modes, source SHA, finite-difference error margins, corrected topology and actual
whole-domain scoring. CHECKPOINT.json and run.log preserve progress and failure
evidence; no failed Newton step occurred. No full100 batches or finite descent
probes were run.

In a writable copy with the same relative read-only source paths:

    python3 -B diagnose.py
    python3 -B audit.py

The parent may choose a different structural/event-forced search based on these
bad positive-curvature cores. These numerical soft modes are not promised escape
directions and the corrected roots should not be mistaken for improvements.
