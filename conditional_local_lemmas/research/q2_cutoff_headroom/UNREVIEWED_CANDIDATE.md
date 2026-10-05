# Cutoff .0228 headroom: stopped feasibility notes

Status: **UNREVIEWED CANDIDATE, NOT A THEOREM OR CERTIFICATE**.
Research stopped immediately on the parent's closeout instruction.
This note preserves only derivations already made before that stop.
No further closure, mathematical execution or independent review is
claimed. All accepted/frozen predecessor sources remain unchanged.

The candidate concerns central lambda=7/50,t=1071883/20000000 and
the Pro actual signed-A two-face norm domain: actual consecutive
positive triangles, circumradii<=R, covered actual unit anchors,
the SAME actual middle norm cap, and both outer squared norms>=.83.
The proposed enlarged cutoff is D<=57/2500=.0228. Radius transfer
to [.11322,.116] below is only a prospective composition with the
separately authored radius paper; this note accepts no new radius
or cutoff theorem and no global covering bound.

## Domain comparison already derived

The central scalar deficit budget remains available up to .025,
as independently restated in PRO_PARAMETER_COROLLARY.md. Under
D<=.0228 it gives the prospective enclosing ranges

    .72<theta_i<.94,
    .7525<s<.9045,
    |Delta|<.151,
    0<=Z<=.0114, 0<=T<=.0228.

Here sqrt(D/2)<.108 and sqrt(D)/2<.0755 justify the individual
theta and s bounds from the frozen peak bracket (.828,.829).
The old whole-slab box was

    .753<=s<=.904,
    |Delta|<=.15,
    omega=Z/2+eta in [-.096,.101625].

Thus old slabs do NOT directly cover the new domain. The new
prospective omega upper bound is .1017, provided |eta|<.096 is
proved as described next.

The paper comparison for the sigma coefficient was: its old exact
guard has a margin greater than .03 above 5/4. Increasing the z
endpoint from .0225 to .0228 costs at most 2(.0003)=.0006, since
the absolute z derivative of 2sin^2(theta)/sin(theta+z) is at most
2 in the relevant acute domain. This would retain coefficient 5/4.
The old sin^2(.14) guard also has enough margin to retain |d_i|<.14.
Then the same c_eta=223532401/90000000 bound yields |eta|<.096,
because the already-derived exact comparison is

    c_eta(.096)^2-.0228
       =126165744/1406250000000>0.

These are new changed-domain claims awaiting independent checking;
none is represented as a replay of an old guard.

## Proposed whole-box transfer of the old slab residual

Clamp (s,omega,Delta) to the OLD full slab box. The displacement
enclosures already derived are

    |ds|<=.0005, |domega|<=.000075, |dDelta|<=.001,
    |dalpha|<=.000825, |dh|<=.000325,
    |dtheta_i|<=.001.

The exact support residual at k0=.205 is

    G_k0=k0 sin^2(theta_i)-N,
    N=R rho sin(alpha)[1-R cos(h)cos(Delta)/u]
                  +R rho cos(h)sin(Delta),
    u=sqrt(1-R^2sin^2(alpha)), rho=u-Rcos(alpha),

with Delta reversed for the other face. A candidate analytic
Lipschitz estimate on the enlarged angular box at R<=.116 is

    |N_alpha|<.09, |N_h|<.02, |N_Delta|<.12.

These estimates were derived on paper from alpha>1, cos(alpha)<.55,
h<.55, |Delta|<.151, u>.99, rho<=1, rho_alpha<=R and R^2<.014.
They have NOT been independently audited or numerically executed.
If valid, the proposed angular residual debit is at most

    .205(.001)+.09(.000825)+.02(.000325)+.12(.001)
        =.00040575.

The separate radius paper's old-box k=.213 residual is greater
than 1559/3000000 throughout R in [.11322,.116]. Subtracting
.00040575 would leave a positive residual. This is a feasible
paper transfer route for the linear coefficient; it is not an
accepted extension of any old slab guard.

## Additional provisional headroom and unresolved closure

With kappa<=.213 and the unchanged conditional 2/3 remainder,
the prospective total tau coefficient is

    .213+(2/3)(.0228)=.2282<.23.

The prospective final angle endpoint replaces .6285 by a value
below .522+.71(.151)=.62921, an increase at most .00071.
The old final norm guard's .0003 margin appears sufficient to
pay a coarse cosine Lipschitz debit 2(.116)(.00071)=.00016472.
The prior scalar derivative estimate also leaves room for the new
sqrt-budget endpoint .151. Neither final margin establishes the
missing support/domain/branch obligations.

The principal exact outstanding issue identified BEFORE stopping:
the old **signed tangent-gradient** sufficient inequalities used
|Delta|<=.15, including coefficients 7/40 and 23/20. They cannot
be copied unchanged when |Delta|<.151. A changed-domain proof of
those signs remains required, in addition to a full independent
audit of the sigma/eta extension, angular Lipschitz bounds, actual
branch, all homotopy denominators, remainder, scalar monotonicity,
and single-budget final closure. This note does not supply that
complete audit or assert D>.0228.

A .023 comparison was only begun before the stop instruction;
it is intentionally not promoted into a theorem or validation plan.
No additional computation or research will be started by this task.

Execution record: mathematical programs/imports/checkers/replays/
optimizers/scans=0. Only source reads, new documentation writes
and SHA256 metadata hashing. Old 48 guards were not rerun.
