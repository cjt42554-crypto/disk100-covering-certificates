# Projection of covering centers into the closed unit disk

Author: `/root/q2_radius_uniform`, with a separate paper identity/coverage check by `/root/q2_radius_uniform/slab_radius_derivative`. Status: new elementary paper lemma awaiting independent static review. Mathematical executions/imports/checker replays: zero. This lemma proves a covering-feasibility normalization and exact limitations; it does not transport a marked fan or prove the optimum value.

## 1. Exact pointwise inequality and covering normalization

For any center c define

Pi(c)=c if |c|≤1, and Pi(c)=c/|c| if |c|>1.

The zero center uses the first case and involves no division. If c=rho e with rho=|c|>1 and |e|=1, then for every target point |z|≤1,

|z−c|²−|z−Pi(c)|²
 =(rho−1)(rho+1−2 z·e)
 ≥(rho−1)²>0.                                      (1)

Indeed z·e≤|z|≤1. For inside/boundary centers Pi is the identity and distances are unchanged.

Thus every point of the closed unit disk covered by an original closed radius-R disk with label i remains covered by the same-radius disk centered at Pi(ci). With an external center and a nonempty old covered set, those old covered points even satisfy strict distance<R after projection. The same statement holds for unequal radii R_i and any target subset K of the closed unit disk. In particular any covering of the entire closed unit disk by at most 100 radius-R disks gives a covering by at most 100 such disks whose centers all lie in the closed unit disk.

The labels may be retained, including coincident projected centers; projection does not establish distinctness or general position. If a formulation uses exactly 100 disks, unchanged labels/duplicates or unused disks retain that count. This reduction does not justify deleting labels in a global graph/count argument without checking its counting convention.

Every old assigned connected boundary arc, **as a fixed point set**, remains covered by its old label after projection. One may keep the old point-set assignment as a legal covering assignment if the subsequent proof permits arbitrary covered assignments. This is stronger and more precise than saying all old arcs disappear, but is weaker than saying a canonical assignment or fan is unchanged.

## 2. Coverage survives; structural assertions need new proofs

Different labels receive different distance reductions in (1), so projection need not preserve nearest-owner order, Voronoi boundaries, equidistant switches, tangencies or a canonical assigned-arc partition. Old covered anchor inequalities survive if the same points and endpoint labels are retained; their equal-distance or touching properties need not survive. A prescribed arc may remain covered while losing its former ownership or adjacency role.

Nor does (1) preserve the original Delaunay fan, consistent orientation, nondegeneracy, circumradius bound, marked incidence or scores. Areas, signed A, x, gamma, sigma, tau and the deficit must be recalculated from the projected physical geometry before a local lemma is applied. Pairwise contraction of a projection, by itself, does not imply a triangle's circumradius bound survives.

An explicit rational example shows the nondegeneracy issue even with a small original circumradius. Let

a=(49/50,1/50), b=(99/100,1/100), c=(51/50,0).

The first two centers lie strictly inside the unit disk, while c projects to (1,0). The original triangle is positively oriented with area 1/10000. Its side lengths are h/sqrt(2), h sqrt(5), h sqrt(5/2), h=1/50, so its circumradius is 5h/2=1/20. After projection a,b,(1,0) are collinear on y=1−x. Thus a perfectly finite circumradius and positive orientation cannot simply be carried over. This example concerns a geometric premise; it is not a full 100-disk-cover counterexample.

For a global argument the safe sequence is: normalize a covering configuration, establish the applicable arc assignment and actual fan for that normalized configuration, then compute and allocate its actual scores. The mean/arc lemma's middle restriction |c1|²≤1+3/2000 is automatic for every normalized actual center, because |c1|≤1. This use eliminates a norm-cap transport step; it does not eliminate the fan, connected arc, anchor and accounting obligations.

The inequality is a target-domain statement. It need not hold for points outside the closed unit disk, and it does not preserve coverage of the entirety of each old radius-R disk in the ambient plane.

## 3. Compactness removes an attainment concern, without finding the value

For a labeled N-center tuple C define its covering radius on the unit disk by

F(C)=max_{|z|≤1} min_i |z−ci|.

The inner minimum and outer maximum exist and are continuous. For any two labeled tuples C,C', the triangle inequality gives

|F(C)−F(C')|≤max_i |ci−ci'|.

Projection ensures that the infimum over all arbitrary centers equals the infimum over tuples in the compact product (closed unit disk)^N. The displayed Lipschitz continuity then implies that this infimum is attained. This establishes existence of an optimal covering configuration with all centers in the closed unit disk, allowing repeats and degeneracies. It supplies neither the numerical optimal radius nor a proof that a chosen triangulation/marked structure describes every minimizer.
