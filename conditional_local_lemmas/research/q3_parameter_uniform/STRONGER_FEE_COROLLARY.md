# Stronger uniform rectangle corollary

This addendum uses the separately accepted actual-A fee 19/500 and its
separate radius extension. Those acceptance records must accompany the
corollary; a successful scalar calculation alone does not establish them.
The original PROOF.md and its fee16027/800000 remain unchanged.

With exactly the triangles, signed anchors, full-apex support and parameter
rectangle of PROOF.md, the following stronger conclusion holds:

    27/200<=lambda<=29/200,
    971883/20000000<=t<=1171883/20000000,
    0<R<=794400/5378209,
    S(lambda,t)<3C(lambda,t)-25473/800000.

In particular it holds for all R<=231/2000=.1155; the uniform fee is
.03184125. This is still a local actual signed-A score theorem, not an
independent-J fee or a whole-cover lower bound.

The proof from Sections 2--4 of PROOF.md applies with the hypothetical
new deficit threshold replaced by E*=25473/800000. Every localization
gap remains larger than E*: the acute endpoint gap exceeds 161/5000,
the obtuse gap exceeds 15109287/70000000, and the convex cut gap is
at least 1/10. In particular 161/5000-E*=287/800000>0. The degenerate
gap exceeds 7/5-47890713/70000000, also larger than E*.

Thus a hypothetical Dnew<=E* again yields L<theta_i<H and z_i<.1,
then the old scalar budget at (lambda0,t0) and the exact sensitivity

    Dnew >= f(Dold),
    f(D)=(199/200)D-(1/40)sqrt(3D/2).

The stronger actual-A theorem now gives Dold>19/500 throughout the
stated radius interval. The radius extension uses its strict envelopes

    K=993/1250, Sfit=5378209/1000000,
    K<2|B_perp|<=R(|B|^2+|b|^2)<R*Sfit<=K.

Hence the endpoint R=K/Sfit=794400/5378209 is included. No new boundary
constant or independent support coefficient is transferred by this step.

The same monotonicity argument for f applies on D>=19/500. Alternatively
sqrt(57/1000)>23/100 makes its secant slope exceed
199/200-(3/80)/(46/100)>0. The strict integer-square comparison

    57/1000 < (191/800)^2

therefore gives

    f(19/500) > (199/200)(19/500)-(1/40)(191/800)
                 =25473/800000=E*,

contradicting Dnew<=E*. This also excludes equality. The original
uniform closure argument is unchanged: any zero-area limiting face is
far below C-E*, so all divisions occur after localization.

Only the existing deficit estimates are used for parameter sensitivity.
No second fee is subtracted from assigned-cell lens, variance, terminal
BI credits, or the global area identity. A whole-cover proof must still
use the same lambda/t on each selected ordinary fan and establish its
actual incidence and all other rows.
