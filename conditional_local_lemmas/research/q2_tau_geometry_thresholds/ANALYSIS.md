# Actual arc-to-norm thresholds and the two-anchor translation branch

This is a new analytic geometry interface for the next radius range. Every
finite decimal below denotes an exact rational. No mathematical program,
import, old checker replay, numerical scan, optimizer, or trig evaluator was
run. All new numerical guards are paper Taylor or rational inequalities.
This note awaits independent static review. It introduces no global radius
exclusion and no claim of the exact optimum. The parent reports a separate
global .11322 certification undergoing portable publication.

The fixed-point input is the accepted Pro proof at
`R0=5661/50000`, `lambda=7/50`, `t=1071883/20000000`, with
`D=2C-S0-S1`. Its low-deficit assumption is `D<=9/400`. Its actual lens
geometry permits negative signed A and arbitrary two-endpoint covered unit
anchors. Neither radial placement nor equal endpoint distances is assumed.
The source files and hashes are recorded in `SOURCE_PINS.json`.

## 1. Arbitrary-center endpoint lemma; both norm branches retained

Let `0<R<1`. A disk centered at an arbitrary vector c covers a connected
unit-circle arc of full angular length L>0. A center c=0 cannot cover a unit
point. For rho=|c|>0, its covered unit-circle set is determined by

    cos(angle(q,c)) >= (1+rho^2-R^2)/(2rho) >0.

Thus the nonempty covered set is a minor connected arc, including its
endpoints, of half-angle at most `asin R<pi/2`. In particular any connected
covered arc has `0<L/2<=asin R`; a wrap across an angle-coordinate cut is
handled by lifting inside this same minor arc. We do not substitute a sum
of disconnected assigned components for L.

Put a=L/2 and let m be the unit midpoint direction of the arc. The two
endpoint vectors have average `cos(a)m`. Averaging their squared distances
to c gives

    1+rho^2-2 cos(a) (c dot m) <= R^2.

Since `c dot m<=rho` and `cos a>0`, necessarily

    1+rho^2-2rho cos a <= R^2.

The discriminant and both roots give the complete necessary interval

    sin a<=R,
    cos a-sqrt(R^2-sin^2 a) <= rho
                                  <= cos a+sqrt(R^2-sin^2 a).       (1)

Both roots are positive because their product is `1-R^2>0`.
The lower and upper radial branches are both permitted by the geometry;
neither branch is discarded by (1). Equality is attained by a center in
direction m at either root: such a disk covers the entire intervening
minor arc, since its endpoint distance is maximal on that arc. Thus the
lower bound is sharp even among arbitrary centers. Direction m was used
only to prove sharpness, not assumed of the given c.

Define

    rho_min(a,R)=cos a-sqrt(R^2-sin^2 a).

For `0<a<asin R`,

    partial_a rho_min = sin a[cos a/sqrt(R^2-sin^2 a)-1] >0,
    partial_R rho_min = -R/sqrt(R^2-sin^2 a) <0.                   (2)

The first strict sign follows from
`cos^2 a-(R^2-sin^2 a)=1-R^2>0`. The function is continuous at
`a=asin R`, where it equals `sqrt(1-R^2)`. Strict full-length input
`L>L0` therefore gives `rho>rho_min(L0/2,R)` whenever the fixed half-angle
is feasible. If it is infeasible the asserted covered arc cannot exist.

## 2. Exact radius and full-arc thresholds for a specified norm floor

Let `u=sqrt(v)` and `a0=L0/2`, with `0<u<cos a0` and `0<a0<pi/2`.
Then all disks covering a connected arc `L>L0` necessarily have `rho>u`
whenever

    R^2 <= 1+v-2sqrt(v) cos(L0/2).                                (3)

Indeed the right side minus `sin^2 a0` equals `(cos a0-u)^2`, so
(3) is exactly the lower-root inequality `rho_min(a0,R)>=u` whenever
the arc is feasible. Define the sharp radius interface

    Rcrit(L0,v)=sqrt(1+v-2sqrt(v)cos(L0/2)).                        (4)

If R exceeds this value and `a0<asin R`, then a radial center at the
lower root for an arc of length slightly greater than L0 still has norm
less than u. Consequently the strict arc premise alone cannot imply the
same floor beyond (4). This is a one-disk interface counterexample,
not a full covering or fan counterexample.

Conversely, in the useful regime `1-u<R<sqrt(1-v)`, the exact critical
full length is

    Lcrit(R,v)=2 acos[(1+v-R^2)/(2sqrt(v))].                        (5)

Here the acos argument lies between u and 1. A strict covered length
`L>Lcrit` implies `rho>u`, and (5) is attained on the inner radial branch.
The outer branch is retained throughout. For other parameter regimes one
must check feasibility and the lower-root branch explicitly, rather than
use (5) without its conditions.

## 3. New paper guards for floor v=.83 and fixed full length .147

The following square identities prove a narrow algebraic root bracket:

    (.91104335)^2=.8299999855792225 <.83,
    (.91104336)^2=.8300000038000896 >.83.                          (6)

Put `a0=.0735`. Its alternating cosine polynomials yield

    .99730009 < cos a0 < .9973001.                               (7)

For transparency, `a0^2=.00540225`,
`a0^4=.0000291843050625`, and

    P4(a0)=1-a0^2/2+a0^4/24=.9973000910127109375.

Moreover `a0^6<.075^6` and `.075^6/720<.00000000025`, so
`P6(a0)>.9973000907627109375>.99730009`.
The usual alternating remainder has the required signs on this interval.
Combining (6) and (7) gives

    .012832731935328 < Rcrit(.147,.83)^2
                                           < .012832770102197. (8)

These are finite rational multiplications:

    1.83-2(.91104336)(.9973001)=.012832731935328,
    1.83-2(.91104335)(.99730009)=.012832770102197.

Since

    (.113281)^2=.012832584961,
    (.113282)^2=.012832811524,

we obtain the certified strict bracket

    .113281 < Rcrit(.147,.83) < .113282.                          (9)

In particular the exact necessary outer norm-square floor .83 follows
from full arcs>.147 at R=.11323, .11325, and throughout `R<=.11328`
when the arcs are feasible. At R=.1133 it does not follow from this arc
premise alone. This supersedes an earlier informal .113277 size estimate;
only the exact expression (4) and paper bracket (9) are results.

## 4. Explicit strengthened arc interfaces

For R=.1133, full arcs>.1471 suffice for the same floor .83. Take
`a=.07355`, for which `a^2=.0054096025` and `a^4<.000029304`.
Thus

    cos a <=P4(a)<1-.00270480125+.000001221=.99729641975
                                                           <.99729642.
    1.83-2(.91104336)(.99729642)=.0128394372144576
                                         >(.1133)^2=.01283689. (10)

For R=.1135, full arcs>.1478 suffice. Take `a=.0739`, with
`a^2=.00546121` and `a^4<.000029904`. Then

    cos a <1-.002730605+.000001246=.997270641<.99727065,
    1.83-2(.91104336)(.99727065)=.012886392389232
                                         >(.1135)^2=.01288225. (11)

For R=.116, full arcs>.1561 suffice. With `a=.07805`,
`a^2=.0060918025`, `a^4<.0000372`, we have

    cos a <1-.00304590125+.00000155=.99695564875<.996956,
    1.83-2(.91104336)(.996956)=.01345971197568
                                             >(.116)^2=.013456. (12)

This final threshold is close to the sharp one in a certified coarse
sense: full length .156 is insufficient for the sole arc-to-.83 interface
at R=.116. At a=.078,

    P4(a)=.996959542294,
    a^6/720<.08^6/720<.0000000004,
    cos a>=P6(a)>.99695954,
    1.83-2(.91104335)(.99695954)=.013453281727882<.013456.

Thus ` .156<Lcrit(.116,.83)<.1561`. All sufficient lengths in this
section are exact convenient rationals; none is called the sharp value.

## 5. Weaker outer floor for the local .0223 radius route

Full arcs>.147 also imply `rho>.91075` throughout `0<R<=.1135`
when feasible. This is a new guard, not a rerun of the old R0 guard.
At `a0=.0735`, the cubic sine and elementary cosine bounds give

    sin a0 >= a0-a0^3/6=.0734338224375>.07343,
    cos a0 >=1-a0^2/2=.997298875.

The exact finite rational guard is

    (.1135)^2-(.07343)^2=.01288225-.0053919649=.0074902851,
    (.086548875)^2=.007490707763765625>.0074902851.

Therefore the square-root term in (1) is below .086548875 and

    rho>.997298875-.086548875=.91075,
    rho^2>(.91075)^2=.8294655625.                                 (13)

This floor is weaker than .83 and cannot be substituted for .83 in the
Pro .0225 theorem without re-proving its final norm contradiction. It
is the appropriate proposed outer-norm interface for the separate local
.0223 radius proof. Its local bound and all radius-dependent payment
constants remain separate obligations.

## 6. Positive translation branch: middle cap can be relaxed

This section uses the Pro proof's own low-D dimensionless localization,
not the local proof's opposite skew-sign convention. For `D<=9/400`,
its full lens error coordinates satisfy

    alpha<=11661/8000<73/50,
    sin alpha>6/7,
    |g|,|d_b|<=P=(15/13)T<=27/1040,
    b=(0,sin alpha)+d_b, B=|b|<=1,
    n=(b_y,-b_x)/B, u=sqrt(1-R^2B^2),
    k=(cos alpha,0)+g,
    c1=+-u n-Rk.                                                (14)

The statement `B<=1` uses actual endpoint coverage and lens convexity,
not any middle norm cap: the squared-norm average of the two actual
anchor-minus-middle offsets equals `|k|^2+B^2<=1`. The two unit-anchor
equations then give both signs in (14).

A new sixth-order cosine paper guard is

    P6(73/50)-11/100=783848711/11250000000000 >0.

Hence `cos alpha>11/100`. Because `B<=1` and
`b_y>6/7-P=6051/7280>0`, we get

    n_x>6051/7280,
    k dot n > (11/100)(6051/7280)-27/1040
             =47661/728000 >.065.                              (15)

For every `R<=.116`, `R^2<=.013456<1-(.993)^2=.013951`, so
`u>.993`. Consequently

    2u(k dot n) > (993/500)(47661/728000)
                  =47327373/364000000 >13/100>.116>=R.          (16)

On the negative branch, the exact norm identity therefore gives

    |c1|^2=u^2+2Ru(k dot n)+R^2|k|^2
           >1-R^2+R(13/100)>1.                                 (17)

Thus, within the same low-D lens and angle domain, the negative branch
is already excluded by the weaker actual-middle cap `|c1|<=1`, for
every `0<R<=.116`. No radial relocation, reflection, new fan selection,
or transport of an assigned arc was used.

There is even an explicit extra margin on `R0<=R<=.116`. Dropping the
nonnegative final term in the negative-branch norm yields

    |c1|^2>1+R(.13-R)
           >=1+R0(.13-.116)=1+.00158508>1+3/2000.                (18)

Hence an actual-middle squared cap `|c1|^2<=1+3/2000` also excludes
the negative branch throughout that interval. A convenient weaker
corollary is `|c1|<=1`.

At fixed R0 this is a proposed strengthening of the Pro .0225 theorem:
its middle norm hypothesis can be replaced by `|c1|^2<=1+3/2000`, while
its two actual outer hypotheses remain `|c0|^2,|c2|^2>=.83`, and all
other original score, apex, lens, and support premises remain unchanged.
The reason is that the middle cap is used in that proof only to select
the positive sign in Sections 3 and 5; (14)-(18) now provide that
selection including T=0. All later Pro formulas and 48 accepted fixed-R0
guards are unchanged and need not be rerun. Independent review must
confirm this premise-use claim before treating the strengthening as
accepted.

For R>R0 this section supplies only the branch-selection lemma. A full
fee theorem still needs the new radius-dependent directional payment,
Taylor remainder, chi envelope, and final outer-norm contradiction. The
parent's separate radius analysis is responsible for those obligations.

The earlier coarser branch derivation using `cos alpha>.1` is preserved
as `history/COARSE_BRANCH_INTERFACE.md`. No old accepted source is edited.

## 7. Bottlenecks and scope map

* A fixed connected full arc>.147 supports squared norm .83 only to
  the sharp radius interface (4), bracketed in (9). This is a real
  one-disk obstruction to that interface at .1133, regardless of how
  much the actual-tau payment is improved.
* The weaker floor .8294655625 retains full arc>.147 through .1135.
  It must be paired with a separately established mean-norm upper
  bound; a theorem with outer floor .83 does not cover it by wording.
* The middle cap `sqrt(1-R^2)` is not forced by this low-deficit branch
  argument. The weaker cap at 1, or the explicit squared cap in (18),
  suffices in the Pro dimensionless localization up to .116.
* New lambda/t rectangles alter the complete support and the low-D
  localization, so the angle and error bounds in (14) would need a new
  proof. This note does not enlarge lambda, t, D, or anchor semantics.
* None of these local interfaces supplies the parent's full-cover
  arc lower bound, admissible actual fan, global graph/profile
  exclusion, or exact optimum. Full-cover integration remains external.

Mathematical execution/import/replay/scan count for this note: 0.
Source reading, text writing, byte hashing, and JSON metadata operations
are recorded separately in `EXECUTION_LEDGER.json`.
