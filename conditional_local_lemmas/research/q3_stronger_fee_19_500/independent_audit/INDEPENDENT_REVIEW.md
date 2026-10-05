# Independent full-source analytic review and one fresh guard replay

Verdict: **PASS for the fixed-radius local actual signed-anchor theorem**

    sum_i(area_i/R^2 + A_i/R^2 -4 lambda x_i^2 +t gamma_i)
                         <3 C_lambda(t)-19/500,

at exactly `R=11321/100000`, `lambda=7/50`,
`t=1071883/20000000`, on the positive-area consistently oriented
consecutive common-apex three-face domain stated in the frozen proof,
and its stated limiting closure. I found no mathematical blocker in
this local claim. The complete newly generated certificate matches
the author's JSON object and serialized bytes exactly.

This conclusion combines analytic review with the new guard execution;
it is not inferred from a numerical PASS alone. The pinned previously
accepted scalar-conjugate certificate remains an explicit arithmetic
premise: I reviewed its entire analytic proof and checker source, but
did not regenerate its root/value enclosures. SHA identity establishes
which accepted dependency is used, not an independent proof of those
saved enclosures. Neither this review nor the new guard run establishes
an independent-J fee, the surrounding global reduction, or a global
radius exclusion.

## Sources and execution boundary

I read all of the new `PROOF.md`, `GUARD_ADDENDUM.md`,
`verify_actual_quadratic_19_500.py`, `certificate.json`, author manifest,
README, attempt ledger and author receipt. I also read the complete
bundled `coupled_edges/verify_support.py` and the full proof, checker
source and saved certificate of `quadratic_lambda14_conjugate`.
For the unchanged-mechanism comparison I read the read-only original
actual-quadratic proof and source review, its degenerate closure
addendum, and inspected the full original-to-new proof diff.

The original-to-new comparison preserves the full-apex reduction,
actual signed-anchor geometry, lifted phases, whole-box Jacobian,
unit-anchor identity and shared-budget mechanism. The substantive
changes are the `.038` deficit, its induced `.998` paid apex bound,
`.86` skew coefficient, `.34/.30` fit errors, `.97` projection,
`.18` final margin guard, and explicit degenerate-closure incorporation.
The dependency path and fixed-radius scope wording were also updated.
I evaluated the changed guards with their changed budget, rather than
transferring the old `.025` acceptance.

The independent execution was exactly

    python3 -I -B independent_audit/run_once.py

using Python 3.12.3, with assertions enabled. The runner imported the
new verifier as a module and called its `verify()` exactly once, thereby
avoiding its author-directory `__main__` certificate write. It imported
the pinned interval primitive definitions; neither historical verifier
function was called. Before and after, it hashed every author-tree file
outside this audit directory and every historical source in the author's
read-only source-pin map. All bytes were unchanged. Every frozen author
manifest entry and every historical pin matched. All new writes belong
to this `independent_audit` directory.

| Frozen source | SHA-256 |
| --- | --- |
| New `PROOF.md` | `5a58cc96f916018885d39888bcd07ca77e768d3c719a695e2d619311fed1610f` |
| `GUARD_ADDENDUM.md` | `e2ce31e43a8ed34b0e8ef006287881b9fc0fa94ad0706dd8f8df344407fca5f7` |
| New verifier | `bd042f4976fbeeb8440da66864814371f82e4fd0e380844cd153dcb7e600f46d` |
| Full new certificate | `b79cc05457486de14e8ed89a8d62bacf596581b152092a3951f9893b3dfa56bb` |
| Interval primitive | `90fd6e3426d8c113bacedd55f9a13f5e47981e9a48d21bfef48183ec4602efec` |
| Accepted conjugate proof | `aaaad3b1629ba98fd44907735b75ae111fecaa45c8d77f078f4e28acf4713ea3` |
| Accepted conjugate checker source | `e95db0a6ead5ace03884ab54743b6238f8fe4cd242b1563617b789fce9cb57f4` |
| Accepted conjugate certificate | `f325d42bfe651599881b7dc652567729a1586389ee385c38060840bc41696bc2` |

The author manifest snapshot SHA is
`b1be5f8251bf97900e481f6c5b4d5fc0ff8e7004697a1b6f8b53f5cfbe5215c6`.
The machine receipt records `new_guard_verify_calls=1`,
`old_verify_calls=0`, `optimizers=0`, `geometry_searches=0`,
`whole_JSONMATCH=true` and `serialized_certificate_bytes_match=true`.
The matching replay certificate has the identical SHA listed above.

## Full one-face domain and scalar premise

The sine rule gives `sin gamma>=x`, hence
`theta<=gamma<=pi-theta`. Writing the two base angles as
`(pi-gamma)/2 +/- d`, their product gives exactly

    area/R^2=x^2 cot(gamma/2)-2x^2 sin(d)^2/sin gamma.

The two radius-R disk inequalities imply
`s^2+h^2+2x|s|<=1-x^2` in outward-height base coordinates.
The determinant convention gives `A/R^2=xh`. Thus `sigma>=0`,
`tau=x(cos theta-h)>=0`, and the actual one-face score is exactly
`F-sigma-tau`, including negative signed anchors.

The displayed positive `F_gamma_gamma` yields a complete convex
endpoint reduction on the feasible apex interval. It makes each
actual face score at most `C` before any local restrictions. The
bundled scalar proof correctly uses strict concavity on `[0,pi/3]`
and decreasing tails on `[pi/3,pi/2]`, rather than assuming concavity
on the whole acute-endpoint domain. I checked its tail derivative
estimates and both second-derivative estimates. Its root-isolation
implementation uses those monotonicity facts, rigorous Taylor
intervals and 1-Lipschitz trigonometric interval expansion. Its saved
certificate, with these exact coefficients, is the accepted premise
used for `C`, `Co` and `theta*` in the new calculation.

The interval backend rounds lower endpoints down and upper endpoints
up using exact fractions. Its multiplication uses all endpoint
products; its division excludes zero; the Machin pi enclosure uses
alternating arctangent remainders. The sine and cosine Taylor errors
are legitimate degree-60 and degree-59 Lagrange bounds respectively,
with the intervening odd/even coefficient zero. New guard calls use
only rational point arguments or independently justified monotone
enclosures. The backend's unrelated old `verify()` is never invoked.

## Complete localization and single budget

If the conclusion fails, every individual deficit is nonnegative and
at most `19/500`. Since each actual score is greater than `Co`, the
convex interpolation forces its acute endpoint `G(theta)` above that
score. Unique-peak monotonicity and the two strictly positive new
endpoint localization margins imply `7/10<theta<24/25` for each face.
Sine concavity and cosine monotonicity on `[2L,2H]` establish
`-G''>4` throughout this window. Integrating from the stationary peak
therefore pays `2(theta-theta*)^2`, including the factor one-half.

For `u=(theta+.1)/2`, the derivative numerator of
`sin theta/sin u` has derivative `-(3/4)sin theta sin u<0`.
Its only possible interior extremum is consequently a maximum.
Both endpoint guards exceed `81/50`. The ensuing bound
`-F_gamma>25172117/20000000>1` is valid all the way to the cut.
The cut lies inside every feasible interval. Convexity then bounds
the entire remaining apex interval by its cut and obtuse-endpoint
values, both below the required high-score threshold. This proves
`z<.1` on the full original domain, not merely near an optimizer.

Combining these integrated bounds with the exact nonnegative losses
gives one simultaneous payment

    2 sum e_i^2+sum z_i+sum sigma_i+sum tau_i <=D<=19/500.

It implies `z_i<=D0` and `gamma_i<H+D0=499/500<pi/2`
before the later acute estimates. No variable is paid twice.

## Actual signed tips, orientation and skew

The lens constraint gives
`|s|<=cos(theta) Delta/x`, where `Delta=cos(theta)-h>=0`.
Thus `|q-a|/R<=Delta/x=tau/x^2` for either sign of `h`.
This is a bound on the actual signed-anchor slack, not an
absolute-anchor substitution.

For the skew convention in the proof, the endpoint vectors relative
to the bisector have lengths `2r cos(gamma/2+d)` and
`2r cos(gamma/2-d)` on directions `m-gamma/2` and `m+gamma/2`.
Their midpoint has exactly the real and imaginary terms stated in
Section 3. The directed-base tangent is `i exp(i(m-d))`; positive fan
orientation places the apex to its left, so the outward normal is
`exp(i(m-d))`. Adding `cos theta` times that normal gives the claimed
`U cos d+i V sin d`, including the sign in `V`.
The consecutive shared rays yield the stated lifted bisector phases.

On the paid acute domain, `1-r<=cot(theta)z`,
`|U-Q|<=(11/4)z`, and the exact V decomposition proves

    -686839/1250000 <=V<=-27/250,

with its lower endpoint strictly above `-.55`. Because `cos d>0`,
`1-cos d<=sin(d)^2`; and the skew formula gives
`sin(d)^2<=625 sigma/512`. The real and imaginary skew corrections
are orthogonal, so the squared bound is

    3081055 sigma/4194304 < (43/50)^2 sigma.

The strict squared guard margin is
`13157649/2621440000`. This is the updated `.86` coefficient
computed at `.038`, not the old `.79` coefficient.

## Whole-box Jacobian and shared error allocation

Direct differentiation gives all endpoint and middle columns in the
proof. The endpoint imaginary component is the difference
`sin(theta_j)sin(a)-Q_j cos(a)/4` (or its sign reversal), so both
its positive and negative ranges must be bounded. The proof uses
both, and the checker confirms their positive upper maximum.
The middle real B column is positive over the full box, while the
other middle columns are differences/sums of positive product
intervals. These facts give exactly

    JB2=76960399529/20000000000,
    Jb2=77583868411/20000000000.

The straight peak-to-actual segment remains in the box, making
integrated Frobenius/operator bounds legitimate everywhere on it.
This does not rely on a Hessian only at the ideal triple.

The common apex translation cancels from B and b. Aggregating the
amplitude and lifted phase errors gives safe `LB=201/50` and
`Lb=529/200`. More explicitly the z coefficients are
`(201/100,201/50,201/100)` for B and
`(201/100,127/100,201/100)` for b. The displayed Lb is conservative
but valid. The anchor coefficients `ct` and `ct/2` are smaller.
The skew weights have squared sums `3/2` and `1/2`. Cauchy combines
the angle spending `2 sum e^2` and skew spending `sum sigma` into

    alphaB2=121336399529/40000000000,
    alphab2=92375868411/40000000000.

For `W=sum z+sum tau`, the same budget therefore gives (15).
The strict `alpha^2-4Llinear^2D0` guards make each continuous
error function decrease with W, so its maximum is at zero linear
spending. The `D=0` endpoint is immediate without differentiation.
The strict squared fit margins are

    6608408949/20000000000000,
    44858500191/20000000000000,

establishing `.34` and `.30` from the same deficit. These are
simultaneous consequences, not two additional budget deductions.

## Unit compatibility, projection and limiting degeneracy

For three unit points, endpoint mean a and physical half-chord h
satisfy `a dot h=0` and `|a|^2+|h|^2=1`. In the perpendicular
direction the displayed scalar identity proves
`|alpha-beta|<=1-alpha beta`, and
`|a-q_1|^2+|h|^2=2(1-alpha beta)`. Dividing by the correct physical
scale gives exactly `2|B_perp_to_b|<=R(|B|^2+|b|^2)`.
The paid b fit implies `|b|>1.43`, so the perpendicular exists.
No arbitrary exclusion of repeated anchors is used.

The ideal size enclosures follow from the accepted cosine bracket.
The ratio bound `.30/1.43` and strict squared projection guard
`3085359/204490000>0` establish a horizontal projection above `.97`.
The reverse triangle inequality then gives a left side greater
than `993/1250`; the right side is less than
`60886704089/100000000000`. Their exact bound difference is

    18553295911/100000000000 >9/50.

This contradiction remains strict at `D=D0`, as required by the
strict theorem statement.

Every zero-area limiting face, including a repeated apex or another
collapsed side, has limiting score at most `1/2+t pi` by the universal
signed-anchor cap and nonpositive quadratic term. This is below
`46790713/70000000<7/10`, whereas
`C-19/500>179/125`. Such a face cannot occur in a hypothesized
high-score triple since each other face has score at most C.
This excludes zero-area closures before any localized division.
Nonzero-area limits retain their geometric parameters and obey the
same strict final contradiction. The closure argument does not
provide a whole-cover generic or topological reduction.

## Scope and remaining external obligations

The accepted objective contains each triangle area, its actual
signed edge anchor, the quadratic base potential and the fixed apex
tilt exactly once. The proof does not consume actual BC-lens losses
or BI/centered credits. A valid global signed-shoelace budget may
replace these same three marked-face/anchor terms once, provided
its eligible consecutive fan and charged incidences are certified
with these same coefficients. Multiple such payments require
disjoint charged face/edge-anchor sets.

The global identity, fan existence/incidence, all remaining profile
rows, nonordinary cases and any monotone extension to smaller radii
remain external. There is no deduction here that replaces signed A
by the independent J allowance. There is no construction or exact
100-circle optimum claim. These scope limitations are explicit in
the frozen source and are not mathematical gaps in its local claim.

A supplemental read-only geometry crosscheck independently confirmed
Sections 3–6 after my full static review; it executed no code and
changed no files. Its observation that Lb is conservative agrees
with the coefficient computation above. It is supporting analysis,
not a replacement for this complete source review or the single
fresh replay receipt.
