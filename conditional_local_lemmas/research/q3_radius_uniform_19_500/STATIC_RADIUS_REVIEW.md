# Independent static review of the 19/500 uniform-radius extension

Verdict: **PASS_STATIC_UNIFORM_RADIUS_LOCAL_ACTUAL_A_FEE_19_500**.

Using the separately accepted stronger normalized premises, the exact
fixed pair `lambda=7/50`, `t=1071883/20000000` supports

```
sum_i(area_i/R^2 + A_i/R^2 -4 lambda x_i^2 +t gamma_i)
                                      <3 C_lambda(t)-19/500
```

for every `0<R<=794400/5378209`, including the upper endpoint, on the
consecutive, consistently oriented, actual common-apex three-face domain
and positive-radius limiting closure stated in `PROOF_ADDENDUM.md`.
In particular this includes all `0<R<=231/2000=.1155`.
I found no remaining mathematical blocker in this local radius extension.
The accepted scalar-conjugate enclosures remain the explicit arithmetic
premise already retained by the stronger independent review; this audit
does not independently regenerate them.

## Reviewed sources and acceptance boundary

I read the final addendum and endpoint source in full, as well as the
entire stronger `PROOF.md`, verifier, generated certificate and guard
addendum. I also read its complete independent analytic review and fresh
replay receipt; and the entire bundled scalar-conjugate proof, checker,
saved certificate and interval primitive source. Finally I read the new
endpoint arithmetic JSON, source-pin map and premise-acceptance record.

The stronger independent analytic review accepts the fixed-radius local
actual signed-anchor theorem, and its one fresh replay matches the
author's complete JSON and serialized certificate bytes. Its accepted
normalized facts are the premises used here. I checked their radius
dependence from the mathematical formulas and complete source, rather
than inferring a range theorem merely from that fixed-radius acceptance.
`PREMISE_STATUS.json` correctly records both the accepted source and the
continuing accepted scalar-certificate premise. The stronger author's
radius metadata remains `11321/100000`.

The reviewed bytes are bound by these SHA-256 values. Extension paths in
the first five rows are relative to this directory; all other paths are
relative to `../q3_stronger_fee_19_500/`.

| Reviewed file | SHA-256 |
| --- | --- |
| `PROOF_ADDENDUM.md` | `eacd376402698410798b1ade9373c786896519a18b7a95c589ba4cb93a480fae` |
| `radius_endpoint_fraction.py` | `7f56e4f3dd808320586ed060a9eb089ec752868a479ea38780cc25cda7783d47` |
| `endpoint_arithmetic.json` | `4b1890a96617817c58481435864678eb867528c90ce7b6448ca94f6e69df85de` |
| `SOURCE_PINS.json` | `c29656471c6dff992ba1f4e87e05e13d7e5408e28d4820e7233e8023de4a891b` |
| `PREMISE_STATUS.json` | `7389a605e8926a2602b2c66ab7395edda8b52fa252874143a515f9dcf398d980` |
| `PROOF.md` | `5a58cc96f916018885d39888bcd07ca77e768d3c719a695e2d619311fed1610f` |
| `verify_actual_quadratic_19_500.py` | `bd042f4976fbeeb8440da66864814371f82e4fd0e380844cd153dcb7e600f46d` |
| `certificate.json` | `b79cc05457486de14e8ed89a8d62bacf596581b152092a3951f9893b3dfa56bb` |
| `GUARD_ADDENDUM.md` | `e2ce31e43a8ed34b0e8ef006287881b9fc0fa94ad0706dd8f8df344407fca5f7` |
| `independent_audit/INDEPENDENT_REVIEW.md` | `9a7b5e434a3b4a7f6989fa8d054afe9bd2554f20d98041cdec7a3ce01480aee1` |
| `independent_audit/execution_receipt.json` | `f19cd7f0b9d249a03e715cf5a44155f2b2de3bf1e2421769628ac41e17ec8ad0` |
| `dependencies/coupled_edges/verify_support.py` | `90fd6e3426d8c113bacedd55f9a13f5e47981e9a48d21bfef48183ec4602efec` |
| `dependencies/quadratic_lambda14_conjugate/PROOF.md` | `aaaad3b1629ba98fd44907735b75ae111fecaa45c8d77f078f4e28acf4713ea3` |
| `dependencies/quadratic_lambda14_conjugate/verify_conjugate.py` | `e95db0a6ead5ace03884ab54743b6238f8fe4cd242b1563617b789fce9cb57f4` |
| `dependencies/quadratic_lambda14_conjugate/certificate.json` | `f325d42bfe651599881b7dc652567729a1586389ee385c38060840bc41696bc2` |

## Radius independence of the normalized premises

The circumradius condition and the normalized half-base
`x=|c_(i+1)-c_i|/(2R)` give `theta<=gamma<=pi-theta` for any positive
`R`. The exact area, skew and signed-anchor formulas scale by `R^2`.
In normalized base coordinates the disk intersection inequality is
`s^2+h^2+2x|s|<=1-x^2`, and the determinant convention gives
`A/R^2=xh`, including negative `h`. Consequently the actual score is
`F-sigma-tau` on the full initial domain for every positive radius.

The scalar support depends on the fixed coefficient pair only. Its
analytic source explicitly starts with arbitrary `R` and eliminates it
by normalization. Its checker and saved scalar certificate have no
geometric radius parameter. Localization, the convex apex cut, paid
acute-apex bound, curvature, and the simultaneous budget

```
2 sum e_i^2+sum z_i+sum sigma_i+sum tau_i <=D<=19/500
```

are therefore radius-independent. The exact normalized tips, outward
normal convention and lifted bisector phases likewise have no hidden
absolute-position or fixed-radius assumption. The common apex cancels
from both `B` and `b`.

The whole-box Jacobians, skew and anchor estimates, Cauchy allocation,
linear-spending optimization, ideal sizes and axis projection all use
only the preceding dimensionless data. The final addendum conservatively
retains the accepted weak interval
`-686839/1250000<=V<=-27/250`, which is sufficient for `|V|<11/20`;
it does not assume an unneeded strict strengthening of those endpoints.
The stronger verifier's only nontrivial use of its geometric `R` is the
final compatibility RHS
`R*((Amax+EB)^2+(Jmax+Eb)^2)`; its other occurrence is radius metadata.
In particular its original `margin>9/50` reporting guard is not asserted
throughout this extension.

The primitive source's unrelated `R=9041/80000` is used by its separate
boundary `verify()` routine. The definitions `I`, `PI`, `sinq`, `cosq`
and their arithmetic helpers do not reference that radius. Reading
those definitions therefore reveals no additional radius restriction.
No primitive or other frozen source was imported or executed here.

## Exact endpoint and strictness

The accepted fit remains, for every hypothesized high-score triple,

```
|B-B*|<17/50, |b-b*|<3/10,
B*=(-A,0), b*=(0,J),
19/25<A<763/1000, 173/100<J<87/50.
```

Hence `|b|>143/100`, which justifies its perpendicular. The accepted
projection guard gives horizontal-axis projection above `97/100`.
The reverse triangle inequality and the size fits then give the strict
envelopes

```
2|B_perp_to_b|>K=2[(19/25)(97/100)-17/50]=993/1250,
|B|^2+|b|^2<S=(1103/1000)^2+(51/25)^2=5378209/1000000.
```

For actual unit anchors, the physical mean `a` and half-chord `h` obey
`a dot h=0`, `|a|^2+|h|^2=1`. Writing the scalar components normal to
`h` as `alpha,beta`, the displayed identity proves
`|alpha-beta|<=1-alpha*beta`, while
`|a-q_1|^2+|h|^2=2(1-alpha*beta)`. Substitution of `a-q_1=RB`,
`h=Rb` and division by positive `R` give exactly

```
2|B_perp_to_b|<=R(|B|^2+|b|^2).
```

Thus `Rmax=K/S=794400/5378209` forces

```
K < 2|B_perp_to_b| <= R(|B|^2+|b|^2) < R*S <= K.
```

The strict fit bounds hold even when `D=19/500`. The high-score assumption
is therefore excluded at the upper-radius endpoint as well as at
loss-threshold equality. It is correct that `K-Rmax*S=0`; the endpoint argument
uses strict actual-fit inequalities, not a positive coarse margin.

The endpoint source and output correctly give

```
K-(231/2000)*S=346433721/2000000000>0,
Rmax-231/2000=346433721/10756418000>0,
11321/100000<231/2000,
19/500-3/100=1/125>0.
```

The resulting endpoint is sufficient for these rational envelopes; no
sharp local validity threshold is proved.

## Signed anchors, degeneracies and zero skew

Negative actual anchor areas remain covered by nonnegative
`tau=x(cos(theta)-h)`, rather than an absolute-anchor substitution.
No ownership rule, assigned arc, bound on the apex, or center norm is
needed. Repeated endpoint anchors are initially allowed and excluded
only in the high-score domain by `|b|>143/100`. The unit-anchor identity
permits all other repetitions without additional assumptions.

The frozen stronger display (12) writes a strict product bound even
at zero skew. The new addendum correctly makes the needed reading
explicit: displacement squared is at most `(43/50)^2 sigma`, and at
`sigma=0` the displacement is zero. This harmless display caveat does
not remove the strict final fits, which follow from the strictly positive
squared-fit coefficient margins. No frozen source was changed to address
it.

For each fixed positive `R`, a zero-area limiting face has score at most
`1/2+t*pi<46790713/70000000<7/10`, whereas `C-19/500>179/125`.
Since the other two faces have score at most `C`, such a limit cannot
be a high-score triple. This excludes collapsed sides before any
localized division. The addendum also correctly distinguishes positive
limiting radius and vanishing normalized area from physical area tending
to zero with `R` tending to zero. It does not extend the theorem to
`R=0` or assign a zero-area bound solely from physical collapse in that
joint limit.

## Execution and scope

This audit used file reads and SHA-256 inspection only. I executed no
Python program, mathematical verifier, trigonometric guard, frozen-module
import, optimizer, scan or geometry search. I wrote only this review;
I did not edit a manifest or any earlier or stronger proof source.

The separate new endpoint source imports only `fractions.Fraction` and
JSON serialization. Its computation consists of the quoted rational
envelopes and comparisons; it contains no frozen-source import, guard
replay or geometric search. I statically reviewed its output, without
executing it again. The owner reports one new endpoint-only invocation,
consistent with the receipt; the stronger fresh replay belongs to the
separate already accepted audit.

The payment is the same actual signed-`A` joint score, with each marked
face and its actual edge anchor charged once. Its physical area payment
is `(19/500)R^2` only within an otherwise valid matching identity.
No BC-lens, BI or centered credit is spent by this extension. It proves
no independent-`J` payment, coefficient rectangle, ordinary-fan existence,
disjoint incidence theorem, boundary/profile budget, nongeneric global
reduction, global radius exclusion, construction or exact optimum.
Those remain separate external obligations, not blockers to the stated
local parameter theorem.
