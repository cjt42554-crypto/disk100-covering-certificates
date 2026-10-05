# Independent static review of fixed-slack directional lens support

Verdict: **PASS_EXACT_FIXED_TAU_DIRECTIONAL_LENS_SUPPORT**.

For `R>0`, `0<x<1`, `y=sqrt(1-x^2)`, and a feasible fixed signed-anchor
slack `tau=x*v` with `0<=v<=2y`, the exact support over the two-endpoint
disk intersection is

```
max (a e+b n) dot (q-q*)/R
 =|a|[sqrt(x^2+2yv-v^2)-x]+b v
```

for every real a and b. The sign of the normal coefficient b is
correctly preserved. Equations (1)--(6), the stated converse,
attainment, zero-slack case and endpoint scope are valid.
I found no mathematical blocker in this one-base lemma.

This is an exact lens support lemma, not a two-face fee. Additional
unit-anchor, fan, arc or global constraints must still be imposed in
an application; exact attainment in the unconstrained fixed-slack lens
does not imply attainment after those additional restrictions.

## Reviewed bytes and clarification

I read the complete original and final `DIRECTIONAL_LENS_LEMMA.md`,
the original and final author manifests, and the revision record.
The original attainment sentence referred to the sign of a without
separately treating a=0. I requested the explicit clarification that
for a!=0 one may take `t=sign(a)T_x(v)`, while for a=0 every feasible
t attains equality. The author made that document-only clarification
and preserved the original source and manifest. The support formula
and all inequalities are unchanged.

The following files are relative to `../q2_actual_anchor_direction/`:

| Reviewed file | SHA-256 |
| --- | --- |
| `DIRECTIONAL_LENS_LEMMA.md` | `2d0df21426259167a02ad874b483255f22f021ff88d56dcce4df4007c9905248` |
| `REVISION_RECORD.md` | `9cfb61e4758677512f7e40495064537bb2a362d1715dcb12a2a1d8c47c948012` |
| `SHA256.json` | `d890fc66b633f72ce67d2bb67245fc744561b6986c531a471fad58c901934e9a` |
| `history/pre_zero_coefficient_DIRECTIONAL_LENS_LEMMA.md` | `43763af68a53506cc895ec838968cbde243b883979be2bcd941cd6bf93add999` |
| `history/pre_zero_coefficient_SHA256.json` | `0cdd0715d7f1bee50f1276a7a059a217e59c7b0869fb2c0fa7b4b6e8bcf815ad` |

## Exact feasible set and signed area

Because `det(e,n)=1`, the directed-base determinant is
`2R^2 x(-y+v)`. Hence the defined actual signed external area is
`A=R^2 x(y-v)` and its slack is exactly `tau=xv`.
No absolute anchor area was substituted. The chosen external tip is
on the -n side, so the correction is exactly `R(t e+v n)`.

Taking the larger of `(t+x)^2` and `(t-x)^2` gives
`t^2+2x|t|+x^2`. Since `x^2+y^2=1`, both endpoint disk inequalities
are therefore equivalent to (1). Its nonnegative tangential part
forces `v(v-2y)<=0`, giving `0<=v<=2y` and nonnegative tau.

Writing w=|t|, (1) is equivalent to

```
(w+x)^2<=x^2+2yv-v^2.
```

The square root is at least x throughout the feasible v interval.
Thus w lies exactly in `[0,T_x(v)]`, and every such t obeys both
original disk inequalities. This proves the converse and shows that
no outer-boundary assumption or relaxed error disk was used.
At either v endpoint, T_x(v)=0. Values v>y give negative signed A
and are still included.

## Exact directional maximum and its bounds

At fixed v the chosen functional is `a t+b v`, with
`-T_x(v)<=t<=T_x(v)`. Maximizing this interval's linear function
gives `|a|T_x(v)+b v`. The corrected attainment statement covers
positive, negative and zero a. It also covers arbitrary signed b
because b v is fixed during this maximization.

Rationalization gives

```
T_x(v)=v(2y-v)/(sqrt(x^2+2yv-v^2)+x).
```

The denominator is positive for x>0. The identity extends to v=0
with both sides zero and also holds at v=2y. Its tangential ratio
is nonnegative and at most `2y/(2x)=y/x`, since the numerator is at
most 2y and the denominator is at least 2x. Multiplication by
`|a|v>=0` therefore proves both inequalities in (4) while leaving
b unchanged. Even if `b+|a|y/x` is negative, the inequality direction
is unchanged: the potentially negative term b v was never replaced
or multiplied by an unknown sign. At tau=0 all directional errors
and both bounds are zero.

The same tangential bound gives

```
sqrt(t^2+v^2)<=v*sqrt(1+y^2/x^2)=v/x=tau/x^2,
```

proving (5). Since the exact normalized squared correction is
`t^2+v^2`, (6) is just (1) rewritten using v=tau/x. These are
different consequences of the same slack, not separate payments.

## Degeneracy, execution and scope

At x=1 the two radius-R endpoint disks are tangent and intersect
only at the midpoint. Thus q=q*, A=0 and tau=0. At x=0 the base
endpoints coincide, A and tau vanish for every anchor in the disk,
so tau cannot bound its motion. The lemma correctly asserts no
division-by-x formulas there. It assumes R>0 throughout.

My mathematical execution count is **zero**: no mathematical program,
Fraction calculation, source import, guard, checker, trigonometric
evaluation, replay, scan or optimizer was run. I used file reads,
hash inspection and handwritten algebra, and wrote only this separate
review directory. The author records zero mathematical executions and
preserved both original files when clarifying a=0. This reviewer made
no author or frozen-source edit.

The lemma is self-contained algebra from the actual endpoint disk
inequalities. It uses no saved scalar enclosure or external numerical
guard as a mathematical premise. A later two-face argument must match
the actual coefficient direction, retain both unit-anchor equations,
use its real middle HALF-arc variance once, and preserve the distinction
between this actual edge-anchor tau and assigned-cell BC lens slack.
No general positive-tau fee, global Q condition or full-cover exclusion
is accepted here.
