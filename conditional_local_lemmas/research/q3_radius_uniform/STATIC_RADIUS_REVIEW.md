# Independent static review of the uniform-radius addendum

**Verdict: PASS for the stated local theorem, including the closed upper
radius endpoint.** The fixed coefficients are `lambda=7/50` and
`t=1071883/20000000`; the normalized three-face fee remains `1/40` for
every real `0<R<=949600/4987489`. This review does not certify a global
covering bound or optimality.

## Reviewed new sources

Both new files were read completely. Their SHA-256 values at review are:

| File | SHA-256 |
| --- | --- |
| `PROOF_ADDENDUM.md` | `4f8ef8ce523091f6dd22dc6ed9064d736eab8d1f4f4cdf2ddbbf9bb8f3b506a3` |
| `radius_endpoint_fraction.py` | `a0c916945229dc79bfdc9e5beea6d933e3c50a2566f923093d7bdacba6324649` |

I also read the complete frozen local proof and checker, its separate
degenerate-closure addendum, the scalar conjugate proof, checker and
certificate, and the complete primitive source. The author's source-pin
manifest and execution receipt are separate from this review. No frozen
source, certificate, or manifest was edited.

## Uniform statement and normalized geometry

The statement uses the same positive `R` consistently in the circumradius
hypothesis, both endpoint coverage inequalities, base normalization,
and area/anchor normalization. The physical anchors remain unit points.
Consequently it quantifies over admitted configurations at each radius;
it does not assert invariance when a fixed configuration is relabeled
with a different radius.

The sine-rule area and skew identities, feasible apex interval, signed
anchor loss and full-apex endpoint support all remove `R` exactly. The
same scalar `C` and `Co` therefore apply. The accepted localized angle
window, curvature and convex cut, single deficit budget, tip corrections,
two-sided Jacobian bounds, and ideal fit have no remaining radius
dependence. All divisions used by the localized argument have positive
denominators.

The signed convention is consistent. For a positively oriented triangle,
the directed base has the apex on its left. Taking positive anchor height
toward its right gives `A/R^2=x*h`; it also gives this identity when
`h<0`. The two disk inequalities imply `tau=x*(cos(theta)-h)>=0`, so
the actual anchor is never replaced by its absolute area. Resolving the
sine-rule ray lengths along the actual bisector gives the stated base
midpoint expression; its outward normal is `exp(i*m-i*delta)`. Adding
the tip height gives the exact `U,V` expression in the addendum. The
common translation `p` cancels in `B,b` without any assumption on `p`.

The corrected skew sentence is valid also when `sigma=0`: its squared
displacement is at most `(79/100)^2*sigma`, and is zero for zero skew
loss. No false strict product inequality at zero is used. The separate
strict rational guards `alphaB2*D0<EB^2` and `alphab2*D0<Eb^2` supply
the strict final fit, including when `D=D0` or all skew losses vanish.

## Strict closed endpoint

The three-unit-point condition retains exactly one physical factor `R`:

```
2|B_perp_to_b| <= R*(|B|^2+|b|^2).
```

The fit gives `|b|>149/100`, so its normal is defined for every positive
`R`; an initially allowed repeated pair `q_0=q_2` cannot occur in the
hypothesized high-score triple. The scalar square identity used to prove
the condition has nonnegative right side, and `1-alpha*beta>=0`, so
taking its square-root comparison is justified.

The unchanged strict envelopes give

```
K=2*((19/25)*(49/50)-27/100)=1187/1250,
S=(1033/1000)^2+(99/50)^2=4987489/1000000.
```

Thus, for every `0<R<=K/S`, the complete contradiction is

```
K < 2|B_perp_to_b| <= R*(|B|^2+|b|^2) < R*S <= K.
```

At `R=K/S=949600/4987489`, the two fit envelopes are still strict,
although the coarse-envelope difference `K-R*S` is zero. The endpoint
is therefore included. The original fixed-radius margin `>19/50` is
correctly omitted from this extension. The addendum correctly describes
this as the endpoint supplied by these particular envelopes, rather
than a sharp threshold for the local fee.

The endpoint source's two simpler-interval comparisons are correct:

```
K-(19/100)*S = 197709/100000000 > 0,
K/S-19/100 = 197709/498748900 > 0.
```

## Radius dependencies and degeneracies

Static inspection confirms that the frozen local checker's `R` is used
only in its final `rhs` and returned radius metadata. Its earlier guard
expressions do not reference `R`. The scalar conjugate checker and its
certificate have no radius parameter. In the primitive source,
`I`, `outward`, `atan`, `PI`, `sinq`, `cosq`, and their rational,
factorial and integer-square-root dependencies do not reference the
unrelated boundary radius `9041/80000`. That radius belongs to the
separate boundary `verify()` routine and its output. This supports the
addendum's dependency audit without executing any imported definition.

The fixed-positive-radius closure argument is valid. As actual area
tends to zero at a fixed `R>0`, normalized area also tends to zero;
the universal anchor cap, nonpositive quadratic term and angle bound
give `limsup score<=1/2+t*pi<7/10<C-D0`. It excludes every zero-area
limiting face, including collapse of a side different from the base.
Positive-area limits retain defined parameters and satisfy the proof.
The addendum explicitly excludes `R=0` and correctly warns that physical
area tending to zero while `R` tends to zero need not mean normalized
area tends to zero. No such varying-radius collapse is silently assigned
the fixed-radius degenerate bound.

## Computational and global scope

The new endpoint source performs only exact standard-library Fraction
arithmetic on the stated rational envelopes, with JSON serialization.
It imports no mathematical project source, performs no trigonometric
calculation, and invokes no frozen verifier, optimizer, scan or geometry
search. I read it but did not execute it. Its output counters describe
the intended endpoint invocation; this review adds no execution receipt.

During this review there were zero checker calls, zero mathematical
module imports, zero optimizer/scan/search calls, and zero new scalar
guard executions. Only file reads, source hashing, and writing this
separate report were performed.

The extension changes no whole-cover boundary constant, profile row,
incidence theorem, topology, or nongeneric reduction. It requires actual
signed anchors, matching coefficients and disjoint charged face/anchor
sets when used in a valid signed area identity. Its normalized fee is
`1/40`; the corresponding physical area payment is `R^2/40`. No local
mathematical blocker was found.
