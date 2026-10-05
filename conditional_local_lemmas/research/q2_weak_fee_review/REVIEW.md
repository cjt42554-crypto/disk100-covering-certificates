# Independent static review of the conditional weak q2 fee

**PASS for the final corrected conditional local fee in Section 2.**
For the two actual
positive-area consecutive consistently oriented faces with covered unit
anchors of either signed side and circumradii at most R, the statement is

```
R=5661/50000, lambda=7/50, t=1071883/20000000,
|p|<=1,
at least one outer center c0 or c2 covers a connected unit-circle
arc with xi=sin(arc_length/2)/R>=3/5,
S0+S1 < 2C-1/10000000.
```

The support C is the accepted complete one-face support at these exact
coefficients. The theorem makes no nearest-owner, reflection or norm
assumption on the other centers. It is a very weak local actual-A fee,
not the requested larger q2-row payment or a global covering exclusion.

I read the complete final `AUTHOR_ANALYSIS.md` for context and independently
checked Section 2 and its Section 1 ideal-coordinate formula. I also read
the complete diffs against both retained correction histories. I did not
audit or accept the new mirrored-skew obstruction in Sections 3--4 or
its targeted-helper guards. The accepted single-face scalar and
area/skew/anchor sources had been read completely in the preceding
independent reviews; their exact source pins and accepted status are
retained here. The assigned-arc norm lemma was independently reviewed
immediately before this task.

This is a static analytic audit. No checker, Fraction program,
mathematical module, optimizer or geometry search was executed by this
reviewer. I wrote only this separate review directory. The root interface,
author analysis and previously frozen files were not edited.

## Complete support and one budget for two faces

The accepted exact single-face decomposition is `F-sigma-tau`, including
negative signed anchors. Full gamma convexity gives every actual score
at most C before localization. Thus each `d_i=C-S_i` is nonnegative.
The accepted endpoint gaps, obtuse gap and convex gamma cut apply
separately to either face; summing the resulting pointwise bounds works
for two faces and uses no three-anchor compatibility argument.

Under `D=d0+d1<=10^(-7)`, those facts give

```
.7<theta_i<.96,
2 sum e_i^2+sum z_i+sum sigma_i+sum tau_i<=D,
e_i=theta_i-theta*, z_i=gamma_i-theta_i>=0.
```

In particular gamma is acute before the radius/tip estimates. The
accepted `.79` skew coefficient, signed anchor error `625/256`, tip
amplitude bound `11/4` and `Qmax=2.54` are legitimate at this tiny
deficit: their earlier valid domain was `D<=1/40`. No fee, fitted
three-anchor error, or joint geometric payment is transferred to q2.

Each zero-area limiting face has score at most `1/2+t*pi<.7`, well
below `C-10^(-7)`. Since the other face is at most C, all such limits
are excluded before any division or phase localization.

## Actual two-anchor error from the same paid budget

Rotate the entire physical coordinates so the nonzero shared ray
`c1-p` is positive real. This preserves unit-anchor norms and `|p|`;
it does not translate the unit circle or force p to the origin.
The actual bisector phases are `-gamma0/2,+gamma1/2`.
The ideal normalized anchor vectors are `(T,-H),(T,H)` with

```
Q=1+2cos(theta*),
T=Q*cos(theta*/2), H=Q*sin(theta*/2).
```

The stated `2.15<T<2.16`, `.94<H<.95` follow directly from the
accepted `.676<y=cos(theta*)<.677`. For example Q is between
2.352 and 2.354, the positive half-angle cosine is between .915
and .916, and the positive half-angle sine is between .4 and .403.
Their corresponding products enclose T and H as claimed.

On the whole theta window the derivative norm of
`Q(theta)*exp(+-i*theta/2)` is at most `2sin(theta)+Q(theta)/2<3`.
The amplitude slack costs `11z/4`; the individual bisector phase
slack costs at most `Qmax*z/2`, giving `4.02z`. The signed anchor
coefficient is smaller than 4.02. Therefore the individual w error
satisfies the author's displayed inequality. Cauchy uses weights
3 and .79 against the single angular/skew budget, giving

```
|w_i-w_i*| <= sqrt(9/2+.79^2)*sqrt(D)+4.02D < eta,
eta=3/4000=.00075.
```

The squared rational comparison is correct:
`(9/2+.79^2)*10^(-7)=.00000051241<.00072^2`, and
`.00072+4.02*10^(-7)<.00075`. Neither loss type is spent again.

## Two unit anchors determine the physical translation up to sign

The two individual errors imply

```
|v-(T,0)|<eta, |b-(0,H)|<eta,
v=(w0+w1)/2, b=(w1-w0)/2,
a=(q0+q1)/2=p+R*v,
a dot b=0, |a|^2+R^2|b|^2=1.
```

These last two identities are exact unit-anchor identities. The b
fit gives `|b|>H-eta>0`, so its normal exists; repeated anchors are
excluded by this derived high-score bound rather than by fiat.

Both `u=sqrt(1-R^2*H^2)` and `|a|` exceed .994. One simple exact
justification uses `R<.114`, `H+eta<.951` and
`1-(.114*.951)^2>.994^2`. The normal n to b with positive x component
exists because `b_y>H-eta>0`. Normalizing b and comparing to `(0,H)`
gives

```
|n-(1,0)| <= 2eta/(H-eta) < (107/50)*eta.
```

The difference of the squared magnitudes gives the precise estimate
for `||a|-u|`. With its positive denominator `|a|+u>1.988`, the stated
coefficient bound `<eta/50` is valid. Since `|a|<=1`, the two estimates
imply `|a-(+-u,0)|<(9/4)eta`.

The negative sign is impossible under the actual `|p|<=1` hypothesis:
using `v_x>T-eta>0`, `R>.11` and `u>.994` gives

```
p_x < -.994+(9/4)*eta-.11*(2.15-eta) < -1.
```

For the positive sign, comparison to `p*=(u-RT,0)` gives
`|p-p*|<(9/4+R)eta<(12/5)eta=.0018`. This is a full two-dimensional
translation bound, not merely a bound on its real component. The
apex norm hypothesis is essential to the sign choice and is retained
in the theorem's statement.

## Correct combined radial and direction error

The outer ray direction differs from its ideal direction by at most
`|e_i|+z_i`. The exact normalized radial side length is
`2r_i*cos(gamma_i/2+-skew_i)` and its ideal value is
`2*cos(theta*/2)<2`. Comparing lengths by trigonometric Lipschitz
bounds and `1-r_i<=(6/5)z_i` gives

```
radial-length error <= |e_i|+(17/5)z_i+2|skew_i|.
```

Comparing the actual and ideal direction vectors costs an additional
`2(|e_i|+z_i)`. Therefore the correct combined vector bound is

```
3|e_i|+(27/5)z_i+2|skew_i|.
```

The earlier 17/5 coefficient for this combined expression omitted
the direction's `2z_i`. The corrected 27/5 coefficient includes that
payment and still leaves the advertised positive rational slack.
The historical source containing the old coefficient is retained
separately; no mathematical rerun is needed for this analytic correction.

The single budget implies `|e_i|<1/4000` and
`sqrt(sigma_i)<1/3000`. The accepted skew identity, `|skew_i|<pi/2`,
and sine concavity give

```
|skew_i| <= (pi/2)*sqrt(625/512)*sqrt(sigma_i)
          <= (9/5)*sqrt(sigma_i) < 3/5000.
```

The coefficient estimates `pi/2<8/5` and
`sqrt(625/512)<111/100` are sufficient. At sigma zero the skew is
exactly zero, so the universally valid multiplied estimate is
nonstrict; it becomes strict when sigma is positive. The final author
source explicitly makes this correction and preserves the previous
display in `history/AUTHOR_ANALYSIS_before_zero_skew_correction.md`.
The next positive bound `<3/5000` is strict for all cases and is all
the later proof uses. Both this correction and the earlier direction
coefficient correction were document-only changes without a
mathematical rerun.

Thus the corrected combined error is strictly less than

```
3/4000+(27/5)*10^(-7)+6/5000=.00195054<1/500.
```

Together with the physical translation estimate, each outer center
is within `.0018+.12/500=.00204<.0021` of its ideal outer center.

## Strict upper/lower norm gap

For the ideal pair, the outer coordinates simplify to
`(u-Rk,+-2Rkx)` with `k=cos(theta*/2)`, and direct squaring yields
exactly

```
rho_ideal^2=1-2Rku+R^2*(2+y-2y^2).
```

The accepted y bracket gives `k>.915` and
`2+y-2y^2<1.764`; the preceding argument gives `u>.994`.
Using `.1132<R<.1133` has the correct signs in the lower negative
term and upper positive term. The stated rational bound is

```
rho_ideal^2 < .81673120996 < .904^2.
```

The outer-center error therefore gives actual outer norms below
`.9061` for both c0 and c2.

If either one covers a connected arc with `xi>=3/5`, the independently
reviewed endpoint norm lemma instead requires

```
rho >= sqrt(1-.36R^2)-.8R
     > .9976-.8*.1133=.90696>.9069.
```

The square-root comparison is strict: at the looser R=.1133,
`1-.36R^2>.9976^2`. A connected unit-circle arc covered by R<1
cannot contain an antipodal pair, so its length is less than pi
and it lies in the small-arc domain of that norm lemma. Actual
endpoint coverage, rather than ownership or formal lens replacement,
is all that lemma needs. The strict .9061/.90696 gap contradicts
the hypothesized high score also at the closed `D=10^(-7)` endpoint.

## Acceptance boundary and remaining scope

The saved full-apex scalar proof/certificate at lambda 7/50 and t0
is an explicit accepted arithmetic premise, independent of R; it
was not regenerated. The already accepted area/skew/anchor formulas,
ratio cut and one-budget estimates are used only in their valid
domain. The q3 joint fee and three-unit-anchor fit are not used.
The independently accepted assigned-arc norm lemma contributes a
necessary geometric condition and consumes no score or lens credit.

This proof does not establish that an arbitrary cover supplies
`xi>=.6` or the needed q2 fan. It does not pay the parent's larger
conditional q2-row target. It does not verify the radial Q premise,
four-cell window, empty fan degree, remaining profiles or completion
to 100 circles. Sections 3--4 of the author analysis and their new
helper arithmetic remain outside this Section 2 verdict.

The author's execution ledger reports helper computations elsewhere
in its q2 work; those must not be mislabeled as zero execution for
that entire project. This static Section 2 review itself made zero
mathematical calls or imports and makes no fresh-runtime claim.
No mathematical blocker remains for the corrected weak-fee mechanism
and its stated hypotheses.
