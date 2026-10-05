# Independent static review of the continuous parameter actual-tau lemma

Reviewer: `/root/q2_uniform_review`. Verdict: **STATIC PASS** for the
closed parameter rectangles and supplementary support statements in
`q2_parameter_uniform_tau/PARAMETER_UNIFORM_LEMMA.md`, SHA256
`6c131c32e6f679f4557f9c1a31b95864a357ce33bf1949fba4c1fae31531fd74`.
This is a new review independent of the parameter author. Mathematical
executions, mathematical imports, checker replays, numerical scans and
optimizers: **0**. Sources and prior reviews were not changed.

The accepted domain is exactly R=5661/50000, two actual consecutive
consistently oriented positive faces of circumradius at most R, actual
covered unit anchors with signed A of either sign, the SAME actual
middle cap |c1|<=sqrt(1-R^2), and both ORIGINAL actual connected outer
full assigned arcs strictly above .147. The support C(lambda,t) is the
matched complete support for the perturbed objective. The review
accepts

    D(lambda,t) > .0223 - .6 |lambda-.14|
                           - .1723 |t-1071883/20000000|

for a<=.0004 and b<=.0002, including the closed rectangle boundary.
Its uniform broad fee is .02202554. On a,b<=.00002 the uniform fee
is .022284554. The latter exceeds the old fixed delta target, with
the source's explicitly stated limitation that a changed global
objective requires its own matched target and budget.

## Matched support and the bootstrap

For fixed x, the gamma curvature of the independent support is
positive; its convex endpoint reduction is unchanged by lambda,t.
Thus both full endpoint branches are required, and evaluating their
new supremum at the central acute peak is always a valid lower
support plane. Direct subtraction yields

    D >= D0 + 4u sum(x_i^2-x*^2)
                - v sum(gamma_i-theta*).

I independently checked the factors and signs: the score's quadratic
term is -4lambda x_i^2 and there are two copies of C, so the central
support plane subtracts 8u x*^2 and adds 2v theta*.

The pinned central cosine interval (.676,.677) gives
.541671<x*^2<.543024. Thus the global base difference is below .544
per face. Positive triangle angles and pi<22/7 give
|gamma_i-theta*|<2.315. The resulting loss 4.352a+4.63b is valid
before any high-score localization. At the broad rectangle's corners
it is .0026668. If D0>=.025, the perturbed deficit is at least
.0223332, already strictly above .0223 and the requested fee.
The only remaining branch is D0<.025, where the normalized central
per-face paid budget is available. This avoids circular use of a
central low-deficit hypothesis inferred from a perturbed objective.

I read the full pinned `actual_anchor_1_40_PROOF.md` and its certificate.
Only its scalar per-face support/localization/budget is used here.
Those formulas are dimensionless; the historical three-face geometry
and its .025 compatibility fee are not transferred to the present
two-face statement. The separately accepted full actual-tau theorem
supplies the geometric premise D0>.0223 at the specified R22.

## One paid central deficit and zero endpoints

The derivative of sin^2 theta has magnitude at most one. Writing
e_i=theta_i-theta*, gamma_i-theta*=e_i+z_i gives
sum|x_i^2-x*^2|<=sum|e_i| and
sum|gamma_i-theta*|<=sum|e_i|+sum z_i. The simultaneous central
budget bounds these by sqrt(D0) and D0, respectively. It follows
that

    D >= f(D0)=(1-b)D0-(4a+b)sqrt(D0).

No extra copy of tau, sigma, variance, or Q is introduced: these
are score-perturbation estimates from the same paid central deficit.
For D0>.0223 and below .025, sqrt(D0)>.1; at the rectangle extremes
f' > 1-.0002-.0018/.2=.9908. Hence f(D0)>f(.0223).
The nonstrict estimate sqrt(.0223)<.15 gives the claimed linear
fee. If a=b=0 the transport loss vanishes and the strict conclusion
comes from D0>.0223. If only one parameter increment vanishes,
the same argument applies without a strict zero-product assertion.

The broad fee arithmetic is
.0223-.00024-.00003446=.02202554. The narrow arithmetic is
.0223-.000012-.000003446=.022284554. I independently checked the
large exact delta comparison by integer paper multiplication:
22283*H-1000000*N equals
2178965591728225095328099000000, as recorded by the author.
Thus delta<.022283 and the narrow margin is greater than .000001554.
This is a comparison to the old fixed target, not a proof that a
perturbed global problem has that same target.

## Full-branch peak stability and residual

On the full branch domains the perturbation magnitude is at most
P=4a+(22/7)b<.00223. The central endpoint localization gaps exceed
.038, while the complete obtuse support gap exceeds .9. Comparing
a possible new maximum with the value at the old peak loses at most
2P, so every new maximum is acute and lies strictly inside (.7,.96).
The old peak monotonicity is used only for this outside-window
exclusion, not as a false concavity statement on (0,pi/2).

On the admitted window, -G_new''>4-8a>=3.9968. There is one new
interior acute peak. At its stationary point
|G0'(theta_new)|<=q=4a+b. Integrating the central derivative with
-G0''>4 proves |theta_new-theta*|<=q/4, including the nonstrict
q=0 endpoint. The uniform maximum is .00045.

The exact central curvature estimate gives
G0(theta)<=C0-2(theta-theta*)^2 on that window. The objective
increment relative to the old-peak value is q-Lipschitz there.
Completing -2r^2+q|r| gives the residual bound q^2/8; no new
maximizer or numerical support value is required. Its largest value
is .000000405. The exact deficit decomposition with +2V agrees
with the support-plane comparison, so there is no factor-of-two
loss in the residual statement. Only V>=0 is used for the fee.

## Scope and evidence

The current source, ledger, local author pin record, and every pinned
dependency are checked by metadata SHA256 and byte count in this
review's SOURCE_PINS.json. That operation does not import or execute
mathematical sources. The complete central full-tau theorem, its
independent audit, mean-envelope audit, and endpoint clarifiers were
read as accepted premises; no predecessor checker was rerun.

The conclusion preserves the geometric premises. It supplies no
radius interval, arc-threshold change, middle-cap relaxation, global
incidence or covering certificate. A conditional local fee cannot be
promoted to a global lower bound without the separate actual-fan
admission and once-only global budget. No unresolved mathematical
defect was found in this frozen parameter lemma.
