# Independent static review of the Pro-domain parameter corollary

Reviewer: `/root/q2_uniform_review`. Verdict: **STATIC PASS** for the
new, separately frozen `PRO_PARAMETER_COROLLARY.md`, SHA256
`43059d525b9e23f4c1639fe6d95cfeea45c5b55fb993daf20aaf6a29b5a12f1b`.
This review is independent of the parameter author. Mathematical
executions, mathematical imports, old/new checker calls, numerical
replays, scans and optimizers: **0**. The earlier one-run Pro replay
is a pinned predecessor receipt; it was not repeated by this review.
The local .0223 parameter theorem and its review remain unchanged.

The accepted domain fixes R=5661/50000. It requires two same actual
consecutive positive nondegenerate faces with circumradii at most R,
actual covered unit anchors on their bases, SAME actual middle squared
norm at most 1-R^2, and both actual outer squared norms at least .83.
It permits arbitrary signed A and actual tau, with no arc, symmetry,
nearest-owner, global Q, incidence, or variance hypothesis. For the
matched complete perturbed support it establishes

    D(lambda,t) > .0225-.6 |lambda-.14|
                           -.1725 |t-1071883/20000000|

on the CLOSED rectangle a<=.0003, b<=.0002. The uniform fee is
.0222855; on a,b<=.0002 it is .0223455. This verdict is conditional
local geometry at R22 and is not a new global radius bound.

## Independently extended .025 scalar budget

Unlike a citation to the Pro theorem's .0225 contradiction range,
the new source rederives its paid per-face implication up to .025.
I checked this derivation in full. The support includes both convex
apex endpoint branches. C0>=G0(5/6)>1.4 and Bo<1.2 imply that
every actual score with central deficit at most .025 exceeds the
complete obtuse branch. Thus its G0 value exceeds its actual score.
The old .035 endpoint gaps and acute full-branch monotonicity force
.7<theta_i<.96. No false global concavity assumption is used.

On the full allowed apex interval,

    L_z=.5 sin^2(theta)/sin^2((theta+z)/2)-t0
        >=1-cos(theta)-t0>.17>0.

The last inequality follows from the fourth-order upper cosine at
.7 being below .77 and t0<.06. Up to z=.1, the elementary ratio
lower bound is 3871/2500: cos(theta/2)>=553/625,
cot(theta/2)<=20/7 and the .05 sine/cosine bounds give exactly that
ratio. Its squared half minus t0 exceeds one. Thus z>=.1 would
already cost more than .1, contradicting a deficit at most .025.
For actual z<.1, L>=z. Integrating -G0''>4 from the unique central
peak and retaining actual nonnegative sigma and tau gives

    2 sum(theta_i-theta*)^2+sum z_i+sum sigma_i+sum tau_i<=D0.

This is the newly proved simultaneous scalar budget at .025. It
depends only on normalized support, circumradius and actual endpoint
coverage; it does not import a local .0223 geometric conclusion.
I also read the separate old full-apex scalar proof and certificate.
They identify the same exact unique maximizer and supply the central
cosine bracket and stronger endpoint/obtuse gaps used later.

## Support-plane bootstrap and single paid transport

Evaluation of the new support at theta* gives the lower plane
C>=C0-4u x*^2+v theta*. Direct subtraction from the SAME actual
fan therefore gives

    D>=D0+4u sum(x_i^2-x*^2)-v sum(gamma_i-theta*).

The central cosine bracket implies |x_i^2-x*^2|<.544 over every
possible base, and the positive apex range implies
|gamma_i-theta*|<2.315. Hence the pre-localization gross loss is
4.352a+4.63b. Its broad maximum is .0022316. If D0>=.025,
the perturbed deficit is at least .0227684, already above .0225.
Only the remaining D0<.025 branch invokes the paid localization.
This establishes the required bootstrap without circularity.

On that branch the Lipschitz bound for sin^2 theta and
gamma_i-theta*=e_i+z_i give the same budget comparison

    D>=f(D0)=(1-b)D0-(4a+b)sqrt(D0).

Cauchy yields sum|e_i|<=sqrt(D0), while sum z_i<=D0, from one
central deficit. These estimate the score perturbation; no additional
copy of tau, skew, variance, Q or the geometric fee is charged.
The derivative of f is positive for d>=.0225 on the new rectangle,
since sqrt(d)>=.15 and q<=.0014. The accepted central theorem
supplies D0>.0225, so f(D0)>f(.0225). Here sqrt(.0225)=.15
exactly; hence the final coefficients are .6 and .1725. At zero
parameters the strictness comes from the central theorem. Zero
slack or a vanishing individual increment does not cause a strict
zero-product error.

The broad fee arithmetic is .0225-.00018-.0000345=.0222855.
The balanced fee arithmetic is .0225-.00012-.0000345=.0223455.
The large cross-product for the old delta target was checked by
paper integer multiplication and is the recorded positive integer
2178965591728225095328099000000. Therefore the margins over
that OLD fixed target are strictly greater than .0000025 and
.0000625, respectively. This does not hold a changed global target
fixed: the source correctly requires a new matched global budget.

## Support movement and composition limitations

The uniform full-branch perturbation is at most
P=4a+(22/7)b<.00183. The pinned common scalar endpoint gaps
exceed .038 and the obtuse gap exceeds .9. Losing 2P still forces
every new support maximum into the acute interval (.7,.96).
There -G_new''>4-8a>=3.9976, so the new peak is unique. Its
stationarity gives |G0'|<=q=4a+b, and the old curvature yields
|theta_new-theta*|<=q/4<=.00035.

The exact support residual follows from
G0(theta)<=C0-2(theta-theta*)^2 and a q-Lipschitz parameter
increment. Completing -2r^2+q|r| gives 0<=V<=q^2/8<=.000000245.
This rederives the supplementary result for this new rectangle;
it is not accepted merely by analogy to the earlier .0223 lemma.

The scalar transport implication is normalized and can be used at
another R ONLY after a central theorem at that R and its actual
geometric hypotheses are separately proved. This review accepts the
present fixed-R22 theorem. It does not silently combine it with the
radius theorem or relax its middle cap; those changes need explicit
composition. The actual outer floor .83 is preserved and is not
inferred from a fixed assigned-arc length at another R. Global fan
admission, profile coverage, once-only allocation and exact optimum
remain external. No unresolved mathematical defect was found.

Metadata identities for the new author source and ledger, its local
pin record, all accepted Pro/central scalar inputs and the unchanged
local parameter source are recorded in SOURCE_PINS.json. Reading and
hashing these sources did not execute or import their mathematics.
