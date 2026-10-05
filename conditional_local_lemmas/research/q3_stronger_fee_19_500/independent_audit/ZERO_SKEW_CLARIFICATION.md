# Zero-skew endpoint clarification; no replay or frozen-source change

This is a targeted independent analytic clarification of the literal
strict inequality at the end of Eq. (12), and the corresponding
multiplied strict expression in the frozen `INDEPENDENT_REVIEW.md`.
Neither frozen file is edited. No checker, mathematical source,
optimizer, search, or new scalar guard calculation was executed for
this clarification.

Let `E_i` be the norm of the complex skew correction from
`exp(i m_i) U_i` to the actual outer tip. The rigorously sufficient
global conclusion on the localized domain is

    E_i^2 <= skew_sq sigma_i <= (43/50)^2 sigma_i,

where `skew_sq=3081055/4194304<(43/50)^2` is the already checked
strict scalar coefficient guard. The second multiplied inequality
is strict if `sigma_i>0`; if `sigma_i=0`, both sides equal zero.
Thus the unqualified literal expression
`skew_sq sigma_i < (43/50)^2 sigma_i` is false at zero skew and
needs this qualification. The same qualification applies to the
corresponding expression in the frozen review.

There is no unbounded or omitted zero case. Eq. (2) states

    sigma_i=2 x_i^2 sin(d_i)^2/sin(gamma_i).

After localization, `x_i>0`, `sin(gamma_i)>0`, and
`|d_i|<pi/2`. Hence `sigma_i=0` implies `sin(d_i)=0` and therefore
`d_i=0`. The exact correction

    U_i(cos(d_i)-1)+i V_i sin(d_i)

is then exactly zero. In particular the non-strict norm bound
`E_i<=(43/50)sqrt(sigma_i)` is valid for every admitted face,
including the zero-skew endpoint.

No subsequent step requires strictness of this multiplied skew
bound. The triangle-inequality/Cauchy shared-budget step uses only
the non-strict per-face bound and gives the non-strict expressions
in Eq. (15). The error-function maximization likewise requires only
that upper bound. Strict final fit errors follow from the distinct
already checked scalar inequalities

    alphaB2 D0 < (17/50)^2,
    alphab2 D0 < (3/10)^2,

whose respective positive margins are

    6608408949/20000000000000,
    44858500191/20000000000000.

Consequently `|B-B*|<=sqrt(alphaB2 D0)<.34` and
`|b-b*|<=sqrt(alphab2 D0)<.30` still hold even when every skew is
zero. The unit-anchor projection and final strict contradiction
therefore retain their stated strength, including at `D=D0`.

The fixed-radius local actual signed-A `.038` acceptance remains
valid with this endpoint notation qualification. No scalar value,
certificate, dependency, or global scope changes. The existing
frozen audit manifest remains a record of the original frozen
review/receipt files; this subsequent clarification is a separate
addendum and is not retroactively inserted into that manifest.
