# Preserved clarification to the mean-envelope comparison

The frozen `MEAN_BUDGET_ENVELOPE.md` SHA is
`4a4e634053549ff56019d93802cfb10da56ae21bbf23b03cc6822b6fc97a0a8c`.
Its equation (4) displays a strict comparison of square roots using
the strictly larger coefficient `47/100`. At `H=B-Delta^2=0`, the
two square roots are both zero, so that coefficient comparison itself
must be non-strict. The uniform intended chain is

    chi0(s,eta+)
       < .523+sqrt[(.6^2/4+1.002^2/(1666/625))H]
       <=.523+sqrt[(47/100)H].

The first inequality is strict at every H, including zero, because
the accepted peak starting angle is strictly less than `.523`.
The second is strict when H is positive and equality when H is zero.
No later bound needs strictness from a zero square-root comparison.
The final envelope, derivative payment and conditional exclusion
are unaffected. Frozen source and its existing pins remain unchanged.

This clarification is a paper correction. Mathematical program
execution, mathematical imports, checker replay and search: 0.
