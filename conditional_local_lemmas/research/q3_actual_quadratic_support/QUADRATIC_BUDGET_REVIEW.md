# Independent algebraic review of the quadratic replacement interface

The root-authored `ROOT_QUADRATIC_BUDGET_INTERFACE.md` passes the symbolic
signed-area and paid-summand review. This review does not establish its
external boundary support, ordinary reduction, global type window or
numerical radius exclusion.

For a CCW center edge, expanding
`A_i=-det(c_(i+1)-c_i,q_i-c_i)/2` gives
`[det(c_i,q_i)+det(q_i,c_(i+1))+det(c_(i+1),c_i)]/2`.
Summing these terms cancels the polygon shoelace determinants and the
opposite determinants in the assigned-cell lens sum. The exact surviving
term is the halfarc sum `pi`. This confirms the anchor sign, including
negative individual anchors.

For the complete ordinary marked-face partition, put
`s_i=area_i/R^2+A_i/R^2-lambda*(l_i/R)^2+t_i*gamma_i`.
Direct substitution gives

```
pi = sum(unmarked actual area)
       +R^2*sum s_i-R^2*sum t_i*gamma_i
       +sum L_i+lambda*sum l_i^2.
```

The remaining quadratic boundary functional is exactly
`Q_actual=sum L_i+lambda*sum l_i^2`. The single edge-anchor loss
`tau_i=x_i sqrt(1-x_i^2)-A_i/R^2` occurs inside `s_i`; it is absent
from `Q_actual` as defined. Replacing three actual scores at the same
coefficients by `3C-1/40` can therefore retain a valid bound for this
`Q_actual` once. Every displacement or variance already used within that
bound stays within it and cannot be charged separately. An implementation
whose boundary functional also consumes edge-anchor loss would require
its own new compatibility audit; the root interface explicitly records
this restriction.

The new local theorem uses `4lambda=14/25`, the base convention
`l=2Rx`, and the exact common angle coefficient
`t=1071883/20000000`. Thus both before and after replacement the identical
quadratic base and apex-potential terms are present. Distinct eligible
interior apexes allocate disjoint marked faces and their unique boundary
base anchors; sharing endpoint centers allocates no extra paid term.
The original q0/BI, unmarked corner and centered terms do not appear in
the new payment. The exact comparison of the local fee with the stated
target is

```
1/40-19931551/1000000000=5068449/1000000000>0.
```

This replaces an allowance for actual signed scores, not the independently
maximized J score. The latter's standalone required fee remains blocked
by the separate explicit local diagram. No anchor-to-lens bridge factor
is needed for the present substitution, and no such bridge is silently
assumed. The quoted reflected boundary support is an external whole-cover
theorem; this algebraic review does not assign the undefined parent
outer-Q parameter to it or transfer any old-radius numerical boundary
constant. The interface correctly preserves these remaining obligations.
