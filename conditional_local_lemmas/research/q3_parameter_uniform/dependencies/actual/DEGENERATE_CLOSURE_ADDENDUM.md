# An explicit bound for every zero-area limiting face

This supplements the closure explanation without changing the frozen
positive-area theorem or its checker. Degenerate limiting faces also
include a repeated apex and base endpoint; one need not classify their
undefined geometric angle directly.

Take a sequence of admitted positive-area faces converging to zero actual
area and pass to a subsequence on which its apex angles converge in
`[0,pi]`. Endpoint disk coverage always gives `0<=x<=1` and
`A/R^2<=x sqrt(1-x^2)<=1/2`. Since the quadratic term is nonpositive
and `t>0`, the limiting actual signed score satisfies

```
limsup score <= 1/2+t*pi
              <1/2+(1071883/20000000)*(22/7)
              =46790713/70000000 <7/10.
```

The accepted conjugate certificate gives `C>147/100`, and hence
`C-1/40>289/200>7/10`. Thus no zero-area limiting face can occur
in a triple whose total score is at least `3C-1/40`: its individual
score would have to be at least `C-1/40`, because every other face is
bounded by `C`. This excludes all zero-area closure cases, including
those beyond the explicitly discussed `x=0` and `theta=pi/2` endpoints.

If a limiting triple has no zero-area face, its geometric angle and
skew parameters are well defined and the frozen proof applies directly.
Its final rational unit-condition contradiction is strict, so it also
excludes equality at the stated threshold. This is a local closure
observation; it supplies no whole-cover nongeneric reduction or global
incidence theorem.
