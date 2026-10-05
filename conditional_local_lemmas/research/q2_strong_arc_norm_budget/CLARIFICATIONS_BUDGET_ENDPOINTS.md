# Budget endpoints and peak constants for the frozen zero-tau proof

The main proof remains byte-identical at SHA
`c0077a0fa8ca16bf56c4c13a93571d048082c11ba60f1da4c4f7a6b3ffa549a1`.
The first clarification remains byte-identical at SHA
`16d4865e47b8878d529edccc119627421319308b8c1dec6e869e191506ab5bf1`.

1. In Section 4, `eta-eta'<v/6` is strict for v>0. Its uniform
   endpoint form is `eta-eta'<=v/6`; at v=0, eta'=eta and both
   sides vanish. Thus the uniform relaxation is

       chi0(s,eta)-v<=chi0(s,eta')-5v/6<=chi0(s,eta').

   The final strict bound in Section 5 comes from the strict
   starting-angle, Cauchy/rational square and cosine inequalities,
   and does not require strictness of this zero-v endpoint.
2. The factor `1+2cos(theta*)<2.354` in Section 5 uses the
   **upper** peak-cosine enclosure `cos(theta*)<.677` from the
   accepted certificate. The sine-halfangle bound `<.403` uses
   its lower enclosure `cos(theta*)>.676`. Both bounds are in
   the bundled, SHA-pinned accepted certificate; they are not
   inferred from the lower bound alone.
3. A sufficient explicit rational guard for `sin(.75)>.68` is
   the lower seventh-degree alternating Taylor polynomial

       .75-(.75)^3/6+(.75)^5/120-(.75)^7/5040>.68.

   The cubic lower polynomial by itself is below .68 and is
   not used for that guard.

This is a document-only clarification. Mathematical program
executions, imports, old checker calls, parameter scans and
optimizer calls remain zero. Neither the frozen main proof, its
original manifest nor the first clarification was edited.
