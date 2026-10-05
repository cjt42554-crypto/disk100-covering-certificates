# Order of the skew localization used in the full actual-tau proof

The frozen full proof imports the accepted paid per-face localization.
Before using `sin^2 d0+sin^2 d1>=2sin^2 eta` to improve the eta
bound, one first uses that same accepted per-face estimate:

    sin^2 d_i <=(625/512)sigma_i
                 <=(625/512)D<.028.

Positive base angles select `|d_i|<pi/2`. Since
`sin(.18)>.17` and `.17^2>.028`, this forces `|d_i|<.18`.
Thus `|eta|<.18<pi/4`, so `cos(2eta)>0` and the identity

    sin^2(m+eta)+sin^2(m-eta)
       =1-cos(2m)cos(2eta)>=2sin^2 eta

is valid in that imported localized range. Only then does the
full proof's single-budget guard improve the conclusion to
`|eta|<.095`. This ordering avoids deriving a small eta from
an identity used outside its positive-cosine range.

The underlying bounds are already explicit in Section 1 of the
accepted zero-tau source. This note makes the imported step's
ordering explicit; no source or earlier manifest is modified.
Mathematical program execution, imports, checker replay, scan and
search: 0.
