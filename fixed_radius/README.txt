Exact fixed-configuration radius and rigorous nonglobal comparison
================================================================
This package certifies the TRUE covering radius of the frozen grid25 rational
configuration: 100 centers, common denominator 10^15, configuration SHA256
84fb2a5dd9869e641f18250b8839ceb44bfe1263fa54bf588cfdb4867fb28a36.
The unique worst point is the rational circumcenter T(24,26,33), with exactly
those three nearest centers. Its radius is the unique positive root of the
primitive irreducible quadratic printed in full in polynomial.txt and
exact_result.json. Its adjacent 10^-60 isolation numerators are
115446954885582804739528389204068749175480606343527731854843,
115446954885582804739528389204068749175480606343527731854844.

Self-contained fixed-radius acceptance, from any working directory:
  python3 -S -O -B ./verify_exact.py
All 138437 conservative critical candidates are regenerated exactly.
At the unique exception, all 100 center distances are rigorously compared
and exactly three ties verified. All 138436 OTHER candidates have a witness
with squared distance <= 6663999696176893/500000000000000000, strictly
below the exception's exact squared radius. The proof covers the WHOLE
CLOSED unit disk and establishes the unique worst point. No sampled-grid
or numerical ranking is trusted. Original critical_candidates.py is pinned
and byte-identical to its saved source. The exact acceptance needs no optimizer or numerical proposal generator.

Self-contained strict nonglobal acceptance:
  python3 -S -O -B ./verify_nonglobal.py
This independently regenerates the fixed proof and the bundled 932-variable
companion's exact root + full closed-disk covering proof. It checks system
semantics, exact Krawczyk inclusion/contraction and rigorous interval
triangulation/boundary coverage. The freshly recomputed companion radius
upper endpoint is strictly below the fixed radius lower endpoint. It NEVER
accepts a stored success flag. The snapshot pins all copied dependencies
and binds to the separately frozen original companion manifest
4bc4d4647c5b5e8e6e7bb01d1ce10f50fae1823139ce0333efffa15f329472e2.
No files from the original project are read or written on either replay.

P is the exact minimal polynomial of this FIXED rational configuration's
radius, NOT a polynomial for the unknown arbitrary-configuration global
optimum. This fixed configuration is rigorously strictly nonglobal because
the companion is a 100-center full closed-disk cover with smaller radius.
No global optimum claim is made for the companion either.

Actual -S -O logs/receipts cover fixed and nonglobal acceptance. Three
meaningful corruptions must be rejected: polynomial constant changed,
exact squared-radius numerator changed, and root box moved two grid units.
The sealed ZIP is extracted, every manifest hash checked, then BOTH
acceptance entry points rerun on the unpacked files.
The companion verifier regenerates rational_midpoint.json and two report
JSON files inside companion_snapshot/; these derived outputs are excluded
from the source manifest/archive, and their stored flags are never inputs.
All previous rational-upper and fixed-radius packages remain frozen.

PUBLIC COPY: this is a static dependency-closed publication of historical
accepted evidence. No new mathematical replay occurred during copying.
Old packing/ranking tools are omitted. Personal launcher locations in
historical receipts are redacted; their actual mathematical outputs remain.
