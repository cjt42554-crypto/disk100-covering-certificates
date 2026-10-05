Independent pure integer/surd strict covering margin
===================================================
V = 1154469548855829/10^16 = 0.1154469548855829.
U = 115446954885583/10^15 = 0.115446954885583.
The unchanged frozen 100-center rational configuration covers the whole
CLOSED unit disk at V, and V<U with exact margin U-V=1/10^16.
Consequently its covering radius and the arbitrary-configuration optimal
radius are strictly below U. This is a separate local strict proof; it
leaves the independent cloud result's non-strict <= statement unchanged.

Run from any working directory, using only the Python standard library:
  python3 -S -O -B /path/to/strict_upper_margin_1154469/verify_margin.py
The source configuration is byte-identical to the frozen grid25 Q15 JSON,
SHA256 84fb2a5dd9869e641f18250b8839ceb44bfe1263fa54bf588cfdb4867fb28a36.
Its old radius_bound is provenance metadata and is ignored for acceptance.
margin_ledger.json supplies the explicit new V; cover_certificate.json
must encode exactly V squared. The checker exhaustively regenerates
138437 candidates (T=128436/B=9900/A=100/E=1) and verifies every witness
with exact integers or one quadratic surd. No candidate count is assumed.
It then recomputes V<U by a pure integer cross-product comparison.

PROOF.txt and savedsource/CLOSED_DOMAIN_PROOF.txt establish the full-domain
extremum reduction. critical_candidates.py is byte-identical to the original
source and cryptographically pinned. savedsource/witness_proposal.json
provides proposed integer covering indices only; the old squared bound or
any stored success flag is never trusted. build_margin.py uses no float,
Decimal, ranking, KKT, optimizer or algebraic companion.

Actual -S -O logs/receipts and three meaningful corruption rejections are
included: truncating the witnesses, decreasing the radius by 2e-16, and
setting the reference bound equal to V (destroying strictness). The ZIP is
unpacked, every manifest hash checked, then the portable entry point is
actually rerun on those files. The package is under 1 MB and has no
acceptance dependency outside its own directory.

This margin uses the previously frozen geometry. It is not a material new
geometry, a proof of optimality or an optimum polynomial. No earlier
certificate, configuration, frozen ZIP or Library item is modified.
