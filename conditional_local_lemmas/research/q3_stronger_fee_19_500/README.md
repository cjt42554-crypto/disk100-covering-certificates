# Actual signed-anchor fan fee 19/500 at a fixed radius

Exact parameters: R=11321/100000, lambda=7/50,
t=1071883/20000000. The proved local score is

    sum_i[area_i/R^2 + A_i/R^2 -4lambda*x_i^2 +t*gamma_i]
                       < 3 C_lambda(t) -19/500,

for three consecutive consistently oriented positive-area faces sharing
an actual apex, each circumradius <=R, with arbitrary common covered
unit-circle anchor q_i on each base. A_i is signed; both sides are
allowed. The closure argument includes degenerate limits.

This is fixed R only. It does not cap the independent J score, prove a
new global lower radius, spend lens/BI credits, prove all fan incidences,
or resolve the exact optimum. A global substitution requires the same
actual signed-A identity and disjoint eligible charged face/anchor sets.

Replay from any copied location with Python 3:

    python3 verify_actual_quadratic_19_500.py

The checker uses exact fractions and outward Taylor intervals. It reads
the bundled accepted scalar-conjugate certificate, verifies coefficient
and SHA pins, and imports only interval/Taylor primitive definitions.
It never calls the accepted conjugate checker or primitive verify().
The author execution count and final frozen SHA manifest are in
ATTEMPT_LEDGER.json and SHA256.json. Independent review is separate.
