PUBLICATION SCOPE: CURATED HISTORICAL RESEARCH RECORD
This is a lean, privacy-edited copy. Original replay/rebuild instructions below describe the complete historical archive, not a ready-to-run public reproduction. Large run snapshots, inventories, caches and some dependency trees are intentionally omitted. No complete replay of these finite runs is claimed here; no program was executed for publication. Read ../README.md and this route's PUBLICATION_SCOPE.json. Original files were retained unchanged.

# Independent review of the five-site projection verifier

**PASS for the corrected 7,173-byte source**, SHA256
`63a9c0f0379796c941dd068a7aa2c714a335143910b9c44300110c8322bdf6b4`.
The review did not execute a full-cover acceptor, optimizer or old proof.
Twelve new arithmetic/input-binding checks passed. The original positive toy
and nearer-blocker execution receipts are reused as original evidence; their
source version predates the two input-binding corrections.

## Corrections and reproducible counterexample

The original 6,798-byte source has reconstructed SHA256
`5d9cb7906d916a43f1808aa4ec1a0067a41f8f377174be44421d880d57c8f528`.
Its snapshots are retained with explicit reconstruction labels. The isolated
shape-guard version has SHA256
`397872f80ca393d06cc5299f52786e123a69b624e0d69d432103a697393c448a`.
Reversing the two exact source edits reproduces both sizes and the separately
observed intermediate SHA. The final root source remains unchanged by review.

The original parser accepted 100 JSON rows with three coordinates each. The
projection used the first two coordinates, while the pinned covering acceptor
flattened all coordinates and paired them into a different set of 150 centers.
The saved `malformed_100_triples.json` gives five intended points on the circle
of radius 3 about the origin and 95 points at (4,4); all intended centers miss
the origin for claimed radius 1. The flattened data nevertheless contains an
origin center, whose radius-1 disk alone covers the unit disk. Both added centers
already lie on the primary circle, so the displacement bound is zero. This was
a real source/geometry binding defect, not a failure of the exact cover acceptor.

The corrected verifier rejects this input before any call to `acceptor.verify`.
An intercepted call counter proves the rejection order without running that
acceptor. The same checks reject duplicate or out-of-range selected indices,
an intervening nearer site, and a coordinate file changing during parser read.

The original reader's digest was of newline-normalized text. A CRLF file could
therefore have a different raw-byte SHA while receiving the original LF digest.
The corrected verifier hashes its byte snapshot, checks unchanged bytes after
parsing, and separately records the reader's normalized-text digest. The CRLF
IO test verifies both distinct values, using identical decoded rows and a
mocked return of the *existing* positive source-cover result. This IO test is
not a new mathematical covering execution or a new accepted cover certificate.

## Mathematical checks

For primary points a,b,c, the determinant `2 cross(b-a,c-a)` is nonzero.
The displayed rational q is the solution to
`2(b-a)·q=|b|²-|a|²` and `2(c-a)·q=|c|²-|a|²`. The verifier independently
checks all three squared distances equal the positive rational r².
The translated-circle arithmetic check confirms both coordinate formulas.

Each added center is defined exactly by
`q + sqrt(r²/|p-q|²) (p-q)`, with the positive square-root branch.
The rational scale square is positive and its product with the ray norm square
equals r². Thus all five projected sites lie on the same nonzero circle.
Their positive scales preserve ray order. The exact half-plane/cross-product
sort rejects equal rays; every successive cyclic cross product must be positive.
For points on a common circle this gives five distinct convex vertices, with q
strictly on the inner side of all five edges. The verifier also proves |q|²<1.
Same-ray collisions and a polygon-boundary q are rejected by edge tests.

The other 95 centers are unchanged. Every exact squared distance to q exceeds
r², so the five projected sites are the nearest sites at their common circle
center. Duplicate unrelated centers cause no issue; the strict gap is checked
against every one of the other 95 labeled entries.

`sqrt_interval(x)` uses integer floor square roots to give rational bounds
lo²≤x≤hi² with nonnegative lo. Zero, exact squares, nonsquares and a large rational
are independently checked. If R=√r² and d=|p-q|, the radial displacement is exactly
|R-d|. From R∈[Rlo,Rhi], d∈[dlo,dhi], the bound
`max(Rhi-dlo, dhi-Rlo, 0)` is valid; an exact rectangle/endpoint check confirms it.
Taking the maximum over the two moved centers also bounds every labeled center's
displacement, since the other 98 centers (three primary plus the other 95) do not move.

The SHA-pinned acceptor is the already accepted integer/surd source. With exactly
two coordinates per row, its finite decimal input and the verifier's Fractions
represent the same rational planar centers. A successful original coverage
acceptance at U implies for every x in the closed unit disk that some original
center has distance ≤U. Its projected counterpart has distance ≤U+D by the
triangle inequality. This proves covering of the *whole closed disk*, without
sampling the disk and without requiring centers themselves to stay inside it.

## Scope and limits

The output specifies one exact algebraic five-site tie configuration. The two
positive square roots generate a field of degree at most four; it may have
smaller degree in special cases. The reported rational U+D and its linear
integer polynomial describe a proved covering **cap**, not the exact covering
radius or an optimal-radius polynomial. Floating approximations are reporting
only and do not enter acceptance.

There is no certified continuous path from the input to the projected configuration,
no first-event classification, KKT/stationary or second-order proof, record
improvement claim, global optimality proof, or claim that arbitrary configurations
must have this geometry. A specific real disk100 projection requires its own
actual source-cover execution and bound output; this independent code review
does not fabricate that execution.

Reproduce the arithmetic and input-order checks only:

```sh
python3 -S -B five_site_projection_certificate/independent_review/verify_review_edges.py
```

`EDGE_CHECK_REPORT.json` records the zero full-cover execution count and the
intercepted/mock gate uses. `dimension_counterexample.json` explains the prior
binding failure with exact arithmetic. Root source and original regressions
were not modified by this review.
