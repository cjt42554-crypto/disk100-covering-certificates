PUBLICATION SCOPE: CURATED HISTORICAL RESEARCH RECORD
This is a lean, privacy-edited copy. Original replay/rebuild instructions below describe the complete historical archive, not a ready-to-run public reproduction. Large run snapshots, inventories, caches and some dependency trees are intentionally omitted. No complete replay of these finite runs is claimed here; no program was executed for publication. Read ../README.md and this route's PUBLICATION_SCOPE.json. Original files were retained unchanged.

CERTIFIED RATIONAL PIECEWISE-LINEAR PATH THROUGH A DELAUNAY FLIP

Independent standard-library replay (no NumPy, SciPy, samples, or optimizers):
  python3 -S certified_event_path/verify_path.py
  python3 -S -O certified_event_path/verify_path.py
Mathematical tamper tests (SHA checking is bypassed to test actual predicates):
  python3 -S -O certified_event_path/negative_tests.py

Rebuild from the read-only historical snapshots (NumPy needed for NPZ decoding
and proposing witnesses; all acceptance uses integer/Fraction/surd predicates):
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    python3 certified_event_path/build_source.py --full

After a rebuild regenerate MANIFEST_SHA256.json before independent replay.
The manifest binds file integrity; mathematical acceptance is recomputed.

Data:
  path.json: ordered 23 rational nodes and 22 exact continuous segment caps.
  nodes/-1: the old certified rational U=0.1155018241946 configuration.
  nodes/00 through nodes/21: exact Q=10^15 rational rounding of accepted nodes.
  nodes/*/upper_certificate.json: one covering index per exhaustive candidate.
  event_geometry.json: every integer coefficient of the quartic event and
    circumcenter/common-radius rational functions; isolating rational box;
    exact interval evidence for convexity, inside-D, strict emptiness, and
    q inside the convex quadrilateral; strict Delaunay triangles on both sides.
  short_segment.json: exact cap for 018_event -> 019_event_after only.
  PROOF.txt: the whole-domain reduction and continuous-path/event arguments.

Certified whole-path cap: 235173307525027/2000000000000000 = 0.1175866537625135.
Certified event short-segment cap: 231185886538821/2000000000000000
  = 0.1155929432694105.

This is a new explicit rational PL path. It certifies neither the original
numerical SLSQP arc nor a true stationary-branch supremum. It is not a proof
of a new record covering upper bound, a lower bound, or global optimality.
