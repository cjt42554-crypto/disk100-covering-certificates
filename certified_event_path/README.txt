CERTIFIED RATIONAL PIECEWISE-LINEAR PATH THROUGH A DELAUNAY FLIP

Independent standard-library replay (no NumPy, SciPy, samples, or optimizers):
  python3 -S certified_event_path/verify_path.py
  python3 -S -O certified_event_path/verify_path.py
Mathematical tamper tests (SHA checking is bypassed to test actual predicates):
  python3 -S -O certified_event_path/negative_tests.py

The numerical discovery snapshots and their builder are omitted from this
clean release. Exact rational nodes, covering witnesses, event data and
the independent acceptance checker are included. Saved PASS receipts
are inherited original evidence; publication made no mathematical replay.

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
