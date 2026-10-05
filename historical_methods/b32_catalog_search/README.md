PUBLICATION SCOPE: CURATED HISTORICAL RESEARCH RECORD
This is a lean, privacy-edited copy. Original replay/rebuild instructions below describe the complete historical archive, not a ready-to-run public reproduction. Large run snapshots, inventories, caches and some dependency trees are intentionally omitted. No complete replay of these finite runs is claimed here; no program was executed for publication. Read ../README.md and this route's PUBLICATION_SCOPE.json. Original files were retained unchanged.

# B32 catalog identity, historical graph audit, and harmonic metric seeds

All work in this directory is an independent graph/seed investigation. No
covering optimizer was run here, no earlier frozen source was changed, and no
global bound or optimality claim is made.

The archived project README and conditional-6049 README were read first. Their
existing machine catalogs contain N26/F106 and N27/F105 cores, rather than the
new B32/P26/F108 catalog reported by the parent. The latter machine package
remains pending. `reconstruct_catalog.py` uses the older supplied containing
regions and independently obtains 16 N68/F108/P26 core classes, then appends a
standard B32 annulus to each. These are independent indices 0..15; completeness
and official family names are not asserted by this reconstruction.

The incumbent's exact saved binary64 source is
`../best147_boundary_transfer/frozen_proposals/grid0025_R0.11544695488558258.npy`,
SHA256 `de82a12038c6f9a705b7a985c7e25c53c8e94f5c506946aba9a6c6f03a02852d`.
Its boundary/reflection/face canonical fingerprint is
`8764ca79e182ece49aaf63c485289f9030d4698480a93a558f58d21e322aa1a7`.
It matches exactly one independently reconstructed class, index 10. The parent
reports official ID `B32_P26_F108_010`; the explicit official-catalog mapping is
pending and is not inferred solely from the matching suffix.

`exact_current_dt.py` independently validates the proposed Delaunay complex
using integers after converting binary64 centers into exact rationals with
common denominator 2^65. All 178 full-DT triangles have positive orientation,
and all 17,266 other-site empty-circle signs pass, with no extra-site ties.
Every circumcenter has a strict unit-disk sign, leaving 166 restricted faces.
The 32 boundary dual segments/rays cross the unit circle, while the other
outside/outside segment and 11 outside hull rays provably miss the disk.
`current_graph_integer_verified.json` records this graph identity certificate;
it does not certify a radius.

The resulting 100-vertex face bijection develops the 68-vertex core into integer
triangular-lattice coordinates. `current_core_integer_development.json` checks
all 166 full-graph faces and all 175 core edges. The core has 108 faces, boundary
length 26, and 42 true core-interior vertices. Its boundary degrees are eight
3s, sixteen 4s and two 5s; turns are eight +60 degrees, sixteen zero turns and
two -60 degrees, totaling 360 degrees. The two concave boundary vertices are
actual center indices 31 and 54; both have marked fan size zero.

The read-only historical audit scans retained coordinate fields from task-3
and task-2, including direct NPY, NPZ fields, stacked paths, numerical JSON and
rational JSON coordinates. A quick Delaunay/profile gate precedes complete
circle-arc validation and canonical face matching. The sealed local
inventory contains 44,995 coordinate fields and 34,868 distinct ordered float64
coordinate arrays; 484 arrays match the B32/F166/68-interior-degree-six profile.
Four reconstructed graph classes occur:

| Independent index | Distinct retained arrays | Joined saved final records |
|---:|---:|---:|
| 8 | 5 | 2 |
| 10 | 426 | 114 |
| 11 | 38 | 16 |
| 13 | 15 | 6 |

All 484 pass the numerical face-boundary/circle-arc agreement, with no singleton
and no graph extraction failure. None of the structural-search 9,840 start or
final records matches these four B32 graph classes. A graph appearing in a
start, stage, path or saved final does not establish a stationary metric or
local minimum. The other 12 graph classes have no match in the scanned retained
inventory, which does not prove they were never explored outside this inventory.
The separate history snapshot records all coordinate/source hashes, role/source
counts, trial accounting, saved final joins and archive-member checks. Its
manifest SHA256 is
`4564114ea32e4853656cc97cba860cb42e81bd64b854f6ccb918ab1e7a109a4d`.
The final ZIP comparison scans 56 archives, 41 unique archive bytes, and finds
all 24,884 archived coordinate hashes already present in the local inventory.
Include `history_inventory/frozen_snapshot_errata.json` when sharing the
snapshot: it marks two early ZIP comparisons as superseded by the final direct
hash check. The immutable snapshot itself was not altered.

The immutable seed handoff is `frozen_handoff/READY.json` with manifest SHA256
`e72c23ba91a7909a65faf5cf2eea1f541bad24bb3324c677f1d2305f36cd952e`.
It includes the 16 graph descriptions, a canonical helper and 96 harmonic
100-center proposals: a radius-.94 outer circle with uniform or core-radial
angles, and core-edge conductance .7, 1 or 1.4. All proposals lie strictly in
the unit disk and have 166 positively oriented target faces. Actual restricted
Delaunay faces equal the target in 93 proposals. Three index-0 uniform-angle
proposals differ by eight faces and are explicitly marked as mismatches; the
three index-0 radial proposals match. These mismatches must not count as target
graph exploration.

`exact_seed_face_audit.py` separately checks every sealed target triangle using
integer orientation, emptiness, inside-disk and center-in-disk signs. Exactly
the same 93 seeds pass strictly, with no additional exact ties or failures among
numerically matched seeds. The three mismatches remain recorded. This is seed
graph validation, not a certified covering bound. The seed audit does not
separately reconstruct full Voronoi dual segments/rays or all circle-arc signs;
those numerical graph checks remain in the sealed ledger. The current-incumbent
audit additionally checks the full 178-face DT and exact dual disk crossings.

The sibling `b32_metric_seed_search/` owns subsequent metric realization and
the authorized actual 100-center covering optimization. It received the sealed
handoff and prioritizes the 12 graph classes absent from retained history, then
the already retained indices 8, 11 and 13. Current index 10 is excluded from
that seed batch. No claim here forces the global optimum into these 16 graphs.

Replay graph identity and seed integer signs without modifying frozen inputs:

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python3 -B b32_catalog_search/exact_current_dt.py
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python3 -B b32_catalog_search/exact_seed_face_audit.py
```

The two scripts update their non-frozen result JSON files in this directory.
Supply a first positional output-JSON path to replay without modifying their
saved reports, including when executing an archived frozen copy.
`reconstruct_catalog.py` and `development_and_seeds.py` are discovery generators;
do not invoke the one-shot `seal_handoff.py` again or overwrite the handoff.
