PUBLICATION SCOPE: CURATED HISTORICAL RESEARCH RECORD
This is a lean, privacy-edited copy. Original replay/rebuild instructions below describe the complete historical archive, not a ready-to-run public reproduction. Large run snapshots, inventories, caches and some dependency trees are intentionally omitted. No complete replay of these finite runs is claimed here; no program was executed for publication. Read ../README.md and this route's PUBLICATION_SCOPE.json. Original files were retained unchanged.

Exact local-to-official B32 catalog identity audit (2026-10-04)

Result: PASS. Local independent_index i uniquely matches the supplied official
B32_P26_F108_{i:03d}, for every i=0,...,15. Each selected witness is the identity
D6 pose, translation (0,0), and outer-cycle direction +1, shift 0. There are 25
valid pose witnesses in total because some graphs have automorphisms; each local
row still has exactly one matching official family.

Current: PASS, independent_index=10, official family B32_P26_F108_010.
Fingerprint: 8764ca79e182ece49aaf63c485289f9030d4698480a93a558f58d21e322aa1a7.
The current-to-local and current-to-official oriented face maps reverse every
triangle, and the current outer cycle is reversed. This orientation belongs to
the explicit vertex bijection and is independently checked. It is not inferred
from the numerical canonicalizer's reflection flag.

Replay (Python3 standard library only):
  python3 verify_official_mapping.py
Run from this directory or give an absolute script path. All required inputs
are under inputs/. The replay has no network calls and does not use numpy,
scipy, optimizers, external projects, configuration files, or archive extraction.
--stage-sources was used once to read selected exact ZIP members and frozen
local sources. Do not use it for ordinary replay; it refuses existing inputs.

Inputs and provenance:
- Official ZIP: disk100_lower_bound_0p11317_20261004.zip, 1550571 bytes,
  SHA256 e7220866aab85d8c457dae49b79b4436bd647179d4311877766b345a2cedc510.
- Local frozen catalog SHA256:
  f3234fd61f521362c69e3442940edc1d381fb5c2e06b1cd12a79e93ab2e57bb4.
- The 18 official member JSON files retain the original ZIP-member bytes.
  SOURCE_SHA256.json records every original member path, byte count and SHA256.
  The official graph index's catalog SHA is verified against the original bytes.
- The 4 local inputs retain their original frozen bytes. current.npy is hashed
  against current_graph_numerical.json's source_sha256; no DT extraction runs.

Checks:
- Pure integer D6 poses: reflection (x,y)->(y,x), rotation (x,y)->(-y,x+y),
  followed by integer translation. For equal finite site sets, translation is
  fixed by the lexicographic minima, so every one of the 12 D6 poses is covered.
- Full 68-site lattice equality and 108 positively oriented unit core triangles.
- Full 100-vertex bijection, 265-edge bijection, and 166-face bijection preserving
  orientation after the explicit reflection compensation and each category:
  108 core_q0 + 26 ring_q0 + 32 q1. q1 apex and fan metadata are preserved.
- Full outer32 and core26 cycles with the same orientation parity as the pose.
- Disk Euler count, opposite internal oriented edges, unique faces, connected
  face traversal, every vertex link, all 68 interior degrees six, official
  index/graph metadata, and independently replayed canonical fingerprints.
- Current full100 and full166 maps, 265 edges, outer32 cycle, and all68 integer
  development records are verified and composed into the official labels.

Deliverables:
- official_catalog_mapping.json: all16 selected actual vertex, face and edge
  bijections; alternate matching poses; rederived graph invariants.
- current_official_identity.json: actual current->official 100-vertex and
  166-face bijections, source coordinate SHA, orientation and family evidence.
- all256_pair_attempts.json: every local/official pair, match counts, and reasons
  for rejected D6 poses or outer shifts. The 240 nonmatching family pairs are
  expected identity-search rejections, not missing or unresolved mappings.
- failure_records.json: execution failures or ambiguous/unmatched rows.
- verification.log and portable_replay.json: original and isolated replay PASS.
- MANIFEST_SHA256.json and READY.json: sealed local file inventory.

Scope: this proves identity of these two supplied 16-member catalogs and the
stored numerical current graph's official family. It does not prove that 16
families exhaust configurations at arbitrary radius, prove optimality, change
any radius bound, or produce a fresh exact geometric-Delaunay certificate.
The official source's conditional completeness claim is not re-proved here.
Source-region and deletion-history metadata need not agree across independent
reconstructions; only the exact developed mesh and attached graph are matched.

No optimizer was run. No configuration or original material was modified.
