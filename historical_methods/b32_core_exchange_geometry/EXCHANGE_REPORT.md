PUBLICATION SCOPE: CURATED HISTORICAL RESEARCH RECORD
This is a lean, privacy-edited copy. Original replay/rebuild instructions below describe the complete historical archive, not a ready-to-run public reproduction. Large run snapshots, inventories, caches and some dependency trees are intentionally omitted. No complete replay of these finite runs is claimed here; no program was executed for publication. Read ../README.md and this route's PUBLICATION_SCOPE.json. Original files were retained unchanged.

# Exact core-site exchange and metric starts

The supplied 16 independent B32 catalog cores were compared by all 12 triangular-lattice D6 elements and all point-difference translations. Two independent enumerations agree on all 256 ordered pairs and all 949 maximum-overlap poses, including every reciprocal inverse pose. Distances count removed sites, so their symmetric differences contain twice as many sites.

All indices here are the zero-based independent catalog indices with catalog SHA256 `f3234fd61f521362c69e3442940edc1d381fb5c2e06b1cd12a79e93ab2e57bb4`. Official parent catalog IDs can be attached by the separately delegated exact graph mapping; they are not inferred from metric proximity or radius.

Baseline index10 has exchange distances `[3,5,4,3,3,3,4,1,1,1,0,1,1,1,1,1]` to indices0..15. There is no distance2 from10. The full matrix, every optimal integer pose, common-site label maps, and changed triangles/edges are in `overlap_analysis/`.

| Directed motif | Target-to-source lattice pose | Removed source site | Added posed target site | Source wall labels |
|---|---|---|---|---|
| 10→14 | identity | label5, (0,9) | target38, (4,9) | 89,90,91,94,95,96 |
| 10→15 | R S + (8,-4) | label0, (0,4) | target38, (4,9) | 68,69,70,89,90,91 |
| 14→15 | identity | label63, (8,0) | target67, (8,5) | 78,79,80,84,85,86 |

Here R(q,r)=(-r,q+r), S(q,r)=(r,q). The outer cyclic maps are respectively shift0/direction+1, shift12/direction−1, and shift0/direction+1. Each motif retains67 core sites,106 core faces,172 core edges,22 core boundary edges,52 cross edges,32 outer-cycle edges and50 annulus faces. The full100 label map is bijective; 14↔15 maps are mutual inverses.

The one-site graph on indices7..15 is K7 on {7,8,9,10,11,14,15} united with K4 on {10,11,12,13}, sharing edge10–11. It has26 edges and diameter2. The larger one-site components are {0}, {1}, {2,3,4}, {5,6} and {7..15}. A path from14 or15 to12 or13 needs two site exchanges within this catalog. Twenty-five of26 edges use six affected outer vertices; 7↔9 uses five because its wall neighborhoods overlap. This finite graph statement is not a global completeness assertion.

Metric proposals preserve the actual coordinates of all67 common core centers. They replace one core center and optionally move only the six outer centers identified by exact incidence differences. The remaining99 or93 actual centers are bitwise equal to their mapped source rows. New points come from exact integer affine coefficients on actual source triangles, or from a labeled67-center similarity fit of independently accepted pinned14/15 frames. Any individual unit-disk projection is recorded. Reflection is carried in face orientation metadata; raw coordinates are not silently reflected to pass a check.

Three immutable handoffs contain19 distinct ordered and label-unordered metric starts:

| Handoff | Starts | Actual initial graphs | Integer graph support |
|---|---:|---|---|
| first_handoff | 6 | 3 index14, 3 index15 | Strict166 target face orientations, empty circles and inside circumcenters; actual labeled graph and circle arcs additionally numerical |
| second_handoff | 3 | Three unexpected B32/F166 graphs, outside the supplied16 all-degree6 classes | Complete178-triangle full DT,166 inside triangles,32 exact dual disk crossings, and all outside dual segments/rays miss the disk |
| third_handoff | 10 | 5 index10, 2 index15, 3 index14 | Same complete integer actual graph and dual disk identity as second_handoff |

The second-handoff fingerprints are `9e2baacbb9f1c1b3e9792d900791e83f6f8f992b73bc26c2f2e3a7bb42ac6429`, `b83a8447d1339bd09d32f29a4a03f9d845e39c1b4aaa1b7cbb0d1d8f33368cde`, and `2336a49d758f6fd5f05545cb05d65af6209f53dd91b3a6dfc987e55622941a7a`. Their intended targets14/10/15 are explicitly mismatches. All three preserve99 source centers and have exact unit-disk/face/empty-circle/dual support. Their interior degree profiles differ from the sixteen all-degree6 B32 graphs.

All19 starts passed the independent numerical oracle using161700 triples,4950 unordered pairs (9900 circle candidates),100 antipodes and four axis points, followed by actual nearest-center distances. This is numerical whole-disk evaluation, not an exact covering-radius certificate. The first6 were independently prepared and optimized by the transport sibling; this geometry writer ran no whole100 covering or proxy optimizer. The other13 were handed off separately, with the unexpected3 first. They remain proposals in this report; later optimization outcomes belong to the sibling's independent execution ledger.

The fifteen novelty probes (the thirteen accepted followups plus two exact-graph rejects) were checked against44995 retained coordinate fields/34868 unique old arrays,174 sealed metric starts/finals, the first6 seeds, and each other. All ordered and lexicographically sorted unordered hashes are distinct. This proves absence of bitwise duplicate input centers in that retained inventory; it does not prove distinct rotational metric basins or global historical absence. The nineteen accepted handoff starts are also mutually distinct by both hashes.

Two affine-only unexpected arrays are preserved as failures in `unexpected_proposal_actual_audit.json`: from10_to14_pose00_affine00 and from10_to15_pose00_affine00. Their numerical full DT contains a disk-exterior triangle with a negative exact empty-circle sign. They are not accepted as exact actual graphs and are excluded from all handoffs. Unselected target mismatches and duplicate trials remain in their original proposal ledgers.

The earlier B32 finite catalog was reconstructed under necessary conditions for a lower radius threshold. Neither these16 graphs nor the historical absence of twelve classes restricts the current best or a global optimum to this catalog. This branch makes no global lower-bound, covering record, or exact radius claim.

The support bundle preserves scripts, matrix results, source frames, every attempted coordinate file, accepted manifests, novelty evidence, and failures. `replay_sealed_handoffs.py` verifies all handoff hashes and numerical fingerprints and replays the thirteen complete integer actual graph audits; give it the support directory and an output JSON path outside that directory. Generation scripts are archived as provenance. Copy them and the documented dependencies to a fresh working tree before rerunning a generation step, so immutable outputs are never overwritten.
