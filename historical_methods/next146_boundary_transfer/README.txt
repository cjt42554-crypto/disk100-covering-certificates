PUBLICATION SCOPE: CURATED HISTORICAL RESEARCH RECORD
This is a lean, privacy-edited copy. Original replay/rebuild instructions below describe the complete historical archive, not a ready-to-run public reproduction. Large run snapshots, inventories, caches and some dependency trees are intentionally omitted. No complete replay of these finite runs is claimed here; no program was executed for publication. Read ../README.md and this route's PUBLICATION_SCOPE.json. Original files were retained unchanged.

New grid25/de82 146-active donor transfer experiment

The prior429 and1467 file trees remain unchanged. Baseline is grid25 SHAde82a12038c6f9a705b7a985c7e25c53c8e94f5c506946aba9a6c6f03a02852d, radius0.11544695488558258. All preparation NNLS/deletion/response statistics belong to this baseline only.

The selected12 donors comprise5 low-debt,4 low-incident-stress and3 mixed-rank centers. All100 donors were deleted and whole-disk scored independently before selection; no nonbinding disk is assumed removable. The core grid covers all32 actual active boundary maxholes, plus3 nearest nonincident active internal maxholes for each of the5 lowest-debt donors. There are399 total planned entries; skips and actual count are explicit in grid.json/RESULT.json. Every start preserves the other99 centers bitwise and checks its hash against all prior465 starts and this grid.

All100 centers are zero-offset epigraph variables. Four BLAS1 workers use separate output files. Acceptance uses the complete numerical whole-disk critical-point oracle, never sampled coverage. Stage inventories, newly exposed binding holes and frozen-objective versus actual-radius gaps are retained in topology_obstacles.json and each stage record. Limits of one or two optimizer rounds are stated and status8 cases remain explicit.

Independent audit covers all100 deletion configurations and every final by every triple circumcenter, every pair-circle crossing, antipodes and axes. This floating replay is independent of Delaunay and is not an exact covering or optimality proof. Exact certification is performed by the parent.

Reproduction: python3 -B source/grid_search.py --prepare-only; python3 -B source/grid_search.py; python3 -B source/audit_finalize.py; python3 -B source/freeze_result.py.
