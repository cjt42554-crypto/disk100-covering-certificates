PUBLICATION SCOPE: CURATED HISTORICAL RESEARCH RECORD
This is a lean, privacy-edited copy. Original replay/rebuild instructions below describe the complete historical archive, not a ready-to-run public reproduction. Large run snapshots, inventories, caches and some dependency trees are intentionally omitted. No complete replay of these finite runs is claimed here; no program was executed for publication. Read ../README.md and this route's PUBLICATION_SCOPE.json. Original files were retained unchanged.

# Exact B32 core overlap within the supplied 16-graph catalog

Status: PASS. Index convention: zero-based `independent_index`.

Catalog SHA256: `f3234fd61f521362c69e3442940edc1d381fb5c2e06b1cd12a79e93ab2e57bb4`

All poses map target coordinates into the unchanged source frame. Axial basis has a 60-degree angle; R(q,r)=(-r,q+r), S(q,r)=(r,q), pose R^k S^s p+t.

For every ordered pair all 12 D6 elements and every point-difference translation were enumerated. Coordinate, face, edge, and boundary comparisons are exact integers.

The site exchange distance is 68 minus maximum coordinate intersection. Common triangles break maximum-site ties first; the next criterion minimizes the symmetric difference of core boundary edge sets. All maximum-site poses are retained.

## Distance matrix

| index | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 0 | 5 | 4 | 3 | 3 | 3 | 4 | 4 | 4 | 3 | 3 | 4 | 3 | 4 | 4 | 4 |
| 1 | 5 | 0 | 4 | 5 | 5 | 5 | 5 | 4 | 4 | 5 | 5 | 5 | 6 | 6 | 5 | 5 |
| 2 | 4 | 4 | 0 | 1 | 1 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 5 | 4 |
| 3 | 3 | 5 | 1 | 0 | 1 | 3 | 4 | 4 | 4 | 3 | 3 | 4 | 3 | 4 | 4 | 4 |
| 4 | 3 | 5 | 1 | 1 | 0 | 3 | 4 | 4 | 4 | 3 | 3 | 4 | 3 | 4 | 4 | 4 |
| 5 | 3 | 5 | 4 | 3 | 3 | 0 | 1 | 4 | 4 | 3 | 3 | 4 | 3 | 4 | 4 | 4 |
| 6 | 4 | 5 | 4 | 4 | 4 | 1 | 0 | 4 | 5 | 4 | 4 | 4 | 4 | 5 | 5 | 4 |
| 7 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 0 | 1 | 1 | 1 | 1 | 2 | 2 | 1 | 1 |
| 8 | 4 | 4 | 4 | 4 | 4 | 4 | 5 | 1 | 0 | 1 | 1 | 1 | 2 | 2 | 1 | 1 |
| 9 | 3 | 5 | 4 | 3 | 3 | 3 | 4 | 1 | 1 | 0 | 1 | 1 | 2 | 2 | 1 | 1 |
| 10 | 3 | 5 | 4 | 3 | 3 | 3 | 4 | 1 | 1 | 1 | 0 | 1 | 1 | 1 | 1 | 1 |
| 11 | 4 | 5 | 4 | 4 | 4 | 4 | 4 | 1 | 1 | 1 | 1 | 0 | 1 | 1 | 1 | 1 |
| 12 | 3 | 6 | 4 | 3 | 3 | 3 | 4 | 2 | 2 | 2 | 1 | 1 | 0 | 1 | 2 | 2 |
| 13 | 4 | 6 | 4 | 4 | 4 | 4 | 5 | 2 | 2 | 2 | 1 | 1 | 1 | 0 | 2 | 2 |
| 14 | 4 | 5 | 5 | 4 | 4 | 4 | 5 | 1 | 1 | 1 | 1 | 1 | 2 | 2 | 0 | 1 |
| 15 | 4 | 5 | 4 | 4 | 4 | 4 | 4 | 1 | 1 | 1 | 1 | 1 | 2 | 2 | 1 | 0 |

## Preferred baseline-10 poses

| target | distance | sites | common triangles | common edges | boundary edge changes | pose (k,s,t) | all max-site poses |
|---|---|---|---|---|---|---|---|
| 0 | 3 | 65 | 103 | 167 | 14 | (4,0,[0, 12]) | 2 |
| 1 | 5 | 63 | 99 | 161 | 22 | (1,1,[7, -3]) | 4 |
| 2 | 4 | 64 | 101 | 164 | 20 | (0,0,[0, 0]) | 4 |
| 3 | 3 | 65 | 103 | 167 | 14 | (4,0,[0, 12]) | 2 |
| 4 | 3 | 65 | 103 | 167 | 14 | (4,0,[0, 12]) | 2 |
| 5 | 3 | 65 | 103 | 167 | 14 | (4,0,[0, 12]) | 2 |
| 6 | 4 | 64 | 101 | 164 | 20 | (1,0,[8, -4]) | 4 |
| 7 | 1 | 67 | 106 | 172 | 8 | (0,0,[0, 0]) | 1 |
| 8 | 1 | 67 | 106 | 172 | 8 | (0,0,[0, 0]) | 2 |
| 9 | 1 | 67 | 106 | 172 | 8 | (0,0,[0, 0]) | 1 |
| 11 | 1 | 67 | 106 | 172 | 8 | (0,0,[0, 0]) | 3 |
| 12 | 1 | 67 | 106 | 172 | 8 | (0,0,[0, 0]) | 2 |
| 13 | 1 | 67 | 106 | 172 | 8 | (0,0,[0, 0]) | 4 |
| 14 | 1 | 67 | 106 | 172 | 8 | (0,0,[0, 0]) | 2 |
| 15 | 1 | 67 | 106 | 172 | 8 | (1,1,[8, -4]) | 2 |

## Exchange interpretation

The separate pair14_to15_exact_exchange_details.json also retains all four maximum-site poses for 14 to 15. Its preferred identity pose removes (8,0) and adds (8,5), retains 67 sites, 106 triangles and 172 edges, and changes 8 boundary edges in symmetric difference.

Each baseline pose includes the complete retained-site label map, the removed baseline coordinates, the added transformed target coordinates, all removed/added core triangles and edges, boundary-edge changes, and retained-site incidence/fan changes. Changed triangles are also grouped into components sharing an edge.

Boundary ring changes refer to the 26-edge core boundary. The 32 outer ring vertices do not have canonical integer coordinates in this input, so their metric geometry is not inferred.

Positive directed edges are computed after normalizing every triangle to counterclockwise orientation, including reflected poses.

Validation checks: 68 distinct vertices, 175 unit edges, 108 unit triangles, 26 core boundary edges per graph; source boundary cycles agree with triangle incidence; distance diagonal/symmetry/triangle inequality; reciprocal common-triangle matrix symmetry.

These results describe only the supplied 16 cores. They do not establish a globally complete catalog, a covering-radius optimum, or any approximate geometric correspondence.

Reproduce from task-3 root:

```sh
python3 b32_core_exchange_geometry/overlap_analysis/enumerate_core_overlap.py
```
