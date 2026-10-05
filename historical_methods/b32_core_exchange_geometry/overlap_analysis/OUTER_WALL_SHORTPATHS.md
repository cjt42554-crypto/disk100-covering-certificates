PUBLICATION SCOPE: CURATED HISTORICAL RESEARCH RECORD
This is a lean, privacy-edited copy. Original replay/rebuild instructions below describe the complete historical archive, not a ready-to-run public reproduction. Large run snapshots, inventories, caches and some dependency trees are intentionally omitted. No complete replay of these finite runs is claimed here; no program was executed for publication. Read ../README.md and this route's PUBLICATION_SCOPE.json. Original files were retained unchanged.

# Integer outer-wall and one-site path checks

Status: PASS. All labels are source/target catalog labels; all calculations use exact axial coordinates and graph incidences.

14→15 identity and inverse 15→14 identity each have a unique best outer32 map: shift 0, direction +1. Each retains 67 core sites, 106 core triangles, 172 core edges, 52 cross edges, 32 outer cycle edges, and 50 annulus faces. The six affected source outer labels are 78,79,80,84,85,86, forming two three-vertex cycle clusters. The two full 100-bijections are exact inverses.

The nine graphs 7..15 are connected with diameter 2, with 26 of 36 possible one-site edges. They are the union of K7 on {7,8,9,10,11,14,15} and K4 on {10,11,12,13}, sharing edge {10,11}. The ten missing edges join {7,8,9,14,15} to {12,13}.

25 of the 26 one-site edges have six affected outer vertices under the best incidence alignment. Edge 7↔9 is the exception: five distinct affected labels because the wall neighborhoods overlap. Its annulus-face symmetric difference is 14 and cross-edge symmetric difference is 10. Routes requiring exactly six distinct affected outer vertices avoid this edge.

Minimum-hop recommendations, with one site and six outer vertices per step:

- From 10 to any of the other eight graphs: one direct step.
- 14→15: one direct step.
- From 14 or 15 to 12 or 13: two steps via 10 or 11.
- For 7→9 under the six-vertex condition: 7→10→9, two steps.
- To visit all nine once starting at 14: 14→15→7→8→9→10→11→12→13, eight steps. Eight is the minimum number required to reach eight other distinct graphs.

The full 16-graph one-site components are {0}, {1}, {2,3,4}, {5,6}, and {7,8,9,10,11,12,13,14,15}. No one-site path to 0..6 exists inside this finite catalog.

Full maps: pair14_to15_outer_wall_maps.json and nine_graph_one_site_outer_edge_maps.json. All shortest paths and per-step wall labels: one_site_shortpath_recommendations.json.

Reproduce:

```sh
python3 b32_core_exchange_geometry/overlap_analysis/analyze_outer_wall_shortpaths.py
```

No metric coordinate matching, optimizer execution, or global catalog completeness inference is used.
