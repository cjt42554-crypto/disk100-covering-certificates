# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Integer-only outer-wall maps and finite-catalog one-site shortest paths.

The supplied outer_alignment.py is executed with runpy, with its main disabled.
This avoids bytecode writes outside this output directory. No metric file,
writerfinal, optimizer, or old seed inventory is read.
"""
from collections import deque
from hashlib import sha256
import itertools
import json
from pathlib import Path
import runpy


HERE=Path(__file__).resolve().parent
ALIGNMENT=HERE.parent/'outer_alignment.py'
CATALOG=HERE.parent.parent/'b32_catalog_search/frozen_handoff/catalog.json'
functions=runpy.run_path(str(ALIGNMENT))
raw=CATALOG.read_bytes()
graphs=json.loads(raw)['graphs']
matrix=json.loads((HERE/'site_exchange_distance_matrix.json').read_text())['site_exchange_distance_matrix']


def write(name,value):
    (HERE/name).write_text(json.dumps(value,indent=2)+'\n')


def full_map_counts(source,target,pose,maps,row):
    common=dict(maps['core_common_target_to_source'])
    full=row['target_to_source_catalog_bijection']
    common.update({k:v for k,v in full.items() if k>=68})
    source_edges=functions['edges'](source['faces'])
    target_edges=functions['edges'](target['faces'])
    assert set(full)==set(full.values())==set(range(100))
    crossed={e for e in target_edges if e[0]<68<=e[1]}
    crossed_shared={tuple(sorted((common[u],common[v]))) for u,v in crossed
                    if u in common and v in common}&source_edges
    outer={e for e in target_edges if e[0]>=68}
    outer_shared={tuple(sorted((full[u],full[v]))) for u,v in outer}&source_edges
    core={e for e in target_edges if e[1]<68}
    core_shared={tuple(sorted((common[u],common[v]))) for u,v in core
                 if u in common and v in common}&source_edges
    q0_source={tuple(sorted(f)) for f in source['ring_q0_triangles']}
    q1_source={tuple(sorted(f)) for f in source['q1_triangles']}
    def shared_faces(faces,source_faces):
        return len({tuple(sorted(common[v] for v in face)) for face in faces
                    if set(face)<=set(common)}&source_faces)
    q0_shared=shared_faces(target['ring_q0_triangles'],q0_source)
    q1_shared=shared_faces(target['q1_triangles'],q1_source)
    assert len(crossed_shared)+len(outer_shared)==row['common_annulus_edges']
    assert q0_shared+q1_shared==row['common_annulus_faces']
    affected=row['affected_source_outer_ids']
    # Wall clusters use adjacency on the source's outer32 label cycle.
    clusters=[]; unvisited=set(affected)
    while unvisited:
        seed=min(unvisited); unvisited.remove(seed); stack=[seed]; cluster=[]
        while stack:
            v=stack.pop(); cluster.append(v)
            neighbors={68+(v-68-1)%32,68+(v-68+1)%32}&unvisited
            unvisited.difference_update(neighbors); stack.extend(neighbors)
        clusters.append(sorted(cluster))
    return {
        'common_core_sites':len(maps['core_common_target_to_source']),
        'common_core_triangles':pose['core_face_overlap'],
        'common_core_edges':len(core_shared),
        'common_cross_edges':len(crossed_shared),
        'common_outer_cycle_edges':len(outer_shared),
        'common_annulus_edges':row['common_annulus_edges'],
        'common_annulus_faces':row['common_annulus_faces'],
        'common_ring_q0_faces':q0_shared,'common_q1_faces':q1_shared,
        'affected_outer_vertex_count':len(affected),
        'affected_outer_cycle_clusters':clusters,
        'full_bijection_size':len(full),'full_bijection_is_permutation_of_0_to_99':True,
    }


def record(i,j,pose,maps,row):
    counts=full_map_counts(graphs[i],graphs[j],pose,maps,row)
    result={'source_index':i,'target_index':j,
            'integer_lattice_pose':pose,
            'removed_source_core_ids':maps['removed_source_core_ids'],
            'added_target_core_ids':maps['added_target_core_ids'],
            'removed_source_axial':maps['removed_source_axial'],
            'added_target_posed_axial':maps['added_target_posed_axial'],
            'common_core_target_to_source_label_map':maps['core_common_target_to_source'],
            'outer_alignment':row,'counts':counts}
    result['target_to_source_catalog_bijection_array']=[row['target_to_source_catalog_bijection'][k] for k in range(100)]
    inverse={v:k for k,v in row['target_to_source_catalog_bijection'].items()}
    result['source_to_target_catalog_bijection_array']=[inverse[k] for k in range(100)]
    return result


identity_rows=[]
for i,j in ((14,15),(15,14)):
    maximum,poses=functions['optimal_poses'](graphs[i],graphs[j])
    identity=next(p for p in poses if not p['reflection'] and p['rotation']==0 and p['translation']==[0,0])
    maps=functions['outer_maps'](graphs[i],graphs[j],identity)
    assert maximum==67 and len(maps['best_outer_alignments'])==1
    row=maps['best_outer_alignments'][0]
    result=record(i,j,identity,maps,row)
    result['all_outer32_alignment_scores']=maps['all_outer_alignment_scores']
    identity_rows.append(result)
forward,backward=identity_rows
assert forward['target_to_source_catalog_bijection_array']==backward['source_to_target_catalog_bijection_array']
assert forward['source_to_target_catalog_bijection_array']==backward['target_to_source_catalog_bijection_array']
metadata={
    'status':'PASS','catalog_sha256':sha256(raw).hexdigest(),
    'outer_alignment_script_sha256':sha256(ALIGNMENT.read_bytes()).hexdigest(),
    'scope':'Exact integer core-coordinate and annulus-incidence comparisons within the supplied 16 graphs only.',
    'index_convention':'Zero-based independent_index; labels 0..67 core, 68..99 source or target outer cycle.',
    'map_convention':'Stored 100-bijection maps target catalog labels into source catalog labels. Common core labels match exact axial coordinates. Each added target site occupies the one removed source label.',
    'outer_convention':'target outer label 68+j maps to source label 68+(shift+direction*j) modulo 32.',
    'affected_outer_definition':'Source outer labels occurring in the symmetric difference of annulus face sets after the full replacement-site bijection.',
    'integer_only':True,
}
write('pair14_to15_outer_wall_maps.json',dict(metadata,
    forward_and_inverse_identity_pose_maps=identity_rows,
    exact_forward_inverse_full_bijection_check=True))


def candidate_rank(z):
    p,m,r=z
    return (-p['core_face_overlap'],-r['common_annulus_edges'],-r['common_annulus_faces'],
            r['annulus_face_symmetric_difference'],r['cross_edge_symmetric_difference'],
            len(r['affected_source_outer_ids']),p['reflection'],p['rotation'],
            tuple(p['translation']),r['outer_shift'])


def edge_record(i,j):
    _,poses=functions['optimal_poses'](graphs[i],graphs[j])
    candidates=[]
    for p in poses:
        maps=functions['outer_maps'](graphs[i],graphs[j],p)
        for row in maps['best_outer_alignments']:
            candidates.append((p,maps,row))
    candidates.sort(key=candidate_rank)
    best=candidates[0]
    result=record(i,j,*best)
    result['six_distinct_outer_vertices_affected']=len(best[2]['affected_source_outer_ids'])==6
    result['maximum_site_pose_outer_optimum_summaries']=[{
        'reflection':p['reflection'],'rotation':p['rotation'],'translation':p['translation'],
        'core_face_overlap':p['core_face_overlap'],'outer_shift':r['outer_shift'],
        'outer_direction':r['outer_direction'],'common_annulus_edges':r['common_annulus_edges'],
        'common_annulus_faces':r['common_annulus_faces'],
        'annulus_face_symmetric_difference':r['annulus_face_symmetric_difference'],
        'cross_edge_symmetric_difference':r['cross_edge_symmetric_difference'],
        'affected_source_outer_ids':r['affected_source_outer_ids']}
        for p,m,r in candidates]
    return result


vertices=list(range(7,16))
directed_edges={}
for i in vertices:
    for j in vertices:
        if matrix[i][j]==1:
            directed_edges[i,j]=edge_record(i,j)
assert len(directed_edges)==52
write('nine_graph_one_site_outer_edge_maps.json',dict(metadata,
    directed_edge_count=len(directed_edges),
    directed_edges=list(directed_edges.values())))
one_adjacency={i:[j for j in vertices if matrix[i][j]==1] for i in vertices}
six_adjacency={i:[j for j in vertices if (i,j) in directed_edges
                 and directed_edges[i,j]['six_distinct_outer_vertices_affected']] for i in vertices}
assert all((j in six_adjacency[i])==(i in six_adjacency[j]) for i in vertices for j in vertices)


def shortest(adjacency,start,end):
    queue=deque([[start]])
    seen={start}
    while queue:
        route=queue.popleft();here=route[-1]
        if here==end:return route
        for next_vertex in adjacency[here]:
            if next_vertex not in seen:
                seen.add(next_vertex);queue.append(route+[next_vertex])
    return None


def path_record(route):
    steps=[]
    for source,target in zip(route,route[1:]):
        e=directed_edges[source,target]
        steps.append({'source_index':source,'target_index':target,
            'site_exchange_distance':1,
            'affected_source_outer_ids':e['outer_alignment']['affected_source_outer_ids'],
            'affected_outer_vertex_count':e['counts']['affected_outer_vertex_count'],
            'outer_shift':e['outer_alignment']['outer_shift'],
            'outer_direction':e['outer_alignment']['outer_direction'],
            'integer_lattice_pose':{k:v for k,v in e['integer_lattice_pose'].items()
                                   if k!='target_posed_axial'},
            'full_map_reference':f'nine_graph_one_site_outer_edge_maps.json directed edge {source}->{target}'})
    return {'route':route,'hop_count':len(route)-1,
            'every_step_one_site_and_six_outer_vertices':all(s['affected_outer_vertex_count']==6 for s in steps),
            'steps':steps}


all_paths=[]
for i in vertices:
    for j in vertices:
        p=shortest(one_adjacency,i,j);six=shortest(six_adjacency,i,j)
        assert p is not None and six is not None
        all_paths.append({'source':i,'target':j,
                          'minimum_one_site_hops':len(p)-1,'one_site_shortest_route':p,
                          'minimum_hops_requiring_six_distinct_outer_vertices_each_step':len(six)-1,
                          'six_outer_vertices_shortest_path':path_record(six)})
missing=[(i,j) for i,j in itertools.combinations(vertices,2) if matrix[i][j]!=1]
assert missing==[(7,12),(7,13),(8,12),(8,13),(9,12),(9,13),(12,14),(12,15),(13,14),(13,15)]
assert len(set(six_adjacency[7])^set(one_adjacency[7]))==1
assert [e for e in directed_edges if not directed_edges[e]['six_distinct_outer_vertices_affected']]==[(7,9),(9,7)]
all_adjacency={i:[j for j in range(16) if matrix[i][j]==1] for i in range(16)}
unvisited=set(range(16));components=[]
while unvisited:
    start=min(unvisited);component=[j for j in range(16) if shortest(all_adjacency,start,j) is not None]
    components.append(component);unvisited.difference_update(component)
assert components==[[0],[1],[2,3,4],[5,6],vertices]
cover=[14,15,7,8,9,10,11,12,13]
assert set(cover)==set(vertices) and len(cover)==9
assert all(b in six_adjacency[a] for a,b in zip(cover,cover[1:]))
recommendations={
    'baseline10_to_any_other_graph_in_7_to_15':[
        path_record([10,j]) for j in vertices if j!=10],
    'stable14_to15':path_record([14,15]),
    'stable14_or15_to12_or13':[
        path_record([i,10,j]) for i in (14,15) for j in (12,13)],
    'equal_minimum_hop_alternatives_via11':[
        path_record([i,11,j]) for i in (14,15) for j in (12,13)],
    'cover_all_nine_once_starting14':dict(path_record(cover),
         minimum_hop_reason='Starting at one of 9 distinct graphs requires at least 8 transitions to visit the other 8. This 8-hop route achieves that bound.'),
    'strict_six_outer_vertices_path7_to9':path_record([7,10,9]),
}
write('one_site_shortpath_recommendations.json',dict(metadata,
    induced_nine_graph={
        'vertices':vertices,'connected':True,'clique':False,
        'undirected_one_site_edge_count':26,'possible_pair_count':36,
        'missing_one_site_edges':missing,'one_site_graph_diameter':2,
        'structure':'Union of the clique {7,8,9,10,11,14,15} and clique {10,11,12,13}, sharing edge {10,11}.',
        'one_site_adjacency':one_adjacency,
        'strict_six_affected_outer_vertex_adjacency':six_adjacency,
        'strict_six_affected_outer_vertex_edge_count':25,
        'strict_six_affected_outer_vertex_graph_diameter':2,
        'one_site_edge_exception_to_six_outer_vertices':{
            'edge':[7,9],'affected_source_outer_ids_in_7_frame':[83,84,85,86,87],
            'affected_vertex_count':5,'annulus_face_symmetric_difference':14,
            'cross_edge_symmetric_difference':10,
            'explanation':'The two local wall neighborhoods overlap; a one-site exchange need not affect six distinct outer vertices.'}},
    all_sixteen_one_site_connected_components=components,
    all_nine_ordered_pair_shortest_paths=all_paths,
    recommendations=recommendations,
    scope_limit='No one-site path exists between the nine-graph component and indices 0..6 within this supplied catalog. This says nothing about cores absent from the catalog.'))

report=['# Integer outer-wall and one-site path checks','',
        'Status: PASS. All labels are source/target catalog labels; all calculations use exact axial coordinates and graph incidences.',
        '', '14→15 identity and inverse 15→14 identity each have a unique best outer32 map: shift 0, direction +1. Each retains 67 core sites, 106 core triangles, 172 core edges, 52 cross edges, 32 outer cycle edges, and 50 annulus faces. The six affected source outer labels are 78,79,80,84,85,86, forming two three-vertex cycle clusters. The two full 100-bijections are exact inverses.',
        '', 'The nine graphs 7..15 are connected with diameter 2, with 26 of 36 possible one-site edges. They are the union of K7 on {7,8,9,10,11,14,15} and K4 on {10,11,12,13}, sharing edge {10,11}. The ten missing edges join {7,8,9,14,15} to {12,13}.',
        '', '25 of the 26 one-site edges have six affected outer vertices under the best incidence alignment. Edge 7↔9 is the exception: five distinct affected labels because the wall neighborhoods overlap. Its annulus-face symmetric difference is 14 and cross-edge symmetric difference is 10. Routes requiring exactly six distinct affected outer vertices avoid this edge.',
        '', 'Minimum-hop recommendations, with one site and six outer vertices per step:',
        '', '- From 10 to any of the other eight graphs: one direct step.',
        '- 14→15: one direct step.',
        '- From 14 or 15 to 12 or 13: two steps via 10 or 11.',
        '- For 7→9 under the six-vertex condition: 7→10→9, two steps.',
        '- To visit all nine once starting at 14: 14→15→7→8→9→10→11→12→13, eight steps. Eight is the minimum number required to reach eight other distinct graphs.',
        '', 'The full 16-graph one-site components are {0}, {1}, {2,3,4}, {5,6}, and {7,8,9,10,11,12,13,14,15}. No one-site path to 0..6 exists inside this finite catalog.',
        '', 'Full maps: pair14_to15_outer_wall_maps.json and nine_graph_one_site_outer_edge_maps.json. All shortest paths and per-step wall labels: one_site_shortpath_recommendations.json.',
        '', 'Reproduce:', '', '```sh',
        'python3 b32_core_exchange_geometry/overlap_analysis/analyze_outer_wall_shortpaths.py',
        '```','',
        'No metric coordinate matching, optimizer execution, or global catalog completeness inference is used.','']
(HERE/'OUTER_WALL_SHORTPATHS.md').write_text('\n'.join(report))
print(json.dumps({'status':'PASS','identity_maps':[{'source':r['source_index'],'target':r['target_index'],
    'counts':r['counts'],'outer_shift':r['outer_alignment']['outer_shift'],
    'outer_direction':r['outer_alignment']['outer_direction']} for r in identity_rows],
    'one_site_edge_count':26,'six_outer_vertex_edge_count':25,
    'cover_route':cover,'connected_components':components}))
