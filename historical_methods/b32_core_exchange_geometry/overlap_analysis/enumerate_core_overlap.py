# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Exact finite-catalog axial-lattice overlap, with no numerical matching.

Only writes beneath --output. Reads core.axial_vertices / triangles, together
with the core boundary cycle for integer consistency checks. Python stdlib only.
"""
import argparse
from collections import Counter
from hashlib import sha256
import itertools
import json
from pathlib import Path


def det(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def point(p, rotation, reflected):
    q, r = p
    if reflected:
        q, r = r, q
    for _ in range(rotation):
        q, r = -r, q+r
    return q, r


def positive_face(face):
    a, b, c = face
    if det(a,b,c) < 0:
        b,c = c,b
    cycle = (a,b,c)
    return min(cycle[i:]+cycle[:i] for i in range(3))


def edge(a,b):
    return tuple(sorted((a,b)))


def geometry(vertices, triangles):
    faces = {}
    positive = {}
    edge_faces = Counter()
    incidence = Counter()
    directed = set()
    for index, labels in enumerate(triangles):
        vv = tuple(vertices[k] for k in labels)
        key = tuple(sorted(vv))
        assert key not in faces
        assert abs(det(*vv)) == 1
        for a,b in itertools.combinations(vv,2):
            dq,dr = a[0]-b[0],a[1]-b[1]
            assert dq*dq+dq*dr+dr*dr == 1
        faces[key] = index
        cycle = positive_face(vv)
        positive[key] = cycle
        for a,b in zip(cycle,cycle[1:]+cycle[:1]):
            edge_faces[edge(a,b)] += 1
            directed.add((a,b))
        incidence.update(vv)
    edges = set(edge_faces)
    boundary = {e for e,n in edge_faces.items() if n == 1}
    assert all(n in (1,2) for n in edge_faces.values())
    boundary_vertices = set(itertools.chain.from_iterable(boundary))
    boundary_directed = {d for d in directed if edge(*d) in boundary}
    profile = {(p,incidence[p]) for p in boundary_vertices}
    return {
        'vertices':set(vertices), 'faces':set(faces), 'face_indices':faces,
        'positive_faces':positive, 'edges':edges, 'directed_edges':directed,
        'boundary_edges':boundary, 'boundary_directed_edges':boundary_directed,
        'boundary_vertices':boundary_vertices, 'boundary_profile':profile,
        'incidence':incidence,
    }


def translated_geometry(graph, rotation, reflected, translation):
    dx,dy = translation
    vertices = [point(p,rotation,reflected) for p in graph['vertices']]
    vertices = [(q+dx,r+dy) for q,r in vertices]
    return vertices,geometry(vertices,graph['triangles'])


def transform_data(rotation,reflected,translation):
    e1,e2 = point((1,0),rotation,reflected),point((0,1),rotation,reflected)
    return {
        'rotation_k':rotation, 'rotation_degrees_ccw':rotation*60,
        'reflection_swap_q_r_first':reflected,
        'linear_matrix':[[e1[0],e2[0]],[e1[1],e2[1]]],
        'translation':list(translation),
        'formula':'target_in_source_frame = R^k S^s target + translation',
    }


def metrics(a,b):
    common = {key:len(a[key]&b[key]) for key in (
        'vertices','faces','edges','directed_edges','boundary_edges',
        'boundary_directed_edges','boundary_vertices','boundary_profile')}
    changes = {key:len(a[key]^b[key]) for key in (
        'boundary_edges','boundary_directed_edges','boundary_vertices',
        'boundary_profile')}
    return {
        'common_sites':common['vertices'],
        'site_exchange_distance':68-common['vertices'],
        'common_core_triangles':common['faces'],
        'common_core_edges':common['edges'],
        'common_positive_directed_core_edges':common['directed_edges'],
        'common_boundary_edges':common['boundary_edges'],
        'common_positive_boundary_edges':common['boundary_directed_edges'],
        'common_boundary_sites':common['boundary_vertices'],
        'common_boundary_incidence_profiles':common['boundary_profile'],
        'boundary_ring_edge_changes':changes['boundary_edges'],
        'positive_boundary_edge_changes':changes['boundary_directed_edges'],
        'boundary_site_role_changes':changes['boundary_vertices'],
        'boundary_incidence_profile_changes':changes['boundary_profile'],
    }


def rank(pose):
    m = pose['metrics']; t = pose['transform']
    return (
        -m['common_core_triangles'],m['boundary_ring_edge_changes'],
        m['boundary_incidence_profile_changes'],-m['common_core_edges'],
        -m['common_positive_directed_core_edges'],
        m['positive_boundary_edge_changes'],m['boundary_site_role_changes'],
        t['reflection_swap_q_r_first'],t['rotation_k'],tuple(t['translation']))


def enumerate_pair(a,b):
    best = -1
    candidates = []
    for reflected in (False,True):
        for rotation in range(6):
            vv = [point(p,rotation,reflected) for p in b['vertices']]
            translations = Counter((p[0]-q[0],p[1]-q[1])
                                   for p in a['vertices'] for q in vv)
            rotation_best = max(translations.values())
            if rotation_best < best:
                continue
            if rotation_best > best:
                best = rotation_best
                candidates = []
            candidates.extend((rotation,reflected,tr)
                              for tr,n in translations.items() if n == best)
    poses = []
    for rotation,reflected,tr in candidates:
        _,bgeo = translated_geometry(b,rotation,reflected,tr)
        m = metrics(a['geometry'],bgeo)
        assert m['common_sites'] == best
        poses.append({'transform':transform_data(rotation,reflected,tr),'metrics':m})
    poses.sort(key=rank)
    for index,p in enumerate(poses):
        p['tie_break_rank'] = index
    return {
        'source_index':a['index'], 'target_index':b['index'],
        'maximum_site_overlap':best, 'site_exchange_distance':68-best,
        'max_site_optimal_pose_count':len(poses),
        'maximum_common_triangles_at_maximum_site_overlap':
            max(p['metrics']['common_core_triangles'] for p in poses),
        'preferred_pose':poses[0], 'all_max_site_optimal_poses':poses,
    }


def pose_details(a,b,pose):
    t = pose['transform']
    bv,bgeo = translated_geometry(b,t['rotation_k'],
        t['reflection_swap_q_r_first'],tuple(t['translation']))
    av = a['vertices']; ageo = a['geometry']
    ai = {p:i for i,p in enumerate(av)}; bi = {p:i for i,p in enumerate(bv)}
    common = ageo['vertices']&bgeo['vertices']
    removed = ageo['vertices']-bgeo['vertices']
    added = bgeo['vertices']-ageo['vertices']
    common_map = [{'baseline_label':ai[p],'target_label':bi[p],
                   'coordinate':list(p)} for p in sorted(common)]

    def face_records(keys,geo,labels,side):
        return [{side+'_face_index':geo['face_indices'][f],
                 side+'_labels':[labels[p] for p in f],
                 'coordinates':[list(p) for p in f],
                 'positive_cycle_coordinates':[list(p) for p in geo['positive_faces'][f]]}
                for f in sorted(keys)]

    def edge_records(keys,labels,side):
        return [{side+'_labels':[labels[p] for p in e],
                 'coordinates':[list(p) for p in e]} for e in sorted(keys)]

    removed_faces = ageo['faces']-bgeo['faces']; added_faces = bgeo['faces']-ageo['faces']
    removed_edges = ageo['edges']-bgeo['edges']; added_edges = bgeo['edges']-ageo['edges']
    removed_ring = ageo['boundary_edges']-bgeo['boundary_edges']
    added_ring = bgeo['boundary_edges']-ageo['boundary_edges']
    role_profiles = []
    for p in sorted((ageo['boundary_vertices']|bgeo['boundary_vertices'])&common):
        ab = p in ageo['boundary_vertices']; bb = p in bgeo['boundary_vertices']
        ac,bc = ageo['incidence'][p],bgeo['incidence'][p]
        if ab != bb or ac != bc:
            role_profiles.append({
                'coordinate':list(p), 'baseline_label':ai[p], 'target_label':bi[p],
                'baseline_on_boundary':ab, 'target_on_boundary':bb,
                'baseline_core_triangle_incidence':ac,
                'target_core_triangle_incidence':bc,
                'baseline_outer_q1_fan_size_if_boundary':4-ac if ab else None,
                'target_outer_q1_fan_size_if_boundary':4-bc if bb else None,
            })
    # Components of changed triangles sharing an edge. Isolated vertex contacts
    # are kept separate, so this is a concrete local face-change motif.
    changed_faces = sorted(removed_faces|added_faces)
    components = []
    unvisited = set(range(len(changed_faces)))
    while unvisited:
        seed=min(unvisited); unvisited.remove(seed); queue=[seed]; component=[]
        while queue:
            current=queue.pop(); component.append(current)
            neighbors=[j for j in sorted(unvisited)
                       if len(set(changed_faces[current])&set(changed_faces[j])) >= 2]
            unvisited.difference_update(neighbors); queue.extend(neighbors)
        cc=[changed_faces[i] for i in component]
        coords=set(itertools.chain.from_iterable(cc))
        components.append({
            'removed_triangle_coordinates':[[list(p) for p in f] for f in cc if f in removed_faces],
            'added_triangle_coordinates':[[list(p) for p in f] for f in cc if f in added_faces],
            'support_coordinates':[list(p) for p in sorted(coords)],
            'removed_sites':[list(p) for p in sorted(coords&removed)],
            'added_sites':[list(p) for p in sorted(coords&added)],
        })
    return {
        'transform':t, 'metrics':pose['metrics'],
        'common_label_map':common_map,
        'target_to_baseline_common_label_map':{str(row['target_label']):row['baseline_label'] for row in common_map},
        'baseline_to_target_common_label_map':{str(row['baseline_label']):row['target_label'] for row in common_map},
        'removed_baseline_sites':[{'baseline_label':ai[p],'coordinate':list(p)} for p in sorted(removed)],
        'added_transformed_target_sites':[{'target_label':bi[p],
             'original_coordinate':list(b['vertices'][bi[p]]),
             'transformed_coordinate':list(p)} for p in sorted(added)],
        'all_transformed_target_coordinates':[list(p) for p in bv],
        'removed_core_faces':face_records(removed_faces,ageo,ai,'baseline'),
        'added_core_faces':face_records(added_faces,bgeo,bi,'target'),
        'removed_core_edges':edge_records(removed_edges,ai,'baseline'),
        'added_core_edges':edge_records(added_edges,bi,'target'),
        'removed_boundary_ring_edges':edge_records(removed_ring,ai,'baseline'),
        'added_boundary_ring_edges':edge_records(added_ring,bi,'target'),
        'changed_retained_site_boundary_profiles':role_profiles,
        'changed_core_face_edge_motif':{
            'removed_face_count':len(removed_faces),'added_face_count':len(added_faces),
            'removed_edge_count':len(removed_edges),'added_edge_count':len(added_edges),
            'removed_faces_using_only_common_sites':sum(set(f)<=common for f in removed_faces),
            'added_faces_using_only_common_sites':sum(set(f)<=common for f in added_faces),
            'removed_edges_using_only_common_sites':sum(set(e)<=common for e in removed_edges),
            'added_edges_using_only_common_sites':sum(set(e)<=common for e in added_edges),
            'edge_connected_changed_triangle_components':components,
        },
    }


def dump(path,value):
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalog',type=Path,default=Path('b32_catalog_search/frozen_handoff/catalog.json'))
    parser.add_argument('--output',type=Path,default=Path('b32_core_exchange_geometry/overlap_analysis'))
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    raw=args.catalog.read_bytes(); data=json.loads(raw)
    graphs=[]; checks=[]
    neighbor_directions=((1,0),(0,1),(1,-1),(-1,0),(0,-1),(-1,1))
    for index,g in enumerate(data['graphs']):
        assert g['independent_index'] == index
        c=g['core']; vv=[tuple(p) for p in c['axial_vertices']]
        geo=geometry(vv,c['triangles'])
        assert len(vv)==len(set(vv))==68
        assert (len(geo['edges']),len(geo['faces']),len(geo['boundary_edges'])) == (175,108,26)
        ring=[vv[k] for k in c['boundary_ids']]
        assert {edge(a,b) for a,b in zip(ring,ring[1:]+ring[:1])} == geo['boundary_edges']
        induced_edges={edge(p,(p[0]+dq,p[1]+dr)) for p in vv
                       for dq,dr in neighbor_directions if (p[0]+dq,p[1]+dr) in geo['vertices']}
        induced_faces={tuple(sorted((p,q,r))) for p in vv
                       for q,r in itertools.combinations(
                           [(p[0]+dq,p[1]+dr) for dq,dr in neighbor_directions
                            if (p[0]+dq,p[1]+dr) in geo['vertices']],2)
                       if (q[0]-r[0])**2+(q[0]-r[0])*(q[1]-r[1])+(q[1]-r[1])**2 == 1}
        checks.append({'index':index,'vertex_count':len(vv),'edge_count':len(geo['edges']),
                       'face_count':len(geo['faces']),'boundary_count':len(geo['boundary_edges']),
                       'all_unit_edges_induced_by_sites_present':geo['edges']==induced_edges,
                       'all_elementary_triangles_induced_by_sites_present':geo['faces']==induced_faces})
        graphs.append({'index':index,'vertices':vv,'triangles':c['triangles'],'geometry':geo})
    assert len(graphs)==16
    pairs=[]; matrix=[[None]*16 for _ in range(16)]
    face_matrix=[[None]*16 for _ in range(16)]
    for i,a in enumerate(graphs):
        for j,b in enumerate(graphs):
            pair=enumerate_pair(a,b); pairs.append(pair)
            matrix[i][j]=pair['site_exchange_distance']
            face_matrix[i][j]=pair['maximum_common_triangles_at_maximum_site_overlap']
        print('completed source',i,flush=True)
    assert all(matrix[i][i]==0 for i in range(16))
    assert all(matrix[i][j]==matrix[j][i] for i in range(16) for j in range(16))
    assert all(face_matrix[i][j]==face_matrix[j][i] for i in range(16) for j in range(16))
    assert all(matrix[i][k]<=matrix[i][j]+matrix[j][k]
               for i in range(16) for j in range(16) for k in range(16))
    metadata={
        'status':'PASS', 'scope':'Exact overlap over the supplied 16 cores only; no global catalog completeness or optimality claim.',
        'catalog_path':str(args.catalog),'catalog_sha256':sha256(raw).hexdigest(),
        'catalog_fingerprints':[g['fingerprint'] for g in data['graphs']],
        'indices':'zero-based independent_index in frozen_handoff/catalog.json',
        'source_target_convention':'Apply each stored pose to target coordinates to express them in the unchanged source frame.',
        'axial_embedding':'(q,r) -> (q+r/2, sqrt(3)*r/2); all computations use integer axial coordinates.',
        'rotation':'R(q,r)=(-r,q+r), 60 degrees counterclockwise',
        'reflection':'S(q,r)=(r,q), applied before R^k',
        'd6_elements':[transform_data(k,s,(0,0)) for s in (False,True) for k in range(6)],
        'translation_enumeration':'For each D6 element, enumerate every source_point - transformed_target_point. Its Counter multiplicity is the exact site intersection count.',
        'completeness_of_translation_search':'Any pose with nonzero site overlap matches a source and transformed target point, so its translation appears among point differences. The maximum is positive.',
        'triangle_matching':'Unordered triples of exact integer axial coordinates. Positive directed core edges are the union of directed edges of counterclockwise normalized core faces.',
        'boundary_ring_change_definition':'Cardinality of symmetric difference of undirected core boundary edge sets, each boundary edge having one core face. Removed and added edges both count.',
        'boundary_profile_definition':'Set of (boundary coordinate, core triangle incidence); changed incidences count as removal plus addition. q1 fan size at a boundary site is 4 minus core triangle incidence.',
        'tie_break_order':[
            'maximum common core triangles','minimum boundary ring edge changes',
            'minimum boundary incidence profile changes','maximum common core edges',
            'maximum common positive directed core edges','minimum positive boundary edge changes',
            'minimum boundary site role changes','prefer no reflection, then lower rotation k, then lexicographic translation'],
        'consistency_checks':checks,
        'matrix_checks':{'diagonal_zero':True,'symmetric':True,'triangle_inequality':True},
    }
    dump(args.output/'site_exchange_distance_matrix.json',dict(metadata,site_exchange_distance_matrix=matrix,
         max_common_triangle_matrix_at_max_site_overlap=face_matrix))
    dump(args.output/'all_pair_max_site_optimal_poses.json',dict(metadata,pairs=pairs))
    baseline=[]
    for j,b in enumerate(graphs):
        if j==10:continue
        pair=pairs[10*16+j]
        details=[pose_details(graphs[10],b,pose) for pose in pair['all_max_site_optimal_poses']]
        baseline.append({'baseline_index':10,'target_index':j,
                         'site_exchange_distance':pair['site_exchange_distance'],
                         'maximum_site_overlap':pair['maximum_site_overlap'],
                         'max_site_optimal_pose_count':len(details),
                         'preferred_pose':details[0],
                         'all_max_site_optimal_poses_with_full_exchange_details':details})
    dump(args.output/'baseline10_exact_exchange_details.json',dict(metadata,baseline_exchanges=baseline))
    special=pairs[14*16+15]
    special_details=[pose_details(graphs[14],graphs[15],pose)
                     for pose in special['all_max_site_optimal_poses']]
    dump(args.output/'pair14_to15_exact_exchange_details.json',dict(metadata,
         source_index=14,target_index=15,
         label_note='Fields named baseline_label refer to source index 14 in this file.',
         site_exchange_distance=special['site_exchange_distance'],
         maximum_site_overlap=special['maximum_site_overlap'],
         maximum_common_triangles_at_maximum_site_overlap=
             special['maximum_common_triangles_at_maximum_site_overlap'],
         preferred_pose=special_details[0],
         all_max_site_optimal_poses_with_full_exchange_details=special_details))
    csv='index,'+','.join(map(str,range(16)))+'\n'
    csv+=''.join(str(i)+','+','.join(map(str,row))+'\n' for i,row in enumerate(matrix))
    (args.output/'site_exchange_distance_matrix.csv').write_text(csv)
    report=['# Exact B32 core overlap within the supplied 16-graph catalog','',
            'Status: PASS. Index convention: zero-based `independent_index`.',
            '',f'Catalog SHA256: `{metadata["catalog_sha256"]}`','',
            'All poses map target coordinates into the unchanged source frame. Axial basis has a 60-degree angle; R(q,r)=(-r,q+r), S(q,r)=(r,q), pose R^k S^s p+t.',
            '', 'For every ordered pair all 12 D6 elements and every point-difference translation were enumerated. Coordinate, face, edge, and boundary comparisons are exact integers.',
            '', 'The site exchange distance is 68 minus maximum coordinate intersection. Common triangles break maximum-site ties first; the next criterion minimizes the symmetric difference of core boundary edge sets. All maximum-site poses are retained.',
            '', '## Distance matrix','', '| index | '+' | '.join(map(str,range(16)))+' |',
            '|---|'+'---|'*16]
    report+=['| '+str(i)+' | '+' | '.join(map(str,row))+' |' for i,row in enumerate(matrix)]
    report+=['','## Preferred baseline-10 poses','',
             '| target | distance | sites | common triangles | common edges | boundary edge changes | pose (k,s,t) | all max-site poses |',
             '|---|---|---|---|---|---|---|---|']
    for e in baseline:
        p=e['preferred_pose']; m=p['metrics']; t=p['transform']
        report.append(f'| {e["target_index"]} | {e["site_exchange_distance"]} | {m["common_sites"]} | {m["common_core_triangles"]} | {m["common_core_edges"]} | {m["boundary_ring_edge_changes"]} | ({t["rotation_k"]},{int(t["reflection_swap_q_r_first"])},{t["translation"]}) | {e["max_site_optimal_pose_count"]} |')
    report+=['','## Exchange interpretation','',
             'The separate pair14_to15_exact_exchange_details.json also retains all four maximum-site poses for 14 to 15. Its preferred identity pose removes (8,0) and adds (8,5), retains 67 sites, 106 triangles and 172 edges, and changes 8 boundary edges in symmetric difference.',
             '',
             'Each baseline pose includes the complete retained-site label map, the removed baseline coordinates, the added transformed target coordinates, all removed/added core triangles and edges, boundary-edge changes, and retained-site incidence/fan changes. Changed triangles are also grouped into components sharing an edge.',
             '', 'Boundary ring changes refer to the 26-edge core boundary. The 32 outer ring vertices do not have canonical integer coordinates in this input, so their metric geometry is not inferred.',
             '', 'Positive directed edges are computed after normalizing every triangle to counterclockwise orientation, including reflected poses.',
             '', 'Validation checks: 68 distinct vertices, 175 unit edges, 108 unit triangles, 26 core boundary edges per graph; source boundary cycles agree with triangle incidence; distance diagonal/symmetry/triangle inequality; reciprocal common-triangle matrix symmetry.',
             '', 'These results describe only the supplied 16 cores. They do not establish a globally complete catalog, a covering-radius optimum, or any approximate geometric correspondence.',
             '', 'Reproduce from task-3 root:', '',
             '```sh','python3 b32_core_exchange_geometry/overlap_analysis/enumerate_core_overlap.py','```','']
    (args.output/'README.md').write_text('\n'.join(report))
    # Check exact stored full maps and loss/gain counts against every baseline pose.
    for e in baseline:
        for d in e['all_max_site_optimal_poses_with_full_exchange_details']:
            assert len(d['common_label_map'])==e['maximum_site_overlap']
            assert len(d['removed_baseline_sites'])==len(d['added_transformed_target_sites'])==e['site_exchange_distance']
            m=d['metrics']
            assert len(d['removed_core_faces'])==len(d['added_core_faces'])==108-m['common_core_triangles']
            assert len(d['removed_core_edges'])==len(d['added_core_edges'])==175-m['common_core_edges']
    summary={'status':'PASS','source_sha256':metadata['catalog_sha256'],
             'ordered_pair_count':len(pairs),
             'all_max_site_pose_count':sum(p['max_site_optimal_pose_count'] for p in pairs),
             'baseline_pose_count':sum(e['max_site_optimal_pose_count'] for e in baseline),
             'distance_matrix':matrix,
             'baseline10_summary':[{'target':e['target_index'],'distance':e['site_exchange_distance'],
                 'pose_count':e['max_site_optimal_pose_count'],
                 'preferred_transform':e['preferred_pose']['transform'],
                 'preferred_metrics':e['preferred_pose']['metrics']} for e in baseline]}
    dump(args.output/'run_summary.json',summary)
    print(json.dumps(summary,ensure_ascii=False),flush=True)


if __name__=='__main__':
    main()
