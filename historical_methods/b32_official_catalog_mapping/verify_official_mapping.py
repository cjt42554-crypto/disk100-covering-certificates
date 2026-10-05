# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Exact integer graph identity audit; stdlib only, no optimizer or DT extraction.

Run once with --stage-sources in the original task workspace, then replay anywhere
with `python3 verify_official_mapping.py`. All writes stay beside this script.
"""
from collections import Counter, defaultdict, deque
from pathlib import Path
import argparse
import hashlib
import json
import sys
import zipfile

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'inputs'
ZIP_SHA = 'e7220866aab85d8c457dae49b79b4436bd647179d4311877766b345a2cedc510'
LOCAL_SHA = 'f3234fd61f521362c69e3442940edc1d381fb5c2e06b1cd12a79e93ab2e57bb4'
PREFIX = 'disk100_lower_bound_0p11317_20261004/'
MEMBERS = {
    'official_core_catalog.json': PREFIX + 'lower_bounds_extension_20261004/B32_P26_F108_core_catalog.json',
    'official_graph_index.json': PREFIX + 'lower_bounds_extension_20261004/B32_P26_F108_graph_index.json',
}
for i in range(16):
    name = f'B32_P26_F108_{i:03d}.json'
    MEMBERS[name] = PREFIX + 'signed_edge_metric_B32_20261004/graphs/' + name
LOCAL = {
    'local_catalog.json': 'b32_catalog_search/frozen_handoff/catalog.json',
    'current_core_integer_development.json': 'b32_catalog_search/frozen_handoff/current_core_integer_development.json',
    'current_graph_numerical.json': 'b32_catalog_search/frozen_handoff/current_graph_numerical.json',
    'current.npy': 'b32_catalog_search/frozen_support/inputs/current.npy',
}
CATEGORIES = ('core_q0_triangles', 'ring_q0_triangles', 'q1_triangles')


def check(ok, reason):
    if not ok:
        raise ValueError(reason)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def save(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, sort_keys=True) + '\n')


def stage():
    check(not INPUT.exists(), 'inputs_already_exist_do_not_overwrite')
    base = ROOT.parent
    archive = base / 'global_lower_sync_11317/disk100_lower_bound_0p11317_20261004.zip'
    raw = archive.read_bytes()
    check(sha(raw) == ZIP_SHA, 'official_archive_sha_mismatch')
    sources = {}
    blobs = {}
    with zipfile.ZipFile(archive) as z:
        counts = Counter(z.namelist())
        for name, member in MEMBERS.items():
            check(counts[member] == 1, f'nonunique_archive_member:{member}')
            data = z.read(member)
            json.loads(data)
            blobs[name] = data
            sources[name] = dict(origin='official_zip_member', member=member,
                                 sha256=sha(data), bytes=len(data))
    for name, relative in LOCAL.items():
        data = (base / relative).read_bytes()
        if name.endswith('.json'):
            json.loads(data)
        blobs[name] = data
        sources[name] = dict(origin='local_frozen_source', relative_path=relative,
                             sha256=sha(data), bytes=len(data))
    check(sources['local_catalog.json']['sha256'] == LOCAL_SHA, 'local_catalog_sha_mismatch')
    INPUT.mkdir()
    for name, data in blobs.items():
        (INPUT / name).write_bytes(data)
    save('SOURCE_SHA256.json', dict(official_archive=dict(relative_path=str(archive.relative_to(base)),
                                                        sha256=ZIP_SHA, bytes=len(raw)), sources=sources))


def load_sources():
    receipt = json.loads((ROOT / 'SOURCE_SHA256.json').read_text())
    check(receipt['official_archive']['sha256'] == ZIP_SHA, 'archive_provenance_changed')
    check(receipt['sources']['local_catalog.json']['sha256'] == LOCAL_SHA, 'local_provenance_changed')
    check(set(receipt['sources']) == set(MEMBERS) | set(LOCAL), 'source_inventory_changed')
    result = {}
    for name, row in receipt['sources'].items():
        raw = (INPUT / name).read_bytes()
        check(sha(raw) == row['sha256'] and len(raw) == row['bytes'], f'source_hash_mismatch:{name}')
        if name.endswith('.json'):
            result[name] = json.loads(raw)
    return result, receipt


def cycle_norm(cycle):
    t = tuple(cycle)
    return min(t[i:] + t[:i] for i in range(len(t)))


def triangle_norm(f):
    check(len(f) == 3 and len(set(f)) == 3, 'degenerate_face')
    return cycle_norm(f)


def faces_of(graph, category=None):
    if category is None:
        return [f for c in CATEGORIES for f in faces_of(graph, c)]
    return [tuple(f['vertices'] if isinstance(f, dict) else f) for f in graph[category]]


def face_set(faces):
    return {tuple(sorted(f)) for f in faces}


def oriented_set(faces):
    return {triangle_norm(f) for f in faces}


def edges_of(faces):
    return {tuple(sorted((f[j], f[(j+1) % 3]))) for f in faces for j in range(3)}


def mapped(faces, bijection, reflection=False):
    out = [tuple(bijection[v] for v in f) for f in faces]
    return [(f[0], f[2], f[1]) for f in out] if reflection else out


def canonical_disk(faces, boundary):
    """Independent stdlib replay of boundary-rooted oriented disk encoding."""
    directed = {}
    vertices = set()
    for f in faces:
        triangle_norm(f)
        vertices.update(f)
        for j in range(3):
            edge = (f[j], f[(j+1) % 3])
            check(edge not in directed, 'duplicate_oriented_edge')
            directed[edge] = f[(j+2) % 3]
    boundary_edges = [e for e in directed if e[::-1] not in directed]
    following = dict(boundary_edges)
    check(len(following) == len(boundary_edges), 'boundary_branches')
    cycle, v = [], min(following)
    while v not in cycle:
        cycle.append(v)
        check(v in following, 'open_boundary')
        v = following[v]
    check(v == cycle[0] and len(cycle) == len(following), 'multiple_boundary_cycles')
    check(set(cycle) == set(boundary), 'boundary_vertex_set_mismatch')
    check(len(vertices)-len(edges_of(faces))+len(faces) == 1, 'non_disk_euler')
    links = defaultdict(lambda: defaultdict(set))
    for f in faces:
        for j in range(3):
            v, a, b = f[j], f[(j+1) % 3], f[(j+2) % 3]
            links[v][a].add(b)
            links[v][b].add(a)
    for v, link in links.items():
        seen = {min(link)}
        todo = list(seen)
        while todo:
            x = todo.pop()
            for y in link[x] - seen:
                seen.add(y)
                todo.append(y)
        check(seen == set(link), 'disconnected_vertex_link')
        ds = Counter(len(n) for n in link.values())
        check((ds[1] == 2 and not (set(ds)-{1, 2})) if v in cycle else set(ds) == {2},
              'invalid_vertex_link')
    candidates = []
    for reflected in (False, True):
        thirds = {(b, a): c for (a, b), c in directed.items()} if reflected else directed
        seeds = [(b, a) for a, b in boundary_edges] if reflected else boundary_edges
        for a, b in seeds:
            c = thirds[a, b]
            labels = {a: 0, b: 1, c: 2}
            todo = deque([(a, b, c)])
            seen, out = set(), []
            while todo:
                f = todo.popleft()
                key = frozenset(f)
                if key in seen:
                    continue
                seen.add(key)
                out.append(tuple(sorted(labels[x] for x in f)))
                for j in range(3):
                    x, y = f[j], f[(j+1) % 3]
                    if (y, x) in thirds:
                        z = thirds[y, x]
                        if z not in labels:
                            labels[z] = len(labels)
                        todo.append((y, x, z))
            check(len(seen) == len(faces) and set(labels) == vertices, 'disconnected_faces')
            candidates.append(tuple(sorted(out)))
    payload = dict(vertices=len(vertices), boundary=len(cycle), faces=min(candidates))
    return sha(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode())


def validate_graph(g, xy):
    allfaces = faces_of(g)
    check(len(xy) == 68 and len(set(xy)) == 68, 'invalid_core_sites')
    check(all(type(x) is int and type(y) is int for x, y in xy), 'noninteger_axial_site')
    check(set(g['interior_ids']) == set(range(68)), 'interior_label_partition')
    check(set(g['outer_boundary_ids']) == set(range(68, 100)), 'outer_label_partition')
    check(len(allfaces) == 166 and len(face_set(allfaces)) == 166, 'full_face_count')
    check(set(v for f in allfaces for v in f) == set(range(100)), 'full_vertex_count')
    edge_counts = Counter(tuple(sorted((f[j], f[(j+1) % 3]))) for f in allfaces for j in range(3))
    check(len(edge_counts) == 265 and set(edge_counts.values()) == {1, 2}, 'full_edge_count')
    outer = g['outer_boundary_ids']
    boundary_edges = {tuple(sorted((outer[j], outer[(j+1) % 32]))) for j in range(32)}
    check({e for e, n in edge_counts.items() if n == 1} == boundary_edges, 'outer_cycle_edges')
    adj = defaultdict(set)
    for a, b in edge_counts:
        adj[a].add(b)
        adj[b].add(a)
    check(all(len(adj[i]) == 6 for i in range(68)), 'full_interior_degrees')
    for category, count, outer_count in zip(CATEGORIES, (108, 26, 32), (0, 1, 2)):
        fs = faces_of(g, category)
        check(len(fs) == count, f'category_face_count:{category}')
        check(all(sum(v >= 68 for v in f) == outer_count for f in fs), f'category_partition:{category}')
    for f in faces_of(g, 'core_q0_triangles'):
        a, b, c = (xy[v] for v in f)
        det = (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
        check(det == 1, 'core_face_integer_orientation')
        check(all((xy[u][0]-xy[v][0])**2+(xy[u][0]-xy[v][0])*(xy[u][1]-xy[v][1])+
                  (xy[u][1]-xy[v][1])**2 == 1 for u, v in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0]))),
              'core_face_nonunit_edge')
    core_edges = Counter(tuple(sorted((f[j], f[(j+1) % 3]))) for f in faces_of(g, 'core_q0_triangles') for j in range(3))
    ring = g['ring_ids']
    check(len(ring) == 26 and len(set(ring)) == 26, 'ring_count')
    check(len(core_edges) == 175 and {e for e, n in core_edges.items() if n == 1} ==
          {tuple(sorted((ring[j], ring[(j+1) % 26]))) for j in range(26)}, 'core_boundary_cycle')
    apex_counts = Counter(next(v for v in f if v < 68) for f in faces_of(g, 'q1_triangles'))
    check([apex_counts[v] for v in ring] == g['q1_fan_sizes_in_ring_order'], 'q1_fan_metadata')
    for f in g['q1_triangles']:
        if isinstance(f, dict):
            check(f['interior_apex'] == next(v for v in f['vertices'] if v < 68), 'q1_apex_metadata')
    counts = dict(vertices=100, edges=265, faces=166, B=32, I=68, N=26,
                  marked_apex_count_K=len(apex_counts),
                  interior_degree_histogram={str(d): n for d, n in sorted(Counter(len(adj[v]) for v in range(68)).items())},
                  max_marked_fan_size=max(apex_counts.values()),
                  outer_degree_histogram={str(d): n for d, n in sorted(Counter(len(adj[v]) for v in range(68, 100)).items())},
                  oriented_edge_check='PASS', all_vertex_links_check='PASS', connected='PASS')
    fingerprint = canonical_disk(allfaces, outer)
    if 'faces' in g:
        check(oriented_set(g['faces']) == oriented_set(allfaces), 'redundant_full_faces')
    if 'fingerprint' in g:
        check(fingerprint == g['fingerprint'], 'canonical_fingerprint')
    if 'verification' in g:
        check(g['verification'] == counts, 'official_verification_metadata_replay')
    return dict(vertices=100, edges=265, faces=166, core_vertices=68, core_edges=175,
                core_faces=108, core_boundary=26, outer_boundary=32, fingerprint=fingerprint,
                marked_apex_count_K=len(apex_counts),
                ring_fan_histogram={str(p): n for p, n in sorted(Counter(apex_counts[v] for v in ring).items())},
                interior_p_degree_histogram={f'{p},{d}': n for (p, d), n in sorted(Counter((apex_counts[v], len(adj[v])) for v in range(68)).items())})


def transform(p, reflection, rotation):
    x, y = p
    if reflection:
        x, y = y, x
    for _ in range(rotation):
        x, y = -y, x+y
    return x, y


def find_matches(local, official):
    A = [tuple(p) for p in local['core']['axial_vertices']]
    B = [tuple(p) for p in official['core_axial_labels']]
    lookup = {p: i for i, p in enumerate(B)}
    failures = Counter()
    results = []
    for reflection in (False, True):
        for rotation in range(6):
            posed = [transform(p, reflection, rotation) for p in A]
            # Lexicographic minimum commutes with translation; hence one possible
            # translation per D6 pose for exact equality of these finite sets.
            p, q = min(posed), min(B)
            dx, dy = q[0]-p[0], q[1]-p[1]
            posed = [(x+dx, y+dy) for x, y in posed]
            if set(posed) != set(B):
                failures['core_site_set_mismatch'] += 1
                continue
            coremap = {i: lookup[p] for i, p in enumerate(posed)}
            if oriented_set(mapped(faces_of(local, CATEGORIES[0]), coremap, reflection)) != oriented_set(faces_of(official, CATEGORIES[0])):
                failures['core_oriented_faces_mismatch'] += 1
                continue
            direction = -1 if reflection else 1
            for shift in range(32):
                bijection = dict(coremap)
                bijection.update({local['outer_boundary_ids'][j]: official['outer_boundary_ids'][(shift+direction*j) % 32] for j in range(32)})
                if not all(oriented_set(mapped(faces_of(local, cat), bijection, reflection)) == oriented_set(faces_of(official, cat)) for cat in CATEGORIES):
                    failures['annulus_oriented_faces_mismatch'] += 1
                    continue
                check(set(bijection) == set(range(100)) and set(bijection.values()) == set(range(100)), 'invalid_full_bijection')
                check({tuple(sorted((bijection[a], bijection[b]))) for a, b in edges_of(faces_of(local))} == edges_of(faces_of(official)), 'edge_bijection_failure')
                outer_mapped = [bijection[v] for v in local['outer_boundary_ids']]
                ring_mapped = [bijection[v] for v in local['ring_ids']]
                check(cycle_norm(outer_mapped) == cycle_norm(official['outer_boundary_ids'][::direction]), 'outer_cycle_bijection_failure')
                check(cycle_norm(ring_mapped) == cycle_norm(official['ring_ids'][::direction]), 'ring_cycle_bijection_failure')
                officialfans = dict(zip(official['ring_ids'], official['q1_fan_sizes_in_ring_order']))
                check(all(officialfans[bijection[v]] == p for v, p in zip(local['ring_ids'], local['q1_fan_sizes_in_ring_order'])), 'mapped_q1_fan_sizes')
                targetfaces = {triangle_norm(f): i for i, f in enumerate(faces_of(official))}
                facemap = [targetfaces[triangle_norm(f)] for f in mapped(faces_of(local), bijection, reflection)]
                check(sorted(facemap) == list(range(166)), 'face_bijection_failure')
                targetedges = {e: i for i, e in enumerate(sorted(edges_of(faces_of(official))))}
                edgemap = [targetedges[tuple(sorted((bijection[a], bijection[b])))] for a, b in sorted(edges_of(faces_of(local)))]
                check(sorted(edgemap) == list(range(265)), 'edge_index_bijection_failure')
                results.append(dict(reflection=reflection, rotation_steps_60deg=rotation, translation=[dx, dy],
                                    outer_shift=shift, outer_direction=direction, local_to_official_vertex_bijection=bijection,
                                    local_to_official_face_index_bijection=facemap,
                                    local_to_official_sorted_edge_index_bijection=edgemap,
                                    checks='PASS:100 vertices,265 edges,166 oriented/category-preserving faces,outer32 cycle,core26 cycle'))
    return results, dict(failures)


def current_audit(docs, chosen, localgraphs, officialgraphs):
    dev, cur = docs['current_core_integer_development.json'], docs['current_graph_numerical.json']
    i = dev['independent_index']
    local = localgraphs[i]
    m = {int(a): b for a, b in dev['all100_vertex_bijection'].items()}
    check(set(m) == set(range(100)) and set(m.values()) == set(range(100)), 'current_bijection_partition')
    check(face_set(mapped(cur['faces'], m)) == face_set(faces_of(local)), 'current_local166_faces')
    orientation = [r for r in (False, True) if oriented_set(mapped(cur['faces'], m, r)) == oriented_set(faces_of(local))]
    check(len(orientation) == 1, 'current_local_orientation')
    check(canonical_disk(cur['faces'], cur['boundary_cycle']) == cur['fingerprint'] == local['fingerprint'], 'current_fingerprint_replay')
    check(docs['current_graph_numerical.json']['source_sha256'] == sha((INPUT / 'current.npy').read_bytes()), 'current_coordinate_source_hash')
    check(cycle_norm([m[v] for v in cur['boundary_cycle']]) == cycle_norm(local['outer_boundary_ids'][::-1 if orientation[0] else 1]), 'current_local_outer_cycle')
    check(set(m[v] for v in dev['boundary_current_order']) == set(local['ring_ids']), 'current_core_boundary_vertices')
    rows = dev['current_core_vertices']
    check(len(rows) == 68 and {row['catalog_vertex'] for row in rows} == set(range(68)), 'current_core_development_labels')
    for row in rows:
        check(m[row['current_vertex']] == row['catalog_vertex'], 'current_core_development_bijection')
        check(row['axial'] == local['core']['axial_vertices'][row['catalog_vertex']], 'current_integer_axial_labels')
    match = chosen[i]
    om = {int(a): b for a, b in match['local_to_official_vertex_bijection'].items()}
    fullmap = {v: om[m[v]] for v in range(100)}
    family = match['official_family_id']
    official = officialgraphs[family]
    reflection = orientation[0] != match['reflection']
    check(oriented_set(mapped(cur['faces'], fullmap, reflection)) == oriented_set(faces_of(official)), 'current_official_oriented166_faces')
    check({tuple(sorted((fullmap[a], fullmap[b]))) for a, b in edges_of(cur['faces'])} == edges_of(faces_of(official)), 'current_official265_edges')
    check(cycle_norm([fullmap[v] for v in cur['boundary_cycle']]) == cycle_norm(official['outer_boundary_ids'][::-1 if reflection else 1]), 'current_official_outer32_cycle')
    targetfaces = {tuple(sorted(f)): j for j, f in enumerate(faces_of(official))}
    facemap = [targetfaces[tuple(sorted(fullmap[v] for v in f))] for f in cur['faces']]
    check(sorted(facemap) == list(range(166)), 'current_official_face_index_bijection')
    return dict(status='PASS', local_independent_index=i, official_family_id=family,
                parent_reported_family_id=cur.get('reported_family_id'), parent_report_confirmed=cur.get('reported_family_id') == family,
                fingerprint=cur['fingerprint'], current_coordinate_file_sha256=cur['source_sha256'],
                current_to_official_vertex_bijection=fullmap, current_to_official_face_index_bijection=facemap,
                current_to_local_orientation_reversed=orientation[0], current_to_official_orientation_reversed=reflection,
                checks='PASS:100 vertices,265 edges,166 oriented faces,outer32 cycle,68 integer-developed core vertices',
                scope='Exact identity of the stored numerical extracted graph and its integer development; no fresh geometric-DT certificate or radius claim')


def main():
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument('--stage-sources', action='store_true')
    opt = args.parse_args()
    if opt.stage_sources:
        stage()
    docs, receipt = load_sources()
    corecat, index = docs['official_core_catalog.json'], docs['official_graph_index.json']
    check(index['core_catalog_sha256'] == receipt['sources']['official_core_catalog.json']['sha256'], 'official_index_core_catalog_hash')
    cores = {c['core_id']: c for c in corecat['cores']}
    officials = {f'B32_P26_F108_{i:03d}': docs[f'B32_P26_F108_{i:03d}.json'] for i in range(16)}
    locals_ = {g['independent_index']: g for g in docs['local_catalog.json']['graphs']}
    check(set(locals_) == set(range(16)) and set(cores) == set(officials), 'catalog_inventory')
    check({r['family_id'] for r in index['families']} == set(officials), 'official_index_inventory')
    official_validation, local_validation = {}, {}
    for family, graph in officials.items():
        check(graph['family_id'] == graph['core_id'] == family, 'official_family_labels')
        core = cores[family]
        check(core['axial_vertices'] == graph['core_axial_labels'], 'official_core_axial_consistency')
        check(oriented_set(core['triangles']) == oriented_set(faces_of(graph, CATEGORIES[0])), 'official_core_face_consistency')
        check(core['boundary_ids'] == graph['ring_ids'], 'official_core_ring_consistency')
        official_validation[family] = validate_graph(graph, [tuple(p) for p in graph['core_axial_labels']])
        check(next(r['verification'] for r in index['families'] if r['family_id'] == family) == graph['verification'], 'official_index_graph_verification_consistency')
    for i, graph in locals_.items():
        check(oriented_set(graph['core']['triangles']) == oriented_set(faces_of(graph, CATEGORIES[0])), 'local_core_face_consistency')
        check(graph['core']['boundary_ids'] == graph['ring_ids'], 'local_core_ring_consistency')
        local_validation[i] = validate_graph(graph, [tuple(p) for p in graph['core']['axial_vertices']])
    attempts, chosen, ambiguities = [], {}, []
    for i, local in sorted(locals_.items()):
        families = []
        for family, official in sorted(officials.items()):
            matches, failures = find_matches(local, official)
            attempts.append(dict(local_independent_index=i, official_family_id=family, full_matches=len(matches), failure_counts=failures))
            if matches:
                families.append(family)
                chosen[i] = dict(local_independent_index=i, official_family_id=family, exact_match_pose_count=len(matches),
                                 **matches[0], alternate_matching_poses=[{k: v for k, v in m.items() if 'bijection' not in k} for m in matches[1:]])
        if len(families) != 1:
            ambiguities.append(dict(local_independent_index=i, matching_families=families))
    save('all256_pair_attempts.json', dict(status='PASS' if not ambiguities else 'FAIL', attempts=attempts, failed_or_ambiguous_local_rows=ambiguities))
    check(not ambiguities and len({m['official_family_id'] for m in chosen.values()}) == 16, 'catalog_not_one_to_one')
    current = current_audit(docs, chosen, locals_, officials)
    save('current_official_identity.json', current)
    result = dict(status='PASS', scope='Identity bijection of these two supplied16-member catalogs only; no global completeness, optimality, radius, or geometric-DT claim',
                  algorithm='Integer D6 reflection(x,y)->(y,x), rotations(x,y)->(-y,x+y), uniquely determined translation; outer cyclic/reflection maps; category-preserving oriented face/edge/cycle checks',
                  official_archive_sha256=ZIP_SHA, local_catalog_sha256=LOCAL_SHA,
                  graph_count=16, pair_attempts=256, unmatched_or_ambiguous_rows=0,
                  local_validation=local_validation, official_validation=official_validation, mappings=list(chosen.values()),
                  current=dict(status=current['status'], local_independent_index=current['local_independent_index'], official_family_id=current['official_family_id'], fingerprint=current['fingerprint']))
    save('official_catalog_mapping.json', result)
    save('failure_records.json', dict(status='PASS', execution_failures=[], failed_or_ambiguous_local_rows=[],
                                      nonmatching_family_pairs=sum(not r['full_matches'] for r in attempts),
                                      note='Rejected poses and nonmatching pairs are recorded in all256_pair_attempts.json; these are normal identity-search rejections.'))
    print(json.dumps(dict(status='PASS', graph_count=16,
                          mappings={i: m['official_family_id'] for i, m in chosen.items()}, current=result['current']), indent=2))


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        save('failure_records.json', dict(status='FAIL', execution_failures=[dict(type=type(exc).__name__, reason=str(exc))]))
        print(f'FAIL {type(exc).__name__}: {exc}', file=sys.stderr)
        raise
