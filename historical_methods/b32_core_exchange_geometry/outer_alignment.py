# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Exact core-lattice pose matching plus outer32 cyclic/reflected label maps.

No radius or coordinate proximity is used for graph alignment.  Faces are
matched through integer axial sites and annulus incidences.
"""
from pathlib import Path
from collections import Counter
import json,hashlib
ROOT=Path(__file__).resolve().parent
CAT=ROOT.parent/'b32_catalog_search/frozen_handoff/catalog.json'


def transform(p,reflection,rotation):
    x,y=p
    if reflection:x,y=y,x
    for _ in range(rotation):x,y=-y,x+y
    return x,y


def edges(faces):return {tuple(sorted((f[j],f[(j+1)%3]))) for f in faces for j in range(3)}


def optimal_poses(base,target):
    B=[tuple(p) for p in base['core']['axial_vertices']];bs=set(B);facesB={frozenset(B[j] for j in f) for f in base['core_q0_triangles']}
    best=0;poses=[]
    for reflection in (False,True):
        for rotation in range(6):
            A=[transform(p,reflection,rotation) for p in target['core']['axial_vertices']]
            shifts=Counter((b[0]-a[0],b[1]-a[1]) for a in A for b in B)
            overlap=max(shifts.values())
            if overlap<best:continue
            if overlap>best:best=overlap;poses=[]
            for (dx,dy),count in sorted(shifts.items()):
                if count!=overlap:continue
                posed=[(x+dx,y+dy) for x,y in A]
                facesT={frozenset(posed[j] for j in f) for f in target['core_q0_triangles']}
                poses.append(dict(reflection=reflection,rotation=rotation,translation=[dx,dy],site_overlap=count,
                                  core_face_overlap=len(facesB&facesT),target_posed_axial=posed))
    return best,poses


def outer_maps(base,target,pose):
    B=[tuple(p) for p in base['core']['axial_vertices']];lookup={p:j for j,p in enumerate(B)}
    A=[tuple(p) for p in pose['target_posed_axial']]
    core_common={j:lookup[p] for j,p in enumerate(A) if p in lookup}
    missingT=[j for j in range(68) if j not in core_common]
    removedB=sorted(set(range(68))-set(core_common.values()),key=lambda j:B[j])
    assigned=dict(core_common);assigned.update(zip(sorted(missingT,key=lambda j:A[j]),removedB))
    annulusB=base['ring_q0_triangles']+base['q1_triangles'];annulusT=target['ring_q0_triangles']+target['q1_triangles']
    bfaces={tuple(sorted(f)) for f in annulusB};bedges={e for e in edges(base['faces']) if e[1]>=68}
    rows=[];direction=-1 if pose['reflection'] else 1
    for shift in range(32):
        outer={68+j:68+(shift+direction*j)%32 for j in range(32)}
        common=dict(core_common);common.update(outer)
        full=dict(assigned);full.update(outer)
        sharedfaces={tuple(sorted(common[v] for v in f)) for f in annulusT if set(f)<=set(common)}&bfaces
        sharededges={tuple(sorted((common[a],common[b]))) for a,b in edges(target['faces']) if b>=68 and a in common and b in common}&bedges
        fullann={tuple(sorted(full[v] for v in f)) for f in annulusT}
        fullcross={tuple(sorted((full[a],full[b]))) for a,b in edges(target['faces']) if a<68<=b}
        bcross={e for e in edges(base['faces']) if e[0]<68<=e[1]}
        changed_annulus_faces=sorted(fullann^bfaces)
        changed_cross_edges=sorted(fullcross^bcross)
        affected_outer=sorted({v for f in changed_annulus_faces for v in f if v>=68})
        rows.append(dict(outer_shift=shift,outer_direction=direction,common_annulus_faces=len(sharedfaces),common_annulus_edges=len(sharededges),
                         annulus_face_symmetric_difference=len(changed_annulus_faces),cross_edge_symmetric_difference=len(changed_cross_edges),
                         affected_source_outer_ids=affected_outer,changed_annulus_faces_in_source_labels=changed_annulus_faces,
                         changed_cross_edges_in_source_labels=changed_cross_edges,target_to_source_catalog_bijection=full))
    # Retained common-site incidences rank ahead of arbitrary replacement labels.
    rows.sort(key=lambda r:(-r['common_annulus_edges'],-r['common_annulus_faces'],r['annulus_face_symmetric_difference'],r['outer_shift']))
    bestscore=tuple(rows[0][k] for k in ('common_annulus_edges','common_annulus_faces','annulus_face_symmetric_difference'))
    optimal=[r for r in rows if tuple(r[k] for k in ('common_annulus_edges','common_annulus_faces','annulus_face_symmetric_difference'))==bestscore]
    return dict(core_common_target_to_source=core_common,removed_source_core_ids=removedB,added_target_core_ids=missingT,
                removed_source_axial=[B[j] for j in removedB],added_target_posed_axial=[A[j] for j in missingT],
                best_outer_alignments=optimal,all_outer_alignment_scores=[{k:v for k,v in r.items() if k not in ('target_to_source_catalog_bijection','changed_annulus_faces_in_source_labels','changed_cross_edges_in_source_labels')} for r in rows])


def main():
    cat=json.loads(CAT.read_text());graphs={g['independent_index']:g for g in cat['graphs']};B=graphs[10];results=[]
    for index,target in sorted(graphs.items()):
        if index==10:continue
        overlap,poses=optimal_poses(B,target)
        rows=[]
        for pose in poses:
            details=outer_maps(B,target,pose)
            rows.append(dict(pose,**details))
        rows.sort(key=lambda r:(-r['core_face_overlap'],-r['best_outer_alignments'][0]['common_annulus_edges'],-r['best_outer_alignments'][0]['common_annulus_faces'],r['reflection'],r['rotation'],r['translation']))
        results.append(dict(source_index=10,target_index=index,site_exchange_distance=68-overlap,maximum_common_sites=overlap,
                            optimal_lattice_pose_count=len(rows),poses=rows))
    out=dict(status='PASS',catalog_sha256=hashlib.sha256(CAT.read_bytes()).hexdigest(),source_index=10,
             scope='Exact D6+translation core site overlap and integer-incidence outer cyclic/reflection maps; no metric or graph-completeness claim',targets=results)
    (ROOT/'baseline10_outer_exchange_maps.json').write_text(json.dumps(out,indent=2)+'\n')
    for r in results:
        p=r['poses'][0];a=p['best_outer_alignments'][0]
        print(r['target_index'],'distance',r['site_exchange_distance'],'corefaces',p['core_face_overlap'],'annulusedges',a['common_annulus_edges'],'annulusfaces',a['common_annulus_faces'],'outershift',a['outer_shift'],'direction',a['outer_direction'],'affectedouter',a['affected_source_outer_ids'],'remove',p['removed_source_axial'],'add',p['added_target_posed_axial'])

if __name__=='__main__':main()
