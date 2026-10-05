# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Twelve reverse/cross stable-frame site+wall proposals; no optimizer."""
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[n]='1'
import json,hashlib,shutil
from pathlib import Path
import numpy as np
from current_site_exchange_proposals import ROOT,GRAPHS,catalog_current,prediction_frames,graph_audit
from outer_alignment import optimal_poses,outer_maps
from coupled_wall_proposals import PIN,PINBASE,metric_orientation,align_frame,exact_faces,prior_hashes,csha,unordered


def main():
    hand=json.loads((PIN/'HANDOFF.json').read_text());frames={}
    for r in hand['records']:
        if r['independent_index'] not in frames:frames[r['independent_index']]=r
    arrays={10:catalog_current()[0]}
    for i,r in frames.items():
        p=PINBASE/r['file'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['file_sha256'];arrays[i]=np.load(p)
    folder=ROOT/'reverse_pair_proposals';folder.mkdir(exist_ok=True)
    old,oldu,nadd=prior_hashes();seen={};seenu={};records=[];pose_details=[]
    for src,dst in ((14,10),(15,10),(14,15),(15,14)):
        S=arrays[src];P=arrays[dst];source_orientation=metric_orientation(S,src);target_orientation=metric_orientation(P,dst)
        overlap,poses=optimal_poses(GRAPHS[src],GRAPHS[dst]);assert overlap==67
        decorated=[dict(pose,**outer_maps(GRAPHS[src],GRAPHS[dst],pose)) for pose in poses]
        decorated.sort(key=lambda r:(-r['core_face_overlap'],-r['best_outer_alignments'][0]['common_annulus_edges'],-r['best_outer_alignments'][0]['common_annulus_faces'],r['reflection'],r['rotation'],r['translation']))
        pose=decorated[0];outer=pose['best_outer_alignments'][0]
        mapping={int(a):int(b) for a,b in outer['target_to_source_catalog_bijection'].items()};inverse={b:a for a,b in mapping.items()};common={int(a):int(b) for a,b in pose['core_common_target_to_source'].items()}
        new=pose['added_target_core_ids'][0];wallsource=outer['affected_source_outer_ids'];wall=[inverse[j] for j in wallsource]
        desired=source_orientation*(-1 if pose['reflection'] else 1)
        Y,alignment=align_frame(S,P,common,desired/target_orientation)
        extrapolations=prediction_frames(S,pose['added_target_posed_axial'][0],limit=8,base_index=src)
        affine_checks=[]
        for frame in extrapolations:
            C=S[np.array([mapping[j] for j in range(100)])].copy();C[new]=frame['point']
            a=graph_audit(C,dst,pose['reflection'],source_orientation)
            affine_checks.append(dict(frame=frame,actual_graph=a))
        good=[r for r in affine_checks if r['actual_graph'].get('target_labelled_faces_match') and r['actual_graph'].get('target_fingerprint_match')]
        # Integer frame distance and target graph validity rank; never radius.
        best=good[0] if good else affine_checks[0];extrapolated=np.array(best['frame']['point'])
        pose_details.append(dict(source_index=src,target_index=dst,chosen_pose=pose,affine_predictions=affine_checks,valid_affine_count=len(good)))
        for alpha,beta in ((0.,0.),(0.,1.),(1.,1.)):
            C=S[np.array([mapping[j] for j in range(100)])].copy();C[new]=(1-alpha)*extrapolated+alpha*Y[new]
            C[wall]=(1-beta)*C[wall]+beta*Y[wall]
            projected=[]
            for j in [new]+wall:
                norm=np.linalg.norm(C[j])
                if norm>=1:C[j]*=(1-1e-8)/norm;projected.append(j)
            assert np.array_equal(C[sorted(common)],S[np.array([common[j] for j in sorted(common)])])
            modified=[j for j in range(100) if not np.array_equal(C[j],S[mapping[j]])]
            path=folder/f'from{src:02d}_to{dst:02d}_new{alpha:g}_wall{beta:g}.npy';np.save(path,C);cs=csha(C);us=unordered(C)
            v=graph_audit(C,dst,pose['reflection'],source_orientation)
            row=dict(source_graph_index=src,target_graph_index=dst,source_pinned_frame=frames[src],target_frame_source=frames.get(dst,{'kind':'current_exact_catalog10_frame'}),
                     target_to_source_catalog_bijection=mapping,source_modified_catalog_labels=[mapping[j] for j in modified],target_modified_labels=modified,
                     retained_common_core_centers_bitwise=67,retained_outer_centers_bitwise=32-len(set(modified)&set(range(68,100))),retained_real_centers_bitwise=100-len(modified),
                     exchanged_source_core=pose['removed_source_core_ids'],added_target_core=pose['added_target_core_ids'],exchanged_source_axial=pose['removed_source_axial'],added_posed_axial=pose['added_target_posed_axial'],
                     pose_reflection=pose['reflection'],pose_rotation=pose['rotation'],pose_translation=pose['translation'],source_metric_orientation=source_orientation,
                     outer_shift=outer['outer_shift'],outer_direction=outer['outer_direction'],wall_source_catalog_labels=wallsource,wall_target_labels=wall,
                     target_frame_alignment=alignment,affine_reference=best['frame'],newpoint_targetframe_fraction=alpha,outer_wall_targetframe_fraction=beta,unitdisk_projected_target_labels=projected,
                     coordinate_file=str(path),file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),coordinate_sha256=cs,unordered_geometry_sha256=us,
                     old_ordered_duplicate=cs in old,old_unordered_duplicate=us in oldu,new_ordered_duplicate_of=seen.get(cs),new_unordered_duplicate_of=seenu.get(us),
                     scope='Reverse/cross1site exchange from pinned14/15; all67 common actual core centers fixed, only motif sixwall may move. No optimizer/cover certificate',**v)
            if v.get('target_labelled_faces_match') and v.get('target_fingerprint_match'):row['integer_target_face_audit']=exact_faces(C,dst,desired)
            row['qualified_target_metric_start']=bool(row.get('target_labelled_faces_match') and row.get('target_fingerprint_match') and row.get('integer_target_face_audit',{}).get('exact_target_face_strict_pass') and not(row['old_ordered_duplicate'] or row['old_unordered_duplicate'] or row['new_ordered_duplicate_of'] or row['new_unordered_duplicate_of']))
            seen.setdefault(cs,path.name);seenu.setdefault(us,path.name);records.append(row)
    out=dict(status='COMPLETE',count=len(records),qualified_target_metric_starts=sum(r['qualified_target_metric_start'] for r in records),records=records,pose_details=pose_details)
    (ROOT/'reverse_pair_proposal_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('records','pose_details')},indent=2))
    for r in records:print(Path(r['coordinate_file']).name,'retained',r['retained_real_centers_bitwise'],'actual',r.get('actual_graph_index'),'target',r.get('target_fingerprint_match'),'strict',r.get('integer_target_face_audit',{}).get('exact_target_face_strict_pass'),'new',r['qualified_target_metric_start'],'R',r.get('numerical_dynamic_whole_disk_radius'))

if __name__=='__main__':main()
