# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Few site-plus-wall metric proposals using pinned, accepted numeric frames.

Preserves all 67 common actual core centers. Only the exchanged site and the six
outer vertices identified by integer face/incidence differences may move.
Neither a covering optimizer nor a metric optimizer is run.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[name]='1'
import json,hashlib,shutil
from pathlib import Path
import numpy as np
from current_site_exchange_proposals import ROOT,GRAPHS,catalog_current,graph_audit
PINBASE=ROOT.parent/'b32_metric_seed_search'
PIN=PINBASE/'final_handoff'


def csha(C):return hashlib.sha256(np.asarray(C,dtype='<f8').tobytes()).hexdigest()
def unordered(C):return csha(C[np.lexsort((C[:,1],C[:,0]))])


def metric_orientation(C,index):
    T=np.array(GRAPHS[index]['core_q0_triangles']);P=C[T]
    u=P[:,1]-P[:,0];v=P[:,2]-P[:,0];d=u[:,0]*v[:,1]-u[:,1]*v[:,0]
    assert np.all(d>0) or np.all(d<0)
    return 1 if np.all(d>0) else -1


def align_frame(source,targetC,common,desired_det):
    ids=sorted(common);A=source[np.array([common[j] for j in ids])];B=targetC[ids]
    am=A.mean(0);bm=B.mean(0);a=A-am;b=B-bm
    U,s,Vt=np.linalg.svd(b.T@a)
    adjust=np.eye(2);adjust[1,1]=desired_det*np.linalg.det(U@Vt)
    rotation=U@adjust@Vt
    scale=float(np.sum((b@rotation)*a)/np.sum(b*b))
    offset=am-scale*(bm@rotation);Y=scale*(targetC@rotation)+offset
    error=np.linalg.norm(Y[ids]-A,axis=1)
    return Y,dict(method='Common67 labeled-center similarity Procrustes; prescribed reflection from exact lattice pose and source metric orientation',
                  matched_core_count=len(ids),rotation_matrix=rotation.tolist(),determinant=float(np.linalg.det(rotation)),scale=scale,translation=offset.tolist(),
                  common_core_residual_rms=float(np.sqrt(np.mean(error*error))),common_core_residual_max=float(error.max()),
                  aligned_full_frame_max_center_norm=float(np.linalg.norm(Y,axis=1).max()))


def prior_hashes():
    ordered=json.loads((PINBASE/'prior_ordered_hashes.json').read_text());uns=json.loads((PINBASE/'prior_unordered_hashes.json').read_text())
    additions=0
    for pattern in ('workers/worker*/finals/*.npy','workers/worker*/starts/*.npy'):
        for path in PINBASE.glob(pattern):
            C=np.load(path);ordered.setdefault(csha(C),[]).append(str(path));uns.setdefault(unordered(C),[]).append(str(path));additions+=1
    return ordered,uns,additions


def exact_faces(C,target,orientation):
    ratios=[[float(x).as_integer_ratio() for x in row] for row in C];Q=max(d for row in ratios for _,d in row)
    pts=[tuple(n*(Q//d) for n,d in row) for row in ratios]
    neg=0;ties=0;inverted=0;outside=0;checks=0
    dot=lambda a,b:a[0]*b[0]+a[1]*b[1]
    for face in GRAPHS[target]['faces']:
        f=list(face)
        if orientation<0:f[1],f[2]=f[2],f[1]
        a,b,c=(pts[j] for j in f);u=(b[0]-a[0],b[1]-a[1]);v=(c[0]-a[0],c[1]-a[1]);D=2*(u[0]*v[1]-u[1]*v[0])
        if D<=0:inverted+=1;continue
        uu=dot(u,u);vv=dot(v,v);pr=(uu*v[1]-vv*u[1],u[0]*vv-v[0]*uu)
        N=(a[0]*D+pr[0],a[1]*D+pr[1]);r2=dot(pr,pr)
        if dot(N,N)>=(Q*D)**2:outside+=1
        for j,p in enumerate(pts):
            if j in f:continue
            w=(p[0]*D-N[0],p[1]*D-N[1]);sign=dot(w,w)-r2;checks+=1
            neg+=sign<0;ties+=sign==0
    center_inside=all(dot(p,p)<Q*Q for p in pts)
    return dict(exact_target_face_strict_pass=not(neg or ties or inverted or outside) and center_inside,
                integer_emptycircle_checks=checks,negative_other_site_emptycircle_signs=neg,extra_site_emptycircle_ties=ties,
                nonpositive_oriented_faces=inverted,nonstrict_inside_circumcenter_faces=outside,all_centers_exact_strictly_inside_unitdisk=center_inside,
                scope='Target166-face integer signs; complete actual graph and circle arcs additionally checked numerically')


def main():
    ready=json.loads((PIN/'READY.json').read_text());assert ready['status']=='SEALED'
    manifest=json.loads((PIN/'MANIFEST_SHA256.json').read_text())
    assert hashlib.sha256((PIN/'MANIFEST_SHA256.json').read_bytes()).hexdigest()=='50e8990b6d02ef2a644fe8795d2a59bcf446d881470a29f0b9a36831e45c49b7'
    assert all(hashlib.sha256((PIN/f).read_bytes()).hexdigest()==s for f,s in manifest.items())
    handoff=json.loads((PIN/'HANDOFF.json').read_text());source,inverse=catalog_current();source_orientation=metric_orientation(source,10)
    exchanges=json.loads((ROOT/'baseline10_outer_exchange_maps.json').read_text());targetrows={r['target_index']:r for r in exchanges['targets']}
    affine=json.loads((ROOT/'site_only_proposal_audit.json').read_text())['records']
    old,oldu,nadd=prior_hashes();folder=ROOT/'coupled_proposals';folder.mkdir(exist_ok=True)
    pininputs=ROOT/'inputs/pinned_frames';pininputs.mkdir(parents=True,exist_ok=True)
    records=[];seen={};seen_unordered={}
    for frame in handoff['records']:
        target=frame['independent_index'];pinpath=PINBASE/frame['file'];assert hashlib.sha256(pinpath.read_bytes()).hexdigest()==frame['file_sha256']
        pinC=np.load(pinpath);assert frame['final_label_to_catalog_label']==list(range(100));assert graph_audit(pinC,target,False,metric_orientation(pinC,target))['target_labelled_faces_match']
        dst=pininputs/pinpath.name;shutil.copy2(pinpath,dst)
        pose=targetrows[target]['poses'][0];outer=pose['best_outer_alignments'][0]
        mapping={int(a):int(b) for a,b in outer['target_to_source_catalog_bijection'].items()};common={int(a):int(b) for a,b in pose['core_common_target_to_source'].items()}
        newtarget=pose['added_target_core_ids'][0];modified_outer_source=outer['affected_source_outer_ids'];inversemap={b:a for a,b in mapping.items()}
        walltarget=[inversemap[j] for j in modified_outer_source]
        desired_orientation=source_orientation*(-1 if pose['reflection'] else 1)
        Y,alignment=align_frame(source,pinC,common,desired_orientation/metric_orientation(pinC,target))
        af=next(r for r in affine if r['target_index']==target and r['pose_number']==0 and r['target_labelled_faces_match'])
        extrapolated=np.array(af['prediction_frame']['point'])
        for alpha in (0.,1.):
            for beta in (0.,0.5,1.):
                C=source[np.array([mapping[j] for j in range(100)])].copy()
                C[newtarget]=(1-alpha)*extrapolated+alpha*Y[newtarget]
                C[walltarget]=(1-beta)*C[walltarget]+beta*Y[walltarget]
                projected=[]
                for j in [newtarget]+walltarget:
                    norm=np.linalg.norm(C[j])
                    if norm>=1:C[j]*=(1-1e-8)/norm;projected.append(j)
                keptcore=sorted(common);assert np.array_equal(C[keptcore],source[np.array([common[j] for j in keptcore])])
                fixedouter=sorted(set(range(68,100))-set(walltarget));assert np.array_equal(C[fixedouter],source[np.array([mapping[j] for j in fixedouter])])
                name=f'from10_to{target:02d}_pin{frame["grid_id"]:04d}_new{alpha:g}_wall{beta:g}.npy';path=folder/name;np.save(path,C)
                cs=csha(C);us=unordered(C)
                changed=[j for j in range(100) if not np.array_equal(C[j],source[mapping[j]])]
                audit=graph_audit(C,target,pose['reflection'],source_orientation)
                row=dict(source_graph_index=10,target_graph_index=target,source_original_label_to_catalog=DEVELOPMENT_MAP if False else None,
                         target_to_source_catalog_bijection=mapping,target_to_source_original_label={j:inverse[mapping[j]] for j in range(100)},
                         source_modified_catalog_labels=[mapping[j] for j in changed],source_modified_original_labels=[inverse[mapping[j]] for j in changed],
                         target_modified_labels=changed,retained_common_core_centers_bitwise=67,retained_outer_centers_bitwise=32-len(set(changed)&set(range(68,100))),
                         retained_real_centers_bitwise=100-len(changed),exchanged_source_core=pose['removed_source_core_ids'],added_target_core=pose['added_target_core_ids'],
                         exchanged_source_axial=pose['removed_source_axial'],added_posed_axial=pose['added_target_posed_axial'],
                         pose_reflection=pose['reflection'],pose_rotation=pose['rotation'],pose_translation=pose['translation'],source_metric_orientation=source_orientation,
                         outer_shift=outer['outer_shift'],outer_direction=outer['outer_direction'],wall_source_catalog_labels=modified_outer_source,wall_target_labels=walltarget,
                         pinframe_file=str(dst),pinframe_source_file=str(pinpath),pinframe_file_sha256=frame['file_sha256'],pinframe_prior_independent_alltriple_acceptance=frame['independent_alltriple_acceptance'],
                         pinframe_alignment=alignment,newpoint_pin_fraction=alpha,outer_wall_pin_fraction=beta,unitdisk_projected_target_labels=projected,
                         affine_reference=af['prediction_frame'],coordinate_file=str(path),file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),coordinate_sha256=cs,unordered_geometry_sha256=us,
                         old_ordered_duplicate=cs in old,old_unordered_duplicate=us in oldu,new_ordered_duplicate_of=seen.get(cs),new_unordered_duplicate_of=seen_unordered.get(us),
                         scope='Exact integer-pose coupled1site+6outer metric proposal; all67 common actual core centers held bitwise. No optimizer, coverage certificate, or new-record assertion',**audit)
                if audit.get('target_labelled_faces_match') and audit.get('target_fingerprint_match'):
                    row['integer_target_face_audit']=exact_faces(C,target,desired_orientation)
                row['qualified_new_target_start']=bool(row.get('target_labelled_faces_match') and row.get('target_fingerprint_match') and row.get('integer_target_face_audit',{}).get('exact_target_face_strict_pass') and not(row['old_ordered_duplicate'] or row['old_unordered_duplicate'] or row['new_ordered_duplicate_of'] or row['new_unordered_duplicate_of']))
                records.append(row);seen.setdefault(cs,name);seen_unordered.setdefault(us,name)
    out=dict(status='COMPLETE',count=len(records),qualified_new_target_starts=sum(r['qualified_new_target_start'] for r in records),
             prior_ordered_hashes_count=len(old),prior_unordered_hashes_count=len(oldu),added_sealed87_final_and_start_arrays=nadd,
             pinned_handoff_manifest_sha256=ready['manifest_sha256'],records=records)
    (ROOT/'coupled_proposal_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='records'},indent=2))
    for r in records:print(Path(r['coordinate_file']).name,'retained',r['retained_real_centers_bitwise'],'actual',r.get('actual_graph_index'),'target',r.get('target_fingerprint_match'),'strict',r.get('integer_target_face_audit',{}).get('exact_target_face_strict_pass'),'new',r['qualified_new_target_start'],'R',r.get('numerical_dynamic_whole_disk_radius'))

if __name__=='__main__':main()
