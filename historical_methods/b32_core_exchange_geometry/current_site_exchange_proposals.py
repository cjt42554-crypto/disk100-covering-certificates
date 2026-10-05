# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""A few 99-center-retained metric proposals from exact lattice exchange maps.

One center is replaced using affine extrapolation from actual incumbent centers.
No covering optimizer runs.  Initial graphs are explicitly audited, and a target
graph mismatch is retained as a mismatch rather than assigned the intended ID.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[name]='1'
import sys,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
FROZEN=ROOT.parent/'b32_catalog_search/frozen_handoff'
sys.dont_write_bytecode=True;sys.path.insert(0,str(FROZEN))
from graph_tools import extract
CAT=json.loads((FROZEN/'catalog.json').read_text())
GRAPHS={g['independent_index']:g for g in CAT['graphs']}
DEVELOPMENT=json.loads((ROOT.parent/'b32_catalog_search/frozen_support/current_core_integer_development.json').read_text())
SOURCE=ROOT.parent/'b32_catalog_search/frozen_support/inputs/current.npy'


def catalog_current():
    raw=np.load(SOURCE);C=np.empty_like(raw)
    m={int(a):int(b) for a,b in DEVELOPMENT['all100_vertex_bijection'].items()}
    for a,b in m.items():C[b]=raw[a]
    assert extract(C)['fingerprint']==GRAPHS[10]['fingerprint']
    return C,{b:a for a,b in m.items()}


def distance_axial(p,q):
    x,y=p[0]-q[0],p[1]-q[1]
    return x*x+x*y+y*y


def prediction_frames(baseC,q,limit=4,base_index=10):
    xy=np.array(GRAPHS[base_index]['core']['axial_vertices'],int);frames=[]
    for f in GRAPHS[base_index]['core_q0_triangles']:
        p=xy[f];a,b,c=p
        # Exact axial unit triangles have determinant one.
        u=b-a;v=c-a;d=int(u[0]*v[1]-u[1]*v[0]);assert d==1
        w=np.array(q)-a
        coeff=[int(w[0]*v[1]-w[1]*v[0]),int(u[0]*w[1]-u[1]*w[0])]
        pred=baseC[f[0]]+coeff[0]*(baseC[f[1]]-baseC[f[0]])+coeff[1]*(baseC[f[2]]-baseC[f[0]])
        rank=int(sum(distance_axial(tuple(p),tuple(q)) for p in xy[f]))
        frames.append(dict(reference_core_triangle=f,integer_affine_coefficients=coeff,rank=rank,point=pred.tolist()))
    frames.sort(key=lambda x:(x['rank'],x['reference_core_triangle']))
    unique=[]
    for frame in frames:
        if any(np.linalg.norm(np.array(frame['point'])-np.array(v['point']))<1e-12 for v in unique):continue
        unique.append(frame)
        if len(unique)==limit:break
    return unique


def numerical_radius(C,actual):
    T=np.array(actual['faces']);P=C[T];a=P[:,0];u=P[:,1]-a;v=P[:,2]-a;d=u[:,0]*v[:,1]-u[:,1]*v[:,0]
    uu=np.sum(u*u,1);vv=np.sum(v*v,1)
    pr=np.stack((uu*v[:,1]-vv*u[:,1],u[:,0]*vv-v[:,0]*uu),1)/(2*d[:,None])
    tmax=np.sum(pr*pr,1).max()
    anchors=np.array(actual['boundary_anchors']);pairs=np.array(actual['boundary_pairs'])
    bmax=np.sum((anchors-C[pairs[:,0]])**2,1).max()
    return float(np.sqrt(max(tmax,bmax)))


def graph_audit(C,target_index,reflection=False,source_metric_orientation=-1):
    graph=GRAPHS[target_index];T=np.array(graph['faces'])
    expected_orientation=source_metric_orientation*(-1 if reflection else 1)
    if expected_orientation<0:T=T[:,[0,2,1]]
    P=C[T];u=P[:,1]-P[:,0];v=P[:,2]-P[:,0];d=u[:,0]*v[:,1]-u[:,1]*v[:,0]
    row=dict(target_min_oriented_double_area=float(d.min()),target_inverted_faces=int(np.sum(d<=0)),
             maximum_center_norm=float(np.linalg.norm(C,axis=1).max()),intended_graph_index=target_index,
             intended_fingerprint=graph['fingerprint'],expected_metric_orientation=expected_orientation)
    try:
        actual=extract(C,True)
        actualindex=next((j for j,g in GRAPHS.items() if g['fingerprint']==actual['fingerprint']),None)
        row.update(actual_graph_index=actualindex,actual_fingerprint=actual['fingerprint'],actual_B=actual['boundary'],actual_F=actual['restricted_faces'],
                   actual_interior_p_degree_histogram=actual['interior_p_degree_histogram'],actual_q_face_histogram=actual['q_face_histogram'],
                   numerical_dynamic_whole_disk_radius=numerical_radius(C,actual),
                   target_fingerprint_match=actual['fingerprint']==graph['fingerprint'],
                   target_labelled_faces_match={tuple(sorted(f)) for f in actual['faces']}=={tuple(sorted(f)) for f in graph['faces']},status='GRAPH_EXTRACTED')
    except Exception as e:row.update(status='GRAPH_EXTRACTION_FAILED',error=f'{type(e).__name__}: {e}',target_fingerprint_match=False)
    return row


def main():
    C,inverse=catalog_current();outfolder=ROOT/'site_only_proposals';outfolder.mkdir(exist_ok=True)
    baseT=np.array(GRAPHS[10]['core_q0_triangles']);baseP=C[baseT]
    baseu=baseP[:,1]-baseP[:,0];basev=baseP[:,2]-baseP[:,0];based=baseu[:,0]*basev[:,1]-baseu[:,1]*basev[:,0]
    source_orientation=-1 if np.all(based<0) else 1
    assert np.all(source_orientation*based>0),'source_core_orientation_inconsistent'
    (ROOT/'inputs').mkdir(exist_ok=True);np.save(ROOT/'inputs/current_in_catalog10_labels.npy',C)
    exchanges=json.loads((ROOT/'baseline10_outer_exchange_maps.json').read_text());records=[]
    for targetrow in exchanges['targets']:
        target=targetrow['target_index']
        if target not in (14,15):continue
        for posei,pose in enumerate(targetrow['poses']):
            assert targetrow['site_exchange_distance']==1
            align=pose['best_outer_alignments'][0];mapping={int(k):int(v) for k,v in align['target_to_source_catalog_bijection'].items()}
            newtarget=pose['added_target_core_ids'][0];donor=mapping[newtarget];q=pose['added_target_posed_axial'][0]
            for fi,frame in enumerate(prediction_frames(C,q)):
                X=C[np.array([mapping[j] for j in range(100)])].copy();X[newtarget]=frame['point']
                name=f'from10_to{target:02d}_pose{posei:02d}_affine{fi:02d}.npy';path=outfolder/name;np.save(path,X)
                kept=[j for j in range(100) if j!=newtarget]
                assert np.array_equal(X[kept],C[np.array([mapping[j] for j in kept])])
                row=dict(source_index=10,target_index=target,source_file=str(SOURCE),source_file_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                         pose_number=posei,reflection=pose['reflection'],rotation=pose['rotation'],translation=pose['translation'],
                         source_metric_orientation=source_orientation,
                         outer_shift=align['outer_shift'],outer_direction=align['outer_direction'],
                         donor_source_catalog=donor,donor_source_original_label=inverse[donor],added_target_core_label=newtarget,
                         removed_source_axial=pose['removed_source_axial'][0],added_posed_axial=q,
                         retained_common_core_centers_exactly=67,retained_outer_centers_exactly=32,retained_real_centers_total=99,
                         core_to_outer_affected_source_catalog=align['affected_source_outer_ids'],prediction_frame=frame,
                         source_donor_to_new_center_distance=float(np.linalg.norm(X[newtarget]-C[donor])),
                         coordinate_file=str(path),file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                         coordinate_sha256=hashlib.sha256(X.astype('<f8').tobytes()).hexdigest(),
                         scope='Single lattice-site exchange with99 actual centers held bitwise; no optimizer or radius certificate')
                row.update(graph_audit(X,target,pose['reflection'],source_orientation));records.append(row)
    out=dict(status='COMPLETE',count=len(records),target_matching_proposals=sum(r.get('target_fingerprint_match',False) for r in records),
             scope='Initial graph audit of small exact-lattice site exchange proposals; mismatches are not target-family realization claims',records=records)
    (ROOT/'site_only_proposal_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='records'},indent=2))
    for r in records:print(Path(r['coordinate_file']).name,r['status'],'actualindex',r.get('actual_graph_index'),'B',r.get('actual_B'),'F',r.get('actual_F'),'target',r.get('target_fingerprint_match'),'inverted',r['target_inverted_faces'],'Rnum',r.get('numerical_dynamic_whole_disk_radius'))

if __name__=='__main__':main()
