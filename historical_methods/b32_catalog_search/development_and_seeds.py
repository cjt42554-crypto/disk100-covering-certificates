# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Exact graph-to-lattice bijection and small harmonic metric seed collection.

All coordinates are numerical seed proposals; this file never runs a coverage
optimizer.  Embedded target faces and actual restricted-DT faces are audited
separately.  Catalog names remain independent indices until parent verification.
"""
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[n]='1'
import json,hashlib
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from scipy.spatial import Delaunay
from graph_tools import extract
ROOT=Path(__file__).resolve().parent


def edge_set(faces):
    return sorted({tuple(sorted((f[j],f[(j+1)%3]))) for f in faces for j in range(3)})


def harmonic(graph, outer, core_weight=1., annulus_weight=1.):
    n=100;boundary=graph['outer_boundary_ids'];inside=graph['interior_ids']
    assert inside==list(range(68)) and boundary==list(range(68,100))
    A=np.zeros((68,68));b=np.zeros((68,2))
    for a,c in edge_set(graph['faces']):
        weight=core_weight if a<68 and c<68 else annulus_weight
        for u,v in ((a,c),(c,a)):
            if u>=68:continue
            A[u,u]+=weight
            if v<68:A[u,v]-=weight
            else:b[u]+=weight*outer[v-68]
    C=np.zeros((100,2));C[:68]=np.linalg.solve(A,b);C[68:]=outer
    return C


def embedding_audit(C,graph):
    T=np.array(graph['faces']);P=C[T];u=P[:,1]-P[:,0];v=P[:,2]-P[:,0]
    det=u[:,0]*v[:,1]-u[:,1]*v[:,0]
    lengths=np.linalg.norm(C[np.array(edge_set(graph['faces']))[:,0]]-C[np.array(edge_set(graph['faces']))[:,1]],axis=1)
    actual=extract(C,True)
    target={tuple(sorted(f)) for f in graph['faces']}
    got={tuple(sorted(f)) for f in actual['faces']}
    # geometry triangles are CCW for all valid generated target orientations.
    a=P[:,0];uu=np.sum(u*u,1);vv=np.sum(v*v,1)
    centers=a+np.stack((uu*v[:,1]-vv*u[:,1],u[:,0]*vv-v[:,0]*uu),1)/(2*det[:,None])
    targetR=np.sqrt(np.sum((centers-a)**2,1))
    actualT=np.array(actual['faces']);actualP=C[actualT];aa=actualP[:,0];au=actualP[:,1]-aa;av=actualP[:,2]-aa
    adet=au[:,0]*av[:,1]-au[:,1]*av[:,0];auu=np.sum(au*au,1);avv=np.sum(av*av,1)
    ap=aa+np.stack((auu*av[:,1]-avv*au[:,1],au[:,0]*avv-av[:,0]*auu),1)/(2*adet[:,None])
    interiorR2=np.sum((ap-aa)**2,1)
    anchor=np.array(actual['boundary_anchors']);pairs=np.array(actual['boundary_pairs'])
    boundaryR2=np.sum((anchor-C[pairs[:,0]])**2,1)
    return dict(target_min_oriented_double_area=float(det.min()),target_inverted_faces=int(np.sum(det<=0)),
                target_max_circumradius=float(targetR.max()),target_outside_circumcenters=int(np.sum(np.sum(centers*centers,1)>1)),
                target_edge_length_min=float(lengths.min()),target_edge_length_max=float(lengths.max()),
                restricted_DT_target_faces_equal=got==target,restricted_DT_face_symmetric_difference=len(got^target),
                actual_fingerprint=actual['fingerprint'],actual_B=actual['boundary'],actual_F=actual['restricted_faces'],
                actual_target_fingerprint_equal=actual['fingerprint']==graph['fingerprint'],max_center_norm=float(np.linalg.norm(C,axis=1).max()),
                actual_dynamic_whole_disk_numerical_radius=float(np.sqrt(max(interiorR2.max(),boundaryR2.max()))))


def development(current,graph):
    can_to_catalog={int(v):int(k) for k,v in graph['vertex_to_canonical'].items()}
    current_to_catalog={int(k):can_to_catalog[int(v)] for k,v in current['vertex_to_canonical'].items()}
    inverse={v:k for k,v in current_to_catalog.items()}
    mapped={tuple(sorted(current_to_catalog[x] for x in f)) for f in current['faces']}
    assert mapped=={tuple(sorted(f)) for f in graph['faces']}
    core=graph['core'];xy=np.array(core['axial_vertices'],int)
    core_edges=edge_set(core['triangles'])
    assert all((xy[a,0]-xy[b,0])**2+(xy[a,0]-xy[b,0])*(xy[a,1]-xy[b,1])+(xy[a,1]-xy[b,1])**2==1 for a,b in core_edges)
    for f in core['triangles']:
        a,b,c=xy[f];assert np.linalg.det(np.array([b-a,c-a]))>.999
    ring=core['boundary_ids'];steps=[tuple(map(int,xy[ring[(j+1)%26]]-xy[v])) for j,v in enumerate(ring)]
    dirs=[(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
    dindices=[dirs.index(s) for s in steps]
    turns=[((dindices[j]-dindices[j-1]+3)%6)-3 for j in range(26)]
    hist=Counter(turns);assert sum(turns)==6 and set(hist)<={-1,0,1}
    deg=Counter()
    for a,b in core_edges:deg[a]+=1;deg[b]+=1
    p=dict(zip(ring,graph['q1_fan_sizes_in_ring_order']))
    rows=[dict(current_vertex=inverse[v],catalog_vertex=v,axial=list(map(int,xy[v])),core_degree=deg[v],p=p.get(v,0),is_core_boundary=v in p) for v in range(68)]
    out=dict(status='PASS',scope='Integer triangular-lattice development of the numerical extracted restricted-DT core via exact face bijection; not an exact geometric-DT certificate',
             independent_index=graph['independent_index'],all100_vertex_bijection=current_to_catalog,all166_faces_bijection_check='PASS',
             core_vertices=68,core_faces=108,core_edges=len(core_edges),core_boundary=26,
             core_interior_vertices=42,core_boundary_degree_histogram=dict(sorted(Counter(deg[v] for v in ring).items())),
             boundary_catalog_order=ring,boundary_current_order=[inverse[v] for v in ring],boundary_step_directions=steps,
             boundary_turn_units_60deg=turns,boundary_turn_histogram=dict(sorted(hist.items())),boundary_total_turn_deg=60*sum(turns),
             current_core_vertices=rows)
    (ROOT/'current_core_integer_development.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


def main():
    data=json.loads((ROOT/'independent_catalog_numerical_index.json').read_text())
    current=json.loads((ROOT/'current_graph_numerical.json').read_text())
    graph=next(g for g in data['graphs'] if g['fingerprint']==current['fingerprint'])
    dev=development(current,graph)
    print('DEVELOPMENT',json.dumps({k:v for k,v in dev.items() if k not in ('all100_vertex_bijection','current_core_vertices')},separators=(',',':')))
    folder=ROOT/'metric_seeds';folder.mkdir(exist_ok=True)
    records=[]
    for g in data['graphs']:
        xy=np.array(g['core']['axial_vertices'],float)
        E=np.column_stack((xy[:,0]+xy[:,1]/2,np.sqrt(3)*xy[:,1]/2));E-=E.mean(axis=0)
        adj=defaultdict(set)
        for a,b in edge_set(g['faces']):adj[a].add(b);adj[b].add(a)
        raw=np.array([np.mean(E[sorted(adj[v]&set(range(68)))],axis=0) for v in range(68,100)])
        angles=np.unwrap(np.arctan2(raw[:,1],raw[:,0]))
        if np.any(np.diff(angles)<=0):raise ValueError('nonmonotone_radial_outer_angles')
        for mode in ('uniform','radial'):
            theta=np.arange(32)*2*np.pi/32+angles[0] if mode=='uniform' else angles
            outer=.94*np.column_stack((np.cos(theta),np.sin(theta)))
            for weight in (0.7,1.,1.4):
                C=harmonic(g,outer,weight,1.)
                name=f'independent_{g["independent_index"]:02d}_{mode}_coreweight{weight:g}.npy'
                path=folder/name;np.save(path,C)
                row=dict(independent_index=g['independent_index'],method='weighted harmonic',boundary_mode=mode,core_edge_weight=weight,
                         target_fingerprint=g['fingerprint'],coordinate_file=str(path),coordinate_sha256=hashlib.sha256(C.astype('<f8').tobytes()).hexdigest(),
                         file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),numerical_only=True)
                try:row.update(embedding_audit(C,g));row['status']='PASS'
                except Exception as e:row.update(status='FAIL',error=f'{type(e).__name__}: {e}')
                records.append(row)
    out=dict(scope='96 harmonic metric proposals over16 independent reconstructed graph classes; no coverage optimization or catalog completeness claim',
             count=len(records),unfolded_inside_unitdisk=sum(r.get('target_inverted_faces',1)==0 and r.get('max_center_norm',2)<1 for r in records),
             matching_actual_DT=sum(r.get('actual_target_fingerprint_equal',False) for r in records),records=records)
    (ROOT/'metric_seed_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print('SEEDS',json.dumps({k:v for k,v in out.items() if k!='records'}))

if __name__=='__main__':main()
