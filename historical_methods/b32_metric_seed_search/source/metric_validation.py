# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Target metric checks and actual restricted-DT graph, all numerical only."""
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[n]='1'
import numpy as np
import graph_tools
import geometry


def validate(C, graph):
    C=np.asarray(C,dtype=float)
    T=np.array(graph['faces'],dtype=int)
    P=C[T];u=P[:,1]-P[:,0];v=P[:,2]-P[:,0]
    det=u[:,0]*v[:,1]-u[:,1]*v[:,0]
    tg,q,_=geometry.triangles(C,T)
    d2=np.sum((q[:,None,:]-C[None,:,:])**2,axis=2)
    own=d2[np.arange(len(T)),T[:,0]]
    for j in range(3): d2[np.arange(len(T)),T[:,j]]=np.inf
    empty=np.min(d2,axis=1)-own
    result=dict(numerical_only=True,target_fingerprint=graph['fingerprint'],
        target_face_count=len(T),target_min_positive_double_area=float(det.min()),
        target_inverted_faces=int(np.sum(det<=0)),
        target_max_circumradius=float(np.sqrt(tg.max())),
        target_outside_circumcenters=int(np.sum(np.sum(q*q,axis=1)>1+1e-12)),
        target_min_circumcenter_circle_margin=float(np.min(1-np.sum(q*q,axis=1))),
        target_min_empty_squared_margin=float(empty.min()),
        target_nonempty_faces=int(np.sum(empty < -2e-11)),
        max_center_norm=float(np.linalg.norm(C,axis=1).max()))
    try:
        actual=graph_tools.extract(C,True)
        actualfaces={tuple(sorted(t))for t in actual['faces']}
        targetfaces={tuple(sorted(t))for t in T}
        expectededges={tuple(sorted((68+i,68+(i+1)%32)))for i in range(32)}
        actualedges={tuple(sorted(e))for e in actual['boundary_pairs']}
        result.update(actual_graph_status='PASS',actual_fingerprint=actual['fingerprint'],
            actual_B=actual['boundary'],actual_F=actual['restricted_faces'],
            actual_target_fingerprint_equal=actual['fingerprint']==graph['fingerprint'],
            actual_target_faces_equal=actualfaces==targetfaces,
            actual_target_face_symmetric_difference=len(actualfaces^targetfaces),
            actual_outer_circle_arc_edges_equal=expectededges==actualedges,
            actual_circle_arc_edge_difference=[list(e)for e in sorted(expectededges^actualedges)],
            actual_inside_outside_circumcenter_margin=actual['inside_outside_circumcenter_margin'],
            actual_graph=actual)
    except Exception as e:
        result.update(actual_graph_status='FAIL',actual_graph_error=str(e),
            actual_target_fingerprint_equal=False,actual_target_faces_equal=False,
            actual_outer_circle_arc_edges_equal=False)
    result['target_metric_realized']=bool(
        result['target_inverted_faces']==0 and result['target_outside_circumcenters']==0
        and result['target_nonempty_faces']==0 and result['actual_target_faces_equal']
        and result['actual_target_fingerprint_equal'] and result['actual_outer_circle_arc_edges_equal'])
    return result
