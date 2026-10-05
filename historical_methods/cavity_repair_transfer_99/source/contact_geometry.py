# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
import engine as run
import numpy as np
from scipy.optimize import linear_sum_assignment
def skeleton_compare(target,source,seeds=144):
    """Explicit99-to99 or99-to100 subset matching, rotation/reflection only."""
    best=None
    for reflection in (False,True):
        V=source.copy()
        if reflection:V[:,1]*=-1
        for initial in np.linspace(0,2*np.pi,seeds,endpoint=False):
            angle=float(initial)
            for it in range(10):
                Q=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]);W=V@Q;rows,cols=linear_sum_assignment(np.sum((target[:,None,:]-W[None,:,:])**2,axis=2));A=target[rows];Z=V[cols];angle=float(np.arctan2(np.sum(A[:,0]*Z[:,1]-A[:,1]*Z[:,0]),np.sum(A*Z)))
            Q=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]);d=np.linalg.norm(target[rows]-V[cols]@Q,axis=1);error=float(np.sum(d*d))
            if best is None or error<best['sum_squared_error']:best=dict(sum_squared_error=error,maximum_binding_center_distance=float(d.max()),rotation_angle=angle,reflection=reflection,target_selected_rows=rows,source_selected_columns=cols,target_unselected_indices=sorted(set(range(len(target)))-set(rows)),coordinate_match_below2e6=bool(d.max()<2e-6))
    best.update(rotation_seeds=seeds,scope='Successful99 coordinate alignment is a geometry witness; failed finite alignment is not a proof of global geometric inequivalence. Radius is not used for matching.');return best

def contact_compare(K,C,wide=False):
    """Explicit binding geometry+real active graph, with inactive completion kept."""
    kd=run.diagnose(K);cd=run.diagnose(C);ku=np.array(kd['binding_indices']);cu=np.array(cd['binding_indices']);kf=np.array(kd['nonbinding_indices'],dtype=int);cf=np.array(cd['nonbinding_indices'],dtype=int)
    if len(ku)!=len(cu) or kd['active_count']!=cd['active_count']:return dict(matches=False,reason='binding_or_active_contact_cardinality_differs')
    normbound=float(np.max(np.abs(np.sort(np.linalg.norm(K[ku],axis=1))-np.sort(np.linalg.norm(C[cu],axis=1)))))
    if normbound>2e-6+1e-12:return dict(matches=False,reason='center_norm_matching_lower_bound_exceeds_geometry_tolerance',norm_matching_lower_bound=normbound)
    if wide:
        sk=skeleton_compare(K[ku],C[cu],144);candidates=[sk]
    else:
        candidates=[]
        for reflection in (False,True):
            V=C[cu].copy()
            if reflection:V[:,1]*=-1
            angle=0.
            for it in range(20):
                Q=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]);rows,cols=linear_sum_assignment(np.sum((K[ku,None,:]-(V@Q)[None,:,:])**2,axis=2));A=K[ku[rows]];Z=V[cols];angle=float(np.arctan2(np.sum(A[:,0]*Z[:,1]-A[:,1]*Z[:,0]),np.sum(A*Z)))
            candidates.append(dict(reflection=reflection,rotation_angle=angle,target_selected_rows=rows,source_selected_columns=cols))
    best=None
    for sk in candidates:
        V=C.copy()
        if sk['reflection']:V[:,1]*=-1
        angle=sk['rotation_angle'];Q=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]);V=V@Q;W=np.empty_like(C);rows=np.array(sk['target_selected_rows']);cols=np.array(sk['source_selected_columns']);W[ku[rows]]=V[cu[cols]]
        if len(kf):
            r,j=linear_sum_assignment(np.sum((K[kf,None,:]-V[cf][None,:,:])**2,axis=2));W[kf[r]]=V[cf[j]]
        wd=run.diagnose(W);d=np.linalg.norm(K-W,axis=1);active=kd['active_signature']==wd['active_signature'];inactive=kd['inactive_signature']==wd['inactive_signature'];match=bool(d[ku].max()<2e-6 and active);result=dict(matches=match,binding_center_count=len(ku),maximum_binding_center_distance=float(d[ku].max()),all100_maximum_center_distance=float(d.max()),active_inventory_matches=active,inactive_inventory_matches=inactive,full100_geometry_match=bool(d.max()<2e-6 and active and inactive),nonbinding_completion_displacements={int(f):float(d[f]) for f in kf},reflection=sk['reflection'],rotation_angle=angle,binding_target_rows=rows,binding_source_columns=cols,search='144-angle rotation/reflection+binding assignment' if wide else 'explicit20-iteration assignment/procrustes near orientation+reflection',scope='A binding contact family requires successful binding geometry and actual active-key match. Nonbinding completion positions and inactive topology remain separate; no disk is treated as deletable.')
        if best is None or result['maximum_binding_center_distance']<best['maximum_binding_center_distance']:best=result
        if match:return result
    return best
