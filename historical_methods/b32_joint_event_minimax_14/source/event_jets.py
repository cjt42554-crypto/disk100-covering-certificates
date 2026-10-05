# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Independent event jets in all200 center coordinates; t derivative is zero.

DT event sign: third-site circle clearance is positive when fourth site is
outside a CCW support triangle's circumcircle. Circle event sign: positive
means its circumcenter lies outside the unit disk. Scaling is fixed at the
source, never a varying normalization whose derivative is silently omitted.
"""
import numpy as np
import branches

def triangle_jet(C,ids):
    ids=np.asarray(ids,dtype=int);P=C[ids];u=P[1]-P[0];v=P[2]-P[0];orientation=float(np.cross(u,v))
    if abs(orientation)<1e-14:raise ValueError('degenerate support triangle')
    r2,q,grad=branches.geom.triangles(C,ids[None,:]);q=q[0];r2=float(r2[0]);j=np.zeros(200);j[(2*ids[:,None]+[0,1]).ravel()]=grad[0].ravel()
    B=P[1:]-P[0];rhs=np.zeros((2,6))
    for k in (1,2):rhs[k-1,:2]=-(P[0]-q);rhs[k-1,2*k:2*k+2]=P[k]-q
    localQ=np.linalg.solve(B,rhs);Q=np.zeros((2,200));Q[:,(2*ids[:,None]+[0,1]).ravel()]=localQ
    oj=np.zeros(200);oj[2*ids[1]:2*ids[1]+2]=[v[1],-v[0]];oj[2*ids[2]:2*ids[2]+2]=[-u[1],u[0]];oj[2*ids[0]:2*ids[0]+2]=[u[1]-v[1],v[0]-u[0]]
    return dict(radius_squared=r2,center=q,radius_jacobian=j,center_jacobian=Q,orientation=orientation,orientation_jacobian=oj)

def event_jet(C,kind,labels,scale=1.):
    labels=list(map(int,labels));tri=labels[:3];t=triangle_jet(C,tri);q=t['center'];Q=t['center_jacobian']
    if kind in ('dt_incircle_zero','DT_COCIRCULAR_QUAD','incircle','dt_determinant','standard_incircle_determinant'):
        assert len(labels)==4;d=C[labels[3]]-q;value=float(d@d-t['radius_squared']);J=-2*d@Q-t['radius_jacobian'];J[2*labels[3]:2*labels[3]+2]+=2*d
        identity='clearance=|fourth-circumcenter|²-radius²; positive outside; incircleDet=-orientation*clearance'
        if kind in ('dt_determinant','standard_incircle_determinant'):
            J=-t['orientation']*J-value*t['orientation_jacobian'];value=-t['orientation']*value
            identity='standard det[x,y,x²+y²,1]=-orientation*circle_clearance; positive inside when supporttriangleCCW'
    elif kind in ('circle_circumcenter','CIRCUMCENTER_UNIT_CIRCLE','circumcenter_unit_circle','circumcenter_unit_polynomial'):
        assert len(labels)==3;value=float(q@q-1.);J=2*q@Q;identity='circumcenter_norm²-1; positive outside closed unit disk'
        if kind=='circumcenter_unit_polynomial':
            J=-4*t['orientation']**2*J-8*t['orientation']*value*t['orientation_jacobian'];value=-4*t['orientation']**2*value
            identity='4orientation²*(1-circumcenter_norm²); positive inside closed unit disk'
    else:raise ValueError('unsupported event kind '+str(kind))
    return dict(raw_value=value,raw_jacobian=J,scaled_value=value/scale,scaled_jacobian=J/scale,normalization_constant=float(scale),gradient_norm=float(np.linalg.norm(J)),definition=identity,support_triangle=tri,fourth_site=labels[3] if len(labels)==4 else None,center=q,radius_squared=t['radius_squared'],orientation=t['orientation'],epigraph_event_row=np.r_[J/scale,0.])

def differential_dependency(C,kind,labels):
    e=event_jet(C,kind,labels);n=e['raw_jacobian'];result=[]
    if len(labels)==4:
        triples=[labels[:3],[labels[0],labels[1],labels[3]],[labels[0],labels[2],labels[3]],[labels[1],labels[2],labels[3]]]
        jets=[triangle_jet(C,t) for t in triples];reference=jets[0]
        for ids,j in zip(triples[1:],jets[1:]):
            difference=j['radius_jacobian']-reference['radius_jacobian'];factor=float(difference@n/(n@n));error=difference-factor*n
            rows=np.vstack((np.r_[reference['radius_jacobian'],-1],np.r_[j['radius_jacobian'],-1],np.r_[n,0]))
            result.append(dict(other_support=ids,event_value=e['raw_value'],radius_value_difference=j['radius_squared']-reference['radius_squared'],gradient_dependency_factor=factor,gradient_dependency_error_norm=float(np.linalg.norm(error)),gradient_difference_norm=float(np.linalg.norm(difference)),augmented_three_row_singular_values=np.linalg.svd(rows,compute_uv=False).tolist(),expected_rank_at_exact_event_at_most=2))
    else:
        reference=triangle_jet(C,labels)
        for pair in ([labels[0],labels[1]],[labels[1],labels[2]],[labels[2],labels[0]]):
            E=np.array([sorted(pair)]);best=None
            for sign in (-1,1):
                radius,q,j=branches.geom.boundary(C,E,np.array([sign]));distance=float(np.linalg.norm(q[0]-reference['center']))
                if best is None or distance<best[0]:best=(distance,float(radius[0]),q[0],j[0],sign)
            distance,value,q,j,sign=best;J=np.zeros(200);J[(2*E[0,:,None]+[0,1]).ravel()]=j.ravel();difference=J-reference['radius_jacobian'];factor=float(difference@n/(n@n));rows=np.vstack((np.r_[reference['radius_jacobian'],-1],np.r_[J,-1],np.r_[n,0]))
            result.append(dict(pair=E[0].tolist(),pair_sign=sign,event_value=e['raw_value'],pair_point_minus_triangle_point=distance,radius_value_difference=value-reference['radius_squared'],gradient_dependency_factor=factor,gradient_dependency_error_norm=float(np.linalg.norm(difference-factor*n)),gradient_difference_norm=float(np.linalg.norm(difference)),augmented_three_row_singular_values=np.linalg.svd(rows,compute_uv=False).tolist(),expected_rank_at_exact_event_at_most=2))
    return result
