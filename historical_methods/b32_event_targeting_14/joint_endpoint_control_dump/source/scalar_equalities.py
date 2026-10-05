# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Two common-primary event equalities, all200 spatial columns plus zero t.

All chart scales are fixed source constants. This adapter computes derivatives
only; no covering objective, nonlinear solve, geometry qualification or gauge.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[k]='1'
import numpy as np
from event_jets import jets

def target_raw_levels(target,progress=None):
    if progress is None:return np.asarray(target['final_raw_power_goals'],float)
    if progress['phase']=='prejoint':return (1-float(progress['f']))*np.asarray(target['source_driver_raw_e0'],float)
    if progress['phase']=='postjoint':return float(progress['epsilon'])*np.asarray(target['post_raw_eta'],float)
    raise ValueError('phase must be prejoint or postjoint')

def values_and_jacobian(x,target,chart='source_normalized',progress=None,with_hessians=False):
    x=np.asarray(x,float);assert x.shape==(201,) and np.isfinite(x).all();C=x[:200].reshape(100,2)
    labels=target['common_primary_triangle'];P=C[labels];a=P[1]-P[0];b=P[2]-P[0]
    if abs(a[0]*b[1]-a[1]*b[0])<1e-20:raise ValueError('common_primary_triangle_degenerate')
    scales={'raw':np.ones(2),'source_normalized':np.asarray(target['source_driver_raw_e0'],float),'gradient_scaled':np.asarray(target['source_driver_gradient_norm'],float)}
    if chart not in scales:raise ValueError('unknown chart')
    den=scales[chart];assert np.all(den>0);goals=target_raw_levels(target,progress);values=[];J=[];Hs=[]
    for i,d in enumerate(target['added_sites_in_event_order']):
        e,g,H,cols=jets(C,'dt_flip',labels+[d]);row=np.zeros(201);row[:200]=g/den[i];values.append((e-goals[i])/den[i]);J.append(row)
        if with_hessians:
            HH=np.zeros((201,201));HH[np.ix_(cols,cols)]=H/den[i];Hs.append(HH)
    result=(np.array(values),np.array(J))
    return result+(np.array(Hs),) if with_hessians else result
