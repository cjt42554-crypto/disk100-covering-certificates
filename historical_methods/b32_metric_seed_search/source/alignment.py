# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
import engine as run
import numpy as np
from scipy.optimize import linear_sum_assignment
def compare(A,B,seeds=144):
    best=None
    for reflection in (False,True):
        V=B.copy()
        if reflection:V[:,1]*=-1
        for initial in np.linspace(0,2*np.pi,seeds,endpoint=False):
            angle=float(initial)
            for it in range(10):
                rot=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]);W=V@rot;_,j=linear_sum_assignment(np.sum((A[:,None,:]-W[None,:,:])**2,axis=2));Z=V[j];angle=float(np.arctan2(np.sum(A[:,0]*Z[:,1]-A[:,1]*Z[:,0]),np.sum(A*Z)))
            rot=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]);W=V[j]@rot;dist=np.linalg.norm(A-W,axis=1);error=float(np.sum(dist*dist))
            if best is None or error<best['sum_squared_error']:
                da=run.diagnose(A);db=run.diagnose(W);best=dict(sum_squared_error=error,maximum_center_distance=float(dist.max()),rotation_angle=angle,reflection=reflection,permutation=j,active_inventory_matches=da['active_signature']==db['active_signature'],inactive_inventory_matches=da['inactive_signature']==db['inactive_signature'])
    best.update(matches=bool(best['maximum_center_distance']<2e-6 and best['active_inventory_matches'] and best['inactive_inventory_matches']),rotation_seeds=seeds,scope='All100 coordinate alignment plus fresh active and inactive inventories. A failed finite search is not proof of inequivalence. Radius does not define a basin.')
    return best
