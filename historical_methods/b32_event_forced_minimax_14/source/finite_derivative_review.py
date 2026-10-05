# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Finite independent checks of the new proxy's t/event/gauge Jacobians.

No proxy/whole minimizer is run and no candidate is claimed. A source local
adjacent-edge circuit is used solely for differentiating the generic functions.
"""
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[n]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from collections import defaultdict
import numpy as np
from proxy_driver import ROOT,inventory,write,sha
import event_jets

def main():
    source=ROOT.parent/'b32_event_targeting_14_review/inputs/pinned14.npy';C=np.load(source);g,J,keys,top=inventory(C);T=top[0];edge=defaultdict(list)
    for f in T:
        for j in range(3):edge[tuple(sorted((int(f[j]),int(f[(j+1)%3]))))].append(int(f[(j+2)%3]))
    labels=None
    for (a,b),other in edge.items():
        if len(other)!=2:continue
        c,d=other
        if np.cross(C[b]-C[a],C[c]-C[a])<0:a,b=b,a
        trial=[a,b,c,d];jet=event_jets.event_jet(C,'dt_determinant',trial)
        if abs(jet['raw_value'])>1e-10 and jet['gradient_norm']>1e-10:labels=trial;break
    assert labels is not None
    jet=event_jets.event_jet(C,'dt_determinant',labels);scale=jet['gradient_norm'];e0=jet['raw_value'];w=np.column_stack((-C[:,1],C[:,0])).ravel();w/=np.linalg.norm(w)
    x=np.r_[C.ravel(),g.max()];A=np.column_stack((-J,np.ones(len(g))));eventJ=np.r_[jet['raw_jacobian']/scale,0];gaugeJ=np.r_[w,0];rows=[]
    for h in (1e-5,1e-6,1e-7):
        fdA=np.empty_like(A);fdE=np.empty(201);fdG=np.empty(201)
        for j in range(201):
            plus=x.copy();minus=x.copy();plus[j]+=h;minus[j]-=h;P=plus[:200].reshape(100,2);M=minus[:200].reshape(100,2)
            fdA[:,j]=((plus[-1]-inventory(P,top)[0])-(minus[-1]-inventory(M,top)[0]))/(2*h)
            fdE[j]=(event_jets.event_jet(P,'dt_determinant',labels)['raw_value']-event_jets.event_jet(M,'dt_determinant',labels)['raw_value'])/(2*h*scale)
            fdG[j]=(w@(plus[:200]-C.ravel())-w@(minus[:200]-C.ravel()))/(2*h)
        rows.append(dict(h=h,cover_epigraph_J_max_error=float(np.max(np.abs(fdA-A))),event_fixed_source_scaled_J_max_error=float(np.max(np.abs(fdE-eventJ))),rotation_J_max_error=float(np.max(np.abs(fdG-gaugeJ)))))
    assert max(r['cover_epigraph_J_max_error'] for r in rows)<1e-6
    assert max(r['event_fixed_source_scaled_J_max_error'] for r in rows)<1e-6
    write(ROOT/'FINITE_DERIVATIVE_REVIEW.json',dict(status='PASS_FINITE_SOURCE_DERIVATIVES',source_file=str(source),source_sha256=sha(source),diagnostic_local_circuit=labels,event_definition=jet['definition'],raw_e0=e0,fixed_source_gradient_scale=scale,objective='t=R^2',epigraph_constraint_definition='t-g_hole>=0',event_level_equation='[e(C)-(1-f)e0]/fixed_scale=0; tderivative0',rotation_slice_only=True,FD_rows=rows,proxy_solves=0,whole_optimizers=0,qualified_target_count=0,scope='Function/Jacobian calibration at one retained source; this local circuit is not a prepared pilot or an attained event. No earlier trajectory or optimizer was replayed.'))
    print('NEW_PROXY_FINITE_DERIVATIVES_PASS',rows,flush=True)
if __name__=='__main__':main()
