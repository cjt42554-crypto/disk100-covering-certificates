# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Actual incoming-cover contacts and rank sensitivity of approximate KKT data."""
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[n]='1'
import sys,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent;sys.dont_write_bytecode=True;sys.path.insert(0,str(ROOT))
from continuation_modes import inventory,extract

def main():
    D=json.loads((ROOT/'continuation_linear_algebra.json').read_text());sources=[];incoming=[]
    for r in D['records']:
        C=np.load(r['source_file']);a=extract(C,True);g,J,keys=inventory(C,np.array(a['faces']),np.array(a['boundary_pairs']),np.array(a['boundary_signs']));slack=g.max()-g
        active=np.flatnonzero(slack<5e-14);A=J[active];lam=np.array(r['NNLS_lambda']);rot=np.column_stack((-C[:,1],C[:,0])).ravel();rot/=np.linalg.norm(rot)
        spectra={};stress=A.T@lam
        for name,M in [('spatial_active_J',A),('spatial_active_J_plus_rotation_gauge',np.vstack((A,rot))),('heuristic_strong_J',A[lam>=1e-7])]:
            sv=np.linalg.svd(M,compute_uv=False)
            spectra[name]=dict(singular_values=sv.tolist(),rank_by_relative_tolerance={str(t):int(sum(sv>t*sv[0])) for t in (1e-12,1e-10,1e-9,1e-8,1e-7,1e-6)})
        sr=dict(index=r['index'],active_count=len(active),squared_tolerance=5e-14,lambda_l2_norm=float(np.linalg.norm(lam)),stationarity_residual_l2=float(np.linalg.norm(stress)),stationarity_residual_inf=float(abs(stress).max()),
            numerical_sigma_min_upper_bound_from_lambda=float(np.linalg.norm(stress)/np.linalg.norm(lam)),sigma_min_bound_is_exact_certificate=False,spectra=spectra,
            exact_stationary_active_J_rank_would_be_at_most=len(active)-1,rotation_gauged_spatial_tangent_dimension_at_exact_stationarity_would_be_at_least=200-len(active),
            warning='Finite nonzero KKT residual can lift one singular value and drop a cheap near-null mode. Positive reduced Hessian after a single numerical rank truncation does not certify the full exact critical cone.')
        sources.append(sr)
        if r['index']==14:
            near=np.flatnonzero((slack>=5e-14)&(slack<1e-8))
            for j in near:
                ir=dict(key=keys[j],squared_radius_slack=float(slack[j]),contact_index=int(j),source14_squared_radius=float(g.max()),mode_directional_predictions=[])
                for m in D['modes']:
                    if m['index']!=14:continue
                    p=Path(m['file']);v=np.load(p)['V'];slope=float(J[j]@v.ravel());alpha=float(slack[j]/slope) if slope>0 else None
                    ir['mode_directional_predictions'].append(dict(mode_file=str(p),directional_squared_radius_slope=slope,linearized_first_cover_contact_alpha=alpha,
                        alpha_is_exact_event_isolation=False))
                incoming.append(ir)
    out=dict(status='PASS_NUMERICAL_DIAGNOSTIC',source_rank_audits=sources,source14_near_inactive_incoming_contacts=incoming,
        source14_incoming_count=len(incoming),existing_mode_handoff_manifest_sha256='67b8e22ff27a45161ccf2d587038aa909a18f7e8a9eb7870fad27c47329feb8b',
        correction='The146-active/tiny14-only models do not control all finite-scale actual covering contacts. Include all8 near-inactive constraints as finite inequalities before treating a direction as a near-cover continuation. No exact zero multiplier or stationary-rank claim is made.',
        no_covering_or_proxy_optimizer=True)
    (ROOT/'incoming_contact_rank_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print('incoming14',len(incoming))
    for r in sources:print('index',r['index'],'sigma_min',r['spectra']['spatial_active_J']['singular_values'][-1],'lambda_bound',r['numerical_sigma_min_upper_bound_from_lambda'],'rank',r['spectra']['spatial_active_J_plus_rotation_gauge']['rank_by_relative_tolerance'])
    for r in incoming:print(r['key'],r['squared_radius_slack'],r['mode_directional_predictions'][0]['linearized_first_cover_contact_alpha'])

if __name__=='__main__':main()
