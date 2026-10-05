# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Freeze the completed grid and compare proposals by actual geometry."""
import grid_search as run
from pathlib import Path
import json,collections,hashlib,shutil
import numpy as np
from scipy.optimize import linear_sum_assignment,nnls
ROOT=run.ROOT
def compare(A,B):
    best=None
    for reflection in (False,True):
        V=B.copy()
        if reflection:V[:,1]*=-1
        for initial in np.linspace(0,2*np.pi,24,endpoint=False):
            angle=float(initial);j=np.arange(len(A))
            for it in range(8):
                rot=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]);W=V@rot;_,j=linear_sum_assignment(np.sum((A[:,None]-W[None,:])**2,axis=2));Z=V[j];angle=float(np.arctan2(np.sum(A[:,0]*Z[:,1]-A[:,1]*Z[:,0]),np.sum(A*Z)))
            rot=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]);W=V[j]@rot;dist=np.linalg.norm(A-W,axis=1);error=float(np.sum(dist*dist))
            if best is None or error<best['sum_squared_error']:
                da=run.diagnose(A);db=run.diagnose(W);best=dict(sum_squared_error=error,maximum_center_distance=float(dist.max()),rotation_angle=angle,reflection=reflection,permutation=j,active_inventory_matches=da['active_signature']==db['active_signature'],inactive_inventory_matches=da['inactive_signature']==db['inactive_signature'])
    best['matches']=bool(best['maximum_center_distance']<2e-6 and best['active_inventory_matches'] and best['inactive_inventory_matches']);best['scope']='Successful all100 coordinate alignment plus fresh active and inactive inventory identity supports a numerical group; failed finite alignment is not a proof of inequivalence. Radius is not a grouping criterion.'
    return best
def main():
    cp=json.loads((ROOT/'checkpoint.json').read_text());rows=[r for r in cp['records'] if 'final_radius' in r];proposals=sorted([r for r in rows if r['final_radius']<run.THRESHOLD and run.REFERENCE-r['final_radius']>run.MATERIAL],key=lambda r:r['grid_id']);groups=[];comparisons=[]
    for row in proposals:
        C=np.load(ROOT/row['final_coordinate_file']);matched=False
        for g in groups:
            rep=g['representative'];cmp=compare(np.load(ROOT/rep['final_coordinate_file']),C);comparisons.append(dict(a=rep['grid_id'],b=row['grid_id'],comparison=cmp))
            if cmp['matches']:g['members'].append(dict(grid_id=row['grid_id'],radius=row['final_radius'],file=row['final_coordinate_file'],sha256=row['final_file_sha256']));matched=True;break
        if not matched:groups.append(dict(representative=row,members=[dict(grid_id=row['grid_id'],radius=row['final_radius'],file=row['final_coordinate_file'],sha256=row['final_file_sha256'])]))
    run.write(ROOT/'proposal_alignment.json',dict(proposal_count=len(proposals),numerical_geometry_groups=len(groups),groups=groups,comparisons=comparisons))
    best=min(rows,key=lambda r:r['final_radius']);C=np.load(ROOT/best['final_coordinate_file']);g,p,J,kk,sites,top=run.inventory(C);ac=np.flatnonzero(g.max()-g<1e-9);M=np.r_[J[ac].T,np.ones((1,len(ac)))];rhs=np.r_[np.zeros(C.size),1.];lam,res=nnls(M,rhs,maxiter=10000)
    diagnostic=dict(branch=f'next146 grid{best["grid_id"]:04d} fresh actual winner; baseline grid25/de82 statistics are separate',input_file=best['final_coordinate_file'],input_file_sha256=run.sha(ROOT/best['final_coordinate_file']),input_coordinate_sha256=run.csha(C),radius=run.true_cover(C)[0],fresh_diagnostic=run.diagnose(C),active_squared_gap_spread=float(np.ptp(g[ac])),minimum_inactive_squared_gap=float(g.max()-max(g[[i for i in range(len(g)) if i not in set(ac)]])),NNLS_stationarity_residual_inf=float(np.max(np.abs(M@lam-rhs))),NNLS_lambda_min=float(lam.min()),NNLS_lambda_max=float(lam.max()),NNLS_lambda_sum=float(lam.sum()),active_keys=[kk[i] for i in ac],NNLS_lambda=lam,scope='Fresh floating diagnostic only. It is not a stationarity, local-minimum, root-isolation, exact full-cover or global-optimality certificate; baseline grid25/de82 lambda/KKT statistics are not transferred to this winner.')
    run.write(ROOT/'fresh_winner_diagnostic.json',diagnostic)
    result=json.loads((ROOT/'RESULT.json').read_text());result['proposal_numerical_geometry_groups']=len(groups);result['proposal_grouping_scope']='All100 alignment plus complete active/inactive inventory, not radius-only';result['fresh_winner_diagnostic']='fresh_winner_diagnostic.json';result['input_baseline_file_sha256']=json.loads((ROOT/'input_checks.json').read_text())['baseline_file_sha256'];result['input_baseline_scope']='grid25/de82 .11544695488558258 baseline; all preparation lambda and response stats belong only to this input';run.write(ROOT/'RESULT.json',result)
    old=ROOT.parent/'rattler_contact_transfer';oldmanifest=json.loads((old/'MANIFEST_SHA256.json').read_text());oldpass=all(run.sha(old/p)==h for p,h in oldmanifest.items());assert oldpass
    run.write(ROOT/'freeze_checks.json',dict(original429_unchanged=oldpass,original1467_unchanged=all(run.sha(ROOT.parent/'best147_boundary_transfer'/p)==h for p,h in json.loads((ROOT.parent/'best147_boundary_transfer/MANIFEST_SHA256.json').read_text()).items()),original_manifest_file_sha256=run.sha(old/'MANIFEST_SHA256.json'),new_baseline_copy_matches_original_grid25=run.sha(ROOT/'inputs/baseline.npy')==run.sha(ROOT.parent/'best147_boundary_transfer/best.npy'),best_proposal_file_sha256=best['final_file_sha256'],best_grid=best['grid_id']))
    manifest={str(p.relative_to(ROOT)):run.sha(p) for p in sorted(ROOT.rglob('*')) if p.is_file() and p.name!='MANIFEST_SHA256.json'};run.write(ROOT/'MANIFEST_SHA256.json',manifest)
    print('FREEZE_DONE',json.dumps(dict(files=len(manifest),groups=len(groups),best_grid=best['grid_id'],best_radius=best['final_radius'],original429_unchanged=oldpass,original1467_unchanged=all(run.sha(ROOT.parent/'best147_boundary_transfer'/p)==h for p,h in json.loads((ROOT.parent/'best147_boundary_transfer/MANIFEST_SHA256.json').read_text()).items()))),flush=True)
if __name__=='__main__':main()
