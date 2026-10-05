# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Explicit restricted-DT/circle target jets and finite quadratic diagnostics.

No covering optimizer or continuation is run. Curvature models are attached
only to the pinned near-stationary frame. Cap frames receive geometric jets.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[k]='1'
import sys,json,hashlib,shutil
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.linalg import null_space
from scipy.spatial import Delaunay
ROOT=Path(__file__).resolve().parent.parent;BASE=ROOT.parent
sys.dont_write_bytecode=True
sys.path.insert(0,str(BASE/'b32_core_exchange_geometry/frozen_support/replay_source'))
from graph_tools import extract
from exact_actual_graph import verify
from event_jets import jets,event,circum,convex_quad

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def orient(a,b,c):return float((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))

def qp(S,c,W,slacks):
    """Small convex jet QP via finite active-set linear algebra, no C objective."""
    invc=np.linalg.solve(S,c);den=float(c@invc)
    if den<1e-18:return None,dict(status='EVENT_GRADIENT_ZERO_IN_MODEL')
    active=[];seen=set();y=-invc/den
    for it in range(60):
        violation=W@y-slacks
        if violation.max(initial=-1)<=2e-10:
            return y,dict(status='PASS_FINITE_JET_QP',iterations=it,active_inequalities=active,maximum_linear_inequality_excess=float(violation.max(initial=0)),quadratic_squared_radius_cost=float(.5*y@S@y),event_linear_residual=float(c@y+1))
        add=int(np.argmax(violation))
        if add not in active:active.append(add)
        state=tuple(active)
        if state in seen:return None,dict(status='ACTIVE_SET_CYCLE',iterations=it)
        seen.add(state);B=np.vstack((c,W[active]));b=np.r_[-1,slacks[active]];K=np.linalg.solve(S,B.T);G=B@K
        z,res,rank,s=np.linalg.lstsq(G,b,rcond=1e-10)
        if rank<len(b) and np.linalg.norm(G@z-b)>1e-8:return None,dict(status='LINEAR_EVENT_INEQUALITY_INCOMPATIBLE',iterations=it)
        mu=-z[1:]
        if len(mu) and mu.min()<-1e-9:
            active.pop(int(np.argmin(mu)));continue
        y=K@z
    return None,dict(status='ACTIVE_SET_ITERATION_LIMIT')

def make_models(C):
    saved=BASE/'b32_weak_contact_continuation/matrices/incoming_analysis.npz';data=np.load(saved);assert np.array_equal(C,data['C'])
    J=data['J'];H=data['H'];strong=data['strong'];weak=data['weak'];nearJ=data['near_J'];nearSlack=data['near_slacks'];lam=data['lambda_values'];rot=np.column_stack((-C[:,1],C[:,0])).ravel();rot/=np.linalg.norm(rot)
    models={}
    for name,A in [('all_positive_146',J),('strong_132_heuristic',J[strong])]:
        aug=np.vstack((np.column_stack((A,-np.ones(len(A)))),np.r_[rot,0]));N=null_space(aug,rcond=1e-9);assert N.shape[1]==(54 if name.startswith('all') else 68)
        S=N[:200].T@H@N[:200];S=(S+S.T)/2
        W=np.column_stack((np.r_[J[weak],nearJ],-np.ones(22)))@N
        slacks=np.r_[np.zeros(len(weak)),nearSlack]
        models[name]=dict(N=N,S=S,W=W,slacks=slacks,minimum_curvature=float(np.linalg.eigvalsh(S)[0]),lambda_is_exact_zero=False)
    np.savez(ROOT/'reference_models.npz',C=C,J=J,H=H,lambda_values=lam,strong=strong,weak=weak,near_J=nearJ,near_slacks=nearSlack,N_all=models['all_positive_146']['N'],N_strong=models['strong_132_heuristic']['N'])
    return models

def model_cost(g,e0,models,R):
    result={};directions={}
    for name,m in models.items():
        c=g@m['N'][:200]/e0;y,stats=qp(m['S'],c,m['W'],m['slacks'])
        stats.update(tangent_dimension=m['N'].shape[1],minimum_model_curvature=m['minimum_curvature'],incoming_inequalities=8,tiny_positive_multiplier_inequalities=14,nearzero_lambda_is_exact_zero=False,
            scope='Numerical quadratic event jet diagnostic only; not a covering optimizer, exact critical cone or event-path certificate')
        if y is not None:
            d=m['N'][:200]@y;t=float(m['N'][-1]@y);stats.update(predicted_epigraph_squared_radius_linear_change=t,predicted_max_center_displacement=float(np.linalg.norm(d.reshape(100,2),axis=1).max()),predicted_radius_change_from_quadratic_cost=stats['quadratic_squared_radius_cost']/(2*R));directions[name]=d
        result[name]=stats
    return result,directions

def main():
    assert not (ROOT/'candidate_jets.npz').exists(),'Use a new stage instead of overwriting target jets'
    (ROOT/'inputs').mkdir(exist_ok=True)
    base=BASE/'b32_metric_seed_search/final_handoff/index14_grid0061.npy';assert sha(base)=='78d5a0ae08b07010da7c52e9634d1a14ecdcf12cf0e04ff5d1895895231452c0'
    parent=BASE/'b32_weak_contact_event_extension';end=json.loads((parent/'RESULT.json').read_text())['endpoints_at_or_below_cap']
    source=[dict(name='pinned14_best61',file=base,file_sha256=sha(base),radius=.11551473632645183,is_stationary_reference_only=True)]
    for r in end:source.append(dict(name=f'cap_curve{r["curve"]:02d}',file=parent/r['file'],file_sha256=r['file_sha256'],radius=r['actual_radius'],is_stationary_reference_only=False))
    records=[];allgrads=[];localh=[];allcols=[];direction_rows=[];source_reports=[];models=None
    for sid,s in enumerate(source):
        p=s['file'];assert sha(p)==s['file_sha256'];C=np.load(p);dst=ROOT/'inputs'/f'{s["name"]}.npy';shutil.copy2(p,dst);s['file']=str(dst.resolve());s['input_original_file']=str(p.resolve())
        a=extract(C,True);integer=verify(C);assert a['fingerprint']=='3040bfc74bafa7eb3652b512fdf1824069df40e28efb56cb17146a807212ede8';source_reports.append(dict(s,integer_actual_graph=integer))
        if sid==0:models=make_models(C)
        full=[];edges=defaultdict(list)
        for tri in Delaunay(C).simplices:
            f=list(map(int,tri))
            if orient(*C[f])<0:f[1],f[2]=f[2],f[1]
            O,r2=circum(C[f].tolist());O=np.asarray(O);j=len(full);full.append(dict(face=f,O=O,r2=r2,inside=float(O@O)<1))
            for k in range(3):edges[tuple(sorted((f[k],f[(k+1)%3])))].append((j,f[(k+2)%3]))
        candidates=[];convex_rejected=outside_rejected=0
        for edge,adj in edges.items():
            if len(adj)!=2:continue
            one,two=(full[x[0]] for x in adj);opp=[adj[0][1],adj[1][1]];ok,values=convex_quad(C,edge,opp)
            if not ok:convex_rejected+=1;continue
            if not (one['inside'] or two['inside']):outside_rejected+=1;continue
            primary=one if one['inside'] else two;other=adj[1][1] if primary is one else adj[0][1];labels=primary['face']+[other]
            O=primary['O'];r2=float(primary['r2']);remaining=[j for j in range(100) if j not in labels];power=np.sum((C[remaining]-O)**2,1)-r2
            tie=C.copy();d=C[other]-O;tie[other]=O+np.sqrt(r2/(d@d))*d;tieok,tievalues=convex_quad(tie,edge,opp)
            candidates.append(dict(kind='dt_flip',labels=labels,edge=list(edge),opposites=opp,adjacent_faces=[one['face'],two['face']],both_adjacent_circumcenters_inside=one['inside'] and two['inside'],core_annulus='core' if max(labels)<68 else 'annulus',source_convex_orientation_margins=values,source_primary_circumcenter=O.tolist(),source_primary_circumcenter_inside_margin=float(1-O@O),source_primary_circumradius_squared=r2,source_minimum_other_site_emptycircle_power=float(power.min()),ideal_radial_tie_convex=tieok,ideal_radial_tie_max_center_norm=float(np.linalg.norm(tie,axis=1).max()),ideal_radial_tie_event_residual=event(tie,'dt_flip',labels),ideal_radial_tie_is_exact_binary64_cocircular=False,
                geometry_eligible=bool(tieok and power.min()>1e-12 and np.linalg.norm(tie,axis=1).max()<1),geometry_eligibility_scope='Source adjacent exact DT faces; convex four-site flip, ideal nearest four-site circle strictly inside unit disk. Ideal radial tie is floating diagnostic, not an exact binary64 tie.'))
        for f in full:
            O=f['O'];r2=float(f['r2']);labels=f['face'];remaining=[j for j in range(100) if j not in labels];power=np.sum((C[remaining]-O)**2,1)-r2
            candidates.append(dict(kind='circumcenter_unit',labels=labels,core_annulus='core' if max(labels)<68 else 'annulus',source_primary_circumcenter=O.tolist(),source_primary_circumcenter_inside_margin=float(1-O@O),source_primary_circumradius_squared=r2,source_minimum_other_site_emptycircle_power=float(power.min()),geometry_eligible=bool(f['inside'] and power.min()>1e-12),geometry_eligibility_scope='Exact source DT face and nearest circle; prospective unit-circle equality and all competing sites remain requirements of the future path.'))
        counts=dict(full_DT_triangles=len(full),internal_full_DT_edges=sum(len(v)==2 for v in edges.values()),convex_quad_rejects=convex_rejected,irrelevant_disk_exterior_quad_rejects=outside_rejected,DT_targets=sum(r['kind']=='dt_flip' for r in candidates),circle_targets=sum(r['kind']=='circumcenter_unit' for r in candidates))
        source_reports[-1]['candidate_counts']=counts
        rot=np.column_stack((-C[:,1],C[:,0])).ravel();rot/=np.linalg.norm(rot)
        for r in candidates:
            cid=len(records);e0,g,H,cols=jets(C,r['kind'],r['labels']);assert abs(e0)>1e-12
            gp=g-rot*(rot@g);eu=-e0*gp/(gp@gp);directions={'rotation_quotient_euclidean':eu};cost={}
            if sid==0:cost,more=model_cost(g,e0,models,s['radius']);directions.update(more)
            for name,d in directions.items():direction_rows.append(dict(candidate_id=cid,model=name,direction=d))
            r.update(candidate_id=cid,source_id=sid,source_name=s['name'],source_file=str(dst.resolve()),source_file_sha256=s['file_sha256'],source_radius=s['radius'],source_is_stationary_reference=sid==0,raw_event_value=e0,normalized_source_event_value=1.,event_function_definition=('Squared distance of fourth site to primary triangle circumcenter minus primary circumradius squared' if r['kind']=='dt_flip' else '1 minus squared norm of triangle circumcenter'),normalized_event_definition='raw_event(C)/raw_event(sourceC); source constant scale, homotopy target1-f',gradient_units='Raw squared-length event per center coordinate; normalized derivatives divide by fixed raw_source_value',jet_index=cid,spatial_dimension=200,rotation_quotient_dimension=199,local_jet_columns=cols.tolist(),full_gradient_norm=float(np.linalg.norm(g)),rotation_derivative=float(g@rot),euclidean_predicted_max_center_displacement=float(np.linalg.norm(eu.reshape(100,2),axis=1).max()),cost_models=cost,
                endpoint_KKT_or_stress_Hessian_applied=False,rank_or_Hessian_stability_claim=False)
            if sid==0:
                valid=[(name,v['quadratic_squared_radius_cost']) for name,v in cost.items() if v['status']=='PASS_FINITE_JET_QP'];r['ranking_squared_radius_cost']=min((v for _,v in valid),default=None);r['ranking_cost_model']=min(valid,key=lambda x:x[1])[0] if valid else None
            records.append(r);allgrads.append(g);hh=np.zeros((8,8));hh[:len(cols),:len(cols)]=H;localh.append(hh);cc=np.full(8,-1);cc[:len(cols)]=cols;allcols.append(cc)
        print(s['name'],counts,'geometry_eligible',sum(r['geometry_eligible'] for r in records if r['source_id']==sid),flush=True)
    np.savez(ROOT/'candidate_jets.npz',raw_gradient=np.array(allgrads),local_hessian=np.array(localh),local_columns=np.array(allcols),raw_source_event=np.array([r['raw_event_value'] for r in records]))
    np.savez(ROOT/'predicted_directions.npz',candidate_id=np.array([r['candidate_id'] for r in direction_rows]),model=np.array([r['model'] for r in direction_rows]),spatial_direction=np.array([r['direction'] for r in direction_rows]))
    sources_clean=[{k:v for k,v in r.items() if k!='integer_actual_graph'} for r in source_reports]
    (ROOT/'SOURCE_PROVENANCE.json').write_text(json.dumps(dict(status='PASS_INPUTS_AND_FULL_INTEGER_GRAPH',sources=source_reports,reference_model_source_file=str((BASE/'b32_weak_contact_continuation/matrices/incoming_analysis.npz').resolve()),reference_model_source_sha256=sha(BASE/'b32_weak_contact_continuation/matrices/incoming_analysis.npz')),indent=2)+'\n')
    (ROOT/'candidate_targets.json').write_text(json.dumps(dict(status='JET_ENUMERATION_COMPLETE',records=records,sources=sources_clean,raw_jets_file='candidate_jets.npz',predicted_directions_file='predicted_directions.npz',nonlinear_cover_optimizer_calls=0,continuation_paths_run=0),indent=2)+'\n')
    budget=dict(status='ACTUAL_ENUMERATED',source_frames=3,event_jets=len(records),reference_models_only_on_pinned_frame=True,cap_KKT_models=0,full_integer_source_checks=3,finite_jet_QP_diagnostics=sum(len(r['cost_models']) for r in records),new_covering_optimizer_calls=0,repeated_old_continuation_frames=0,independent_new_whole_disk_calls=0,maximum_planned_targets=6,first_delivery_target_count=2,planned_new_whole_disk_probe_calls_first_delivery=6,BLAS_threads=1,max_compute_workers=2)
    (ROOT/'BUDGET.json').write_text(json.dumps(budget,indent=2)+'\n');print(json.dumps(budget,indent=2))

if __name__=='__main__':main()
