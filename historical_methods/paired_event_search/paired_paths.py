# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Numerical codimension-two incircle forcing, without released offsets.

Initial predictors use all 48 tangent dimensions and a complete 2x2 mixed
second-order equality response. Finite optimization uses all 199 fixed-slice
coordinates and the current actual hole inventory. Numeric event acceptance
requires both disk-interior nearest empty convex circles, not only two raw
polynomial zeros. No exact whole-path/global-optimum claims.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import argparse,json,time,hashlib,itertools,traceback
import numpy as np
from scipy.optimize import minimize,lsq_linear
from scipy.linalg import null_space
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent
sys.path.insert(0,str(BASE/'event_forced_paths'))
import forced_paths as solo
geometry=solo.geometry

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):solo.write_json(Path(p),v)


def prepare():
    m=np.load(ROOT/'inputs/matrices.npz');gm=np.load(ROOT/'inputs/full48_geometry.npz');j=np.load(ROOT/'inputs/full48_event_jets.npz')
    C=m['C'];Q=m['Q'];w=m['w'];t0=float(m['t']);Z=gm['Z'];G=gm['G'];A=m['A'];H=m['H'];B=gm['B_eq']
    assert np.linalg.matrix_rank(A)==151 and np.linalg.matrix_rank(np.r_[A,w[None,:]])==152 and np.linalg.matrix_rank(B)==153 and Z.shape==(200,48)
    assert np.linalg.eigvalsh(G).min()>0
    ids=list(map(str,j['event_ids']));ix={e:i for i,e in enumerate(ids)}
    records=json.loads((ROOT/'inputs/ranked_full48_events.json').read_text());ed={r['event']['id']:r['event'] for r in records}
    # Normalized event coordinates q=e/e0, with initial q=1 and forcing
    # q=1-p. This normalizes primitive determinants of different scale.
    normalized=j['tangent_gradients']/j['values'][:,None];Gi=np.linalg.solve(G,normalized.T);full_gram=normalized@Gi
    pairs=[]
    for e1,e2 in itertools.combinations(ids,2):
        if 'e0206' in (e1,e2):continue
        idx=[ix[e1],ix[e2]];gram=full_gram[np.ix_(idx,idx)];eig=np.linalg.eigvalsh(gram)
        if not np.all(np.isfinite(gram)) or eig[0]<=eig[-1]*1e-9:continue
        inverse=np.linalg.inv(gram);p=np.ones(2);cost=.5*p@inverse@p;rho=gram[0,1]/np.sqrt(gram[0,0]*gram[1,1]);shared=len(set(ed[e1]['sites'])&set(ed[e2]['sites']))
        pairs.append(dict(ids=[e1,e2],events=[ed[e1],ed[e2]],gram=gram,gram_eigenvalues=eig,gram_rank=2,gram_condition=float(eig[-1]/eig[0]),normalized_cross_coupling=float(rho),shared_sites=shared,union_sites=sorted(set(ed[e1]['sites'])|set(ed[e2]['sites'])),predicted_squared_barrier=float(cost),predicted_radius_barrier=float(cost/(2*np.sqrt(t0))),local_only=True))
    pairs.sort(key=lambda r:r['predicted_squared_barrier']);selected=[];chosen=set()
    def take(r,family):
        key=tuple(sorted(r['ids']))
        if key not in chosen:selected.append(dict(r,family=family));chosen.add(key)
    for ids2 in [('e0090','e0093'),('e0090','e0094'),('e0093','e0094')]:
        take(next(r for r in pairs if set(r['ids'])==set(ids2)),'weak_three_shared_sites_five_site_cavity')
    used=set()
    for r in pairs:
        if r['shared_sites']==3 and abs(r['normalized_cross_coupling'])>.25 and not set(r['ids'])&used and not set(r['ids'])&{'e0090','e0093','e0094'}:
            take(r,'other_five_site_cavity');used.update(r['ids'])
        if sum(r['family']=='other_five_site_cavity' for r in selected)>=3:break
    for r in pairs:
        if r['shared_sites']==2 and abs(r['normalized_cross_coupling'])>.25:take(r,'two_shared_sites_positive_or_negative_coupling')
        if sum(r['family']=='two_shared_sites_positive_or_negative_coupling' for r in selected)>=3:break
    for r in pairs:
        if r['shared_sites']<2 and abs(r['normalized_cross_coupling'])>.25:take(r,'different_local_cavities_strong_coupling')
        if len(selected)>=12:break
    predictors={}
    for i,r in enumerate(selected):
        idx=[ix[e] for e in r['ids']];aa=normalized[idx].T;inverse=np.linalg.inv(np.asarray(r['gram']));V=-Z@np.linalg.solve(G,aa)@inverse
        W2=np.zeros((2,2,201));rawW2=np.zeros_like(W2);hole_b=np.zeros((2,2,len(m['active'])));event_b=np.zeros((2,2,2));eps=1e-26
        derivatives=[]
        for a in range(2):
            _,aj=solo.evaluate(C.astype(complex)+1j*eps*V[:,a].reshape(100,2),(m['T'],m['E'],m['sg']));derivatives.append(np.imag(aj[m['active']])/eps)
        egrads=np.array([j['gradients'][k]/j['values'][k] for k in idx]);ehess=np.array([j['hessians'][k]/j['values'][k] for k in idx])
        for a in range(2):
            for b in range(a,2):
                hb=.5*(derivatives[a]@V[:,b]+derivatives[b]@V[:,a]);raw=np.linalg.lstsq(B,np.r_[-hb,0.],rcond=1e-11)[0]
                eb=egrads@raw[:200]+np.einsum('i,nij,j->n',V[:,a],ehess,V[:,b]);adjusted=raw.copy();adjusted[:200]+=V@eb
                hole_b[a,b]=hole_b[b,a]=hb;event_b[a,b]=event_b[b,a]=eb;rawW2[a,b]=rawW2[b,a]=raw;W2[a,b]=W2[b,a]=adjusted
        eqres=max(float(np.max(np.abs(B@W2[a,b]-np.r_[-hole_b[a,b],0.]))) for a in range(2) for b in range(2))
        eventres=max(float(np.max(np.abs(egrads@W2[a,b,:200]+np.einsum('i,nij,j->n',V[:,a],ehess,V[:,b])))) for a in range(2) for b in range(2))
        firstres=float(np.max(np.abs(egrads@V+np.eye(2))))
        mixed_energy=V.T@H@V
        r.update(pair_index=i,first_event_response_residual=firstres,second_contact_response_residual=eqres,second_event_response_residual=eventres,mixed_t_objective_curvature=mixed_energy,mixed_curvature_gram_inverse_error=float(np.max(np.abs(mixed_energy-inverse))),maximum_joint_linear_center_shift=float(np.linalg.norm(V.sum(axis=1).reshape(100,2),axis=1).max()),mixed_response_scope='Full geometric contact/event second derivatives with all cross terms; not a third-derivative stationary KKT response')
        if max(eqres,eventres,firstres)>3e-10:raise ValueError('Mixed pair predictor verification failed')
        predictors[i]=dict(V=V,W2=W2,rawW2=rawW2,hole_b=hole_b,event_b=event_b,gram=np.asarray(r['gram']),event_gradients=egrads,event_hessians=ehess)
        np.savez_compressed(ROOT/f'predictor_pair{i:02d}.npz',**predictors[i],C0=C,Z=Z,G=G,H=H,B_eq=B,ids=np.array(r['ids']))
    allfiles={str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'inputs').iterdir())}
    R,detail=solo.true_cover(C);checks=dict(input_sha256=allfiles,baseline_true_radius=R,baseline_score=detail,rank_A=151,rank_A_gauge=152,rank_B_eq=153,initial_tangent_dimension=48,finite_slice_dimension=199,G_min_eigenvalue=float(np.linalg.eigvalsh(G).min()),pair_gram_candidates_screened=len(pairs),geometric_optimization_pairs=len(selected),excluded_known_single_event='e0206 and the .11551610349101131 single-event basin',geometry_budgets=dict(max_true_radius=.125,max_center_shift=.3),one_CPU_BLAS_thread=True,selected=selected,scope='Linear/quadratic pair ranking is not optimization of those candidate configurations; at most twelve finite pair paths')
    write(ROOT/'input_checks.json',checks);write(ROOT/'selected_pairs.json',selected);write(ROOT/'pair_gram_ranking.json',pairs)
    np.savez_compressed(ROOT/'normalized_full48_gram.npz',ids=np.array(ids),normalized_tangent_gradients=normalized,gram=full_gram,G=G,Z=Z)
    return C,Q,w,t0,selected,predictors,checks


class PairProblem(solo.ForcedProblem):
    def __init__(self,C0,Q,w,t0,events,targets,seed,trust):
        super().__init__(C0,Q,w,t0,events[0],targets[0],seed,trust)
        self.events=events;self.targets=np.asarray(targets)
    def eq(self,x):
        C=self.unpack(x)[0]
        return np.array([solo.incircle(C,e)/e['base_value_corrected'] for e in self.events])-self.targets
    def eqjac(self,x):
        C=self.unpack(x)[0]
        return np.array([np.r_[solo.event_jet(C,e)[1]@self.Q/e['base_value_corrected'],0.] for e in self.events])


def correct(C0,Q,w,t0,events,targets,guess,previous,trust=.055,maxiter=280):
    C=guess.copy();rounds=[]
    for outer in range(9):
        problem=PairProblem(C0,Q,w,t0,events,targets,previous,trust);problem.top=geometry.topology(C)
        R,_=geometry.covering(C);x0=np.r_[Q.T@(C-C0).ravel(),min(.125**2/t0,R*R/t0+1e-10)]
        before=geometry.top_key(problem.top);tic=time.monotonic()
        res=minimize(lambda x:x[-1],x0,jac=lambda x:problem.objgrad,method='SLSQP',constraints=[dict(type='ineq',fun=problem.con,jac=problem.jac),dict(type='eq',fun=problem.eq,jac=problem.eqjac)],bounds=[(None,None)]*199+[(.008/t0,.125**2/t0)],options=dict(maxiter=maxiter,ftol=2e-13,disp=False))
        C,t=problem.unpack(res.x)
        if not np.isfinite(C).all():return C,False,dict(reason='nonfinite_coordinates',rounds=rounds)
        R,actual=geometry.covering(C);eq=np.abs(problem.eq(res.x));mincon=float(np.min(problem.con(res.x)));changed=before!=geometry.top_key(actual)
        missing=set(solo.keys((actual[0],actual[1],actual[2])))-set(solo.keys((problem.top[0],problem.top[1],problem.top[2])))
        rounds.append(dict(round=outer,success=bool(res.success),status=int(res.status),message=str(res.message),iterations=int(res.nit),objective_t=t,actual_radius=R,event_relative_residuals=eq,min_scaled_constraint=mincon,topology_changed=changed,missing_actual_keys=sorted(missing),seconds=time.monotonic()-tic))
        feasible=bool(eq.max()<3e-9 and mincon>-2e-8 and R*R<=t+3e-10 and R<=.125+1e-10 and np.linalg.norm(C,axis=1).max()<=1+1e-9 and np.linalg.norm(C-C0,axis=1).max()<=.3+1e-9)
        if feasible and not missing and (res.success or res.status==8):return C,True,dict(reason='feasible_corrected',rounds=rounds,final_t=t,optimizer_success=bool(res.success))
        if not changed and not missing:return C,False,dict(reason='SLSQP_or_feasibility_failure',rounds=rounds,final_t=t)
    return C,False,dict(reason='actual_topology_refresh_limit',rounds=rounds,final_t=t)


def pair_validation(C,events):
    vals=[solo.validate_event(C,e) for e in events];shared=set(events[0]['sites'])&set(events[1]['sites']);union=sorted(set(events[0]['sites'])|set(events[1]['sites']))
    out=dict(individual_events=vals,both_actual_inside_convex_empty=all(v['valid'] for v in vals),shared_sites=sorted(shared),union_sites=union)
    if len(shared)==3 and all('common_circle_center' in v for v in vals):
        q=np.asarray(vals[0]['common_circle_center']);q2=np.asarray(vals[1]['common_circle_center']);ds=np.sum((C-q)**2,axis=1)
        spread=float(np.ptp(ds[union]));center_diff=float(np.linalg.norm(q-q2));r2=float(ds[union[0]])
        valid=bool(out['both_actual_inside_convex_empty'] and center_diff<2e-8 and spread<3e-9 and ds.min()>=r2-2e-9 and q@q<=1+2e-9)
        out.update(five_site_cavity_cocircular_valid=valid,common_center_distance=center_diff,five_site_squared_distance_spread=spread,cavity_point_count=len(union))
        out['valid']=valid
    else:out['valid']=out['both_actual_inside_convex_empty']
    out['scope']='Numerical actual event/cavity acceptance; not an exact algebraic cocircularity certificate'
    return out


def stationarity(C,Q,events):
    top=geometry.topology(C);g,A,kk=solo.all_holes(C,top);ac=np.flatnonzero(g.max()-g<8e-10);Aa=A[ac]
    eg=np.array([solo.event_jet(C,e)[1]/e['base_value_corrected'] for e in events])
    M=np.zeros((200,len(ac)+2));M[:199,:len(ac)]=Q.T@Aa.T;M[:199,-2:]=Q.T@eg.T;M[-1,:len(ac)]=1
    rhs=np.r_[np.zeros(199),1.];ls=lsq_linear(M,rhs,bounds=(np.r_[np.zeros(len(ac)),-np.inf,-np.inf],np.full(len(ac)+2,np.inf)),tol=1e-13,lsmr_tol=1e-13,max_iter=500)
    mu=ls.x[:-2];nu=ls.x[-2:];B=np.r_[Aa@Q,eg@Q];N=null_space(B,rcond=1e-9);Beq=np.r_[np.c_[Aa@Q,-np.ones(len(ac))],np.c_[eg@Q,np.zeros(2)]];sv=np.linalg.svd(Beq,compute_uv=False)
    allmu=np.zeros(len(g));allmu[ac]=mu;H,asym=solo.hessian(C,top,allmu)
    for value,e in zip(nu,events):H+=value*solo.event_jet(C,e)[2]/e['base_value_corrected']
    ev=np.linalg.eigvalsh(N.T@Q.T@H@Q@N) if N.shape[1] else np.empty(0)
    return dict(active_count=len(ac),active_keys=[kk[i] for i in ac],zero_stress_keys=[kk[ac[i]] for i in np.flatnonzero(mu<1e-10)],multipliers=mu,event_multipliers=nu,stationarity_residual_inf=float(np.max(np.abs(M@ls.x-rhs))),epigraph_event_rows=len(Beq),epigraph_event_rank=int(np.linalg.matrix_rank(Beq,tol=1e-9)),epigraph_event_smallest_singular_value=float(sv[-1]),tangent_dimension=N.shape[1],tangent_H_eigenvalues=ev,strict_positive_equality_tangent_curvature=bool(ev.size and ev.min()>1e-9),H_asymmetry=asym,curvature_scope='Numerical active equality subspace diagnostic, no local-minimum/critical-cone certificate')


_KNOWN_SKELETON=None
def known_basin(C):
    global _KNOWN_SKELETON
    if _KNOWN_SKELETON is None:
        path=BASE/'event_forced_paths/path03_e0206/restored_event_after.npy';K=np.load(path)
        top=geometry.topology(K);g,_,kk=solo.all_holes(K,top);active=np.flatnonzero(g.max()-g<1e-9)
        used=sorted({int(z) for i in active for z in (kk[i][1:4] if kk[i][0]=='t' else kk[i][1:3])})
        _KNOWN_SKELETON=(K,used)
    K,used=_KNOWN_SKELETON
    # A disk that touches no binding holes may drift without changing the
    # known covering skeleton. Fitting all centers caused a false negative
    # for the existing 99-binding-center basin with unused center 65.
    angle=np.arctan2(np.sum(C[used,0]*K[used,1]-C[used,1]*K[used,0]),np.sum(C[used]*K[used]));rot=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    shifted=K@rot;distances=np.linalg.norm(C-shifted,axis=1);distance=float(distances.max());skeleton_distance=float(distances[used].max());R=geometry.covering(C)[0]
    return dict(match=bool(skeleton_distance<2e-4 and abs(R-.11551610349101131)<1e-8),coordinate_distance_after_rotation=distance,binding_skeleton_distance_after_rotation=skeleton_distance,binding_center_count=len(used),unbound_center_indices=sorted(set(range(100))-set(used)),radius=R,reference_radius=.11551610349101131,scope='Numerical known binding-skeleton/plateau detection, not identity of all center coordinates')


def restore(C,events,side,folder,label,C0):
    known=known_basin(C)
    if known['match']:
        name=f'restored_{label}_known_basin_skipped.npy';np.save(folder/name,C);R,detail=solo.true_cover(C)
        return dict(skipped=True,reason='known_single_event_basin_not_optimized',known_basin=known,coordinate_file=name,coordinate_sha256=solo.center_sha(C),true_radius=R,true_cover=detail,rounds=[],fully_unconstrained_stage_valid=False)
    best=None;bestR=np.inf;rounds=[];C=C.copy()
    for outer in range(7):
        top=geometry.topology(C);model=geometry.Model(C,top);x0=np.r_[C.ravel(),geometry.covering(C)[0]+1e-10];margin=2e-7 if outer==0 else 0.
        def sidecon(x):return np.array([side*solo.incircle(x[:-1].reshape(100,2),e)/e['base_value_corrected']-margin for e in events])
        def sidejac(x):return np.array([np.r_[side*solo.event_jet(x[:-1].reshape(100,2),e)[1]/e['base_value_corrected'],0.] for e in events])
        constraints=[dict(type='ineq',fun=model.con,jac=model.jac)]
        if outer<2:constraints.append(dict(type='ineq',fun=sidecon,jac=sidejac))
        tic=time.monotonic();res=minimize(lambda x:x[-1],x0,jac=lambda x:np.r_[np.zeros(200),1.],method='SLSQP',constraints=constraints,bounds=[(max(-1,float(z)-.085),min(1,float(z)+.085)) for z in C.ravel()]+[(.09,.125)],options=dict(maxiter=450,ftol=2e-13,disp=False))
        C=res.x[:-1].reshape(100,2);R,nt=geometry.covering(C);changed=geometry.top_key(top)!=geometry.top_key(nt);valid=bool(np.isfinite(C).all() and np.linalg.norm(C,axis=1).max()<=1+1e-9 and np.linalg.norm(C-C0,axis=1).max()<=.3 and (outer>=2 or sidecon(res.x).min()>-2e-9))
        stage='retain_both_diagonals' if outer==0 else ('zero_side_margins' if outer==1 else 'all_actual_holes_no_event_constraints');name=f'restoration_{label}_{outer:02d}_{stage}.npy';np.save(folder/name,C);Rtrue,details=solo.true_cover(C);kb=known_basin(C)
        rounds.append(dict(round=outer,stage=stage,coordinate_file=name,coordinate_sha256=solo.center_sha(C),file_sha256=sha(folder/name),true_radius=Rtrue,true_cover=details,event_signed_margins=sidecon(res.x),event_constraints_present=outer<2,valid=valid,success=bool(res.success),status=int(res.status),message=str(res.message),iterations=int(res.nit),topology_changed=changed,known_single_basin=kb,seconds=time.monotonic()-tic))
        if outer>=2 and valid and Rtrue<bestR:best=C.copy();bestR=Rtrue
        if kb['match']:
            rounds[-1]['remaining_optimization_skipped_known_single_basin']=True;break
        if outer>=2 and res.success and not changed:break
    if best is None:best=C.copy()
    name=f'restored_{label}.npy';np.save(folder/name,best);R,detail=solo.true_cover(best)
    return dict(skipped=False,coordinate_file=name,coordinate_sha256=solo.center_sha(best),file_sha256=sha(folder/name),true_radius=R,true_cover=detail,rounds=rounds,fully_unconstrained_stage_valid=any(r['valid'] and not r['event_constraints_present'] for r in rounds),final_event_constraints_present=False,offsets=0,exact_binary64_signs=[solo.exact_predicates(best,e)['incircle_sign'] for e in events],known_single_basin=known_basin(best))


def run_pair(r,pred,C0,Q,w,t0):
    index=r['pair_index'];folder=ROOT/f'pair{index:02d}_{r["ids"][0]}_{r["ids"][1]}';folder.mkdir(exist_ok=True)
    events=r['events'];rows=[];accepted=[];failure=None;event_validation=None;restores=[];all_events=json.loads((ROOT/'inputs/events.json').read_text())
    def save(C,fraction,phase,log=None):
        R,detail=solo.true_cover(C);top=geometry.topology(C);g,_,kk=solo.all_holes(C,top);values=np.array([solo.incircle(C,e)/e['base_value_corrected'] for e in events]);margins,inventory=solo.event_inventory(C,all_events);name=f'{len(rows):03d}_{phase}.npz'
        np.savez_compressed(folder/name,C=C,fraction=fraction,T=top[0],E=top[1],sg=top[2],F=top[3],actual_hole_values=g,normalized_event_values=values,all_baseline_event_margins=margins)
        previous=set(map(tuple,rows[-1]['actual_hole_keys'])) if rows else set(kk)
        row=dict(phase=phase,fraction=fraction,snapshot=name,file_sha256=sha(folder/name),coordinate_sha256=solo.center_sha(C),true_radius=R,true_cover=detail,event_values_normalized=values,event_target=1-fraction,exact_binary64_predicates=[solo.exact_predicates(C,e) for e in events],fixed_gauge_residual=float(w@(C-C0).ravel()),max_center_shift=float(np.linalg.norm(C-C0,axis=1).max()),max_center_norm=float(np.linalg.norm(C,axis=1).max()),actual_hole_keys=kk,actual_inventory_added=sorted(set(kk)-previous),actual_inventory_removed=sorted(previous-set(kk)),complete_event_inventory=inventory,corrector=log,known_single_basin=known_basin(C),numerical_only=True)
        rows.append(row);write(folder/'steps.json',rows);return row
    schedule=[.0001,.0003,.001,.003,.01,.025,.055,.10,.17,.26,.37,.50,.63,.75,.85,.93,.975,.995,1.,1.001,1.006,1.02]
    if (folder/'steps.json').exists():
        rows=json.loads((folder/'steps.json').read_text())
        for row in rows:
            if row['phase']=='failed_attempt':continue
            C=np.load(folder/row['snapshot'])['C']
            if sha(folder/row['snapshot'])!=row['file_sha256'] or solo.center_sha(C)!=row['coordinate_sha256']:raise ValueError('Resume snapshot hash mismatch')
            accepted.append((row['fraction'],C.copy()))
        last,C=accepted[-1];pending=[x for x in schedule if x>last+1e-12];attempts=len(rows)-1
    else:
        C=C0.copy();save(C,0.,'start');accepted.append((0.,C.copy()));last=0.;pending=list(schedule);attempts=0
    while pending and attempts<75:
        fraction=pending.pop(0);df=fraction-last
        if last==0:
            pp=np.full(2,fraction);guess=C0+(pred['V']@pp+.5*np.einsum('i,j,ijk->k',pp,pp,pred['W2'][:,:,:200])).reshape(100,2)
        elif len(accepted)>1:
            prev,Cprev=accepted[-2];guess=C+df*(C-Cprev)/(last-prev)
        else:guess=C+(df*pred['V'].sum(axis=1)).reshape(100,2)
        shift=np.linalg.norm(guess-C,axis=1).max()
        if shift>.043:guess=C+(guess-C)*(.043/shift)
        np.savez_compressed(folder/f'attempt{attempts:03d}_seed.npz',C=guess,previous=C,fraction=fraction,previous_fraction=last)
        trial,ok,log=correct(C0,Q,w,t0,events,np.full(2,1-fraction),guess,C);attempts+=1
        if not ok:
            save(trial,fraction,'failed_attempt',log)
            if df>.0008 and attempts<75:
                pending.insert(0,fraction);pending.insert(0,(last+fraction)/2);continue
            failure=dict(reason=log['reason'],fraction=fraction,previous_fraction=last,attempts=attempts);break
        row=save(trial,fraction,'double_event' if fraction==1 else ('after_both_events' if fraction>1 else 'corrected'),log);row['stationary_diagnostic']=stationarity(trial,Q,events);write(folder/'steps.json',rows)
        C=trial;last=fraction;accepted.append((fraction,C.copy()))
        write(ROOT/'checkpoint.json',dict(status='running',current_pair=index,current_ids=r['ids'],current_fraction=last,last_snapshot=str(folder/row['snapshot']),completed_pairs=index,driver_sha256=sha(Path(__file__)),updated_unix=time.time()))
        if fraction==1:
            event_validation=pair_validation(C,events);row['pair_event_validation']=event_validation;write(folder/'steps.json',rows)
            if not event_validation['valid']:
                failure=dict(reason='simultaneous_primitive_roots_are_not_both_actual_inside_events',fraction=fraction,validation=event_validation);break
        if row['known_single_basin']['match']:
            failure=dict(reason='known_single_event_basin_not_reoptimized',fraction=fraction);break
    before=[row for row in rows if row['phase']=='corrected' and row['fraction']<1];after=[row for row in rows if row['phase']=='after_both_events']
    if event_validation and event_validation['valid'] and before and after:
        br=before[-1];ar=after[-1];event_validation.update(exact_binary64_sign_brackets=[[br['exact_binary64_predicates'][i]['incircle_sign'],ar['exact_binary64_predicates'][i]['incircle_sign']] for i in range(2)],before_snapshot=br['snapshot'],after_snapshot=ar['snapshot'])
        sources=[(br,'before_both_events',1),(ar,'after_both_events',-1)]
    else:
        row=[row for row in rows if row['phase']!='failed_attempt'][-1];sources=[(row,'last_accepted',1 if row['fraction']<=1 else -1)]
    for row,label,side in sources:restores.append(restore(np.load(folder/row['snapshot'])['C'],events,side,folder,label,C0))
    acceptedrows=[row for row in rows if row['phase']!='failed_attempt'];result=dict(pair=r,steps=rows,stop=failure,event_validation=event_validation,restorations=restores,accepted_nodes=len(acceptedrows),failed_attempt_nodes=len(rows)-len(acceptedrows),corrector_attempts=attempts,last_accepted_fraction=last,sampled_peak_true_radius=max(row['true_radius'] for row in acceptedrows),best_fully_restored_true_radius=min([rr['true_radius'] for rr in restores if rr.get('fully_unconstrained_stage_valid')],default=None),peak_scope='Maximum over saved accepted nodes, not an interval supremum',scope='Numerical two-event/cavity construction only; no exact path or global optimum proof',offsets=0,one_CPU_BLAS_thread=True)
    write(folder/'result.json',result);print(json.dumps(solo.clean(dict(pair=index,ids=r['ids'],accepted=result['accepted_nodes'],failed=result['failed_attempt_nodes'],last_fraction=last,event_valid=event_validation and event_validation['valid'],peak=result['sampled_peak_true_radius'],restored=[rr['true_radius'] for rr in restores],stop=failure))),flush=True)
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--pairs',type=int,default=12);ap.add_argument('--prepare-only',action='store_true');args=ap.parse_args()
    if not 0<=args.pairs<=12:raise ValueError('At most twelve pair paths authorized')
    tic=time.monotonic();C0,Q,w,t0,selected,predictors,checks=prepare()
    if args.prepare_only:return
    results=[]
    for r in selected[:args.pairs]:
        folder=ROOT/f'pair{r["pair_index"]:02d}_{r["ids"][0]}_{r["ids"][1]}';resultfile=folder/'result.json'
        if resultfile.exists():results.append(json.loads(resultfile.read_text()));continue
        try:result=run_pair(r,predictors[r['pair_index']],C0,Q,w,t0)
        except Exception:
            folder.mkdir(exist_ok=True);result=dict(pair=r,exception=traceback.format_exc());write(folder/'failure.json',result);print(json.dumps(solo.clean(result)),flush=True)
        results.append(result);fully=[rr['best_fully_restored_true_radius'] for rr in results if rr.get('best_fully_restored_true_radius') is not None]
        write(ROOT/'summary.json',dict(input_checks=checks,pairs=results,completed_pairs=len(results),best_fully_restored_true_radius=min(fully,default=None),seconds_this_run=time.monotonic()-tic,scope='Numerical full48 coupled event predictions and full199 forced paths; no global optimum claim'))
        write(ROOT/'checkpoint.json',dict(status='completed' if len(results)==args.pairs else 'running',completed_pairs=len(results),requested_pairs=args.pairs,current_pair=None,driver_sha256=sha(Path(__file__)),updated_unix=time.time()))
    print('DONE',len(results),'seconds',time.monotonic()-tic,flush=True)


if __name__=='__main__':main()
