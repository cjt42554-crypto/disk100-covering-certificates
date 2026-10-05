# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""203-variable physical five-site circle homotopy, numerical research only."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,hashlib,time,itertools,argparse,traceback
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import ConvexHull,Delaunay
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent
sys.path.insert(0,str(BASE/'event_forced_paths'))
import forced_paths as f
geometry=f.geometry

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):f.write_json(Path(p),v)

class Jet:
    def __init__(self,v,g=None,H=None):
        self.v=v;self.g=np.zeros(10) if g is None else np.asarray(g);self.H=np.zeros((10,10)) if H is None else np.asarray(H)
    def __add__(self,b):
        if not isinstance(b,Jet):b=Jet(b)
        return Jet(self.v+b.v,self.g+b.g,self.H+b.H)
    __radd__=__add__
    def __neg__(self):return Jet(-self.v,-self.g,-self.H)
    def __sub__(self,b):return self+-b
    def __rsub__(self,b):return -self+b
    def __mul__(self,b):
        if not isinstance(b,Jet):b=Jet(b)
        return Jet(self.v*b.v,self.g*b.v+b.g*self.v,self.H*b.v+b.H*self.v+np.outer(self.g,b.g)+np.outer(b.g,self.g))
    __rmul__=__mul__
    def reciprocal(self):
        if abs(self.v)<1e-13:raise ValueError('degenerate common triangle')
        return Jet(1/self.v,-self.g/self.v**2,-self.H/self.v**2+2*np.outer(self.g,self.g)/self.v**3)
    def __truediv__(self,b):return self*(b.reciprocal() if isinstance(b,Jet) else 1/b)

def area(C,t):
    a,b,c=C[np.array(t)];u=b-a;v=c-a
    return u[0]*v[1]-u[1]*v[0]

def area_gradient(C,t):
    a,b,c=C[np.array(t)];out=np.zeros(200)
    out[np.array([2*i+k for i in t for k in (0,1)])]=[b[1]-c[1],c[0]-b[0],c[1]-a[1],a[0]-c[0],a[1]-b[1],b[0]-a[0]]
    return out

def physical_jets(C,case):
    sites=case['common']+case['extras'];eye=np.eye(10);pts=[[Jet(C[i,k],eye[2*j+k]) for k in (0,1)] for j,i in enumerate(sites)]
    a,b,c=pts[:3];u=[b[k]-a[k] for k in (0,1)];v=[c[k]-a[k] for k in (0,1)];d=u[0]*v[1]-u[1]*v[0]
    uu=u[0]*u[0]+u[1]*u[1];vv=v[0]*v[0]+v[1]*v[1]
    q=[a[0]+(uu*v[1]-vv*u[1])/(2*d),a[1]+(u[0]*vv-v[0]*uu)/(2*d)]
    radius=(a[0]-q[0])*(a[0]-q[0])+(a[1]-q[1])*(a[1]-q[1])
    cols=np.array([2*i+k for i in sites for k in (0,1)]);values=[];grads=[];hess=[]
    for point in pts[3:]:
        gap=(point[0]-q[0])*(point[0]-q[0])+(point[1]-q[1])*(point[1]-q[1])-radius
        g=np.zeros(200);g[cols]=gap.g;H=np.zeros((200,200));H[np.ix_(cols,cols)]=gap.H
        values.append(gap.v);grads.append(g);hess.append(H)
    return np.array(values),np.array(grads),np.array(hess),np.array([j.v for j in q]),float(radius.v)

def cycle_spec(C,cycle):
    ts=[];initial=[]
    for k in range(5):
        a,b=cycle[k],cycle[(k+1)%5]
        for c in cycle:
            if c in (a,b):continue
            ts.append([a,b,c]);initial.append(float(area(C,[a,b,c])))
    if min(initial)<=0:raise ValueError('cyclic hull orientation not positive')
    return ts,np.array(initial)

def prepare():
    m=np.load(ROOT/'inputs/matrices.npz');gm=np.load(ROOT/'inputs/full48_geometry.npz');C=m['C'];Q=m['Q'];w=m['w'];t0=float(m['t']);Z=gm['Z'];G=gm['G'];H=m['H'];B=gm['B_eq']
    assert Z.shape==(200,48) and np.linalg.matrix_rank(B)==153 and np.linalg.eigvalsh(G).min()>0
    es={r['event']['id']:r['event'] for r in json.loads((ROOT/'inputs/ranked_full48_events.json').read_text())};previous=json.loads((BASE/'paired_event_search/selected_pairs.json').read_text());excluded={frozenset(r['ids']) for r in previous};weak=frozenset(['e0093','e0094']);candidates=[];rejections=[]
    for e1,e2 in itertools.combinations(es,2):
        common=sorted(set(es[e1]['sites'])&set(es[e2]['sites']));union=sorted(set(es[e1]['sites'])|set(es[e2]['sites']))
        if len(common)!=3:continue
        ids=[e1,e2]
        if 65 in union:rejections.append(dict(ids=ids,reason='separate_65_agent_scope'));continue
        if frozenset(ids) in excluded and frozenset(ids)!=weak:rejections.append(dict(ids=ids,reason='old_primitive_pair_not_repeated'));continue
        hull=ConvexHull(C[union]).vertices
        if len(hull)!=5:rejections.append(dict(ids=ids,reason='initial_hull_not_five'));continue
        if area(C,common)<0:common[1],common[2]=common[2],common[1]
        extras=sorted(set(union)-set(common));cycle=[union[int(k)] for k in hull];tri,initial=cycle_spec(C,cycle)
        case=dict(ids=ids,events=[es[e1],es[e2]],common=common,extras=extras,union=union,cycle=cycle,cycle_triangles=tri,cycle_initial_areas=initial,common_initial_area=float(area(C,common)))
        gaps,grad,eh,q,rho=physical_jets(C,case);ds=np.sum((C-q)**2,axis=1)-rho;outside=[i for i in range(100) if i not in common]
        if q@q>=.999**2 or gaps.min()<=1e-10 or ds[outside].min()<-1e-10:rejections.append(dict(ids=ids,reason='initial_not_positive_inside_empty_circle'));continue
        validators=[f.validate_event(C,e) for e in case['events']]
        if not all(v.get('inside_unit_disk') and v.get('convex') and v.get('nearest_squared_gap',-1)>=-1e-10 for v in validators):rejections.append(dict(ids=ids,reason='initial_individual_event_not_actual_inside'));continue
        normalized=grad/gaps[:,None];a=Z.T@normalized.T;gi=np.linalg.solve(G,a);gram=a.T@gi;eig=np.linalg.eigvalsh(gram)
        if eig[0]<=eig[-1]*1e-9:rejections.append(dict(ids=ids,reason='physical_gap_gram_rank'));continue
        inverse=np.linalg.inv(gram);V=-Z@gi@inverse;cost=float(.5*np.ones(2)@inverse@np.ones(2))
        case.update(gaps=gaps,q0=q,rho0_squared=rho,physical_gap_gram=gram,gram_eigenvalues=eig,gram_condition=float(eig[-1]/eig[0]),cross_coupling=float(gram[0,1]/np.sqrt(gram[0,0]*gram[1,1])),predicted_squared_barrier=cost,predicted_radius_barrier=cost/(2*np.sqrt(t0)),maximum_joint_linear_shift=float(np.linalg.norm(V.sum(axis=1).reshape(100,2),axis=1).max()),initial_empty_min_nonunion_gap=float(ds[[i for i in range(100) if i not in union]].min()),initial_individual_event_checks=validators)
        candidates.append((case,grad,eh,V))
    candidates.sort(key=lambda z:z[0]['predicted_squared_barrier']);selected=[next(z for z in candidates if frozenset(z[0]['ids'])==weak)];used={frozenset(selected[0][0]['union'])}
    for z in candidates:
        if frozenset(z[0]['union']) in used:continue
        selected.append(z);used.add(frozenset(z[0]['union']))
        if len(selected)>=8:break
    outputs=[];preds={}
    for index,(case,grad,eh,V) in enumerate(selected):
        case['path_index']=index;gaps=np.asarray(case['gaps']);ng=grad/gaps[:,None];nh=eh/gaps[:,None,None];eps=1e-26;der=[]
        for a in range(2):der.append(np.imag(f.evaluate(C.astype(complex)+1j*eps*V[:,a].reshape(100,2),(m['T'],m['E'],m['sg']))[1][m['active']])/eps)
        W2=np.zeros((2,2,201));q0=np.asarray(case['q0']);common=case['common'];rel=C[common]-q0;Maux=np.c_[-2*rel/t0,-np.ones(3)];X=np.zeros((203,2));X[:199]=Q.T@V
        for a in range(2):X[200:,a]=np.linalg.solve(Maux,-2*np.sum(rel*V[:,a].reshape(100,2)[common],axis=1)/t0)
        XX=np.zeros((2,2,203));contact_res=[];event_res=[]
        for a in range(2):
            for b in range(a,2):
                hb=.5*(der[a]@V[:,b]+der[b]@V[:,a]);raw=np.linalg.lstsq(B,np.r_[-hb,0.],rcond=1e-11)[0];eb=ng@raw[:200]+np.einsum('i,nij,j->n',V[:,a],nh,V[:,b]);raw[:200]+=V@eb;W2[a,b]=W2[b,a]=raw
                xp=np.zeros(203);xp[:199]=Q.T@raw[:200];xp[199]=raw[-1]/t0
                va=V[:,a].reshape(100,2)[common]-X[200:202,a];vb=V[:,b].reshape(100,2)[common]-X[200:202,b]
                rhs=-(2*np.sum(rel*raw[:200].reshape(100,2)[common],axis=1)+2*np.sum(va*vb,axis=1))/t0
                xp[200:]=np.linalg.solve(Maux,rhs);XX[a,b]=XX[b,a]=xp
                contact_res.append(float(np.max(np.abs(B@raw-np.r_[-hb,0.]))));event_res.append(float(np.max(np.abs(ng@raw[:200]+np.einsum('i,nij,j->n',V[:,a],nh,V[:,b])))))
        x0=np.r_[np.zeros(199),1.,q0,case['rho0_squared']/t0];prob=PhysicalProblem(C,Q,w,t0,case,0.,C)
        eqJ=prob.eqjac(x0);activeJ=np.zeros((152,203));activeJ[:,:199]=m['A']@Q/t0;activeJ[:,199]=-1;aug=np.r_[activeJ,eqJ]
        firstcontrol=np.zeros((5,2));firstcontrol[3:]=np.diag(gaps/t0);firstres=float(np.max(np.abs(eqJ@X+firstcontrol)))
        mixedres=[]
        unionordered=common+case['extras']
        for a in range(2):
            for b in range(2):
                ca=(Q@X[:199,a]).reshape(100,2)[unionordered]-X[200:202,a];cb=(Q@X[:199,b]).reshape(100,2)[unionordered]-X[200:202,b]
                mixedres.append(float(np.max(np.abs(eqJ@XX[a,b]+2*np.sum(ca*cb,axis=1)/t0))))
        case.update(initial_aux_equalities_residual=float(np.max(np.abs(prob.eq(x0)))),initial_physical_J_rank=int(np.linalg.matrix_rank(eqJ,tol=1e-9)),initial_contact_plus_physical_J_rank=int(np.linalg.matrix_rank(aug,tol=1e-9)),initial_contact_plus_physical_rows=len(aug),first_physical_response_residual=firstres,mixed_physical_response_residual=max(mixedres),mixed_contact_response_residual=max(contact_res),mixed_gap_response_residual=max(event_res),mixed_curvature_gram_inverse_error=float(np.max(np.abs(V.T@H@V-np.linalg.inv(np.asarray(case['physical_gap_gram']))))),mixed_scope='Complete second-order geometric equality response; not stationary third-derivative KKT sensitivity')
        if max(firstres,max(mixedres),max(contact_res),max(event_res))>1e-9 or case['initial_physical_J_rank']!=5 or case['initial_contact_plus_physical_J_rank']!=157:raise ValueError('Physical local predictor/rank verification failure')
        preds[index]=dict(X=X,XX=XX,V=V,W2=W2,x0=x0,physical_gap_gradient=grad,physical_gap_hessian=eh)
        np.savez_compressed(ROOT/f'predictor{index:02d}.npz',**preds[index],Z=Z,G=G,gram=np.asarray(case['physical_gap_gram']),common=np.array(common),extras=np.array(case['extras']))
        outputs.append(case)
    checks=dict(input_sha256={str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'inputs').iterdir())},initial_tangent_dimension=48,finite_variable_count=203,coordinate_layout='0:199 fixed-rotation center displacement;199=t/t0;200:202=q;202=rho_squared/t0',rank_A=151,rank_A_gauge=152,rank_B_eq=153,G_min_eigenvalue=float(np.linalg.eigvalsh(G).min()),physical_candidates=len(candidates),selected_paths=len(outputs),geometry_budgets=dict(max_true_radius=.125,max_center_shift=.3,max_q_norm=.999,common_area_fraction=.08,cycle_edge_area_fraction=.05,min_rho=.03),one_CPU_BLAS_thread=True,source_actual_coverage_centers=100,rattler_scope='99 binding centers are a family fingerprint only; disk65 remains in all empty-circle and complete-cover checks and is never deleted',selected=outputs)
    write(ROOT/'input_checks.json',checks);write(ROOT/'selected_cavities.json',outputs);write(ROOT/'screened_cavities.json',[z[0] for z in candidates]);write(ROOT/'screen_rejections.json',rejections)
    return C,Q,w,t0,outputs,preds,checks


class PhysicalProblem:
    def __init__(self,C0,Q,w,t0,case,s,seed,trust=.055):
        self.C0=C0;self.Q=Q;self.w=w;self.t0=t0;self.case=case;self.s=s;self.nonunion=np.array([i for i in range(100) if i not in case['union']]);self.last=None
        self.core=f.ForcedProblem(C0,Q,w,t0,case['events'][0],0.,seed,trust);self.objgrad=np.zeros(203);self.objgrad[199]=1
    def unpack(self,x):return self.C0+(self.Q@x[:199]).reshape(100,2),x[199]*self.t0,x[200:202],x[202]*self.t0
    def eq(self,x):
        C,t,q,rho=self.unpack(x);sites=self.case['common']+self.case['extras'];offset=np.r_[np.zeros(3),(1-self.s)*np.asarray(self.case['gaps'])]
        return (np.sum((C[sites]-q)**2,axis=1)-rho-offset)/self.t0
    def eqjac(self,x):
        C,t,q,rho=self.unpack(x);sites=self.case['common']+self.case['extras'];rel=C[sites]-q;J=np.zeros((5,203));A=np.zeros((5,200));A[np.arange(5)[:,None],2*np.array(sites)[:,None]+np.array([0,1])]=2*rel/self.t0
        J[:,:199]=A@self.Q;J[:,200:202]=-2*rel/self.t0;J[:,202]=-1
        return J
    def ev(self,x):
        if self.last is not None and np.array_equal(self.last,x):return self.vals,self.J
        C,t,q,rho=self.unpack(x);corev,coreJ=self.core.evaluate(x[:200]);vals=[corev];Js=[np.c_[coreJ,np.zeros((len(corev),3))]]
        rel=C[self.nonunion]-q;v=(np.sum(rel*rel,axis=1)-rho)/self.t0;J=np.zeros((len(v),203));A=np.zeros((len(v),200));A[np.arange(len(v))[:,None],2*self.nonunion[:,None]+np.array([0,1])]=2*rel/self.t0;J[:,:199]=A@self.Q;J[:,200:202]=-2*rel/self.t0;J[:,202]=-1;vals.append(v);Js.append(J)
        J=np.zeros((1,203));J[0,200:202]=-2*q;vals.append(np.array([.999**2-q@q]));Js.append(J)
        tris=[self.case['common']]+self.case['cycle_triangles'];initial=np.r_[self.case['common_initial_area'],np.asarray(self.case['cycle_initial_areas'])];floors=np.r_[.08,np.full(15,.05)]
        vals.append(np.array([area(C,t)/a-f for t,a,f in zip(tris,initial,floors)]));J=np.zeros((16,203))
        for i,(tri,a0) in enumerate(zip(tris,initial)):J[i,:199]=area_gradient(C,tri)@self.Q/a0
        Js.append(J);self.last=x.copy();self.vals=np.concatenate(vals);self.J=np.concatenate(Js);return self.vals,self.J
    def con(self,x):return self.ev(x)[0]
    def jac(self,x):return self.ev(x)[1]


def seed_state(C,C0,Q,t0,case):
    g,q,_=geometry.triangles(C,np.array([case['common']]));R,_=geometry.covering(C)
    return np.r_[Q.T@(C-C0).ravel(),min(.125**2/t0,R*R/t0+1e-10),q[0],g[0]/t0]

def correct(C0,Q,w,t0,case,s,guess,previous,maxiter=300):
    x=guess.copy();logs=[]
    for outer in range(9):
        prob=PhysicalProblem(C0,Q,w,t0,case,s,previous);C=prob.unpack(x)[0];prob.core.top=geometry.topology(C);before=geometry.top_key(prob.core.top);tic=time.monotonic()
        bounds=[(None,None)]*203;bounds[199]=(.008/t0,.125**2/t0);bounds[202]=(.03**2/t0,.125**2/t0)
        res=minimize(lambda x:x[199],x,jac=lambda x:prob.objgrad,method='SLSQP',constraints=[dict(type='ineq',fun=prob.con,jac=prob.jac),dict(type='eq',fun=prob.eq,jac=prob.eqjac)],bounds=bounds,options=dict(maxiter=maxiter,ftol=2e-13,disp=False))
        x=res.x;C,t,q,rho=prob.unpack(x)
        if not np.isfinite(x).all():return x,False,dict(reason='nonfinite_state',rounds=logs)
        R,top=geometry.covering(C);eq=float(np.max(np.abs(prob.eq(x))));con=float(np.min(prob.con(x)));missing=set(f.keys((top[0],top[1],top[2])))-set(f.keys((prob.core.top[0],prob.core.top[1],prob.core.top[2])))
        logs.append(dict(round=outer,success=bool(res.success),status=int(res.status),message=str(res.message),iterations=int(res.nit),actual_radius=R,objective_t=t,physical_eq_scaled_residual=eq,minimum_scaled_inequality=con,q_norm=float(np.linalg.norm(q)),rho_squared=rho,topology_changed=before!=geometry.top_key(top),missing_actual_keys=sorted(missing),seconds=time.monotonic()-tic))
        feasible=bool(eq<3e-9 and con>-2e-8 and R*R<=t+3e-10 and R<=.125+1e-10)
        if feasible and not missing and (res.success or res.status==8):return x,True,dict(reason='feasible_corrected',rounds=logs)
        if not missing and before==geometry.top_key(top):return x,False,dict(reason='SLSQP_or_physical_feasibility_failure',rounds=logs)
    return x,False,dict(reason='actual_topology_refresh_limit',rounds=logs)

def validate(C,q,rho,case):
    ds=np.sum((C-q)**2,axis=1);union=case['union'];common=case['common'];cyc=np.array([area(C,t)/a for t,a in zip(case['cycle_triangles'],case['cycle_initial_areas'])]);spread=float(np.ptp(ds[union]));gap=float(ds.min()-rho)
    events=[f.validate_event(C,e) for e in case['events']];valid=bool(spread<3e-9 and gap>=-2e-9 and np.linalg.norm(q)<=.999+1e-9 and area(C,common)>=.08*case['common_initial_area']-1e-10 and cyc.min()>=.05-1e-8 and all(e['valid'] for e in events))
    return dict(valid=valid,q=q,rho_squared=rho,q_norm=float(np.linalg.norm(q)),five_squared_distance_spread=spread,all100_nearest_empty_gap=gap,minimum_cycle_area_ratio=float(cyc.min()),common_area_ratio=float(area(C,common)/case['common_initial_area']),individual_actual_events=events,scope='Numerical physical five-point circle, not an exact algebraic event/path certificate')

_FAMILY=None
def known_family(C):
    global _FAMILY
    if _FAMILY is None:
        K=np.load(BASE/'event_forced_paths/path03_e0206/restored_event_after.npy');gt=geometry.topology(K);g,A,kk=f.all_holes(K,gt);ac=np.flatnonzero(g.max()-g<1e-9);used=sorted({z for i in ac for z in (kk[i][1:4] if kk[i][0]=='t' else kk[i][1:3])});_FAMILY=(K,used,set(kk[i] for i in ac),f.true_cover(K)[0])
    K,used,oldkeys,R0=_FAMILY;ang=np.arctan2(np.sum(C[used,0]*K[used,1]-C[used,1]*K[used,0]),np.sum(C[used]*K[used]));rot=np.array([[np.cos(ang),-np.sin(ang)],[np.sin(ang),np.cos(ang)]]);distance=float(np.linalg.norm(C[used]-(K@rot)[used],axis=1).max());g,A,kk=f.all_holes(C,geometry.topology(C));currentkeys=set(kk[i] for i in np.flatnonzero(g.max()-g<1e-9));R,detail=f.true_cover(C);samekeys=currentkeys==oldkeys
    return dict(match=bool(distance<2e-4 and samekeys and abs(R-R0)<1e-8),binding_skeleton_distance=distance,active_hole_keys_match=samekeys,actual_all100_oracle_radius=R,reference_radius=R0,binding_center_count=len(used),unbound_center_indices=sorted(set(range(100))-set(used)),scope='Same-family recurrence fingerprint only; all100 oracle and active holes checked. An unbound disk remains essential to inactive-hole coverage and is never deleted')

def restore(C,case,side,folder,label,C0):
    family=known_family(C)
    if family['match']:
        name=f'restored_{label}_known_family.npy';np.save(folder/name,C);return dict(coordinate_file=name,true_radius=family['actual_all100_oracle_radius'],known_family=family,skipped=True,fully_unconstrained_valid=False,rounds=[])
    stages=[];best=None;bestR=np.inf
    for outer in range(7):
        top=geometry.topology(C);model=geometry.Model(C,top);R,_=geometry.covering(C);x0=np.r_[C.ravel(),R+1e-10];margin=2e-7 if outer==0 else 0.
        def sidefun(x):return np.array([side*f.incircle(x[:-1].reshape(100,2),e)/e['base_value_corrected']-margin for e in case['events']])
        def sidejac(x):return np.array([np.r_[side*f.event_jet(x[:-1].reshape(100,2),e)[1]/e['base_value_corrected'],0.] for e in case['events']])
        cons=[dict(type='ineq',fun=model.con,jac=model.jac)]
        if outer<2:cons.append(dict(type='ineq',fun=sidefun,jac=sidejac))
        res=minimize(lambda x:x[-1],x0,jac=lambda x:np.r_[np.zeros(200),1.],method='SLSQP',constraints=cons,bounds=[(max(-1,float(c)-.08),min(1,float(c)+.08)) for c in C.ravel()]+[(.09,.125)],options=dict(maxiter=450,ftol=2e-13,disp=False));C=res.x[:-1].reshape(100,2);R,detail=f.true_cover(C);changed=geometry.top_key(top)!=geometry.top_key(geometry.topology(C));family=known_family(C);valid=bool(np.linalg.norm(C,axis=1).max()<=1+1e-9 and np.linalg.norm(C-C0,axis=1).max()<=.3 and R<=.125+1e-10 and (outer>=2 or sidefun(res.x).min()>-2e-9));name=f'restoration_{label}_{outer:02d}.npy';np.save(folder/name,C)
        stages.append(dict(stage=outer,event_side_constraints=outer<2,physical_circle_constraints=False,coordinate_file=name,coordinate_sha256=f.center_sha(C),true_radius=R,true_cover=detail,valid=valid,success=bool(res.success),status=int(res.status),message=str(res.message),iterations=int(res.nit),known_family=family))
        if outer>=2 and valid and R<bestR:best=C.copy();bestR=R
        if family['match']:stages[-1]['stop_known_family']=True;break
        if outer>=2 and res.success and not changed:break
    if best is None:best=C.copy()
    name=f'restored_{label}.npy';np.save(folder/name,best);R,detail=f.true_cover(best)
    return dict(coordinate_file=name,coordinate_sha256=f.center_sha(best),file_sha256=sha(folder/name),true_radius=R,true_cover=detail,rounds=stages,skipped=False,fully_unconstrained_valid=any(s['valid'] and not s['event_side_constraints'] for s in stages),offsets=0,physical_circle_constraints=False,final_side_constraints=False,exact_binary64_event_signs=[f.exact_predicates(best,e)['incircle_sign'] for e in case['events']])

def run(case,pred,C0,Q,w,t0):
    index=case['path_index'];folder=ROOT/f'path{index:02d}_{case["ids"][0]}_{case["ids"][1]}';folder.mkdir(exist_ok=True);rows=[];accepted=[];failure=None;event=None;restorations=[];inventory=json.loads((ROOT/'inputs/events.json').read_text())
    def save(x,s,phase,log=None):
        prob=PhysicalProblem(C0,Q,w,t0,case,s,C0,trust=.4);C,t,q,rho=prob.unpack(x);R,detail=f.true_cover(C);top=geometry.topology(C);g,A,kk=f.all_holes(C,top);margins,inv=f.event_inventory(C,inventory);name=f'{len(rows):03d}_{phase}.npz';np.savez_compressed(folder/name,x=x,C=C,t=t,q=q,rho_squared=rho,s=s,T=top[0],E=top[1],sg=top[2],F=top[3],actual_hole_values=g,all_baseline_event_margins=margins)
        prob.core.top=top
        physical_rank=int(np.linalg.matrix_rank(prob.eqjac(x),tol=1e-9));row=dict(phase=phase,s=s,snapshot=name,coordinate_sha256=f.center_sha(C),file_sha256=sha(folder/name),true_radius=R,true_cover=detail,objective_t=t,q=q,rho_squared=rho,physical_eq_scaled_residual=float(np.max(np.abs(prob.eq(x)))),minimum_physical_inequality=float(prob.con(x)[len(g)+300:].min()),common_area_ratio=float(area(C,case['common'])/case['common_initial_area']),minimum_cycle_area_ratio=min(float(area(C,tri)/a) for tri,a in zip(case['cycle_triangles'],case['cycle_initial_areas'])),max_q_norm=float(np.linalg.norm(q)),max_center_shift=float(np.linalg.norm(C-C0,axis=1).max()),max_center_norm=float(np.linalg.norm(C,axis=1).max()),fixed_gauge_residual=float(w@(C-C0).ravel()),physical_J_rank=physical_rank,actual_hole_keys=kk,event_inventory=inv,exact_binary64_predicates=[f.exact_predicates(C,e) for e in case['events']],corrector=log)
        rows.append(row);write(folder/'steps.json',rows);return row
    schedule=[.0001,.0003,.001,.003,.01,.025,.055,.1,.17,.26,.37,.5,.63,.75,.85,.93,.975,.995,1.,1.001,1.005,1.02]
    if (folder/'steps.json').exists():
        rows=json.loads((folder/'steps.json').read_text())
        for row in rows:
            if row['phase']=='failed_attempt':continue
            arr=np.load(folder/row['snapshot']);assert sha(folder/row['snapshot'])==row['file_sha256'];accepted.append((row['s'],arr['x']))
        last,x=accepted[-1];pending=[s for s in schedule if s>last+1e-12];attempts=len(rows)-1
    else:x=pred['x0'];last=0.;save(x,0.,'start');accepted.append((0.,x.copy()));pending=list(schedule);attempts=0
    while pending and attempts<75:
        s=pending.pop(0);ds=s-last
        if last==0:pp=np.full(2,s);guess=pred['x0']+pred['X']@pp+.5*np.einsum('i,j,ijk->k',pp,pp,pred['XX'])
        elif len(accepted)>1:prev,xprev=accepted[-2];guess=x+ds*(x-xprev)/(last-prev)
        else:guess=x+ds*pred['X'].sum(axis=1)
        Cprev=C0+(Q@x[:199]).reshape(100,2);Cg=C0+(Q@guess[:199]).reshape(100,2);shift=np.linalg.norm(Cg-Cprev,axis=1).max()
        if shift>.043:Cg=Cprev+(Cg-Cprev)*(.043/shift);guess=seed_state(Cg,C0,Q,t0,case)
        np.savez_compressed(folder/f'attempt{attempts:03d}_seed.npz',x=guess,s=s,previous_x=x,previous_s=last)
        trial,ok,log=correct(C0,Q,w,t0,case,s,guess,Cprev);attempts+=1
        if not ok:
            save(trial,s,'failed_attempt',log)
            if ds>.0008 and attempts<75:pending.insert(0,s);pending.insert(0,(last+s)/2);continue
            failure=dict(reason=log['reason'],s=s,previous_s=last);break
        row=save(trial,s,'five_site_circle' if s==1 else ('after_circle' if s>1 else 'corrected'),log);x=trial;last=s;accepted.append((s,x.copy()))
        C,t,q,rho=PhysicalProblem(C0,Q,w,t0,case,s,C0).unpack(x)
        write(ROOT/'checkpoint.json',dict(status='running',current_path=index,current_ids=case['ids'],s=s,last_snapshot=str(folder/row['snapshot']),completed_paths=index,driver_sha256=sha(Path(__file__))))
        if s==1:
            event=validate(C,q,rho,case);row['physical_five_circle_validation']=event;write(folder/'steps.json',rows)
            if not event['valid']:failure=dict(reason='physical_circle_event_validation_failure',s=s,validation=event);break
        family=known_family(C)
        if family['match']:failure=dict(reason='known_family_with_actual_oracle_and_inventory',s=s,family=family);break
    before=[r for r in rows if r['phase']=='corrected'];after=[r for r in rows if r['phase']=='after_circle'];sources=[]
    if event and event['valid'] and before and after:
        event['exact_binary64_sign_brackets']=[[before[-1]['exact_binary64_predicates'][i]['incircle_sign'],after[-1]['exact_binary64_predicates'][i]['incircle_sign']] for i in range(2)];sources=[(before[-1],'before_circle',1),(after[-1],'after_circle',-1)]
    else:sources=[([r for r in rows if r['phase']!='failed_attempt'][-1],'last_accepted',1)]
    for row,label,side in sources:restorations.append(restore(np.load(folder/row['snapshot'])['C'],case,side,folder,label,C0))
    validrows=[r for r in rows if r['phase']!='failed_attempt'];result=dict(case=case,steps=rows,event_validation=event,restorations=restorations,stop=failure,accepted_nodes=len(validrows),failed_attempt_nodes=len(rows)-len(validrows),last_s=last,sampled_peak_true_radius=max(r['true_radius'] for r in validrows),best_fully_restored_radius=min([r['true_radius'] for r in restorations if r['fully_unconstrained_valid']],default=None),one_CPU_BLAS_thread=True,source_actual_coverage_centers=100,scope='Numerical physical five-site homotopy, sampled nodes only; no algebraic/full-path/global optimum certificate')
    write(folder/'result.json',result);print(json.dumps(f.clean(dict(path=index,ids=case['ids'],last_s=last,accepted=result['accepted_nodes'],failed=result['failed_attempt_nodes'],event_valid=event and event['valid'],peak=result['sampled_peak_true_radius'],restored=[r['true_radius'] for r in restorations],stop=failure))),flush=True)
    if any(r['true_radius']<.1155018241945471-1e-10 for r in restorations):print('MATERIAL_NEW_UPPER',index,[r['true_radius'] for r in restorations],flush=True)
    return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--paths',type=int,default=8);ap.add_argument('--prepare-only',action='store_true');args=ap.parse_args();assert 0<=args.paths<=8
    tic=time.monotonic();C,Q,w,t0,cases,preds,checks=prepare()
    if args.prepare_only:return
    results=[]
    for case in cases[:args.paths]:
        folder=ROOT/f'path{case["path_index"]:02d}_{case["ids"][0]}_{case["ids"][1]}';out=folder/'result.json'
        if out.exists():results.append(json.loads(out.read_text()));continue
        try:result=run(case,preds[case['path_index']],C,Q,w,t0)
        except Exception:folder.mkdir(exist_ok=True);result=dict(case=case,exception=traceback.format_exc());write(folder/'failure.json',result);print(json.dumps(f.clean(result)),flush=True)
        results.append(result);write(ROOT/'summary.json',dict(input_checks=checks,paths=results,completed_paths=len(results),seconds=time.monotonic()-tic,scope='Physical five-circle numerical searches only'))
        write(ROOT/'checkpoint.json',dict(status='completed' if len(results)==args.paths else 'running',completed_paths=len(results),requested_paths=args.paths,current_path=None,driver_sha256=sha(Path(__file__))))
    print('DONE',len(results),'seconds',time.monotonic()-tic,flush=True)

if __name__=='__main__':main()
