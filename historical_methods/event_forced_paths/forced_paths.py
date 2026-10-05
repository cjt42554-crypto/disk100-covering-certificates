# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Full-48-dimensional local prediction, actual-event-constrained continuation.

All reported paths/corrected minima are numerical. A sampled path peak is not
an interval supremum certificate, and a local minimum is not a global optimum.
The finite corrector optimizes ALL 199 fixed-rotation-slice coordinates, with
the current actual-hole inventory and one genuine incircle equality. No
released-contact offsets, source soft modes, or fixed original active keys.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json, hashlib, time, argparse, traceback
from fractions import Fraction
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize, lsq_linear
from scipy.spatial import Delaunay

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent
sys.path.insert(0,str(BASE))
from pro_core_audit.guarded_core import tangent_event_predictor as guarded_tangent_predictor
sys.path.insert(0,str(BASE/'event_release_diagnostic'))
from derive import evaluate, keys
from true_cover import true_cover
sys.path.insert(0,str(BASE/'latest_mixed_modes'))
import geometry_optimizer as geometry


def clean(v):
    if isinstance(v,np.ndarray):return clean(v.tolist())
    if isinstance(v,(np.integer,np.bool_)):return v.item()
    if isinstance(v,(np.floating,float)):return float(v) if np.isfinite(v) else None
    if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    return v

def write_json(path,value):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(clean(value),indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    tmp.replace(path)

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def center_sha(C):return hashlib.sha256(np.asarray(C,dtype='<f8').tobytes()).hexdigest()


class Jet:
    def __init__(self,v,g=None,H=None):
        self.v=v;self.g=np.zeros(8) if g is None else np.asarray(g);self.H=np.zeros((8,8)) if H is None else np.asarray(H)
    def __add__(self,b):
        if not isinstance(b,Jet):b=Jet(b)
        return Jet(self.v+b.v,self.g+b.g,self.H+b.H)
    __radd__=__add__
    def __neg__(self):return Jet(-self.v,-self.g,-self.H)
    def __sub__(self,b):return self+(-b if isinstance(b,Jet) else -b)
    def __rsub__(self,b):return -self+b
    def __mul__(self,b):
        if not isinstance(b,Jet):b=Jet(b)
        return Jet(self.v*b.v,self.g*b.v+b.g*self.v,self.H*b.v+b.H*self.v+np.outer(self.g,b.g)+np.outer(b.g,self.g))
    __rmul__=__mul__


def det3(rows):
    a,b,c=rows
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])

def incircle(C,event):
    sites=event['sites'];p=C[np.asarray(sites[:3])]-C[sites[3]]
    return -event['base_sign']*det3([[x,y,x*x+y*y] for x,y in p])

def event_jet(C,event):
    pts=[];eye=np.eye(8)
    for i in range(4):pts.append([Jet(C[event['sites'][i],j],eye[2*i+j]) for j in (0,1)])
    rows=[]
    for pt in pts[:3]:
        x,y=pt[0]-pts[3][0],pt[1]-pts[3][1]
        rows.append([x,y,x*x+y*y])
    out=-event['base_sign']*det3(rows)
    cols=np.array([2*i+k for i in event['sites'] for k in (0,1)])
    grad=np.zeros(200);grad[cols]=out.g
    H=np.zeros((200,200));H[np.ix_(cols,cols)]=out.H
    return float(out.v),grad,H

def exact_predicates(C,event):
    # Exact rational predicates of the SAVED binary64 coordinates, not a proof
    # that the numerically corrected event point has incircle exactly zero.
    p=[[Fraction(float(z)) for z in C[i]] for i in event['sites']]
    rows=[]
    for a in p[:3]:
        x,y=a[0]-p[3][0],a[1]-p[3][1];rows.append([x,y,x*x+y*y])
    e=-int(event['base_sign'])*det3(rows)
    def orient(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    values=[orient(p[0],p[1],p[2]),orient(p[0],p[1],p[3]),orient(p[2],p[3],p[0]),orient(p[2],p[3],p[1])]
    sign=lambda x:int(x>0)-int(x<0)
    return dict(incircle_sign=sign(e),incircle_exactly_zero=(e==0),orientation_signs=[sign(x) for x in values],convex=(values[0]*values[1]<0 and values[2]*values[3]<0),arithmetic='exact fractions of saved binary64 coordinates')

def validate_event(C,event,tolerance=2e-9):
    a,b,c=C[np.array(event['sites'][:3])];u=b-a;v=c-a
    D=u[0]*v[1]-u[1]*v[0]
    predicates=exact_predicates(C,event)
    if abs(D)<1e-12:return dict(valid=False,reason='degenerate_base_triangle',exact_predicates=predicates)
    q=a+np.array([(u@u)*v[1]-(v@v)*u[1],u[0]*(v@v)-v[0]*(u@u)])/(2*D)
    distances=np.sum((C-q)**2,axis=1);r2=float(distances[event['sites'][0]])
    four=distances[np.array(event['sites'])];tie=float(np.ptp(four))
    norm=float(np.linalg.norm(q));empty=float(distances.min()-r2)
    inside=norm<=1+tolerance
    convex=predicates['convex']
    # Base orientation must stay the original positive orientation, as well
    # as having a convex quadrilateral and a nearest empty circle.
    base_orientation=predicates['orientation_signs'][0]
    valid=bool(convex and base_orientation==event['base_orientation_sign'] and inside and empty>=-tolerance and tie<tolerance)
    return dict(valid=valid,inside_unit_disk=inside,common_circle_center=q,common_circle_center_norm=norm,common_circle_squared_radius=r2,common_circle_radius=np.sqrt(r2),four_site_squared_spread=tie,nearest_squared_gap=empty,nearest_center=int(np.argmin(distances)),convex=convex,base_orientation_preserved=base_orientation==event['base_orientation_sign'],exact_predicates=predicates,validation_tolerance=tolerance)


def all_holes(C,top):
    T,E,sg,F=top;g,A=evaluate(C,(T,E,sg))
    if len(F):
        fg,fj=geometry.single_boundary(C,F);fa=np.zeros((len(F),200),dtype=C.dtype)
        fa[np.arange(len(F))[:,None],2*F[:,None]+np.array([0,1])]=fj
        g=np.r_[g,fg];A=np.r_[A,fa]
    kk=keys((T,E,sg))+[('s',int(i)) for i in F]
    return g,A,kk


class ForcedProblem:
    def __init__(self,C0,Q,w,t0,event,target,seed,trust=.045):
        self.C0=C0;self.Q=Q;self.w=w;self.t0=t0;self.event=event;self.target=target
        self.scale=event['base_value_corrected'];self.seed=seed.copy();self.top=geometry.topology(seed)
        self.last=None;self.trust=trust
        self.objgrad=np.r_[np.zeros(199),1.]
    def unpack(self,x):return self.C0+(self.Q@x[:199]).reshape(100,2),x[-1]*self.t0
    def evaluate(self,x):
        if self.last is not None and np.array_equal(x,self.last):return self.v,self.j
        C,t=self.unpack(x);g,A,_=all_holes(C,self.top)
        norms=np.sum(C*C,axis=1);disp=C-self.C0
        delta=C-self.seed
        vals=np.r_[(t-g)/self.t0,1-norms,(.3**2-np.sum(disp*disp,axis=1))/.09,(self.trust**2-np.sum(delta*delta,axis=1))/self.trust**2]
        J=np.zeros((len(vals),200));J[:len(g),:199]=-A@self.Q/self.t0;J[:len(g),199]=1
        for block,arr,scale in ((0,C,1.),(1,disp,.09),(2,delta,self.trust**2)):
            blockJ=np.zeros((100,200));blockJ[np.arange(100)[:,None],2*np.arange(100)[:,None]+np.array([0,1])]=-2*arr/scale
            J[len(g)+block*100:len(g)+(block+1)*100,:199]=blockJ@self.Q
        self.last=x.copy();self.v=vals;self.j=J
        return vals,J
    def con(self,x):return self.evaluate(x)[0]
    def jac(self,x):return self.evaluate(x)[1]
    def eq(self,x):return np.array([incircle(self.unpack(x)[0],self.event)/self.scale-self.target])
    def eqjac(self,x):return np.r_[event_jet(self.unpack(x)[0],self.event)[1]@self.Q/self.scale,0.][None,:]


def correct(C0,Q,w,t0,event,target,guess,previous,maxiter=230,trust=.045):
    C=guess.copy();logs=[];result=None
    for outer in range(8):
        problem=ForcedProblem(C0,Q,w,t0,event,target,previous,trust)
        problem.top=geometry.topology(C)
        r,_=geometry.covering(C);x0=np.r_[Q.T@(C-C0).ravel(),r*r/t0+1e-10]
        before=geometry.top_key(problem.top);tic=time.monotonic()
        result=minimize(lambda x:x[-1],x0,jac=lambda x:problem.objgrad,method='SLSQP',constraints=[dict(type='ineq',fun=problem.con,jac=problem.jac),dict(type='eq',fun=problem.eq,jac=problem.eqjac)],bounds=[(None,None)]*199+[(.008/t0,.13**2/t0)],options=dict(maxiter=maxiter,ftol=2e-13,disp=False))
        C,t=problem.unpack(result.x)
        if not np.isfinite(C).all():return C,False,dict(reason='nonfinite_coordinates',rounds=logs)
        rr,actual=geometry.covering(C);eq=float(abs(problem.eq(result.x)[0]));mincon=float(np.min(problem.con(result.x)));changed=before!=geometry.top_key(actual)
        missing=set(keys((actual[0],actual[1],actual[2])))-set(keys((problem.top[0],problem.top[1],problem.top[2])))
        row=dict(round=outer,success=bool(result.success),status=int(result.status),message=str(result.message),iterations=int(result.nit),objective_t=t,actual_radius=rr,event_relative_residual=eq,min_scaled_constraint=mincon,topology_changed=changed,missing_actual_keys=sorted(missing),seconds=time.monotonic()-tic)
        logs.append(row)
        feasible=bool(eq<2e-9 and mincon>-2e-8 and rr*rr<=t+3e-10 and np.max(np.linalg.norm(C,axis=1))<=1+1e-9 and np.max(np.linalg.norm(C-C0,axis=1))<=.3+1e-9)
        # Updating at a flip removes obsolete contacts and admits both new
        # triangles when Qhull chooses the post-event diagonal.
        if feasible and not missing and (result.success or result.status==8):
            return C,True,dict(reason='feasible_corrected',rounds=logs,optimizer_success=bool(result.success),final_t=t)
        if not changed and not missing:return C,False,dict(reason='SLSQP_or_feasibility_failure',rounds=logs,final_t=t)
    return C,False,dict(reason='actual_topology_refresh_limit',rounds=logs,final_t=t)


def hessian(C,top,mu):
    H=np.empty((200,200));eps=1e-27
    for j in range(200):
        z=C.astype(complex);z.ravel()[j]+=1j*eps
        H[:,j]=np.imag(all_holes(z,top)[1].T@mu)/eps
    return (H+H.T)/2,float(np.max(np.abs(H-H.T)))

def stationary_diagnostic(C,Q,event):
    top=geometry.topology(C);g,A,kk=all_holes(C,top);T=float(np.max(g));active=np.flatnonzero(T-g<8e-10)
    Aa=A[active];eg=event_jet(C,event)[1]/event['base_value_corrected']
    # Stationarity for min t with g<=t and normalized event equality:
    # Aa.T mu + nu eg + eta w=0, sum(mu)=1. nu is unrestricted.
    M=np.zeros((200,len(active)+1));M[:199,:len(active)]=Q.T@Aa.T;M[:199,-1]=Q.T@eg;M[199,:len(active)]=1
    rhs=np.r_[np.zeros(199),1.]
    ls=lsq_linear(M,rhs,bounds=(np.r_[np.zeros(len(active)),-np.inf],np.full(len(active)+1,np.inf)),tol=1e-13,lsmr_tol=1e-13,max_iter=500)
    mu=ls.x[:-1];nu=float(ls.x[-1]);res=float(np.max(np.abs(M@ls.x-rhs)))
    B=np.r_[Aa@Q,(eg@Q)[None,:]];sv=np.linalg.svd(B,compute_uv=False);rank=int(np.linalg.matrix_rank(B,tol=1e-9));N=null_space(B,rcond=1e-9)
    Beq=np.r_[np.c_[Aa@Q,-np.ones(len(Aa))],np.r_[eg@Q,0.][None,:]]
    esv=np.linalg.svd(Beq,compute_uv=False);eqrank=int(np.linalg.matrix_rank(Beq,tol=1e-9))
    fm=np.zeros(len(g));fm[active]=mu;H,asym=hessian(C,top,fm);H+=nu*event_jet(C,event)[2]/event['base_value_corrected']
    eigen=np.linalg.eigvalsh(N.T@Q.T@H@Q@N) if N.shape[1] else np.empty(0)
    return dict(active_count=len(active),active_keys=[kk[i] for i in active],positive_inequality_multiplier_count=int(np.sum(mu>1e-10)),zero_stress_keys=[kk[active[i]] for i in np.flatnonzero(mu<1e-10)],minimum_multiplier=float(mu.min(initial=1)),event_multiplier=nu,stationarity_residual_inf=res,equality_gradient_rank=rank,epigraph_event_equality_rows=len(Beq),epigraph_event_equality_rank=eqrank,epigraph_event_smallest_singular_value=float(esv[-1]),epigraph_event_rank_deficiency=len(Beq)-eqrank,rank_scope='At exact incircle tie, the two neighboring circumradius differences are proportional to the event equality; redundant contact rows are recorded, not treated as independent equations',tangent_dimension=N.shape[1],smallest_singular_value=float(sv.min(initial=1)),tangent_H_eigenvalues=eigen,positive_tangent_curvature=bool(eigen.size and eigen.min()>1e-9),H_complex_step_asymmetry=asym,active_tolerance_squared=8e-10,numerical_only=True)


def event_inventory(C,events):
    margins=[];crossed=[];circlecross=[];boundarycross=[];near=[];pairpoints={};dtgeometry=[]
    for ev in events:
        if ev['kind']=='dt_flip':
            val=incircle(C,ev)
            if val<-1e-12:
                crossed.append(ev['id']);a,b,c=C[np.array(ev['sites'][:3])];u=b-a;v=c-a;d=u[0]*v[1]-u[1]*v[0]
                if abs(d)>1e-12:
                    q=a+np.array([(u@u)*v[1]-(v@v)*u[1],u[0]*(v@v)-v[0]*(u@u)])/(2*d);r2=float(np.sum((q-a)**2));nearest=float(np.sum((C-q)**2,axis=1).min())
                    dtgeometry.append(dict(id=ev['id'],base_triangle_circumcenter_norm=float(np.linalg.norm(q)),base_triangle_inside_unit_disk=bool(q@q<=1+1e-10),base_triangle_nearest_empty=bool(nearest>=r2-1e-10),scope='Geometry at the sampled state, not an exact reconstructed first-crossing point'))
        elif ev['kind']=='circle_crossing':
            a,b,c=C[np.array(ev['sites'])];u=b-a;z=c-a;d=u[0]*z[1]-u[1]*z[0]
            k=2*d*a+np.array([(u@u)*z[1]-(z@z)*u[1],u[0]*(z@z)-z[0]*(u@u)])
            val=ev['base_sign']*(4*d*d-k@k)
            if val<-1e-12:circlecross.append(ev['id'])
        else:
            key=tuple(ev['pair'])+(ev['sign'],)
            if key not in pairpoints:
                i,j=ev['pair'];a,b=C[[i,j]];u=b-a;ss=u@u;h=(b@b-a@a)/2;disc=ss-h*h
                pairpoints[key]=(h*u+ev['sign']*np.sqrt(disc)*np.array([-u[1],u[0]]))/ss if disc>=0 and ss>1e-20 else None
            q=pairpoints[key]
            if q is None:val=np.nan
            else:
                a=C[ev['pair'][0]];b=C[ev['competitor']];val=ev['base_sign']*((q-b)@(q-b)-(q-a)@(q-a))
            if np.isfinite(val) and val<-1e-12:boundarycross.append(ev['id'])
        margins.append(val)
        if np.isfinite(val) and abs(val)<1e-12:near.append(ev['id'])
    return np.array(margins),dict(DT_events_crossed=crossed,circle_events_crossed=circlecross,boundary_switch_events_crossed=boundarycross,near_zero_events=near,crossed_DT_geometry=dtgeometry,no_event=None if not(crossed or circlecross or boundarycross) else False,scope='All baseline primitive events checked at every node; current all-actual-hole topology independently refreshed')


def restore(C,event,side,folder,label,C0):
    # All original actual-hole inequalities, zero offsets. Preserve the side
    # of the target diagonal to avoid immediate return through the event.
    r0,_=geometry.covering(C);seed=C.copy();best=C.copy();bestR=r0;rounds=[];best_free=None;best_free_R=np.inf
    staged=[]
    for outer in range(7):
        top=geometry.topology(C);model=geometry.Model(C,top);x=np.r_[C.ravel(),geometry.covering(C)[0]+1e-10]
        eventscale=event['base_value_corrected'];margin=2e-7 if outer==0 else 0.
        def sidecon(x):return side*incircle(x[:-1].reshape(100,2),event)/eventscale-margin
        def sidejac(x):return np.r_[side*event_jet(x[:-1].reshape(100,2),event)[1]/eventscale,0.]
        constraints=[dict(type='ineq',fun=model.con,jac=model.jac)]
        if outer<2:constraints.append(dict(type='ineq',fun=sidecon,jac=sidejac))
        res=minimize(lambda x:x[-1],x,jac=lambda x:np.r_[np.zeros(200),1.],method='SLSQP',constraints=constraints,bounds=[(max(-1,float(z)-.08),min(1,float(z)+.08)) for z in C.ravel()]+[(.09,.13)],options=dict(maxiter=450,ftol=2e-13,disp=False))
        C=res.x[:-1].reshape(100,2);rr,nt=geometry.covering(C);changed=geometry.top_key(top)!=geometry.top_key(nt)
        valid=bool(np.isfinite(C).all() and np.linalg.norm(C,axis=1).max()<=1+1e-9 and (outer>=2 or sidecon(res.x)>-1e-9) and np.linalg.norm(C-C0,axis=1).max()<=.3)
        stage='retain_diagonal' if outer==0 else ('zero_side_margin' if outer==1 else 'all_actual_holes_no_event_constraint')
        sname=f'restoration_{label}_{outer:02d}_{stage}.npy';np.save(folder/sname,C);stageR,stageDetail=true_cover(C)
        rounds.append(dict(round=outer,stage=stage,success=bool(res.success),status=int(res.status),message=str(res.message),iterations=int(res.nit),actual_radius=rr,independent_true_radius=stageR,independent_true_cover=stageDetail,event_signed_relative_margin=sidecon(res.x),event_constraint_present=outer<2,coordinate_file=sname,coordinate_sha256=center_sha(C),topology_changed=changed,valid=valid))
        if valid and rr<bestR:best=C.copy();bestR=rr
        if outer>=2 and valid and rr<best_free_R:best_free=C.copy();best_free_R=rr
        if outer>=2 and not changed and res.success:break
    # The headline restoration must have passed the completely original
    # problem, after removing the temporary diagonal condition.
    if best_free is not None:best=best_free
    R,detail=true_cover(best);name=f'restored_{label}.npy';np.save(folder/name,best)
    return dict(source_phase=label,side=side,coordinate_file=name,coordinate_sha256=center_sha(best),file_sha256=sha(folder/name),true_radius=R,true_cover=detail,rounds=rounds,restoration_offsets=0,final_event_constraint_present=False,fully_unconstrained_stage_valid=best_free is not None,retained_event_diagonal_side=bool(side*incircle(best,event)>0),max_center_shift_from_baseline=float(np.linalg.norm(best-C0,axis=1).max()),numerical_only=True)


def prepare():
    m=np.load(ROOT/'inputs/matrices.npz');C=m['C'];Q=m['Q'];w=m['w'];t0=float(m['t']);A=m['A'];H=m['H']
    Z=null_space(np.r_[A,w[None,:]],rcond=1e-10);G=Z.T@H@Z;values=np.linalg.eigvalsh(G)
    B=np.zeros((153,201));B[:152,:200]=A;B[:152,200]=-1;B[-1,:200]=w
    assert Z.shape==(200,48) and np.linalg.matrix_rank(A)==151 and np.linalg.matrix_rank(np.r_[A,w[None,:]])==152 and np.linalg.matrix_rank(B)==153
    assert values.min()>0
    events=json.loads((ROOT/'inputs/events.json').read_text());ranked=[];jets=[]
    for ev in events:
        if ev['kind']!='dt_flip':continue
        event=dict(ev);e,grad,EH=event_jet(C,event);event['base_value_corrected']=e
        event['base_orientation_sign']=exact_predicates(C,event)['orientation_signs'][0]
        validation=validate_event(C,event)
        # Baseline circle need only be an inside nearest witness; the fourth
        # site's distance is greater until the forced event.
        if not(validation.get('inside_unit_disk') and validation.get('convex') and validation.get('nearest_squared_gap',-1)>=-1e-10):continue
        a=Z.T@grad;b=np.linalg.solve(G,a);den=float(a@b)
        if den<=1e-26:continue
        v,pred_guarded=guarded_tangent_predictor(Z,H,grad,e,contact_gradients=A,gauge=w)
        check_v=-e*Z@b/den
        if np.max(np.abs(v-check_v))>2e-12:raise ValueError('Guarded Pro tangent predictor mismatch')
        eps=1e-26;_,aj=evaluate(C.astype(complex)+1j*eps*v.reshape(100,2),(m['T'],m['E'],m['sg']))
        bholes=np.imag(aj[m['active']])@v/eps
        second=np.linalg.lstsq(B,np.r_[-bholes,0.],rcond=1e-11)[0];kappa=float(second[-1]);pred=e*e/(2*den)
        record=dict(event=event,predicted_squared_barrier=pred,predicted_radius_barrier=pred/(2*np.sqrt(t0)),event_gradient_tangent_norm=float(np.linalg.norm(a)),denominator=den,maximum_linear_center_shift=float(np.linalg.norm(v.reshape(100,2),axis=1).max()),maximum_second_derivative_center_shift=float(np.linalg.norm(second[:200].reshape(100,2),axis=1).max()),second_equality_residual_inf=float(np.max(np.abs(B@second-np.r_[-bholes,0.]))),kappa=kappa,energy_kappa=float(v@H@v),initial_inactive_squared_gap=float(m['inactive_slack'].min()),local_only=True)
        ranked.append(record);jets.append(dict(id=event['id'],event=event,value=e,gradient=grad,hessian=EH,tangent_gradient=a,linear_center_direction=v,center_second_derivative=second[:200],t_second_derivative=kappa))
    ranked.sort(key=lambda z:z['predicted_squared_barrier'])
    first=['e0090','e0093','e0094'];selected=[next(r for r in ranked if r['event']['id']==e) for e in first]
    for record in ranked:
        if record['event']['id'] not in first:selected.append(record)
        if len(selected)>=12:break
    np.savez_compressed(ROOT/'full48_geometry.npz',C=C,Q=Q,w=w,A=A,H=H,Z=Z,G=G,G_eigenvalues=values,B_eq=B,t0=t0)
    np.savez_compressed(ROOT/'full48_event_jets.npz',event_ids=np.array([j['id'] for j in jets]),values=np.array([j['value'] for j in jets]),gradients=np.array([j['gradient'] for j in jets]),hessians=np.array([j['hessian'] for j in jets]),tangent_gradients=np.array([j['tangent_gradient'] for j in jets]),linear_center_directions=np.array([j['linear_center_direction'] for j in jets]),center_second_derivatives=np.array([j['center_second_derivative'] for j in jets]),t_second_derivatives=np.array([j['t_second_derivative'] for j in jets]))
    R,details=true_cover(C)
    manifest=dict(input_files={str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'inputs').iterdir()) if p.is_file()},baseline_coordinate_sha256=center_sha(C),baseline_true_radius=R,independent_baseline_score=details,rank_A=151,rank_A_gauge=152,rank_B_eq=153,full_tangent_dimension=48,finite_slice_dimension=199,positive_G_min_eigenvalue=float(values.min()),positive_G_max_eigenvalue=float(values.max()),gauge_Q_orthogonality=float(np.max(np.abs(w@Q))),tangent_residual_inf=float(np.max(np.abs(np.r_[A,w[None,:]]@Z))),max_paths=12,one_CPU_BLAS_thread=True,geometry_budgets=dict(true_radius=.13,max_center_shift=.3),method='full48 quadratic t-objective event predictor, full199 current actual-hole epigraph minimization with exact-form incircle equality, dynamic inventory refresh, independent true_cover scoring',local_prediction_scope='Taylor optimum only for sufficiently small event progress on a fixed regular geometric stratum; no finite-event or global optimality guarantee',sample_peak_scope='maximum over saved numerical nodes; no certified whole-segment supremum',selected=selected)
    write_json(ROOT/'input_checks.json',manifest);write_json(ROOT/'ranked_full48_events.json',ranked)
    return C,Q,w,t0,Z,G,selected,{j['id']:j for j in jets},manifest


def run_path(index,record,jet,C0,Q,w,t0):
    folder=ROOT/f'path{index:02d}_{record["event"]["id"]}';folder.mkdir(exist_ok=True)
    event=record['event'];v=jet['linear_center_direction'];w2=jet['center_second_derivative'];rows=[];accepted=[];restores=[];failure=None
    if not(folder/'predictor.npz').exists():
        np.savez_compressed(folder/'predictor.npz',C0=C0,Z=np.load(ROOT/'full48_geometry.npz')['Z'],linear_center_direction=v,center_second_derivative=w2,event_gradient=jet['gradient'],event_hessian=jet['hessian'],tangent_gradient=jet['tangent_gradient'],t_second_derivative=jet['t_second_derivative'])
    def save(C,fraction,phase,correction=None):
        R,detail=true_cover(C);top=geometry.topology(C);dt={tuple(sorted(map(int,t))) for t in Delaunay(C).simplices}
        baselineDT={tuple(sorted(map(int,t))) for t in Delaunay(C0).simplices}
        name=f'{len(rows):03d}_{phase}.npz';g,_,kk=all_holes(C,top)
        raw_events=json.loads((ROOT/'inputs/events.json').read_text());margins,other=event_inventory(C,raw_events)
        np.savez_compressed(folder/name,C=C,fraction=fraction,T=top[0],E=top[1],sg=top[2],F=top[3],actual_hole_g=g,event_value=incircle(C,event),all_baseline_event_margins=margins)
        curpred=exact_predicates(C,event)
        # Complete baseline event inventory, checked at every saved step.
        dtcross=other['DT_events_crossed'];circlecross=other['circle_events_crossed']
        previous_keys=set(map(tuple,rows[-1]['actual_hole_keys'])) if rows else set(kk)
        row=dict(phase=phase,event_fraction=fraction,event_target_relative=1-fraction,snapshot=name,file_sha256=sha(folder/name),coordinate_sha256=center_sha(C),true_radius=R,independent_true_cover=detail,event_value=incircle(C,event),event_relative=incircle(C,event)/event['base_value_corrected'],fixed_gauge_residual=float(w@(C-C0).ravel()),max_center_shift=float(np.linalg.norm(C-C0,axis=1).max()),max_center_norm=float(np.linalg.norm(C,axis=1).max()),actual_hole_count=len(g),actual_hole_keys=kk,actual_hole_inventory_added=sorted(set(kk)-previous_keys),actual_hole_inventory_removed=sorted(previous_keys-set(kk)),added_DT_triangles=sorted(dt-baselineDT),removed_DT_triangles=sorted(baselineDT-dt),other_baseline_DT_events_crossed=dtcross,other_baseline_circle_events_crossed=circlecross,complete_baseline_event_inventory=other,exact_predicates=curpred,corrector=correction,numerical_only=True)
        rows.append(row);write_json(folder/'steps.json',rows)
        return row
    schedule=[.0001,.0003,.001,.003,.01,.03,.07,.13,.22,.34,.48,.62,.75,.85,.93,.975,.995,1.,1.002,1.01,1.035]
    if (folder/'steps.json').exists():
        rows=json.loads((folder/'steps.json').read_text())
        for row in rows:
            if row['phase']=='failed_attempt':continue
            state=np.load(folder/row['snapshot'])['C']
            if center_sha(state)!=row['coordinate_sha256'] or sha(folder/row['snapshot'])!=row['file_sha256']:
                raise ValueError('Resume snapshot hash mismatch')
            if row['true_radius']>.13 or row['max_center_shift']>.3:continue
            accepted.append((row['event_fraction'],state.copy()))
        if not accepted:raise ValueError('Checkpoint has no accepted finite state')
        last,C=accepted[-1];pending=[x for x in schedule if x>last+1e-12];attempts=len(rows)-1
        write_json(folder/'resume_receipt.json',dict(accepted_nodes_preserved=len(accepted),saved_attempt_nodes_preserved=len(rows)-1,resumed_fraction=last,optimization_paths_repeated=0,driver_sha256=sha(Path(__file__))))
    else:
        C=C0.copy();save(C,0.,'start');accepted.append((0.,C.copy()))
        last=0.;pending=list(schedule);attempts=0
    while pending and attempts<70:
        fraction=pending.pop(0);df=fraction-last
        if last==0:
            guess=C0+(fraction*v+.5*fraction*fraction*w2).reshape(100,2)
        elif len(accepted)>=2:
            fp,Cp=accepted[-2];slope=(C-Cp)/(last-fp);guess=C+df*slope
        else:guess=C+(df*v).reshape(100,2)
        # Prediction is only a seed; the actual equalities and all current
        # holes determine every accepted finite state.
        prediction_shift=np.linalg.norm(guess-C,axis=1).max()
        if prediction_shift>.035:guess=C+(guess-C)*(.035/prediction_shift)
        np.savez_compressed(folder/f'attempt{attempts:03d}_seed.npz',C=guess,previous=C,fraction=fraction,previous_fraction=last)
        trial,ok,log=correct(C0,Q,w,t0,event,1-fraction,guess,C)
        attempts+=1
        if not ok:
            save(trial,fraction,'failed_attempt',log)
            if df>.0008 and attempts<70:
                midpoint=(last+fraction)/2;pending.insert(0,fraction);pending.insert(0,midpoint);continue
            failure=dict(reason=log['reason'],fraction=fraction,previous_fraction=last,attempts=attempts);break
        row=save(trial,fraction,'event' if fraction==1 else ('event_after' if fraction>1 else 'corrected'),log)
        row['stationary_diagnostic']=stationary_diagnostic(trial,Q,event);write_json(folder/'steps.json',rows)
        if row['true_radius']>.13 or row['max_center_shift']>.3 or row['max_center_norm']>1+1e-9:
            failure=dict(reason='geometry_budget',fraction=fraction);break
        C=trial;last=fraction;accepted.append((fraction,C.copy()))
        write_json(ROOT/'checkpoint.json',dict(current_path=index,current_event=event['id'],current_fraction=last,last_snapshot=str(folder/row['snapshot']),completed_paths=index,status='running',updated_unix=time.time(),driver_sha256=sha(Path(__file__))))
    event_states=[row for row in rows if row['phase']=='event'];event_validation=None
    if event_states:
        Ce=np.load(folder/event_states[-1]['snapshot'])['C'];event_validation=validate_event(Ce,event)
        before=[row for row in rows if row['phase']=='corrected' and row['event_fraction']<1]
        after=[row for row in rows if row['phase']=='event_after']
        if before and after and event_validation['valid']:
            br=before[-1];ar=after[-1]
            for row,label,side in ((br,'event_before',1),(ar,'event_after',-1)):
                Cr=np.load(folder/row['snapshot'])['C'];restores.append(restore(Cr,event,side,folder,label,C0))
            event_validation.update(numerical_crossing=True,before_snapshot=br['snapshot'],after_snapshot=ar['snapshot'],exact_binary64_sign_bracket=[br['exact_predicates']['incircle_sign'],ar['exact_predicates']['incircle_sign']],new_diagonal=event['sites'][2:])
    if not restores:
        row=max((r for r in rows if r['phase']!='failed_attempt'),key=lambda r:r['event_fraction']);Cr=np.load(folder/row['snapshot'])['C'];side=1 if row['event_relative']>=0 else -1
        restores.append(restore(Cr,event,side,folder,'last_accepted',C0))
    result=dict(path_index=index,event=event,initial_prediction=record,steps=rows,accepted_steps=sum(row['phase']!='failed_attempt' for row in rows),corrector_attempts=attempts,stop=failure,event_validation=event_validation,restorations=restores,sampled_peak_true_radius=max(r['true_radius'] for r in rows if r['phase']!='failed_attempt'),sampled_peak_scope='Saved accepted nodes only, not a certified whole-path supremum',best_restored_true_radius=min(r['true_radius'] for r in restores),baseline_true_radius=true_cover(C0)[0],numerical_only=True,global_optimum_claim=False)
    write_json(folder/'result.json',result)
    print(json.dumps(clean(dict(path_index=index,event=event['id'],accepted_steps=result['accepted_steps'],attempts=attempts,last_fraction=last,stop=failure,event_valid=event_validation and event_validation['valid'],peak=result['sampled_peak_true_radius'],restored=[r['true_radius'] for r in restores])),ensure_ascii=False),flush=True)
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--paths',type=int,default=12);ap.add_argument('--prepare-only',action='store_true');args=ap.parse_args()
    if not(0<=args.paths<=12):raise ValueError('First round is limited to twelve paths')
    old=json.loads((ROOT/'checkpoint.json').read_text()) if (ROOT/'checkpoint.json').exists() else None
    # Resume completed outputs, and never silently overwrite accepted results.
    tic=time.monotonic();C0,Q,w,t0,Z,G,selected,jets,manifest=prepare()
    if args.prepare_only:return
    results=[]
    for index,record in enumerate(selected[:args.paths]):
        resultfile=ROOT/f'path{index:02d}_{record["event"]["id"]}'/'result.json'
        if resultfile.exists():results.append(json.loads(resultfile.read_text()));continue
        try:result=run_path(index,record,jets[record['event']['id']],C0,Q,w,t0)
        except Exception:
            failure=dict(path_index=index,event=record['event']['id'],exception=traceback.format_exc())
            write_json(ROOT/f'path{index:02d}_{record["event"]["id"]}'/'failure.json',failure);print(json.dumps(failure),flush=True);result=failure
        results.append(result)
        summary=dict(input_checks=manifest,paths=results,completed_paths=len(results),elapsed_this_run_seconds=time.monotonic()-tic,scope='Full48 local predictors and full199 numerical event-forced paths, sampled coverage, no exact event/path/global optimum certificate',best_restored_true_radius=min([r.get('best_restored_true_radius',1.) for r in results]))
        write_json(ROOT/'summary.json',summary);write_json(ROOT/'checkpoint.json',dict(completed_paths=len(results),current_path=None,status='running' if len(results)<args.paths else 'completed',requested_paths=args.paths,updated_unix=time.time(),driver_sha256=sha(Path(__file__))))
    print('DONE',len(results),'seconds',time.monotonic()-tic,flush=True)


if __name__=='__main__':main()
