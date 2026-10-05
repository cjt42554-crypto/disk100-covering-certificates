# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Exact first DT/circle/hull predicate roots along a binary64 endpoint chord.

All input endpoints are converted to integers with one dyadic denominator.
Integer polynomial Sturm chains certify absence of any earlier predicate root.
This is a straight chord certificate, not a certificate for a curved Newton
path, a radius, or unit-disk center necessity in the original covering problem.
"""
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[n]='1'
import sys,json,hashlib,math
from pathlib import Path
from fractions import Fraction
from collections import defaultdict
import numpy as np
from scipy.spatial import Delaunay
ROOT=Path(__file__).resolve().parent;sys.dont_write_bytecode=True;sys.path.insert(0,str(ROOT.parent))
from exact_actual_graph import verify

def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0:a.pop()
    return a
def add(a,b):return trim([(a[j] if j<len(a) else 0)+(b[j] if j<len(b) else 0) for j in range(max(len(a),len(b)))])
def neg(a):return [-x for x in a]
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for j,x in enumerate(a):
        for k,y in enumerate(b):out[j+k]+=x*y
    return trim(out)
def cross(u,v):return sub(mul(u[0],v[1]),mul(u[1],v[0]))
def dot(u,v):return add(mul(u[0],v[0]),mul(u[1],v[1]))
def psub(u,v):return (sub(u[0],v[0]),sub(u[1],v[1]))
def primitive(a):
    a=trim(a);g=math.gcd(*a)
    return [x//g for x in a] if g else [0]
def remainder(a,b):
    a=trim(a);b=trim(b)
    while a!=[0] and len(a)>=len(b):
        k=len(a)-len(b);lead=a[-1];c=b[-1];a=[abs(c)*x for x in a]
        for j,x in enumerate(b):a[j+k]-=lead*(1 if c>0 else -1)*x
        a=primitive(a)
    return a
def sturm(p):
    p=primitive(p);d=primitive([j*p[j] for j in range(1,len(p))] or [0]);seq=[p,d]
    while seq[-1]!=[0]:
        r=primitive(neg(remainder(seq[-2],seq[-1])))
        if r==[0]:break
        seq.append(r)
    return seq if d!=[0] else [p]
def signat(p,q):
    a,b=q.numerator,q.denominator;value=p[-1];power=b
    for c in p[-2::-1]:value=value*a+c*power;power*=b
    return (value>0)-(value<0)
def variation(seq,q):
    s=[signat(p,q) for p in seq];s=[x for x in s if x]
    return sum(a!=b for a,b in zip(s,s[1:]))
def isolate(p,bits=60):
    seq=sturm(p);v0=variation(seq,Fraction(0));v1=variation(seq,Fraction(1));count=v0-v1
    if count<=0:return None
    left=Fraction(0);right=Fraction(1)
    for _ in range(bits):
        mid=(left+right)/2
        if v0-variation(seq,mid)>0:right=mid
        else:left=mid
    return dict(left=[str(left.numerator),str(left.denominator)],right=[str(right.numerator),str(right.denominator)],
        roots_in_whole_segment=count,roots_before_left=v0-variation(seq,left),roots_in_bracket=variation(seq,left)-variation(seq,right),
        exact_sign_at_left=signat(p,left),exact_sign_at_right=signat(p,right),sturm_chain_sha256=hashlib.sha256(json.dumps(seq,separators=(',',':')).encode()).hexdigest())

def audit(A,B):
    A=np.asarray(A,float);B=np.asarray(B,float);assert A.shape==B.shape==(100,2)
    source=verify(A);assert source['status']=='PASS'
    rats=[[[float(v).as_integer_ratio() for v in row] for row in C] for C in (A,B)];Q=max(d for C in rats for row in C for n,d in row)
    ints=[[[n*(Q//d) for n,d in row] for row in C] for C in rats]
    P=[tuple([ints[0][j][k],ints[1][j][k]-ints[0][j][k]] for k in (0,1)) for j in range(100)]
    T=[];edges=defaultdict(list);candidates=[];initialzero=[]
    def append(kind,labels,p):
        p=primitive(p)
        if len(p)==1:return
        row=dict(kind=kind,labels=list(map(int,labels)),integer_polynomial_coefficients=list(map(str,p)),degree=len(p)-1)
        if p[0]==0:initialzero.append(row)
        else:candidates.append((row,p))
    for tri in Delaunay(A).simplices:
        f=list(map(int,tri));a,b,c=(P[j] for j in f);u=psub(b,a);v=psub(c,a);orient=cross(u,v)
        if orient[0]<0:f[1],f[2]=f[2],f[1];a,b,c=(P[j] for j in f);u=psub(b,a);v=psub(c,a);orient=cross(u,v)
        assert orient[0]>0;D=mul([2],orient);uu=dot(u,u);vv=dot(v,v);pr=(sub(mul(uu,v[1]),mul(vv,u[1])),sub(mul(u[0],vv),mul(v[0],uu)))
        N=(add(mul(a[0],D),pr[0]),add(mul(a[1],D),pr[1]));side=sub(mul([Q*Q],mul(D,D)),dot(N,N));append('CIRCUMCENTER_UNIT_CIRCLE',f,side);append('TRIANGLE_ORIENTATION_ZERO',f,orient)
        j=len(T);T.append(f)
        for k in range(3):edges[tuple(sorted((f[k],f[(k+1)%3])))].append((j,f[(k+2)%3]))
    for e,adj in edges.items():
        if len(adj)!=2:continue
        fi,d=adj[1];f=T[adj[0][0]];a,b,c=(psub(P[j],P[d]) for j in f);aa=dot(a,a);bb=dot(b,b);cc=dot(c,c)
        det=add(sub(mul(a[0],sub(mul(b[1],cc),mul(bb,c[1]))),mul(a[1],sub(mul(b[0],cc),mul(bb,c[0])))),mul(aa,cross(b,c)))
        assert det[0]<0;append('DT_COCIRCULAR_QUAD',f+[d],det)
    following={a:b for f in T for a,b in zip(f,f[1:]+f[:1]) if len(edges[tuple(sorted((a,b)))])==1};hull=[];j=min(following)
    while j not in hull:hull.append(j);j=following[j]
    for k,j in enumerate(hull):
        f=[hull[k-1],j,hull[(k+1)%len(hull)]];p=cross(psub(P[f[1]],P[f[0]]),psub(P[f[2]],P[f[0]]));assert p[0]>=0;append('FULL_HULL_COLLINEARITY',f,p)
    roots=[]
    for row,p in candidates:
        isolation=isolate(p)
        if isolation:roots.append(dict(row,root_isolation=isolation))
    if roots:
        frac=lambda pair:Fraction(int(pair[0]),int(pair[1]))
        first=min(frac(r['root_isolation']['right']) for r in roots)
        group=[r for r in roots if frac(r['root_isolation']['left'])<=first]
        lower=min(frac(r['root_isolation']['left']) for r in group)
        assert all(variation(sturm(p),Fraction(0))-variation(sturm(p),lower)==0 for _,p in candidates)
        bracket=dict(left=[str(lower.numerator),str(lower.denominator)],right=[str(first.numerator),str(first.denominator)])
    else:group=[];bracket=None
    return dict(status='PASS_EXACT_CHORD_FIRST_EVENT' if not initialzero else 'INITIAL_CELL_BOUNDARY_UNKNOWN',source_full_integer_graph_status=source['status'],common_denominator=str(Q),
        affine_coordinate_definition='Integer P(t)=A_integer+t*(B_integer-A_integer), divided by Q; no intermediate floating rounding in predicates',
        polynomial_predicates_checked=len(candidates),initial_zero_predicates=initialzero,predicate_roots_in_segment=len(roots),first_event_cluster=group,first_event_bracket=bracket,
        no_any_DT_circle_hull_predicate_root_before_first_bracket=bool(not initialzero),all_predicate_roots=roots,
        scope='Exact integer polynomial/Sturm first physical event of this endpoint chord only; no whole curved path, unitdisk-center necessity, radius, completeness or optimality assertion')

if __name__=='__main__':
    pa,pb,out=map(Path,sys.argv[1:]);d=audit(np.load(pa),np.load(pb));d.update(endpoint_A_file=str(pa.resolve()),endpoint_B_file=str(pb.resolve()),endpoint_A_file_sha256=hashlib.sha256(pa.read_bytes()).hexdigest(),endpoint_B_file_sha256=hashlib.sha256(pb.read_bytes()).hexdigest());out.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({k:v for k,v in d.items() if k not in ('first_event_cluster','all_predicate_roots')},indent=2))
