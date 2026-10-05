"""Exact polynomial / interval helpers; standard library only, no float acceptance."""
from fractions import Fraction as F
from itertools import permutations
from math import gcd, isqrt

def require(ok, msg):
    if not ok:
        raise ValueError(msg)

def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:p.pop()
    return p

def add(a,b):
    return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0)
                 for i in range(max(len(a),len(b)))])

def neg(a):return [-x for x in a]
def sub(a,b):return add(a,neg(b))
def scale(a,b):return trim([x*b for x in a])

def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return trim(out)

def sq(a):return mul(a,a)
def peval(p,x):
    y=0
    for a in reversed(p):y=y*x+a
    return y

def divrem(a,b):
    a=[F(x) for x in trim(a)];b=[F(x) for x in trim(b)]
    require(b!=[0], 'polynomial division by zero')
    q=[F(0)]*max(1,len(a)-len(b)+1)
    while a!=[0] and len(a)>=len(b):
        j=len(a)-len(b);c=a[-1]/b[-1];q[j]=c
        for i,v in enumerate(b):a[i+j]-=c*v
        a=trim(a)
    return trim(q),trim(a)

def derivative(p):return trim([i*p[i] for i in range(1,len(p))] or [0])
def primitive(p):
    require(all(type(x) is int for x in p),'noninteger coefficient')
    g=0
    for x in p:g=gcd(g,abs(x))
    require(g>0,'zero polynomial')
    return [x//g for x in trim(p)]

def sturm(p):
    seq=[trim(p),derivative(p)]
    require(seq[1]!=[0],'constant event polynomial')
    while True:
        _,r=divrem(seq[-2],seq[-1])
        if r==[0]:break
        r=neg(r)
        # Positive scaling reduces Fraction size without changing signs.
        factor=abs(r[-1]);seq.append([x/factor for x in r])
    return seq

def variations(seq,x):
    signs=[]
    for p in seq:
        y=peval(p,x)
        if y:signs.append((y>0)-(y<0))
    return sum(a!=b for a,b in zip(signs,signs[1:]))

def root_count(p,lo,hi):
    require(lo<hi and peval(p,lo)!=0 and peval(p,hi)!=0,'bad root count endpoints')
    s=sturm(p)
    return variations(s,lo)-variations(s,hi)

def iadd(a,b):return a[0]+b[0],a[1]+b[1]
def imul(a,b):
    v=[x*y for x in a for y in b]
    return min(v),max(v)
def interval(p,lo,hi):
    out=(F(0),F(0))
    for a in reversed(p):out=iadd(imul(out,(lo,hi)),(F(a),F(a)))
    return out

def coords(A,B):
    return [[[a,b-a] for a,b in zip(c,d)] for c,d in zip(A,B)]

def orient(a,b,c):
    u=[sub(x,y) for x,y in zip(b,a)];v=[sub(x,y) for x,y in zip(c,a)]
    return sub(mul(u[0],v[1]),mul(u[1],v[0]))

def norm2(a):return add(sq(a[0]),sq(a[1]))

def determinant(rows):
    n=len(rows);out=[0]
    for perm in permutations(range(n)):
        inversions=sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))
        term=[1]
        for i,j in enumerate(perm):term=mul(term,rows[i][j])
        out=add(out,scale(term,-1 if inversions%2 else 1))
    return out

def incircle(points):
    d=points[3]
    rows=[]
    for a in points[:3]:
        v=[sub(x,y) for x,y in zip(a,d)]
        rows.append([v[0],v[1],norm2(v)])
    return determinant(rows)

def circumcenter(a,b,c,Q):
    u=[sub(x,y) for x,y in zip(b,a)];v=[sub(x,y) for x,y in zip(c,a)]
    det=orient(a,b,c);u2=norm2(u);v2=norm2(v)
    nx=add(scale(mul(det,a[0]),2),sub(mul(u2,v[1]),mul(v2,u[1])))
    ny=add(scale(mul(det,a[1]),2),sub(mul(u[0],v2),mul(v[0],u2)))
    return nx,ny,scale(det,2*Q)

def distance_difference(a,e,nx,ny,den,Q):
    delta=[sub(x,y) for x,y in zip(e,a)]
    dot=add(mul(nx,delta[0]),mul(ny,delta[1]))
    return sub(mul(den,sub(norm2(e),norm2(a))),scale(dot,2*Q))

def ceil_sqrt(n):
    a=isqrt(n)
    return a+(a*a<n)

def displacement(A,B,Q):
    largest=max(sum((x-y)**2 for x,y in zip(a,b)) for a,b in zip(A,B))
    return F(ceil_sqrt(largest),Q),largest

def cap(U,V,D):
    # Exact maximum of min(U+s D, V+(1-s)D), s in [0,1].
    require(U>0 and V>0 and D>=0,'bad radius/displacement')
    if D==0:return min(U,V)
    if U>V+D:return V+D
    if V>U+D:return U+D
    return (U+V+D)/2

def triangle_witness(CP,s,triangle,Q):
    """Exact strict empty inside-disk Delaunay triangle at a rational s."""
    C=[[peval(p,s) for p in c] for c in CP]
    pts=[[[x] for x in C[i]] for i in triangle]
    nx,ny,den=circumcenter(*pts,Q);nx,ny,den=nx[0],ny[0],den[0]
    require(den!=0,'degenerate bracket triangle')
    q=(nx/den,ny/den)
    require(q[0]*q[0]+q[1]*q[1]<1,'bracket triangle circumcenter outside disk')
    distances=[(q[0]-c[0]/Q)**2+(q[1]-c[1]/Q)**2 for c in C]
    radius=distances[triangle[0]]
    require(all(distances[i]==radius for i in triangle),'bracket circumcenter identity failed')
    gaps=[(v-radius,i) for i,v in enumerate(distances) if i not in triangle]
    gap,i=min(gaps);require(gap>0,'bracket triangle not strictly empty')
    pair=lambda x:[x.numerator,x.denominator]
    return {'triangle':list(triangle),'q':[pair(x) for x in q],'radius_squared':pair(radius),
            'minimum_other_site_squared_gap':pair(gap),'minimum_gap_site':i,'strict_other_site_count':len(gaps)}
