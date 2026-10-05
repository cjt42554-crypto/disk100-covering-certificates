# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Two-positive-radical real arithmetic. No inverse or numerical optimizer."""
from fractions import Fraction as F
import itertools,functools,math
def sgn(x):return (x>0)-(x<0)
def root_sign(a,b,k):
    if not a:return sgn(b)
    if not b:return sgn(a)
    if sgn(a)==sgn(b):return sgn(a)
    return sgn(a)*sgn(a*a-b*b*k)
class Field:
    def __init__(self,k,l):self.k=F(k);self.l=F(l);assert self.k>0 and self.l>0
    def q(self,x):return (F(x),F(0),F(0),F(0))
    def add(self,a,b):return tuple(x+y for x,y in zip(a,b))
    def neg(self,a):return tuple(-x for x in a)
    def sub(self,a,b):return self.add(a,self.neg(b))
    def scale(self,a,b):return tuple(x*b for x in a)
    def mul(self,a,b):
        out=[F(0)]*4
        for i,x in enumerate(a):
            for j,y in enumerate(b):
                out[i^j]+=x*y*(self.k if i&j&1 else 1)*(self.l if i&j&2 else 1)
        return tuple(out)
    def sign(self,a):
        # (a0+a1 sqrt(k)) + (a2+a3 sqrt(k)) sqrt(l).
        u=root_sign(a[0],a[1],self.k);v=root_sign(a[2],a[3],self.k)
        if not u:return v
        if not v:return u
        if u==v:return u
        u2=(a[0]*a[0]+a[1]*a[1]*self.k,2*a[0]*a[1])
        v2=(a[2]*a[2]+a[3]*a[3]*self.k,2*a[2]*a[3])
        return u*root_sign(u2[0]-self.l*v2[0],u2[1]-self.l*v2[1],self.k)
    def fl(self,a):return float(a[0])+float(a[1])*math.sqrt(float(self.k))+float(a[2])*math.sqrt(float(self.l))+float(a[3])*math.sqrt(float(self.k*self.l))
    def enc(self,a):return [dict(numerator=str(x.numerator),denominator=str(x.denominator)) for x in a]
    def dot(self,a,b):return self.add(self.mul(a[0],b[0]),self.mul(a[1],b[1]))
    def cross(self,a,b):return self.sub(self.mul(a[0],b[1]),self.mul(a[1],b[0]))
    def orient(self,P,face):
        a,b,c=[P[v] for v in face];u=[self.sub(b[j],a[j]) for j in (0,1)];v=[self.sub(c[j],a[j]) for j in (0,1)];return self.cross(u,v)
    def height_det_coeff(self,P,quad,added):
        out=[]
        for site in added:
            if site not in quad:out.append(self.q(0));continue
            j=quad.index(site);other=quad[:j]+quad[j+1:];out.append(self.scale(self.orient(P,other),(-1)**(j+3)))
        return tuple(out)
    def height_det(self,P,quad,H):
        total=self.q(0)
        for j,site in enumerate(quad):total=self.add(total,self.scale(self.orient(P,quad[:j]+quad[j+1:]),(-1)**(j+3)*H[site]))
        return total
def triangulations(cycle):
    # Independent recursive polygon triangulation, not pairs-of-diagonals.
    def rec(poly):
        if len(poly)<3:return [()]
        result=[]
        for j in range(1,len(poly)-1):
            for left in rec(poly[:j+1]):
                for right in rec(poly[j:]):result.append(left+right+((poly[0],poly[j],poly[-1]),))
        return result
    edges={tuple(sorted((cycle[j],cycle[(j+1)%5]))) for j in range(5)};out=[]
    for faces in rec(tuple(cycle)):
        all_edges={tuple(sorted(e)) for f in faces for e in itertools.combinations(f,2)};out.append(dict(faces=[list(f) for f in faces],diagonals=sorted(all_edges-edges)))
    assert len(out)==5 and len({tuple(t['diagonals']) for t in out})==5;return out
def cone_rays(field,rows):
    candidates=[]
    for a,b in rows:
        for sign in (-1,1):
            ray=(field.scale(b,-sign),field.scale(a,sign))
            if all(field.sign(field.dot(row,ray))>=0 for row in rows):
                if not any(field.sign(field.cross(ray,old))==0 and field.sign(field.dot(ray,old))>0 for old in candidates):candidates.append(ray)
    assert len(candidates)==2
    inner=tuple(field.add(a,b) for a,b in zip(*candidates));assert all(field.sign(field.dot(row,inner))>0 for row in rows)
    return candidates,inner
def circle(P):
    a,b,c=P;u=[b[i]-a[i] for i in (0,1)];v=[c[i]-a[i] for i in (0,1)];d=2*(u[0]*v[1]-u[1]*v[0]);assert d
    U=sum(x*x for x in u);V=sum(x*x for x in v);w=((U*v[1]-V*u[1])/d,(u[0]*V-v[0]*U)/d);return tuple(a[i]+w[i] for i in (0,1)),sum(x*x for x in w)
def positions(source,shared,added,q,r2,levels):
    v={i:tuple(source[i][j]-q[j] for j in (0,1)) for i in added};norm={i:sum(x*x for x in v[i]) for i in added};k=[(r2+levels[i])/norm[i] for i in added];field=Field(*k);P={i:tuple(field.q(x) for x in c) for i,c in enumerate(source)}
    for idx,i in enumerate(added):P[i]=tuple(field.add(field.q(q[j]),tuple(v[i][j] if t==(1<<idx) else F(0) for t in range(4))) for j in (0,1))
    return field,P,k
def circle_inside(field,P,face):
    a,b,c=[P[v] for v in face];u=[field.sub(b[j],a[j]) for j in (0,1)];v=[field.sub(c[j],a[j]) for j in (0,1)];A=field.cross(u,v);U=field.dot(u,u);V=field.dot(v,v)
    Nx=field.add(field.scale(field.mul(A,a[0]),2),field.sub(field.mul(U,v[1]),field.mul(V,u[1])))
    Ny=field.add(field.scale(field.mul(A,a[1]),2),field.sub(field.mul(u[0],V),field.mul(v[0],U)))
    return field.sign(field.sub(field.scale(field.mul(A,A),4),field.add(field.mul(Nx,Nx),field.mul(Ny,Ny))))
