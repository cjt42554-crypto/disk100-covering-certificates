# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Full spatial derivatives of physical DT and unit-circle event functions."""
import numpy as np

class Jet:
    def __init__(self,value,g,H):self.v=float(value);self.g=np.asarray(g,float);self.H=np.asarray(H,float)
    @classmethod
    def variable(cls,value,k,n):
        g=np.zeros(n);g[k]=1;return cls(value,g,np.zeros((n,n)))
    def lift(self,x):return x if isinstance(x,Jet) else Jet(x,np.zeros_like(self.g),np.zeros_like(self.H))
    def __add__(self,x):
        x=self.lift(x);return Jet(self.v+x.v,self.g+x.g,self.H+x.H)
    __radd__=__add__
    def __neg__(self):return Jet(-self.v,-self.g,-self.H)
    def __sub__(self,x):return self+-self.lift(x)
    def __rsub__(self,x):return self.lift(x)+-self
    def __mul__(self,x):
        x=self.lift(x);return Jet(self.v*x.v,self.g*x.v+x.g*self.v,self.H*x.v+x.H*self.v+np.outer(self.g,x.g)+np.outer(x.g,self.g))
    __rmul__=__mul__
    def inverse(self):return Jet(1/self.v,-self.g/self.v**2,2*np.outer(self.g,self.g)/self.v**3-self.H/self.v**2)
    def __truediv__(self,x):return self*self.lift(x).inverse()
    def __rtruediv__(self,x):return self.lift(x)*self.inverse()

def dot(a,b):return sum(x*y for x,y in zip(a,b))
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def sub(a,b):return [x-y for x,y in zip(a,b)]

def circum(P):
    a,b,c=P;u=sub(b,a);v=sub(c,a);d=2*cross(u,v);uu=dot(u,u);vv=dot(v,v)
    pr=[(uu*v[1]-vv*u[1])/d,(u[0]*vv-v[0]*uu)/d]
    return [a[j]+pr[j] for j in (0,1)],dot(pr,pr)

def expression(P,kind):
    O,r2=circum(P[:3])
    if kind=='triangle_radius':return r2
    if kind=='dt_flip':return dot(sub(P[3],O),sub(P[3],O))-r2
    if kind=='circumcenter_unit':return 1-dot(O,O)
    raise ValueError(kind)

def event(C,kind,labels):return float(expression(np.asarray(C)[labels].tolist(),kind))

def jets(C,kind,labels):
    local=np.asarray(C)[labels].ravel();n=len(local);v=[Jet.variable(x,j,n) for j,x in enumerate(local)];P=[v[2*j:2*j+2] for j in range(len(labels))]
    x=expression(P,kind);cols=np.array([2*j+k for j in labels for k in (0,1)]);g=np.zeros(200);g[cols]=x.g
    return x.v,g,x.H,cols

def convex_quad(C,edge,opposites):
    a,b=C[list(edge)];c,d=C[list(opposites)];o=lambda a,b,c:float((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))
    v=[o(a,b,c),o(a,b,d),o(c,d,a),o(c,d,b)]
    return bool(v[0]*v[1]<0 and v[2]*v[3]<0),v

def finite_check(C,kind,labels,g,H,cols):
    d=np.arange(1,len(cols)+1,dtype=float);d=np.sin(d);d/=np.linalg.norm(d);h=2e-6;A=C.copy();B=C.copy();A.ravel()[cols]+=h*d;B.ravel()[cols]-=h*d
    fa=event(A,kind,labels);fb=event(B,kind,labels);f=event(C,kind,labels)
    ga=jets(A,kind,labels)[1][cols];gb=jets(B,kind,labels)[1][cols]
    return dict(step=h,gradient_directional_absolute_error=float(abs((fa-fb)/(2*h)-g[cols]@d)),hessian_vector_max_error=float(abs((ga-gb)/(2*h)-H@d).max()),second_difference_error=float(abs((fa-2*f+fb)/h**2-d@H@d)))
