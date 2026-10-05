# A two-face actual signed-anchor counterexample with actual arc cut above .62

The new fixed algebraic point below has two actual marked-face scores
summing to `2C-D`, where

```
R=5661/50000=.11322,
lambda=7/50, t=1071883/20000000,
D=sin(theta*)/36 < 2228272803/100000000000=.02228272803.
```

All three *actual assigned boundary arcs* have
`xi=sin(full_arc_angle/2)/R>3/5`; the two outer arcs even have
`xi>31/50=.62`. These are actual disk/unit-circle footprints, not a
base-length proxy. Every assigned center is the inner radial endpoint
for its complete short cap, so its assigned BC lens slack is zero.
Both actual signed edge-anchor slacks are also zero; these are separate
terms and are not charged again. Thus an arc cut of `.6`, or even `.62`,
alone cannot yield the parent's requested two-face fee `.02228272803`.

This is a local counterexample, not a full unit-circle or unit-disk
cover. It does not verify the proposed global `Q<=.35754`, remaining
arc restrictions, ordinary reduction, global incidence, or an optimum.
No optimum of the outer arc cut over all admissible two-face diagrams
is claimed: this is one exact witness below the specified loss budget.

## 1. Previously known zero-loss scope

The parent already has the exact unskewed two-face saturation family:
two copies of the scalar maximizing triangle can be translated so
their outer lens tips are unit anchors. Their total actual signed-A
loss is zero, but their outer full footprints have assigned half-sine
about `.5514`, while the middle footprint is larger. That example
handles weak cuts such as `.5`; it does not handle an all-cell `.6`
cut. No old family, old checker or old numerical search was replayed
here. The new point uses a mirrored nonzero skew to enlarge the actual
outer footprints at a controlled quadratic cost.

## 2. Exact new four-center diagram

Let `theta*` be the unique maximizing acute endpoint of the already
accepted full-apex scalar support at the above exact `lambda,t`. Write

```
y=cos(theta*), x=sin(theta*)=sqrt(1-y^2),
c=cos(theta*/2)=sqrt((1+y)/2),
s=sin(theta*/2)=sqrt((1-y)/2),
Q=1+2y, V=1-2y.
```

The accepted enclosure has `169/250<y<677/1000`. The scalar support is

```
G(theta)=sin(theta)+sin(2theta)-4lambda*sin(theta)^2+t*theta,
C=G(theta*)=max_{0<=theta<=pi/2} max(G(theta),Bo(theta)),
Bo(theta)=sin(theta)-4lambda*sin(theta)^2+t*(pi-theta).
```

For an algebraic definition, `y` is the selected root in that interval
of

```
(4y^2+y-2+t)^2=(28/25)^2*y^2*(1-y^2),
4y^2+y-2+t>0.
```

The equation is the exact acute stationary condition. The existing
accepted peak certificate is read as a premise, not recomputed.

Set one fixed skew, with no parameter scan,

```
b=sin(epsilon)=1/12, a=cos(epsilon)=sqrt(143)/12,
k=c*a-s*b, l=s*a+c*b,
T=Q*c*a-V*s*b, H=Q*s*a+V*c*b, M=V*c*a-Q*s*b,
u=sqrt(1-R^2*H^2).
```

Define the actual common apex, centers and unit anchors by

```
p=(u-R*T,0),
Lout=2R(c*a+s*b), Lmid=2R(c*a-s*b),
C0=p+Lout*(y,-x)=(u-R*k,-Lout*x),
C1=p+(Lmid,0)=(u+R*M,0),
C2=(C0_x,-C0_y),
q0=(u,-R*H), q1=(u,R*H).
```

The triangles are `(p,C0,C1)` and `(p,C1,C2)`, with positive fan
orientation. Their actual common ray is `pC1`; each actual apex angle
is exactly `theta*`, and their base-angle skews are respectively
`-epsilon` and `+epsilon`. The positive ray lengths are certified
before any geometric division.

## 3. Exact circumradius, covered unit anchors and loss

Because `a^2+b^2=1` and `x=2sc`, the cosine law gives both base lengths

```
|C1-C0|^2=|C2-C1|^2=4R^2 x^2.
```

For example, with normalized ray lengths
`Lout/R=2(ca+sb)`, `Lmid/R=2(ca-sb)`, their squared base distance is

```
8[c^2 a^2(1-y)+s^2 b^2(1+y)]=4(1-y^2).
```

The sine rule therefore gives circumradius
`(2Rx)/(2sin(theta*))=R` in each triangle.

The anchors are unit points by the definition of `u`. Direct subtraction
gives

```
q0-C0=R(k,l), q0-C1=R(-M,-H),
k^2+l^2=1, M^2+H^2=1.
```

The second identity follows from
`V^2 c^2+Q^2 s^2=1` and rotation by `epsilon`. The reflected equations
hold for `q1`. Thus each actual anchor is covered with equality by the
two endpoint radius-R disks of its base.

For the directed base `C0 -> C1`,

```
det(C1-C0,q0-C0)/R^2=-kH+lM=-2xy<0.
```

The apex is on the positive side of this directed base, while the
anchor is on the opposite side. Hence the actual signed anchor is

```
A0/R^2=A1/R^2=xy,
tau_i=x*sqrt(1-x^2)-A_i/R^2=0.
```

These anchors are the true outer lens tips; no independent or absolute
anchor majorant has replaced them.

The actual triangle area, by the shared-ray determinant, is

```
area_i/R^2=2x(c^2 a^2-s^2 b^2)=x(1+y)-2x b^2.
```

Consequently each exact actual same-coefficient score is

```
area_i/R^2+A_i/R^2-4lambda*x^2+t*theta*
                                      =C-2x b^2.
```

The pair loss is therefore exactly `4x b^2=x/36`, strictly below the
requested fee. The base normalized half-length `x` is also above `.6`,
but that fact is not used as a substitute for the actual assigned-arc
calculation next.

## 4. Three exact contiguous actual assigned arcs

Let `rho=|C0|=|C2|` and `rho1=|C1|`. The new interval certificate
establishes

```
1-R<rho,rho1<1,
rho^2,rho1^2<1-R^2,
0<p_x<1-R.
```

Thus all three boundary disks have proper, nonempty, connected short
footprints on the unit circle. The apex disk misses the whole unit
circle: `|q-p|>=1-|p|>R` for every unit point `q`.

For any center `z` of norm `rho>0`, its unit-circle footprint is exactly

```
q dot z >= (1+rho^2-R^2)/2.
```

It is the cap centered on the ray through `z`, with halfangle `beta`
satisfying

```
cos(beta)=(1+rho^2-R^2)/(2rho)>0.
```

For the lower outer center `C0`, `q0` is an endpoint because its
distance is exactly `R`. Moreover

```
det(C0,q0)=R(l*u+R*k*H)>0.
```

The relevant points and centers lie in the right semicircle. This
positive determinant puts `q0` at the *upper* endpoint of the short
cap of `C0`. Its lower endpoint is the exact reflected unit point

```
w0=2[(C0 dot q0)/rho^2] C0-q0.
```

Reflection preserves unit norm and distance from `C0`. The certificate
also verifies `w0_x>0` and `w0_y<q0_y<0`, so the cap has the intended
counterclockwise order `w0 -> q0`, with no angle-wrap ambiguity.

The center `C1` lies on the positive horizontal ray. Its two exact
distance-R endpoints are `q0,q1`. If
`phi=asin(RH)>0`, its entire short footprint is therefore
`[-phi,+phi]`, counterclockwise from `q0` to `q1`. By reflection the
cap of `C2` runs counterclockwise from `q1` to `w2`, where
`w2=(w0_x,-w0_y)`. Hence choosing these three complete footprints as
assigned cells gives the genuine consecutive arc sequence

```
w0 --[C0]--> q0 --[C1]--> q1 --[C2]--> w2.
```

Their interiors are disjoint. At an interior point of any one assigned
cap, that center is at distance less than `R`, each of the other two
boundary centers is at distance greater than `R`, and the apex is
also farther than `R`. Thus the assigned center is the unique nearest
site among these four local sites. At `q0,q1` the two incident boundary
centers tie at distance `R`, as required by an actual switch. This is
a local nearest-site statement; it asserts nothing about interference
from unspecified sites in a proposed complete covering.

Write the *full* assigned arc angles as `alpha_0,alpha_1,alpha_2`.
Their actual half-sines are

```
xi_0=xi_2=sin(alpha_0/2)/R
           =sin(beta)/R=(l*u+R*k*H)/rho>31/50>.6,
xi_1=sin(alpha_1/2)/R=sin(phi)/R=H>.6.
```

This calculation uses actual footprint endpoints, rather than the
fan's base half-length or its apex half-angle.

The four local centers also have a nondegenerate Delaunay fan: the
first circumcenter is `p+R(k,-l)`, and the reflected one serves the
second triangle. The fourth site's squared distance from the first
circumcenter, minus `R^2`, is

```
R^2[(Lout/R)^2-2(Lout/R)(yk-xl)]>0.
```

Thus the other outer site does not intrude into either marked
circumcircle. This is again only the four-site local check, not a
proof of eligibility in every full covering profile.

## 5. Inner radial cap endpoints and zero BC lens slack

For each complete assigned cap, the norm bound
`rho_i^2<1-R^2` gives

```
cos(beta_i)-rho_i=(1-R^2-rho_i^2)/(2rho_i)>0.
```

The quadratic endpoint-distance equation therefore selects the inner
radial branch exactly:

```
rho_i=cos(beta_i)-sqrt(R^2-sin(beta_i)^2).
```

Here `beta_0=beta_2=beta`, `beta_1=phi`, and the entire arc is centered
on the site's radial ray. Its actual signed BC lens term is

```
L_i=beta_i-rho_i*sin(beta_i)
   =beta_i-sin(beta_i)*cos(beta_i)
                 +sin(beta_i)*sqrt(R^2-sin(beta_i)^2).
```

This equals the usual inner-radial endpoint allowance for that same
actual assigned halfangle. Hence the assigned BC lens slack is zero
at each of the three cells. Separately, `tau_0=tau_1=0` is the actual
edge-anchor slack of Section 3. The pair's entire score deficit is
the skew cost `x/36`; zero lens slack is not a second payment or an
added fee. No global quadratic boundary functional is computed here.

## 6. One new exact check and the limits of the result

`verify_new_skew_point.py` performs exactly one new fixed-point
standalone Fraction/radical enclosure. It reads the already accepted
peak-cosine certificate, verifies its SHA and exact coefficient labels,
and never imports or executes project mathematics. The interval
operations are written here from exact rational endpoint arithmetic;
square roots are enclosed with integer square roots on a fixed rational
grid. It invokes no trigonometry, root isolation, parameter scan,
optimizer, old search, or previous two/three-face checker. All computed
intervals and strict margins are preserved in `certificate.json`.

The perfect zero-loss family was not rerun. Before the unique new
execution, the parent suggested the final fixed skew `sin(epsilon)=1/12`
in place of the initially contemplated `1/16`; neither alternative
was scanned or numerically trialed. Only the final `1/12` point is
executed and certified here.

The result disproves a universal local fee `.02228272803` under these
two-face hypotheses even when every actual assigned arc half-sine is
above `.62`, each center has norm below one in the reflected range,
the cap's radial endpoint allowance is saturated, and the four local
sites do not interfere with the assigned interior arcs or the two
circumcircles. It does not address the unproved global `Q<=.35754`
condition or stronger conditions involving other sites and the complete
boundary. A full covering must separately prove its topology, profiles,
other assignments, total boundary functional and unit-disk coverage.

The later proposed *joint* interface changes the two-face score by also
subtracting `K*(a_mid-pi/33)^2`, where `a_mid=phi` is the actual middle
assigned arc halfangle and `K=(33/20)*(33/32)/R^2`. That changes this
point's deficit to `D+K*(phi-pi/33)^2`. The present certificate evaluates
only the actual signed-A score and does **not** refute the joint target.
It neither evaluates that joint penalty nor proves its global allocation;
any such variance must be allocated once in the complete matching
budget. The separate proposed full-arc minimum `.134` does hold for
these three arcs, since `alpha_i=2asin(R*xi_i)>2R*(3/5)=.135864>.134`.
The proposed global `Q<.35754` and other joint-domain conditions remain
unproved here.
