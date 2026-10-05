# q2_equality_check execution ledger

Owner: child agent `/root/stronger_fee/q2_equality_check`.

This is a recovery record copied from the existing conversation/tool transcript.
Its creation performed file writing only. No command below was rerun to create
this record. No mathematical execution occurred after the parent's instruction
to stop further arithmetic execution.

## Honest execution accounting

- Five mathematical command invocations were attempted.
- One invocation used the unavailable `python` binary and exited 127 before any
  Python mathematical code ran.
- Four invocations used `python3`, imported only Python standard-library modules
  (`math` or `fractions`), and completed with exit 0.
- Each command shown below was attempted exactly once. Command 01 was one new
  exploratory comparison list of seven prescribed angles; it is not accurately
  described as new exploration or scan zero. Command 02 was a subsequent coarse
  fixed-angle bound, command 03 a fixed skew-candidate numerical evaluation, and
  command 04 its fixed exact radical/score guard. No optimizer, historical search,
  historical checker, or historical mathematical source import ran.
- Earlier read-only commands displayed existing proof/certificate text. They are
  not mathematical executions and are not included in the five-command count.
- The commands used the working directory
  `[REDACTED_PRIVATE_LOCATION_0b212b71483f]/task-2` and wrote no files.
- No absolute start/end timestamps were captured by the child. The functions
  result headers rounded the duration of each call to approximately 0.3 seconds.
  The exact shell receipt retained verbatim for command 04 reported
  `wall_time_seconds=0.030096`, `chunk_id=8c0276`, `exit_code=0`.
  Exact shell timing/IDs for commands 00--03 were not separately retained in the
  child's messages; this ledger does not invent them.

## Correction of prior wording

The earlier phrases "static evaluations", "static geometry arithmetic", and
"fresh exact Fraction-only static guards" described fixed-input evaluations,
but could incorrectly suggest paper-only derivation or zero mathematical
execution. The precise corrected statement is: **four new Python mathematical
code executions occurred, including one seven-angle exploratory list**. The
new Fraction guard is an executed guard, not an unexecuted analytical assertion.
"No old
checker/import/optimizer/search" means no historical checker or source import;
standard-library `math` and `fractions` imports did occur. Any interpretation of
the earlier wording as "static-only/no execution" is withdrawn.

## Command 00: absent-python attempt

One attempt, exit 127. The interpreter was unavailable; zero Python calculations
ran. Exact command/source:

```bash
python - <<'PY'
import math
R=5661/50000
lam=.14
t=1071883/20000000
for th in [.8,math.asin(.75),.88,.89,.9,.91,math.acos(.6765)]:
 y=math.cos(th); x=math.sin(th); k=math.cos(th/2); s=math.sin(th/2); Q=1+2*y
 u=math.sqrt(1-R*R*Q*Q*s*s); p=u-R*Q*k
 rho2=1-2*u*R*k+R*R*(2+y-2*y*y)
 rho=math.sqrt(rho2); z=(1+rho2-R*R)/(2*rho); foot=math.sqrt(max(0,1-z*z))/R
 G=x+2*x*y-.56*x*x+t*th
 print('theta',th,'G',G,'two deficit vs C1473500630353',2*(1.473500630352807-G),'rho',rho,'xfoot',foot,'p',p)
rholow=math.sqrt(1-.36*R*R)-.8*R
print('rho threshold',rholow)
PY
```

Combined shell output (stderr; stdout empty):

```text
/bin/bash: line 1: python: command not found
```

## Command 01: fixed equal-angle evaluations

One attempt, exit 0. This is successful mathematical execution 1 of 4 and one
seven-angle exploratory comparison list. The list was prescribed in source;
there was no optimizer, but it was new exploratory comparison work.
Exact command/source:

```bash
python3 - <<'PY'
import math
R=5661/50000
t=1071883/20000000
for th in [.8,math.asin(.75),.88,.89,.9,.91,math.acos(.6765)]:
 y=math.cos(th); x=math.sin(th); k=math.cos(th/2); s=math.sin(th/2); Q=1+2*y
 u=math.sqrt(1-R*R*Q*Q*s*s); p=u-R*Q*k
 rho2=1-2*u*R*k+R*R*(2+y-2*y*y)
 rho=math.sqrt(rho2); z=(1+rho2-R*R)/(2*rho); foot=math.sqrt(max(0,1-z*z))/R
 G=x+2*x*y-.56*x*x+t*th
 print('theta',th,'G',G,'two deficit vs C1473500630353',2*(1.473500630352807-G),'rho',rho,'xfoot',foot,'p',p)
rholow=math.sqrt(1-.36*R*R)-.8*R
print('rho threshold',rholow)
PY
```

Exact stdout:

```text
theta 0.8 G 1.471629147696667 two deficit vs C1473500630353 0.003742965312280244 rho 0.9026203203450096 xfoot 0.5367034658386571 p 0.7448254044686682
theta 0.848062078981481 G 1.472607907919467 two deficit vs C1473500630353 0.001785444866680308 rho 0.9042819022964594 xfoot 0.5613797396029769 p 0.7544239532728848
theta 0.88 G 1.4673945331465812 two deficit vs C1473500630353 0.012212194412451804 rho 0.9054066291237607 xfoot 0.5772004705547136 p 0.7609967878137138
theta 0.89 G 1.464816467384409 two deficit vs C1473500630353 0.017368325936796047 rho 0.9057617713227971 xfoot 0.5820595370041737 p 0.7630846343080592
theta 0.9 G 1.4617926889916142 two deficit vs C1473500630353 0.023415882722385817 rho 0.9061182289913956 xfoot 0.5868736284629122 p 0.7651861246169395
theta 0.91 G 1.4583264982791333 two deficit vs C1473500630353 0.030348264147347592 rho 0.9064759415610689 xfoot 0.5916428214016215 p 0.7673009130137068
theta 0.827796706577295 G 1.4734998335950271 two deficit vs C1473500630353 1.5935155599500206e-06 rho 0.90357644543918 xfoot 0.5511022541594963 p 0.7503322647995635
rho threshold 0.9071139535306547
```

The last angle `acos(.6765)` is merely a fixed nearby value, not the exact peak.
The saved peak certificate read subsequently places the true peak near
`.8283836004345923`.

## Command 02: fresh coarse rational bounds at theta=.9

One attempt, exit 0. This is successful mathematical execution 2 of 4.
Exact command/source:

```bash
python3 - <<'PY'
from fractions import Fraction as F
R=F(5661,50000)
Gup=F(7834,10000)+F(9739,10000)-F(14,25)*F(7833,10000)**2+F(1071883,20000000)*F(9,10)
rho2up=1-2*R*F(993,1000)*F(9,10)+R*R*F(1849,1000)
print('G(.9) rational upper',Gup,float(Gup))
print('2*(1.4735-Gup)',2*(F(14735,10000)-Gup),float(2*(F(14735,10000)-Gup)))
print('rho2 upper',rho2up,float(rho2up),'versus .9064²',float(F(9064,10000)**2),'guard',rho2up < F(9064,10000)**2)
print('threshold lower .9071 guard',1-F(9,25)*R*R > (F(9071,10000)+F(4,5)*R)**2)
# Elementary coarse u bound, using y=cos(.9) between .6216 and .6217
fup=F(22434,10000)**2*F(3784,20000)
print('u>.993 guard',1-R*R*fup > F(993,1000)**2)
PY
```

Exact stdout:

```text
G(.9) rational upper 7309708783/5000000000 1.4619417566
2*(1.4735-Gup) 57791217/2500000000 0.0231164868
rho2 upper 2053331186929/2500000000000 0.8213324747716 versus .9064² 0.82156096 guard True
threshold lower .9071 guard True
u>.993 guard True
```

The source assumes the stated elementary sine/cosine enclosure at `.9`; it did
not independently compute their Taylor enclosure in this command.

## Command 03: fixed skew candidate numerical evaluation

One attempt, exit 0. This is successful mathematical execution 3 of 4.
Exact command/source:

```bash
python3 - <<'PY'
import math
R=5661/50000;y=17/25;x=math.sqrt(336)/25;g=math.acos(y);k=math.sqrt((1+y)/2);l=math.sqrt((1-y)/2);cd=899/901;sd=60/901;Q=1+2*y;W=2*y-1;T=Q*cd*k+W*sd*l;H=Q*cd*l-W*sd*k;u=math.sqrt(1-R*R*H*H);p=u-R*T;a=k*cd-l*sd;b=k*cd+l*sd
c0=(u-R*a,-2*R*b*x);c1=u-R*(W*k*cd+Q*l*sd)
rho=math.hypot(*c0);z=(1+rho*rho-R*R)/(2*rho);xi=math.sqrt(1-z*z)/R
C=1.473500630352807;t=1071883/20000000;G=x+2*x*y-.56*x*x+t*g;D=2*(C-G)+4*x*sd*sd
qminus=(u,-R*H);L=2*R*x
print('g',g,'d',2*math.atan(1/30),'T',T,'H',H,'p',p,'u',u)
print('c0',c0,'rho',rho,'c1',c1,'xi',xi,'max ell',2*math.acos(z),'D',D,'score',2*C-D)
print('dist q c0',math.dist(qminus,c0),'q c1',math.dist(qminus,(c1,0)),'baseL',math.dist(c0,(c1,0)),L)
print('areaR2',a*b*math.sin(g)*2,'isosceles',x*x/math.tan(g/2),'sigma',2*x*sd*sd)
PY
```

Exact stdout:

```text
g 0.823033692134976 d 0.06664199175649439 T 2.167763795215563 H 0.9199326004414993 p 0.7491268843407948 u 0.9945611012351009
c0 (0.8940394455315654, -0.15625241901332373) rho 0.9075909588652281 c1 0.9501701957478659 xi 0.6061274094584276 max ell 0.1373594503712214 D 0.013138477350580838 score 2.9338627833550333
dist q c0 0.11321999999999997 q c1 0.11322 baseL 0.16602855045852807 0.16602855045852807
areaR2 1.225293364914986 isosceles 1.2317963468041297 sigma 0.006502981889144173
```

These are floating numerical observations. They are not an exact global proof,
and do not certify external fan incidence or neighbor windows.

## Command 04: fresh exact Fraction radical/score guard

One attempt, exit 0. This is successful mathematical execution 4 of 4.
Receipt: `chunk_id=8c0276`, `wall_time_seconds=0.030096`, `exit_code=0`.
Exact command/source:

```bash
python3 - <<'PY'
from fractions import Fraction as F
r=F(5661,50000);sl=F(4582575,1000000);sh=F(4582576,1000000)
assert sl*sl<21<sh*sh
alo=(899*sl-120)/4505;ahi=(899*sh-120)/4505
blo=(899*sl+120)/4505;bhi=(899*sh+120)/4505
xl=4*sl/25;xh=4*sh/25
hl=(106082-540*sh)/112625;hh=(106082-540*sl)/112625
ul=F(99456,100000);uh=F(994562,1000000)
assert ul*ul<1-r*r*hh*hh and 1-r*r*hl*hl<uh*uh
rho2lo=(ul-r*ahi)**2+(2*r*blo*xl)**2
rho2hi=(uh-r*alo)**2+(2*r*bhi*xh)**2
assert F(9075,10000)**2<rho2lo and rho2hi<F(9077,10000)**2
# Direct lower alternating Taylor through negative x^14/14! for cos(.823)
z=F(823,1000);fac=1;power=F(1);co=F(1)
for n in range(1,8):
 fac*=(2*n-1)*(2*n);power*=z*z;co+=(-1 if n%2 else 1)*power/fac
assert co>F(17,25)
# Thus gamma=acos(.68)>.823. Static exact score upper from C<1.473501.
t=F(1071883,20000000);Glo=xl*F(59,25)-F(14,25)*F(336,625)+t*F(823,1000)
Dup=2*(F(1473501,1000000)-Glo)+4*xh*F(60,901)**2
assert Dup<F(14,1000)
print('Exact new static radical guards: sqrt21 interval, u interval(.99456,.994562), rho(.9075,.9077), cos(.823)>.68, pair deficit<.014 PASS')
print('deficit rational upper:',Dup,'decimal:',float(Dup))
PY
```

Exact stdout:

```text
Exact new static radical guards: sqrt21 interval, u interval(.99456,.994562), rho(.9075,.9077), cos(.823)>.68, pair deficit<.014 PASS
deficit rational upper: 106697877127691/8118010000000000 decimal: 0.013143353744044538
```

The word "static" inside this historical stdout is retained literally. It does
not negate that this was an actual mathematical code execution. The bound
`C<1.473501` was taken from the read-only saved conjugate certificate, whose
checker was not rerun.

## Local mathematical scope

The independent geometric derivation confirms the parent's proposed two-face
skew coordinates and actual outer-tip anchors. The numerical footprint and exact
coarse rho/score guards address that one local candidate. They do not constitute
a full 100-disk covering, an exact global optimum, a global lower-bound proof,
or a certificate of all assigned arcs and external neighbor-window hypotheses.
