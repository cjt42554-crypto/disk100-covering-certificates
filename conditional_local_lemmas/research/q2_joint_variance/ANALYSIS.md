# One fixed local obstruction to the q2 actual-anchor plus middle-variance target

At the exact coefficients

    R=5661/50000=.11322,
    lambda=7/50,
    t=1071883/20000000,

the single prescribed point `cos gamma=17/25`,
`tan(skew/2)=1/40` gives

    Dactual=2C-S2 approximately .007455478118006579,
    middle variance approximately .012938880787747476,
    Djoint approximately .020394358905754055,

strictly below

    delta=178523206195158271774904671901
          /8011730250000000000000000000000
          approximately .02228272802808835.

All three actual full covered arcs in this local chain exceed `.134`.
Their middle half-arc, rather than the full arc, is used in the
variance term. The generated rational certificate has strictly
positive delta margin, approximately `.001888369122334294`.

This is a local obstruction to the joint inequality from the listed
geometry and single variance term alone. It does **not** refute a
lemma that additionally assumes the unverified global `Q<.35754`;
it does not certify degree-seven empty-circle incidence, the
four-cell neighbor window, a complete high-M profile, exactly 100
circles, or a full disk cover. The auxiliary global hypotheses must
be verified separately before this configuration could be declared
admitted to that complete theorem.

## Exact definitions and actual score

Use `C=C_lambda(t)` from the already accepted complete full-apex
conjugate certificate at these exact lambda,t. That scalar support
is independent of R. Its checker is read-only and is not called.
Let

    y=17/25,
    x=sin gamma=4sqrt(21)/25,
    gamma=acos(17/25),
    k=cos(gamma/2)=sqrt(21)/5,
    l=sin(gamma/2)=2/5,
    cd=cos d=1599/1601,
    sd=sin d=80/1601,
    Qtheta=1+2y=59/25,
    W=2y-1=9/25.

`Qtheta` is only an acute-tip coefficient. It is not the parent's
surrounding global outer-Q quantity. Define

    aa=k cd-l sd=(1599sqrt(21)-160)/8005,
    bb=k cd+l sd=(1599sqrt(21)+160)/8005,
    Ttip=Qtheta cd k+W sd l
        =(94341sqrt(21)+1440)/200125,
    Htip=Qtheta cd l-W sd k
        =(188682-720sqrt(21))/200125,
    v=W k cd+Qtheta l sd
        =(14391sqrt(21)+9440)/200125,
    u=sqrt(1-R^2 Htip^2).

The algebraic centers and actual unit anchors are

    p=(u-R Ttip,0),
    c0=(u-R aa,-2R bb x),
    c1=(u-R v,0),
    c2=(u-R aa,2R bb x),
    q0=(u,-R Htip), q1=(u,R Htip).

Both faces `(p,c0,c1)` and `(p,c1,c2)` are positive-area and
consecutive, with common actual apex p and actual apex angle gamma.
The left radial side lengths are `2R bb=2R cos(gamma/2-d)` and
`2R aa=2R cos(gamma/2+d)`. The right triangle is its mirror.
The sine rule then gives circumradius exactly R in both faces.
Both normalized bases are `x=sin gamma`.

In the left actual bisector frame, the base midpoint is

    R[(1+y)cd-i(1-y)sd].

The outward base normal is `exp(id)`. Adding the tip height Ry
gives `R[Qtheta cd+iW sd]`. Multiplication by
`exp(-i gamma/2)` gives `R(Ttip-iHtip)=q0-p`.
Reflection gives the right outer tip q1. These are actual covered
outward lens tips with `|q_i-c_i|=|q_i-c_(i+1)|=R`.
Their unit norms follow exactly from u. Thus no formal-J replacement
or anchor-to-lens allowance is involved: the actual signed anchor
slack tau is zero and `A_i/R^2=x y`.

The actual area of either face, divided by R squared, is

    2 aa bb x=x(1+y)-2x sd^2.

Therefore, with actual signed A and actual apex angle,

    S2=sum_i[area_i/R^2+A_i/R^2
              -lambda(base_i/R)^2+t gamma_i]
       =2G(gamma)-4x sd^2,

where

    G(gamma)=x(1+2y)-4lambda x^2+t gamma.

The exact actual-score loss is

    Dactual=2(C-G(gamma))+4x sd^2.            (1)

The actual signed terms are included once. The parent variance
term is not presumed already spent elsewhere.

## Actual assigned full arcs and the middle half-arc

The verified center norms are approximately

    |p|=.7490840734820991,
    |c0|=|c2|=.9065107956046333,
    |c1|=.9518345327131940.

All are strictly less than one, and `|p|+R<1`; p covers no
unit-circle point. Let `rho=|c0|` and define

    z=(1+rho^2-R^2)/(2rho).

The full unit-circle footprint length of the c0 ball is
`ell_outer=2acos z`, and its assigned-arc parameter is

    xi_outer=sqrt(1-z^2)/R approximately .5921042566615620.

This point meets the new full-arc condition `.134`. It does not
meet the earlier stronger `xi>.6` condition, which is a different
domain. The exact cut guard in the certificate is

    cos(.067)-z >0,

with lower enclosure

    59194770402467456912058248611
       /10000000000000000000000000000000000.

Since acos is decreasing, this proves `ell_outer>.134` without
confusing base x and assigned xi.

Let

    w0=2[(q0 dot c0)/rho^2]c0-q0,

and let w2 be its mirror. Reflection through the radial line of c0
preserves unit norm and distance to c0, so w0 is its true second
footprint endpoint. The certificate checks `det(c0,q0)>0`, with
a strict lower enclosure about `.06077071055`; q0 is therefore
the counterclockwise outer-cap endpoint. Because c1 is positive
real, q0 is its clockwise endpoint and q1 its counterclockwise
endpoint. Reflection gives the right cap order. The actual connected
assigned arcs are

    c0: counterclockwise [w0,q0],
    c1: counterclockwise [q0,q1],
    c2: counterclockwise [q1,w2].

They are full ball/unit-circle footprints with disjoint interiors.
At every interior point the assigned center has distance <R, while
the other two have distance >R. Thus these are genuine canonical
nearest arcs among the listed centers, not artificial labels; q0
and q1 are genuine equal-distance switches. This argument makes
no assertion about interference from any unlisted global centers.

The actual middle endpoints are
`q0=(u,-R Htip)`, `q1=(u,R Htip)`, with u>0. Hence their true
middle full arc has length `2a_mid`, where

    a_mid=asin(R Htip)
          approximately .1050727413810770042.      (2)

This is the positive connected arc across the positive real axis.
Its full length exceeds `.134` by about `.07614548`. Formula (2)
is the **half** arc used below, not the full length or the apex angle.

## The once-only middle variance payment

The requested modified score is

    S2_modified=S2-K(a_mid-pi/33)^2,
    K=[(33/20)(33/32)]/R^2=472656250/3560769.

Therefore its loss is exactly

    Djoint=2C-S2_modified
          =Dactual+K(a_mid-pi/33)^2.              (3)

The half-arc deviation is positive, approximately
`.009872963999568118`. The unique author guard evaluates (1)--(3)
with outward rational intervals, using the exact actual middle
half-arc (2). Its certificate gives

    Dactual in an interval near .007455478118006579,
    K(a_mid-pi/33)^2 near .012938880787747476,
    Djoint near .020394358905754055.

More precisely, its joint upper bound is

    127464743160962844484011818069449
        /6250000000000000000000000000000000.

The strictly positive lower bound on `delta-Djoint` is

    94418456116714671779720566799411
        /50000000000000000000000000000000000.

Thus the exact joint target is obstructed on this explicitly
described local domain. The variance has been used once in (3);
it must not also be subtracted as an independent local credit in
another row or a surrounding global Cauchy budget.

This diagnosis leaves open whether the extra surrounding global
outer-Q condition, degree-seven empty fan, four-cell window or
full-profile incidence excludes this configuration. They cannot
be inferred merely from the low local loss or the actual cap chain.

## Execution, replay and audit scope

Only the first prespecified point `tan(d/2)=1/40` was evaluated.
It passed both the actual full-arc test and the joint-loss test,
so the fallback prespecified `1/39` point was not run. There was
one author mathematical invocation, one fixed point and one PASS;
no author failure, interval sweep, optimizer or old checker call.

`verify_fixed_joint.py` imports only the source-pinned interval,
Taylor and inverse-sine primitive definitions from the bundled
accepted primitive. It reads the bundled accepted conjugate
certificate, verifies exact coefficients and four dependency SHA
pins, and does not call either old mathematical `verify()` function.
Its source snapshot, complete output, timestamps and command are
retained in `attempts/author_001_denom40/`. This accounting is
separate from the earlier helper's four q2 calls and seven-angle
exploratory list in `../q2_actual_anchor/helper_receipts/`.

For an independent replay of this chosen point only:

    python3 verify_fixed_joint.py 40

The current proof is awaiting independent full-source static review.
A replay is an additional mathematical execution and must be
recorded; it is not requested by this note. The new author guard
and certificate do not certify the external global obligations.
