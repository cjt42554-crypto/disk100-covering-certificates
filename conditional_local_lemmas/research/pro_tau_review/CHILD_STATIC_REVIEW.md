# Independent derivative and remainder review

Reviewer identity: `/root/pro_tau_review/tau_hessian_static`.
Outcome: STATIC PASS. Mathematical executions/imports: zero. Sources unchanged.

The reviewer independently checked Sections 3–6 and the relevant `verify.py` slab reduction. No missing second-order term or mathematical error was found. This report records the reviewer’s completed derivation rather than a new numerical replay.

Put C=cos(alpha), S=sin(alpha), d=(e1−e0)/2, and g=(e0+e1)/2. Directed bases give c0=c1−2Rx0T0 and c2=c1+2Rx1T1. Product-to-sum gives

`−(C,0)−x0 T0+x1 T1=−cos(h)(cos(Delta),sin(Delta))=E*`,

`x0 T0+x1 T1=(−sin(h)sin(Delta), S+sin(h)cos(Delta))=F`.

Thus E=E*−g, the source mean formula (5.5), and the source minimum formula (5.6) follow exactly. At zero slack,

`M*=1+R²−2R(u0 cos(h)−RS sin(h))cos(Delta)`,

which equals (3.3), because u0=cos(a*) and RS=sin(a*).

At b*=(0,S),

`Df(b*)[d]=(−R² S d_y/u0, −u0 d_x/S)`.

The linearization is

`−2R² b*·d+2R Df(b*)[d]·E*−2R f(b*)·g−2R² E*·g`.

Substitution reproduces all four coefficients in (6.1), including both y signs. The tangent signs (6.2) select ell0=−cot(theta0)h0 and ell1=+cot(theta1)h1 in the relaxed lens support. Hence

`kappa0=(−S L0x−C L0y)/x0²`,

`kappa1=(−S L1x+C L1y)/x1²`.

Direct simplification gives exactly (6.3). The full-slab enclosure in the checker is valid. Qlo intentionally drops 1/u0≥1 from its positive second term, giving a safe lower bound.

For the Hessian let r=|b|, e=b/r, t=Jccw e, and w=a e+v t. Differentiation gives

`D²f(b)[w,w]=−Jcw[(R²a²/u³+v²/(u r²))e+2av/(u r²)t]`.

Using `sqrt(v⁴+4a²v²)≤(2/sqrt(3))(a²+v²)` proves (6.7). The path bounds also give the paper bound

`R²/u³+2/(sqrt(3)u r²)<13/970+175/102<7/4`.

With rf=f(b)−f(b*)−Df(b*)d, the complete nonlinear remainder is

`−R²|d|²+R²|g|²+2R[rf·E*−(f(b)−f(b*))·g]`.

Dropping the favorable negative term and applying Taylor's factor one half gives exactly (6.8). Its coefficient obeys the paper bound

`[R²+(7/4)R+(200/83)R](15/13)²<109629/169000<13/20`.

The reviewer did not evaluate the eight numerical slab margins. Their continuous-domain reduction and interval-enclosure directions were checked. The later isolated root replay supplies the numerical guard evidence; this child remains a zero-execution static review.
