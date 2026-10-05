# Supplemental independent static check

Reviewer: child task `/root/stronger_fee_review/chi_static_check`.
The full frozen current SMALL_TAU_STRONG_ARC.md was read separately.
No mathematical execution, source import, checker, replay, search,
optimizer or file write was performed by this supplemental reviewer.

The reviewer reported no additional substantive gap after the separately
recorded zero-product endpoint correction. Its independent hand checks:

    P10(31/20)-1/50
      =29520213502495199/37158912000000000000>0;
    .945-.945^3/6=.8043485625>.8;
    epsilon<=1/10240=.00009765625<.0001;
    .013*(1.001)^2=.013026013<.013951=1-(.993)^2;
    .013*2.001=.026013<.014*1.986=.027804;
    2.628*(625/512)=3285/1024=3.2080078125<3.21;
    (3285/1024)/12500=.000256640625<.0002568;
    (.9108)^2=.82955664>.8294023684;
    .9108+.0002568=.9110568<.911088875.

The final fixed norm gap is `.000032075`. It separately confirmed the
negative common branch exceeds the actual middle cap by at least
`.0022644`, more than the translation error. The positive common branch
obeys the strict middle cap. A single rigid translation preserves the
actual fan and makes both ideal tips unit; the score changes exactly
by T. The reusable zero-tau norm bound precedes its outer-arc premise,
so no assigned-arc hypothesis is moved to the auxiliary configuration.

The supplementary scope is `D<=.0223 => T>.00008`, hence `D>.00008`
under the same actual-middle-cap and original strong-arc hypotheses.
It does not establish `.022282728...` or any global exclusion. Its
reported literal `<3.21T` qualification is strict only for T>0; the
full review explicitly reads the author's preserved zero-endpoint
clarification and excludes T=0 before that comparison.

These are static hand calculations supplied by the supplemental
reviewer. The primary review also gives a separate coarse rational
cosine guard and did not evaluate the polynomial with a program.
