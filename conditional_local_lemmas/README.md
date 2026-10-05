# Conditional local geometry, counterexamples and corrections

These standalone local results do not strengthen the universal .11323 bound
or solve the exact global problem. Read MATHEMATICAL_CLOSEOUT.md and each
named proof for the complete hypotheses and once-only budget obligations.

| Result | Essential accepted hypotheses and scope |
| --- | --- |
| Actual signed-A q3 fee >.025 | One actual shared apex with three positive oriented faces, circumradius <=R, covered unit anchors and the exact core parameters. Negative signed anchors and zero-area closures remain included. |
| Strengthened q3 fee >.038 | Same actual geometry and stated strong support, fixed central parameters; the reviewed radius extension is R<=794400/5378209. See q3_stronger_fee_19_500 and q3_radius_uniform_19_500. |
| Wide q3 fee >.0205 | Matched support, lambda in [.126,.154], t within .014 of 1071883/20000000 and R<=794400/5378209. See q3_parameter_broad; the failed 7/388 estimate is preserved and corrected to 35/1944. |
| Full actual-tau q2 fee >.0223 | Fixed R=.11322; the SAME actual middle norm squared <=1-R^2; both original connected outer full arcs >.147, positive actual fan, circumradius <=R and covered anchors. No tau=0 assumption. |
| Pro norm-domain q2 fee >.0225 | Fixed R=.11322; SAME actual middle squared norm <=1-R^2 and each outer squared norm >=.83, with the full actual-fan/anchor hypotheses of its proof. The retained 48 exact guards and eight whole slabs are historical accepted evidence. |
| Mean-arc q2 fee >.0225 | SAME actual middle squared norm <=1.0015; each of two different actual outer arcs is connected. Total arc >=.294 through R=.11364, or >=.3018 through R=.116. Full covered-fan assumptions still apply. The mean argument does not assume an individual outer norm floor. |
| Radius/parameter fee transfer | On the preceding matched actual domain, fee >.0225-.6 abs(Delta lambda)-.1725 abs(Delta t), with abs(Delta lambda)<=.0003 and abs(Delta t)<=.0002. This is an independently reviewed analytic argument, not a fresh executable PASS. |

Portable standard-library entries, from this package, with assertions enabled:

    python3 -B -S research/quadratic_lambda14_conjugate/verify_conjugate.py
    python3 -B -S research/q3_actual_quadratic_support/verify_actual_quadratic.py
    python3 -B -S research/q3_stronger_fee_19_500/verify_actual_quadratic_19_500.py
    python3 -B -S research/q3_parameter_broad/verify_fraction_guards.py
    python3 -B -S research/q3_parameter_uniform/verify_parameter_uniform.py
    python3 -B -S research/q3_parameter_uniform/verify_stronger_corollary_fraction.py
    python3 -B -S research/pro_tau_local_replay_20261005/verify.py
    python3 -B -S research/q3_unequal_explicit_family/verify_explicit.py
    python3 -B -S research/q2_counterexample/verify_new_skew_point.py
    python3 -B -S research/q2_joint_variance/verify_fixed_joint.py

Continuous q2/full actual-tau extensions also include paper proofs and paired
reviews. Analytic admission is not replaced by a saved numeric success flag.
q2_uniform_fresh_guard is an unreviewed, unexecuted, cancelled draft;
q2_cutoff_headroom/.0228 remains unreviewed. Historical failed and superseded
checkers are not accepted entries. The scalar conjugate sin(2a) coefficient
was corrected from -8 to -4; the safe guard is -2sqrt(3)+4lambda<0.

The explicit 37-site J example proves assigned unit-circle coverage and local
premises, not whole-disk coverage or a 100-site counterexample. The weak-arc
and joint-variance q2 counterexamples do not refute the stronger full-arc/norm
domains. Projection preserves covered point sets but can destroy a chosen
fan's circumradius, orientation or ownership. Global fan admission, complete
incidence/profile classification and once-only allocation remain external.

PUBLIC_SOURCE_PATH_MAP.json maps the original Pro six-source names to the
single byte-identical mathematical copy. Existing SHA/source-pin sidecars and
execution reports denote ORIGINAL historical evidence, not runs of a renamed
or adapted public tree. MANIFEST_SHA256 and PUBLIC_COPY_RECEIPT bind the actual
public bytes. Only counterexample PEAK_CERT was changed from a personal path
to the package-relative actual-peak certificate; all mathematical function
ASTs and its exact source-data pin remain unchanged. Metadata redactions never
change rational constants or exact witnesses. No mathematical import, replay,
optimizer or scorer ran during public curation. Original private files remain.
