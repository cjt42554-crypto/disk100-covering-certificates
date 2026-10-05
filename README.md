# Covering the unit disk with 100 equal disks

This repository preserves rigorous certificates, explicit constructions and selected failed routes for covering the closed unit disk with 100 congruent disks. The exact global optimum remains open.

Let

$$
r_D(100)=\inf_{c_1,\ldots,c_{100}\in\mathbb R^2}
\max_{x\in D}\min_i\lVert x-c_i\rVert,
\qquad D=\{x:\lVert x\rVert\le1\}.
$$

The accepted unconditional interval is

$$
\boxed{\frac{11323}{100000}<r_D(100)\le
\frac{1154469548855829}{10^{16}}}
$$

or **`.11323 < r_D(100) <= .1154469548855829`**. Its gap is `.0022169548855829`. No matching global lower bound, optimal-radius polynomial or globally optimal configuration is claimed.

| Material | Entry point | Scope |
| --- | --- | --- |
| Final global bounds | [Certificate package](certified_bounds_11323/README.md) | Unrestricted strict lower bound and explicit rational upper construction; complete acceptance inputs and code |
| Independent upper margin | [Exact margin proof](strict_upper_margin_1154469/PROOF.txt) | Full closed-disk coverage at `.1154469548855829` and its exact comparison with `.115446954885583` |
| Algebraic construction | [Algebraic package](algebraic/README.md) | Certified algebraic root and whole-disk cover for the stated configuration |
| Strict local optimum | [Local certificate](strict_local/README.md) | Strict local optimality for the true covering objective modulo simultaneous rotations; no global exclusion |
| Fixed rational configuration | [Exact fixed-radius package](fixed_radius/README.txt) | Fixed-configuration radius polynomial and exact proof that this configuration is not globally optimal |
| Conditional geometry | [Local lemmas and corrections](conditional_local_lemmas/README.md) | The stated actual-fan, norm and parameter hypotheses remain necessary |
| New guarded endpoints and finite comparison | [Constructive certificates](constructive/README.md) | Exact endpoint coverage, centerwise transport and finite saved-reference Hall exclusion |
| Continuous event path | [Rational piecewise-linear path](certified_event_path/README.txt) | 23 rational nodes, 22 continuous segments and a certified Delaunay event; no record improvement |
| Earlier unsuccessful methods | [Historical method records](historical_methods/README.md) | Selected mathematical evidence and bounded-run outcomes, with omitted discovery data identified |

Read [FAILURE_ROUTES.md](FAILURE_ROUTES.md) before interpreting a local theorem, conditional fee, fixed-configuration polynomial or finite search as a global statement.

## Verification

Python 3 is required. Integrity and Python syntax can be checked without running any mathematical certificate:

```sh
python3 -B verify_integrity.py
```

The complete global certificate package supplies its own run instructions. Its exact lower-bound chain can take several minutes. Use the documented flags: several checkers rely on assertions and must run with assertions enabled. The independent integer/surd upper-margin checker explicitly supports `-O`.

```sh
python3 -S -B certified_bounds_11323/run_certificates.py
python3 -S -B strict_upper_margin_1154469/verify_margin.py
python3 -S -B certified_event_path/verify_path.py
```

Package READMEs give additional acceptance commands and dependency details. Historical discovery scripts and numerical summaries are distinguished from runnable exact acceptance code.

## Publication identity and evidence

This is a fresh repository containing clean copies, with no inherited Git history. Private conversation links, session identifiers, browser receipts, personal contacts and local user paths are excluded or removed. Original evidence and backups were preserved separately.

`MANIFEST_SHA256.json` binds this public release. Package correspondence maps identify original and public file hashes, omitted files and transformations using project-relative names. Mathematical inputs and acceptance predicates are preserved; any portability adaptation is stated by its package.

Saved mathematical PASS results record earlier executions. This publication phase performed static integrity, syntax, dependency and privacy checks; it made **zero new mathematical replays, optimizations or searches**. A transformed public receipt has a new public hash and does not pretend to be the original byte identity or a fresh execution.

Research was stopped on 2026-10-05. The original search budget had 27 of 96 completed runs and 69 unstarted runs. Publication does not authorize restarting them.
