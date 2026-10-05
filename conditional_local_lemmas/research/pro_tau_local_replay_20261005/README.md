# Two actual anchors: nonzero tau lemma

Run `python verify.py` with Python 3 (standard library only).

Files:
- `PROOF.md`: complete geometric and analytic proof, including all hypotheses,
  phase conventions, full-support localization, exact translation formulas,
  signed lens support, the eight continuous slabs, and the new two-face gap.
- `RESULT_FOR_HANDOFF.md`: concise Chinese interface for the parent proof.
- `verify.py`: exact scalar guard checker, no optimizer or sampling.
- `verification.json`: exact rational bounds from the actual run.
- `verification.log`: actual verification output.

Main theorem, at the fixed parameters in the proof:
if the common middle center satisfies |c1|^2 <= 1-R^2 and both outer centers
satisfy |cj|^2 >= 83/100, then the actual signed-A two-face deficit is >9/400.
No middle-halfarc variance is charged.

This is NOT a global certificate for the exact optimum of 100 disk covering.
