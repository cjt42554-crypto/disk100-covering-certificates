# Static clarifications for the frozen zero-tau proof

The proof file SHA remains
`c0077a0fa8ca16bf56c4c13a93571d048082c11ba60f1da4c4f7a6b3ffa549a1`.
This note adds two explicit premises without altering those bytes.

1. The actual base angles of a positive-area triangle are
   `(pi-gamma_i)/2 +- d_i`, both strictly positive. Therefore
   `|d_i|<(pi-gamma_i)/2<pi/2`. This supplies the monotone sine
   branch used when `sin^2 d_i<.028` implies `|d_i|<.18`.
2. The bound `beta<.114` in Section 2 follows from
   `0<R sin alpha<=R=.11322` and the explicit lower Taylor bound

       sin(.114)>=.114-(.114)^3/6>.11322.

   Sine is increasing on this acute interval, so
   `beta=asin(R sin alpha)<.114`.

Writing this clarification involved no mathematical program execution,
import, checker, scan or optimizer. The author's frozen proof and
its original manifest are unchanged. A reviewer may accept the
original proof together with this separately pinned clarification.
