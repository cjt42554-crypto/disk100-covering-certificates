# Isolated algebraic covering configuration

The exact 932-variable KKT system certifies one isolated local real root and
100-center coverage of the whole closed unit disk. Its radius lies in the
stored exact rational interval near .1154469548855822016397984. This package
does not claim global optimality. Strict local minimality is certified by the
separate strict_local package, modulo simultaneous rotation.

Portable standard-library acceptance:

    python3 -B -S stdioverify.py
    python3 -B -S -O test_negative.py

No NumPy, SciPy, NPZ, optimizer or network is needed on the acceptance path.
Generation branches in the unchanged source can use optional numerical
packages, but they are not called by these acceptance commands. Derived
rational_midpoint and verification JSON may be written inside this package.
Stored reports are historical diagnostics and never acceptance inputs.
LOCAL_PROOF.md contains the root and closed-disk covering reduction.

This publication is a static copy of already accepted sources and data.
No new mathematical replay was run during public copying.
