"""Portable standard-library acceptance entry point.

Print a freshly recomputed JSON report. Exit 0 means acceptance; exit 1 means
rejection. No stored success flag is used. Run: python3 -S stdioverify.py
"""
import json
import sys

if len(sys.argv) != 1:
    print(json.dumps({'verified': False, 'error': 'This acceptance entry point takes no arguments.'}))
    raise SystemExit(1)

try:
    from exact_algebraic_cover import main
    main()
except Exception as exc:
    print(json.dumps({'verified': False, 'error': type(exc).__name__ + ': ' + str(exc)}))
    raise SystemExit(1)
