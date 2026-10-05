#!/usr/bin/env python3
"""Verify this release's bytes and syntax without mathematical execution."""
from pathlib import Path
import ast
import hashlib
import json

ROOT = Path(__file__).resolve().parent

def main():
    manifest = json.loads((ROOT / 'MANIFEST_SHA256.json').read_text())
    if not isinstance(manifest, dict) or not manifest:
        raise ValueError('Empty or malformed publication manifest')
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
              if p.is_file() and '.git' not in p.relative_to(ROOT).parts
              and '__pycache__' not in p.parts
              and p != ROOT / 'MANIFEST_SHA256.json'}
    if actual != set(manifest):
        raise ValueError('Publication file set differs from manifest')
    total = 0
    python_files = 0
    for name, digest in sorted(manifest.items()):
        rel = Path(name)
        if rel.is_absolute() or '..' in rel.parts:
            raise ValueError('Unsafe manifest path')
        data = (ROOT / rel).read_bytes()
        if hashlib.sha256(data).hexdigest() != digest:
            raise ValueError('SHA256 mismatch: ' + name)
        total += len(data)
        if name.endswith('.py'):
            ast.parse(data, filename=name)
            python_files += 1
    print(json.dumps({'status': 'PASS_STATIC_INTEGRITY_AND_SYNTAX',
                      'files': len(manifest), 'bytes': total,
                      'python_files': python_files,
                      'mathematical_executions': 0}, sort_keys=True))

if __name__ == '__main__':
    main()
