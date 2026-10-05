#!/usr/bin/env python3
"""Exactly one fresh explicit-family verify() call, with frozen-byte checks."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys

if not __debug__:
    raise RuntimeError('Assertions must remain enabled')
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
PINS_RAW = (HERE / 'source_pins.json').read_bytes()
PINS = json.loads(PINS_RAW)

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def snapshot():
    actual = {}
    for rel, expected in PINS.items():
        digest = sha((RESEARCH / rel).read_bytes())
        if digest != expected:
            raise ValueError('Frozen source identity: ' + rel)
        actual[rel] = digest
    return actual

before = snapshot()
env = os.environ.copy()
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    env[key] = '1'
env['PYTHONDONTWRITEBYTECODE'] = '1'
launch = r"""
import json, runpy, sys
from pathlib import Path
if not __debug__:
    raise RuntimeError('Assertions must remain enabled')
sys.dont_write_bytecode = True
p = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(p.parent))
namespace = runpy.run_path(str(p), run_name='independent_explicit_reader')
actual = namespace['verify']()
print(json.dumps(actual, separators=(',', ':')))
"""
target = HERE.parent / 'verify_explicit.py'
process = subprocess.run([sys.executable, '-I', '-B', '-c', launch, str(target)],
                         cwd=HERE, env=env, capture_output=True, text=True)
after = snapshot()
if process.returncode:
    raise RuntimeError(process.stderr + process.stdout)
if process.stderr:
    raise RuntimeError('Unexpected stderr: ' + process.stderr)
actual = json.loads(process.stdout)
saved_path = HERE.parent / 'certificate.json'
saved = json.loads(saved_path.read_text())
if actual != saved:
    raise ValueError('Complete saved mathematical JSON mismatch')
if before != after or (HERE / 'source_pins.json').read_bytes() != PINS_RAW:
    raise ValueError('Source bytes changed during replay')
fresh = json.dumps(actual, indent=2) + '\n'
(HERE / 'fresh_certificate.json').write_text(fresh)
result = {
    'status': 'PASS_UNIQUE_FRESH_EXPLICIT_FAMILY_VERIFY_FULL_JSONMATCH',
    'independent_verify_calls': 1,
    'fresh_isolated_verifier_processes': 1,
    'complete_saved_JSON_match': True,
    'exit_status': process.returncode,
    'stderr': process.stderr,
    'assertions_enabled': True,
    'bytecode_writing_disabled': True,
    'BLAS_threads': 1,
    'target_certificate_sha256': sha(saved_path.read_bytes()),
    'fresh_certificate_sha256': sha(fresh.encode()),
    'source_history_and_dependency_files_checked_before_and_after': len(PINS),
    'source_sha256_before': before,
    'source_sha256_after': after,
    'preceding_root_checker_attempts': 4,
    'preceding_root_plumbing_failures': 2,
    'preceding_root_successful_checks': 2,
    'combined_root_and_this_independent_checker_attempts': 5,
    'conjugate_verify_calls': 0,
    'primitive_verify_calls': 0,
    'latest_tiny_fee_verify_calls': 0,
    'optimizer_calls': 0,
    'catalogue_calls': 0,
    'parameter_campaign_calls': 0,
    'additional_candidate_calls': 0,
    'old_batch_calls': 0,
    'scope': 'One algebraic 37-site local obstruction with complete q3,d7 star, 33 connected covered assigned arcs longer than .092 and full unit-circle coverage. Not a whole-unit-disk or 100-site cover, global profile, outer-Q or global radius counterexample.',
}
(HERE / 'replay_result.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
print(json.dumps(result, indent=2, sort_keys=True))
