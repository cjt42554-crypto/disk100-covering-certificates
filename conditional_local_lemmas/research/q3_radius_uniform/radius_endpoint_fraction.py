#!/usr/bin/env python3
"""New endpoint arithmetic only; this is not a replay of the local verifier.

No frozen mathematical source is imported or executed. No scalar/trig fit
guard is rerun. The rational envelopes are premises from the frozen proof.
Prints a machine-readable receipt to stdout and does not modify any file.
"""
from fractions import Fraction as Q
import json


def pair(q):
    return [q.numerator, q.denominator]


def need(condition, message):
    if not condition:
        raise ValueError(message)


K = 2*(Q(19, 25)*Q(49, 50)-Q(27, 100))
S = (Q(763, 1000)+Q(27, 100))**2+(Q(87, 50)+Q(6, 25))**2
RMAX = K/S
R_SIMPLE = Q(19, 100)
R_ORIGINAL = Q(11321, 100000)
need(K == Q(1187, 1250) and S == Q(4987489, 1000000),
     'Exact rational envelopes from frozen equation (18)')
need(RMAX == Q(949600, 4987489) and RMAX*S == K,
     'Exact endpoint, included using the strict geometric envelopes')
need(K-R_SIMPLE*S == Q(197709, 100000000) > 0,
     'The entire requested interval through 19/100 has positive margin')
need(RMAX-R_SIMPLE == Q(197709, 498748900) > 0,
     'Exact endpoint extends past the requested interval')
need(Q(0) < R_ORIGINAL < R_SIMPLE < RMAX,
     'Original radius and simpler interval lie inside the extension')
print(json.dumps({
    'status': 'PASS_NEW_EXACT_RADIUS_ENDPOINT_ARITHMETIC',
    'lambda': pair(Q(7, 50)),
    'tilt': pair(Q(1071883, 20000000)),
    'same_local_normalized_fee': pair(Q(1, 40)),
    'unit_lhs_strict_lower_K': pair(K),
    'unit_sum_norms_squared_strict_upper_S': pair(S),
    'radius_interval': {'lower': pair(Q(0)), 'lower_included': False,
                        'upper': pair(RMAX), 'upper_included': True},
    'K_minus_Rmax_S': pair(K-RMAX*S),
    'requested_simple_upper': pair(R_SIMPLE),
    'K_minus_requested_upper_times_S': pair(K-R_SIMPLE*S),
    'Rmax_minus_requested_upper': pair(RMAX-R_SIMPLE),
    'endpoint_logic': 'K < actual_lhs <= R*actual_sum < R*S <= K; '
                      'strict fit proves the closed endpoint despite zero envelope margin',
    'new_fraction_endpoint_calls_this_receipt': 1,
    'frozen_source_imports_or_executions': 0,
    'local_proof_checker_calls': 0,
    'conjugate_checker_calls': 0,
    'primitive_verify_calls': 0,
    'new_trig_guard_calls': 0,
    'optimizers_or_scans_or_geometry_searches': 0,
    'scope': 'Local same-coefficient actual signed-anchor three-face inequality; '
             'conditional on the unchanged frozen normalized fit and the proved '
             'unitpoint condition. No global radius exclusion or optimum.'
}, indent=2))
