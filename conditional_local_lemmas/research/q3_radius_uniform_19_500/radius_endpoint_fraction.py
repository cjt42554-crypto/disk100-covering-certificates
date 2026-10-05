#!/usr/bin/env python3
"""New exact endpoint arithmetic for the .038 local fit; no verifier replay.

The stronger theorem's normalized rational envelopes are premises. Full
independent acceptance of that theorem is recorded separately. No frozen
source is imported or executed, and no trig or fit guard is recalculated.
This prints JSON to stdout and modifies no file.
"""
from fractions import Fraction as Q
import json


def pair(q):
    return [q.numerator, q.denominator]


def need(condition, message):
    if not condition:
        raise ValueError(message)


K = 2*(Q(19, 25)*Q(97, 100)-Q(17, 50))
S = (Q(763, 1000)+Q(17, 50))**2+(Q(87, 50)+Q(3, 10))**2
RMAX = K/S
R_SIMPLE = Q(231, 2000)
R_ORIGINAL = Q(11321, 100000)
FEE = Q(19, 500)
need(K == Q(993, 1250) and S == Q(5378209, 1000000),
     'Exact rational envelopes from stronger equation (18)')
need(RMAX == Q(794400, 5378209) and RMAX*S == K,
     'Exact endpoint; strict geometric fits include the endpoint')
need(K-R_SIMPLE*S == Q(346433721, 2000000000) > 0,
     'Requested .1155 upper radius has positive coarse-envelope margin')
need(RMAX-R_SIMPLE == Q(346433721, 10756418000) > 0,
     'Exact endpoint extends past .1155')
need(Q(0) < R_ORIGINAL < R_SIMPLE < RMAX,
     'Original and requested radii are inside the extension')
need(FEE-Q(3, 100) == Q(1, 125) > 0,
     'The normalized fee strictly exceeds .03')
print(json.dumps({
    'status': 'PASS_NEW_EXACT_RADIUS_ENDPOINT_ARITHMETIC_19_500',
    'premise': 'Frozen stronger normalized fit; independent theorem '
               'acceptance is tracked separately in PREMISE_STATUS.json',
    'lambda': pair(Q(7, 50)), 'tilt': pair(Q(1071883, 20000000)),
    'same_local_normalized_fee': pair(FEE),
    'fee_minus_3_100': pair(FEE-Q(3, 100)),
    'unit_lhs_strict_lower_K': pair(K),
    'unit_sum_norms_squared_strict_upper_S': pair(S),
    'radius_interval': {'lower': pair(Q(0)), 'lower_included': False,
                        'upper': pair(RMAX), 'upper_included': True},
    'K_minus_Rmax_S': pair(K-RMAX*S),
    'requested_simple_upper': pair(R_SIMPLE),
    'K_minus_requested_upper_times_S': pair(K-R_SIMPLE*S),
    'Rmax_minus_requested_upper': pair(RMAX-R_SIMPLE),
    'endpoint_logic': 'K < actual_lhs <= R*actual_sum < R*S <= K; '
                      'strict fit includes the closed upper endpoint',
    'new_fraction_endpoint_calls_this_receipt': 1,
    'frozen_source_imports_or_executions': 0,
    'stronger_or_earlier_or_conjugate_checker_calls': 0,
    'primitive_verify_calls': 0, 'new_trig_or_fit_guard_calls': 0,
    'optimizers_or_scans_or_geometry_searches': 0,
    'scope': 'Local actual signed-anchor consecutive three-face inequality '
             'at the exact fixed coefficient pair. No coefficient rectangle, '
             'global radius exclusion or optimum.'
}, indent=2))
