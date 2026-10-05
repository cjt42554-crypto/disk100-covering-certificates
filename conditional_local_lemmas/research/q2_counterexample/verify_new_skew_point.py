#!/usr/bin/env python3
"""One new fixed algebraic skew witness; no imported project mathematics.

Reads the already accepted peak-cosine enclosure without rerunning it.
All new arithmetic is standalone exact Fraction interval arithmetic and
integer-square-root enclosures. No trig evaluation, optimizer, scan,
conjugate checker or two/three-face checker is called.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
PEAK_CERT = HERE.parent / 'q3_actual_quadratic_support/certificate.json'
PEAK_PIN = 'c92f9879d82e80ca134eb93379d1e842189f836c5e5603d954361c95b7f286da'
SCALE = 10**28


class Box:
    def __init__(self, lo, hi=None):
        if isinstance(lo, Box) and hi is None:
            self.lo, self.hi = lo.lo, lo.hi
        else:
            self.lo, self.hi = F(lo), F(lo if hi is None else hi)
        if self.lo > self.hi:
            raise ValueError('Reversed rational interval')

    def __add__(self, other):
        other = Box(other)
        return Box(self.lo+other.lo, self.hi+other.hi)
    __radd__ = __add__

    def __neg__(self):
        return Box(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -Box(other)

    def __rsub__(self, other):
        return Box(other) + -self

    def __mul__(self, other):
        other = Box(other)
        values = [self.lo*other.lo, self.lo*other.hi,
                  self.hi*other.lo, self.hi*other.hi]
        return Box(min(values), max(values))
    __rmul__ = __mul__

    def __truediv__(self, other):
        other = Box(other)
        if other.lo <= 0:
            raise ValueError('This witness divides only by positive intervals')
        return self * Box(1/other.hi, 1/other.lo)

    def square(self):
        values = [self.lo*self.lo, self.hi*self.hi]
        return Box(0 if self.lo <= 0 <= self.hi else min(values), max(values))

    def sqrt(self):
        if self.lo < 0:
            raise ValueError('Negative square-root interval')
        lower = isqrt(self.lo.numerator*SCALE*SCALE//self.lo.denominator)
        upper = isqrt(self.hi.numerator*SCALE*SCALE//self.hi.denominator)+1
        lo, hi = F(lower, SCALE), F(upper, SCALE)
        if not (lo*lo <= self.lo and hi*hi > self.hi):
            raise ValueError('Integer square-root enclosure failed')
        return Box(lo, hi)

    def enc(self):
        return [[self.lo.numerator, self.lo.denominator],
                [self.hi.numerator, self.hi.denominator]]


def ratio(q):
    return [q.numerator, q.denominator]


def need(condition, label):
    if not condition:
        raise ValueError(label)


def verify_new_point():
    frozen = PEAK_CERT.read_bytes()
    need(hashlib.sha256(frozen).hexdigest() == PEAK_PIN, 'Accepted peak certificate pin')
    saved = json.loads(frozen)
    need(saved['status'] == 'PASS_SAME_QUADRATIC_ACTUAL_SIGNED_ANCHOR_FAN_FEE_1_40'
         and saved['lambda'] == [7, 50]
         and saved['tilt'] == [1071883, 20000000],
         'Exact accepted scalar-peak coefficients; prior radius metadata is not transferred')
    rows = saved['peak_cosine']
    y = Box(F(*rows[0]), F(*rows[1]))
    need(F(169, 250) < y.lo <= y.hi < F(677, 1000), 'Accepted ideal peak bracket')
    R = F(5661, 50000)
    target = F(2228272803, 100000000000)
    skew_sin = F(1, 12)
    skew_cos = Box(143).sqrt()/12
    x = (1-y.square()).sqrt()
    c = ((1+y)/2).sqrt()
    s = ((1-y)/2).sqrt()
    Q = 1+2*y
    V = 1-2*y
    k = c*skew_cos-s*skew_sin
    l = s*skew_cos+c*skew_sin
    T = Q*c*skew_cos-V*s*skew_sin
    H = Q*s*skew_cos+V*c*skew_sin
    M = V*c*skew_cos-Q*s*skew_sin
    outer_ray = 2*(c*skew_cos+s*skew_sin)
    middle_ray = 2*k
    u = (1-R*R*H.square()).sqrt()
    p = u-R*T
    cx = u-R*k
    cy = -R*outer_ray*x
    middle = u+R*M
    rho2 = cx.square()+cy.square()
    rho = rho2.sqrt()
    xi_outer = (l*u+R*k*H)/rho
    xi_middle = H
    deficit = x/36
    dot = cx*u+cy*(-R*H)
    cap_cos = dot/rho
    reflected_x = 2*dot*cx/rho2-u
    reflected_y = 2*dot*cy/rho2+R*H
    delaunay_gap = outer_ray.square()-2*outer_ray*(y*k-x*l)
    need(x.lo > F(3,5) and k.lo > 0 and l.lo > 0 and outer_ray.lo > 0
         and middle_ray.lo > 0 and H.lo > 0 and M.hi < 0,
         'Positive rays and tips; actual positive-area fans have the stated skew')
    need(F(0) < p.lo <= p.hi < 1-R,
         'Apex disk misses every point of the unit circle')
    need(1-R < rho.lo <= rho.hi < 1 and 1-R < middle.lo <= middle.hi < 1,
         'Three proper nonempty boundary footprints; center norms below one')
    need(rho2.hi < 1-R*R and middle.square().hi < 1-R*R,
         'All three cap centers are the inner radial endpoint, within reflected norm range')
    need(xi_outer.lo > F(31,50) > F(3,5) and xi_middle.lo > F(3,5)
         and xi_outer.hi < 1 and xi_middle.hi < 1,
         'Actual assigned full-cap half-sines are strictly above .62 and .6')
    need(2*R*F(3,5) > F(67,500),
         'All actual full cap angles exceed .134 using asin(z)>z')
    need(deficit.hi < target and deficit.lo > 0,
         'Strict total same-score loss below the requested two-face fee')
    need(cap_cos.lo > rho.hi > 0 and u.lo > middle.hi,
         'Inner radial branch, giving zero assigned BC lens slack')
    need(reflected_x.lo > 0 and reflected_y.hi < (-R*H).lo < 0
         and cx.lo > 0 and cy.hi < 0 and middle.lo > cx.hi > p.hi,
         'Outer reflected endpoint precedes q0; all caps are ordered in the right semicircle')
    need(delaunay_gap.lo > 0,
         'The fourth local center lies strictly outside each marked triangle circumcircle')
    need(PEAK_CERT.read_bytes() == frozen, 'Accepted source unchanged after new arithmetic')
    return {
        'status': 'PASS_NEW_FIXED_Q2_MIRRORED_SKEW_ACTUAL_ARCS_COUNTEREXAMPLE',
        'R': ratio(R), 'lambda': [7,50], 'tilt': [1071883,20000000],
        'skew_sine': ratio(skew_sin), 'skew_cosine': skew_cos.enc(),
        'peak_cosine': y.enc(), 'base_normalized_half_length': x.enc(),
        'half_angle_cosine': c.enc(), 'half_angle_sine': s.enc(),
        'normalized_outer_ray_length': outer_ray.enc(),
        'normalized_middle_ray_length': middle_ray.enc(),
        'T': T.enc(), 'H': H.enc(), 'M': M.enc(),
        'unit_anchor_x': u.enc(), 'unit_anchor_abs_y': (R*H).enc(),
        'apex_x': p.enc(), 'outer_center_x': cx.enc(), 'outer_center_lower_y': cy.enc(),
        'middle_center_radius': middle.enc(), 'outer_center_radius': rho.enc(),
        'actual_assigned_arc_half_sine_outer': xi_outer.enc(),
        'actual_assigned_arc_half_sine_middle': xi_middle.enc(),
        'actual_arc_half_sine_outer_minus_31_50': (xi_outer-F(31,50)).enc(),
        'all_full_assigned_arc_angles_strict_lower': ratio(2*R*F(3,5)),
        'full_arc_threshold_134': ratio(F(67,500)),
        'total_same_score_loss_formula': '4 sin(theta*) (1/12)^2 = sin(theta*)/36',
        'total_same_score_loss': deficit.enc(),
        'requested_fee': ratio(target), 'requested_fee_minus_actual_loss': (target-deficit).enc(),
        'outer_cap_cosine': cap_cos.enc(),
        'outer_cap_cosine_minus_center_radius': (cap_cos-rho).enc(),
        'middle_cap_cosine_minus_center_radius': (u-middle).enc(),
        'outer_left_endpoint_x': reflected_x.enc(), 'outer_left_endpoint_y': reflected_y.enc(),
        'delaunay_other_center_squared_distance_minus_R2_over_R2': delaunay_gap.enc(),
        'assigned_BC_lens_slack': 'exactly zero separately at each full cap inner radial endpoint',
        'actual_edge_anchor_tau': 'exactly zero separately on both marked bases',
        'accepted_peak_certificate_path': str(PEAK_CERT), 'accepted_peak_certificate_sha256': PEAK_PIN,
        'new_fixed_point_calls_this_receipt': 1,
        'old_or_three_face_checker_calls': 0, 'frozen_math_module_imports': 0,
        'trig_calls_or_root_isolation_replays': 0, 'optimizers_or_scans': 0,
        'scope': 'One local four-center diagram and three exact contiguous actual assigned '
                 'unit-circle cap arcs. No full unit-circle/disk cover, complete global Q '
                 'bound, type/profile reduction, or exact optimum. Refutes the signed-A '
                 'fee alone; does not refute the later joint score with a middle-arc '
                 'variance penalty, which is not evaluated here.'
    }


if __name__ == '__main__':
    result = verify_new_point()
    HERE.joinpath('certificate.json').write_text(json.dumps(result, indent=2)+'\n')
    print('PASS: new fixed q2 mirrored-skew witness, sin(epsilon)=1/12.')
    print('Actual assigned outer half-sine > .62; middle > .6; total loss < .02228272803.')
    print('One new standalone Fraction/radical enclosure; old checker/import/search calls 0.')
