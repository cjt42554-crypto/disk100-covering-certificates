#!/usr/bin/env python3
"""Exact diagnostics of separately chosen mirrored-skew q2 points.

Each invocation evaluates one explicitly recorded half-angle denominator.
There is no sweep, optimizer or old checker invocation. The surrounding
global outer-Q condition is not evaluated by this local diagnostic.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys

if not __debug__:
    raise RuntimeError('Assertions must remain enabled')
HERE = Path(__file__).resolve().parent
DEP = HERE / 'dependencies'
PRIMITIVE = DEP / 'coupled_edges/verify_support.py'
PINS = {
    'coupled_edges/verify_support.py': '90fd6e3426d8c113bacedd55f9a13f5e47981e9a48d21bfef48183ec4602efec',
    'quadratic_lambda14_conjugate/PROOF.md': 'aaaad3b1629ba98fd44907735b75ae111fecaa45c8d77f078f4e28acf4713ea3',
    'quadratic_lambda14_conjugate/verify_conjugate.py': 'e95db0a6ead5ace03884ab54743b6238f8fe4cd242b1563617b789fce9cb57f4',
    'quadratic_lambda14_conjugate/certificate.json': 'f325d42bfe651599881b7dc652567729a1586389ee385c38060840bc41696bc2',
}
FROZEN = {name: (DEP/name).read_bytes() for name in PINS}
for name, pin in PINS.items():
    if hashlib.sha256(FROZEN[name]).hexdigest() != pin:
        raise ValueError('Accepted dependency pin changed: '+name)
sys.dont_write_bytecode = True
sys.path.insert(0, str(PRIMITIVE.parent))
from verify_support import I, PI, asin, cosq

R = F(5661, 50000)
LAM = F(7, 50)
T = F(1071883, 20000000)
DELTA = F(178523206195158271774904671901,
          8011730250000000000000000000000)
FULLARC_CUT = F(67, 500)


def need(ok, msg):
    if not ok:
        raise ValueError(msg)


def rat(q):
    return [q.numerator, q.denominator]


def contains_zero(v):
    return v.a <= 0 <= v.b


def verify(m):
    need(m in (40,39), 'Only the explicitly prespecified points 1/40, 1/39 are admitted')
    accepted = json.loads(FROZEN['quadratic_lambda14_conjugate/certificate.json'])
    need(accepted['status'] == 'PASS_FIXED_LAMBDA14_FULL_APEX_CONJUGATE'
         and accepted['lambda_fraction'] == rat(LAM)
         and accepted['tilt_fraction'] == rat(T), 'Exact scalar coefficients match')
    C = I(F(*accepted['conjugate'][0]), F(*accepted['conjugate'][1]))
    root21 = I(21).sqrt()
    y, x = F(17,25), 4*root21/25
    k, l = root21/5, F(2,5)
    cd, sd = F(m*m-1,m*m+1), F(2*m,m*m+1)
    need(cd*cd+sd*sd == 1, 'Exact skew Pythagorean identity')
    Q, W = 1+2*y, 2*y-1
    aa, bb = k*cd-l*sd, k*cd+l*sd
    tipT, tipH = Q*cd*k+W*sd*l, Q*cd*l-W*sd*k
    u = (1-R*R*tipH*tipH).sqrt()
    p = u-R*tipT
    c0x, c0y = u-R*aa, -2*R*bb*x
    c1 = u-R*(W*k*cd+Q*l*sd)
    rho2 = c0x*c0x+c0y*c0y
    rho = rho2.sqrt()
    q0x, q0y = u, -R*tipH
    qnorm = q0x*q0x+q0y*q0y-1
    q_c0 = (q0x-c0x)*(q0x-c0x)+(q0y-c0y)*(q0y-c0y)-R*R
    q_c1 = (q0x-c1)*(q0x-c1)+q0y*q0y-R*R
    need(all(contains_zero(v) for v in (qnorm,q_c0,q_c1)),
         'Interval consequences of symbolic exact unit/tip identities')
    need(0 < p.a <= p.b < 1-R and 0 < c1.a <= c1.b < 1
         and 0 < rho.a <= rho.b < 1,
         'All centers are inside unit disk; p covers no unit point')
    det0 = c0x*q0y-c0y*q0x
    need(det0.a > 0 and tipH.a > 0 and c1.a > 0,
         'q0 is outer-cap upper endpoint and middle-cap lower endpoint')
    z = (1+rho2-R*R)/(2*rho)
    need(0 < z.a <= z.b < 1, 'Nondegenerate outward covered cap')
    xi = (1-z*z).sqrt()/R
    outer_cut_margin = cosq(FULLARC_CUT/2)-z
    amid = asin(R*tipH)
    gamma = asin(x)
    G = x*Q-4*LAM*F(336,625)+T*gamma
    actual_loss = 2*(C-G)+4*x*sd*sd
    K = F(33,20)*F(33,32)/(R*R)
    halfarc_deviation = amid-PI/33
    need(halfarc_deviation.a > 0, 'Middle half-arc is above the comparison mean')
    variance = K*halfarc_deviation*halfarc_deviation
    joint_loss = actual_loss+variance
    delta_margin = DELTA-joint_loss
    middle_cut_margin = 2*amid-FULLARC_CUT
    arc_pass = outer_cut_margin.a > 0 and middle_cut_margin.a > 0
    joint_below = delta_margin.a > 0
    need(all((DEP/name).read_bytes()==data for name,data in FROZEN.items()),
         'Accepted dependencies remain byte-identical')
    return {
        'status': ('PASS_FIXED_POINT_LOCAL_JOINT_OBSTRUCTION'
                   if arc_pass and joint_below else 'PASS_FIXED_POINT_DIAGNOSTIC_ONLY'),
        'radius':rat(R), 'lambda':rat(LAM), 'tilt':rat(T),
        'cos_gamma':rat(y), 'sin_gamma':x.enc(), 'gamma':gamma.enc(),
        'tan_skew_half':rat(F(1,m)), 'cos_skew':rat(cd), 'sin_skew':rat(sd),
        'accepted_C':C.enc(), 'tip_T':tipT.enc(), 'tip_H':tipH.enc(),
        'p_x':p.enc(), 'outer_center_x':c0x.enc(), 'outer_center_y':c0y.enc(),
        'middle_center_x':c1.enc(), 'outer_rho':rho.enc(),
        'outer_arc_half_sine_xi':xi.enc(),
        'full_assigned_arc_cut':rat(FULLARC_CUT),
        'outer_fullarc_cut_cosine_margin':outer_cut_margin.enc(),
        'middle_actual_halfarc':amid.enc(), 'middle_fullarc_cut_margin':middle_cut_margin.enc(),
        'q0_endpoint_orientation_determinant':det0.enc(),
        'actual_score_loss_2C_minus_S2':actual_loss.enc(),
        'middle_variance_K':rat(K), 'middle_halfarc_deviation':halfarc_deviation.enc(),
        'middle_variance_loss':variance.enc(),
        'joint_loss':joint_loss.enc(), 'requested_delta':rat(DELTA),
        'requested_delta_minus_joint_loss':delta_margin.enc(),
        'all_three_full_arcs_above_cut':arc_pass,
        'joint_loss_strictly_below_requested_delta':joint_below,
        'dependency_sha256':PINS,
        'fixed_points_evaluated_this_invocation':1,
        'old_checker_calls':0, 'old_source_verify_calls':0, 'optimizer_calls':0,
        'global_outer_Q_verified':False,
        'scope':'One fixed algebraic local two-face configuration. Actual signed anchors '
                'and actual middle half-arc; joint variance paid here once. Three full '
                'covered cap arcs are consecutive and nearest among listed centers. '
                'Global Q<.35754, complete 100-cover, degree-seven empty fan, '
                'four-cell window and whole-profile incidence are not verified.'
    }


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Supply one explicitly prespecified denominator: 40 or 39')
    denominator = int(sys.argv[1])
    result = verify(denominator)
    (HERE/f'certificate_denom{denominator}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
    for key in ('outer_rho','outer_arc_half_sine_xi','middle_actual_halfarc',
                'actual_score_loss_2C_minus_S2','middle_variance_loss','joint_loss',
                'requested_delta_minus_joint_loss'):
        row=result[key]
        print(key, '['+str(float(F(*row[0])))+', '+str(float(F(*row[1])))+']')
    print('all three actual full arcs >.134:',result['all_three_full_arcs_above_cut'])
    print('global outer Q verified:',result['global_outer_Q_verified'])
