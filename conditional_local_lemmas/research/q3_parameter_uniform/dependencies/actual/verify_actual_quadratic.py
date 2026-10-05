#!/usr/bin/env python3
"""Small exact guards for the same-quadratic actual signed-anchor fan cap.

The accepted scalar-conjugate certificate is read and source-pinned. No old
checker, geometric search, optimizer, or covering classification is run.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import sys

if not __debug__:
    raise RuntimeError('Pinned primitive assertions must remain enabled.')
HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent
BACKEND = RESEARCH / 'coupled_edges/verify_support.py'
PIN = '90fd6e3426d8c113bacedd55f9a13f5e47981e9a48d21bfef48183ec4602efec'
CONJ = RESEARCH / 'quadratic_lambda14_conjugate'
CONJ_PINS = {
    'PROOF.md': 'aaaad3b1629ba98fd44907735b75ae111fecaa45c8d77f078f4e28acf4713ea3',
    'verify_conjugate.py': 'e95db0a6ead5ace03884ab54743b6238f8fe4cd242b1563617b789fce9cb57f4',
    'certificate.json': 'f325d42bfe651599881b7dc652567729a1586389ee385c38060840bc41696bc2',
}
if hashlib.sha256(BACKEND.read_bytes()).hexdigest() != PIN:
    raise ValueError('Frozen interval primitive changed.')
sys.path.insert(0, str(BACKEND.parent))
from verify_support import I, PI, sinq, cosq

R = Q(11321, 100000)
LAM = Q(7, 50)
TILT = Q(1071883, 20000000)
D0 = Q(1, 40)
L, H, ZCUT = Q(7, 10), Q(24, 25), Q(1, 10)


def need(v, message):
    if not v:
        raise ValueError(message)


def qr(x):
    return [x.numerator, x.denominator]


def enclosure(rows):
    return I(Q(*rows[0]), Q(*rows[1]))


def verify():
    frozen = {BACKEND: BACKEND.read_bytes()}
    need(hashlib.sha256(frozen[BACKEND]).hexdigest() == PIN, 'Primitive pin')
    for name, pin in CONJ_PINS.items():
        path = CONJ / name
        frozen[path] = path.read_bytes()
        need(hashlib.sha256(frozen[path]).hexdigest() == pin, 'Accepted conjugate pin: ' + name)
    accepted = json.loads(frozen[CONJ / 'certificate.json'])
    need(accepted['status'] == 'PASS_FIXED_LAMBDA14_FULL_APEX_CONJUGATE'
         and accepted['lambda_fraction'] == qr(LAM)
         and accepted['tilt_fraction'] == qr(TILT)
         and accepted['primitive_sha256'] == PIN,
         'Accepted full-domain conjugate has exactly the present coefficients')
    C = enclosure(accepted['conjugate'])
    Co = enclosure(accepted['obtuse_branch']['value'])
    theta = enclosure(accepted['acute_branch']['peak_angle'])
    y = I(cosq(theta.b).a, cosq(theta.a).b)
    need(Q(169, 250) < y.a <= y.b < Q(677, 1000)
         and L < theta.a <= theta.b < H < PI.a/3,
         'Ideal peak bracket and localized angle window')
    need((C-Co).a > Q(9, 10) > D0,
         'Complete obtuse endpoint is excluded by the same deficit')
    G = lambda a: sinq(a) + sinq(2*a) - 4*LAM*sinq(a)*sinq(a) + TILT*a
    localization = [C-G(L)-D0, C-G(H)-D0]
    need(all(v.a > 0 for v in localization), 'Full high-score acute localization')
    curve = sinq(L) + 4*I(min(sinq(2*L).a, sinq(2*H).a)) + 8*LAM*cosq(2*H)
    need(0 < 2*L < 2*H < PI.a and curve.a > 4,
         'Full-window curvature exceeds four, paying 2 sum e^2')
    ratios = [sinq(a)/sinq((a+ZCUT)/2) for a in (L, H)]
    need(all(v.a > Q(81, 50) for v in ratios),
         'Both endpoints of the ratio with no interior minimum exceed 81/50')
    slope = Q(1, 2)*Q(81, 50)**2-TILT
    need(slope > 1 and D0 < ZCUT and 2*H+ZCUT < PI.a,
         'Convex cut rejects every z>=.1 and pays z with coefficient one')
    S0, C0, S1, C1 = Q(16,25), Q(153,200), Q(41,50), Q(573,1000)
    need(sinq(L).a > S0 and cosq(L).b < C0
         and sinq(H).b < S1 and cosq(H).a > C1,
         'Rational whole-window sine/cosine guards')
    COT, UGAM = Q(6,5), Q(11,4)
    need(C0/S0 < COT and cosq(L/2).b/sinq(L/2).a < UGAM,
         'Circumradius and tip derivative bounds')
    QMAX, QMIN, VMAX = Q(127,50), Q(107,50), Q(11,20)
    need(1+2*C0 < QMAX and 1+2*C1 > QMIN,
         'Whole-window acute tip radius bounds')
    Vlo = 1-2*C0-COT*(1-C1)*D0
    Vhi = 1-2*C1+D0
    need(-VMAX < Vlo < Vhi < VMAX, 'Whole paid-domain tip V bound')
    cs, ct, skew = Q(625,512), Q(625,256), Q(79,100)
    need(cs == 1/(2*S0*S0) and ct == 1/(S0*S0),
         'Exact skew and signed lens-tip coefficient normalization')
    skew_sq = cs*(VMAX*VMAX + QMAX*QMAX*cs*D0)
    need(skew_sq < skew*skew, 'Paid complex-tip skew norm bound')
    jt = S1*C0+QMAX*S1/4
    upos = S1*S1-QMIN*C1/4
    uneg = QMAX*C0/4-S0*S0
    ju = max(upos, uneg)
    jv = 2*S1-QMIN*S0/2
    jw = (QMAX*C0-QMIN*C1)/4
    jvb = (QMAX*S1-QMIN*S0)/4
    jwb = QMAX*C0/2
    need(upos > uneg >= 0 and 2*S0-QMAX*S1/2 > 0,
         'Endpoint imaginary absolute bound includes both signs; middle real remains positive')
    need(min(jt, ju, jv, jw, jvb, jwb) > 0, 'Positive full-box Jacobian bounds')
    JB2 = 2*(jt*jt+ju*ju)+jv*jv+jw*jw
    Jb2 = 2*(jt*jt+ju*ju)+jvb*jvb+jwb*jwb
    alphaB2 = JB2/2+Q(3,2)*skew*skew
    alphab2 = Jb2/2+Q(1,2)*skew*skew
    LB, Lb = Q(201,50), Q(529,200)
    need(LB == UGAM+QMAX/2 and Lb == (UGAM+QMAX)/2
         and LB > ct and Lb > ct/2, 'One-budget linear correction coefficients')
    need(4*LB*LB*D0 < alphaB2 and 4*Lb*Lb*D0 < alphab2,
         'Shared quadratic/linear error maxima occur at zero linear spending')
    EB, Eb = Q(27,100), Q(6,25)
    need(alphaB2*D0 < EB*EB and alphab2*D0 < Eb*Eb,
         'Strict translation-free B and half-chord error bounds')
    yy0, yy1 = Q(169,250), Q(677,1000)
    Amin, Amax, Jmin, Jmax = Q(19,25), Q(763,1000), Q(173,100), Q(87,50)
    need(1+yy1-2*yy1*yy1 > Amin and 1+yy0-2*yy0*yy0 < Amax
         and (1+2*yy0)**2*(1-yy1**2) > Jmin**2
         and (1+2*yy1)**2*(1-yy0**2) < Jmax**2,
         'Uniform ideal sagitta and half-chord bounds')
    axis = Q(49,50)
    need(1-(Eb/(Jmin-Eb))**2 > axis*axis,
         'Actual half-chord perpendicular retains ideal-axis projection >.98')
    lhs = 2*(Amin*axis-EB)
    rhs = R*((Amax+EB)**2+(Jmax+Eb)**2)
    margin = lhs-rhs
    need(lhs == Q(1187,1250) and margin > Q(19,50),
         'Strict actual three-unitpoint contradiction margin exceeds .38')
    need(all(path.read_bytes() == data for path, data in frozen.items()),
         'All imported old sources and accepted certificate remain byte-identical')
    return {
        'status': 'PASS_SAME_QUADRATIC_ACTUAL_SIGNED_ANCHOR_FAN_FEE_1_40',
        'radius': qr(R), 'lambda': qr(LAM), 'tilt': qr(TILT),
        'same_score_fee': qr(D0), 'peak_angle': theta.enc(), 'peak_cosine': y.enc(),
        'accepted_conjugate': C.enc(), 'accepted_obtuse': Co.enc(),
        'localization_margins': [v.enc() for v in localization],
        'curvature_strict_lower_enclosure': curve.enc(),
        'convex_cut': qr(ZCUT), 'cut_ratio_endpoint_enclosures': [v.enc() for v in ratios],
        'paid_gamma_slope_strict_lower': qr(slope),
        'single_paid_budget': '2 sum e^2 + sum z + sum sigma + sum tau <= D <= 1/40',
        'tip_V_bounds': [qr(Vlo), qr(Vhi)], 'tip_skew_squared_upper': qr(skew_sq),
        'tip_skew_coefficient_upper': qr(skew),
        'endpoint_imaginary_two_sided_bounds': [qr(uneg), qr(upos)],
        'B_jacobian_frobenius_squared_upper': qr(JB2),
        'b_jacobian_frobenius_squared_upper': qr(Jb2),
        'B_shared_quadratic_coefficient_squared': qr(alphaB2),
        'b_shared_quadratic_coefficient_squared': qr(alphab2),
        'B_linear_coefficient': qr(LB), 'b_linear_coefficient': qr(Lb),
        'B_error_strict_upper': qr(EB), 'b_error_strict_upper': qr(Eb),
        'unit_lhs_strict_lower': qr(lhs), 'unit_rhs_strict_upper': qr(rhs),
        'unit_compatibility_margin_strict_lower': qr(margin),
        'primitive_sha256': PIN, 'accepted_conjugate_inputs': CONJ_PINS,
        'old_batches_run': 0, 'optimizers_run': 0, 'geometry_searches_run': 0,
        'scope': 'Consecutive actual common-apex three-face fan, circumradii <= R, '
                 'arbitrary common covered unit anchors of either signed side. '
                 'sum(area/R^2 + signedA/R^2 -.56 x^2 + t gamma) <3C-1/40. '
                 'Old BC lens and centered credits are not charged by this theorem; '
                 'global signed-shoelace substitution, fan incidence and other '
                 'profiles remain external.'
    }


if __name__ == '__main__':
    result = verify()
    HERE.joinpath('certificate.json').write_text(json.dumps(result, indent=2)+'\n')
    print('PASS: same-quadratic actual signed-anchor fan fee 1/40=.025.')
    print('One small exact guard calculation; accepted conjugate read only; old batches/searches 0.')
