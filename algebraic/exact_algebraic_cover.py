"""Rational interval proof of whole-disk coverage by the certified KKT root.

The finite critical-point reduction is stated in LOCAL_PROOF.md. This verifier
checks a complete strict Delaunay triangulation, every internal Voronoi vertex,
and a cyclic chain of boundary arcs covered by their supporting centers.
All acceptance tests use Fraction endpoints, with no floating comparisons.
--generate uses SciPy solely to propose triangles/hull/order, fully rechecked.
"""
import os
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import json, sys, time, math, hashlib
from fractions import Fraction as Q
from pathlib import Path
ROOT = Path(__file__).resolve().parent

class I:
    __slots__ = ('lo', 'hi')

    def __init__(self, a, b=None):
        self.lo = Q(a)
        self.hi = Q(a if b is None else b)
        if not self.lo <= self.hi:
            raise ValueError('Certificate acceptance condition failed at line 19')

    def __add__(self, b):
        if not isinstance(b, I):
            b = I(b)
        return I(self.lo + b.lo, self.hi + b.hi)
    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, b):
        return self + -b if isinstance(b, I) else self + I(-b)

    def __rsub__(self, b):
        return I(b) + -self

    def __mul__(self, b):
        if not isinstance(b, I):
            b = I(b)
        v = [self.lo * b.lo, self.lo * b.hi, self.hi * b.lo, self.hi * b.hi]
        return I(min(v), max(v))
    __rmul__ = __mul__

    def __truediv__(self, b):
        if not isinstance(b, I):
            b = I(b)
        if not not b.lo <= 0 <= b.hi:
            raise ValueError('Certificate acceptance condition failed at line 33')
        return self * I(1 / b.hi, 1 / b.lo)

    def sq(self):
        if self.lo <= 0 <= self.hi:
            return I(0, max(self.lo * self.lo, self.hi * self.hi))
        return I(min(self.lo * self.lo, self.hi * self.hi), max(self.lo * self.lo, self.hi * self.hi))

def sub(a, b):
    return [a[0] - b[0], a[1] - b[1]]

def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]

def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]

def norm(a):
    return a[0].sq() + a[1].sq()

def dist(a, b):
    return norm(sub(a, b))

def circumcenter(a, b, c):
    u = sub(b, a)
    v = sub(c, a)
    det = cross(u, v)
    uu = norm(u)
    vv = norm(v)
    rel = [(uu * v[1] - vv * u[1]) / (2 * det), (u[0] * vv - v[0] * uu) / (2 * det)]
    return ([a[0] + rel[0], a[1] + rel[1]], norm(rel), det)

def decimal(q):
    from decimal import Decimal, localcontext, ROUND_FLOOR
    with localcontext() as c:
        c.prec = 75
        c.rounding = ROUND_FLOOR
        return str(Decimal(q.numerator) / Decimal(q.denominator))

def quadrant(q):
    if not (q[0].lo * q[0].hi > 0 and q[1].lo * q[1].hi > 0):
        raise ValueError('Certificate acceptance condition failed at line 52')
    x = q[0].lo > 0
    y = q[1].lo > 0
    return 0 if x and y else 1 if y else 2 if not x else 3

def generate(system, mid):
    import numpy as np
    from scipy.spatial import Delaunay, ConvexHull
    nc = system['n_centers']
    ncoord = 2 * nc
    xx = np.array([float(Q(int(v), mid['denominator'])) for v in mid['numerators']])
    C = xx[:ncoord].reshape(nc, 2)
    dt = Delaunay(C)
    tri = dt.simplices.copy()
    for row in tri:
        if np.cross(C[row[1]] - C[row[0]], C[row[2]] - C[row[0]]) < 0:
            row[1], row[2] = (row[2], row[1])
    nint = len(system['internal_active_triangles'])
    nb = len(system['boundary_active_pairs'])
    q = xx[ncoord + 2 * nint:ncoord + 2 * (nint + nb)].reshape(nb, 2)
    order = np.argsort(np.mod(np.arctan2(q[:, 1], q[:, 0]), 2 * np.pi))
    proposal = {'triangles': tri.tolist(), 'hull': ConvexHull(C).vertices.tolist(), 'boundary_order': order.tolist()}
    (ROOT / 'algebraic_cover_topology.json').write_text(json.dumps(proposal, indent=2))

def check_system_semantics(system):
    """Bind the root equations to the geometrical claims used by this proof."""
    tris = system['internal_active_triangles']
    pairs = system['boundary_active_pairs']
    nc = system['n_centers']
    nt = len(tris)
    nb = len(pairs)
    ncoord = 2 * nc
    if not nc == 100:
        raise ValueError('Certificate acceptance condition failed at line 73')
    nbase = ncoord + 2 * (nt + nb) + 1
    tidx = nbase - 1
    used = set((ci for row in tris + pairs for ci in row))
    unused = sorted(set(range(nc)) - used)
    ncons = 3 * nt + 3 * nb + 1 + 2 * len(unused)
    if not (system['n_base'] == nbase and system['t_index'] == tidx):
        raise ValueError('Certificate acceptance condition failed at line 78')
    if not (system['n_constraints'] == ncons and system['n_variables'] == nbase + ncons):
        raise ValueError('Certificate acceptance condition failed at line 79')
    expected = []

    def distance(ci, hi):
        if not (type(ci) is int and 0 <= ci < nc):
            raise ValueError('Certificate acceptance condition failed at line 82')
        terms = []
        for dim in (0, 1):
            xx = 2 * ci + dim
            yy = ncoord + 2 * hi + dim
            terms.extend([[1, [xx, xx]], [1, [yy, yy]], [-2, [xx, yy]]])
        return terms + [[-1, [tidx]]]
    for hi, tri in enumerate(tris):
        if not (len(tri) == 3 and len(set(tri)) == 3):
            raise ValueError('Certificate acceptance condition failed at line 89')
        for ci in tri:
            expected.append(distance(ci, hi))
    for k, pair in enumerate(pairs):
        if not (len(pair) == 2 and len(set(pair)) == 2):
            raise ValueError('Certificate acceptance condition failed at line 92')
        hi = nt + k
        for ci in pair:
            expected.append(distance(ci, hi))
        qi = ncoord + 2 * hi
        expected.append([[1, [qi, qi]], [1, [qi + 1, qi + 1]], [-1, []]])
    source = ROOT / 'source_configuration.json'
    raw = source.read_bytes()
    if not hashlib.sha256(raw).hexdigest() == system['source_sha256']:
        raise ValueError('Certificate acceptance condition failed at line 97')
    cfg = json.loads(raw)
    gc = system['gauge_center']
    gd = cfg['coordinate_denominator']
    if not (type(gc) is int and gc in used and (len(cfg['centers_numerators']) == nc)):
        raise ValueError('Certificate acceptance condition failed at line 99')
    gnum = cfg['centers_numerators'][gc][1]
    expected.append([[1, [2 * gc + 1]], [[-gnum, gd], []]])
    anchors = []
    for ci in unused:
        for dim in (0, 1):
            anum = cfg['centers_numerators'][ci][dim]
            anchors.append({'center': ci, 'dimension': dim, 'value': [anum, gd]})
            expected.append([[1, [2 * ci + dim]], [[-anum, gd], []]])
    if not system['unused_center_anchors'] == anchors:
        raise ValueError('Certificate acceptance condition failed at line 108')
    station = [[] for _ in range(nbase)]
    station[tidx].append([1, []])
    for i, terms in enumerate(expected):
        lm = nbase + i
        for coeff, vs in terms:
            for k, col in enumerate(vs):
                station[col].append([coeff, [lm] + vs[:k] + vs[k + 1:]])
    expected += station
    if not expected == system['polynomials']:
        raise ValueError('System must have the exact advertised geometry and KKT derivatives')

def main():
    started = time.time()
    system = json.loads((ROOT / 'polynomial_system.json').read_text())
    check_system_semantics(system)
    from exact_krawczyk import verify
    rootcert = verify(generate_requested=False, emit=False)
    if not rootcert['verified']:
        raise ValueError('Certificate acceptance condition failed at line 122')
    mid = json.loads((ROOT / 'rational_midpoint.json').read_text())
    rho = Q(*rootcert['box_radius'])
    D = mid['denominator']
    x = [I(Q(int(v), D) - rho, Q(int(v), D) + rho) for v in mid['numerators']]
    nc = system['n_centers']
    ncoord = 2 * nc
    nbase = system['n_base']
    nt = len(system['internal_active_triangles'])
    nb = len(system['boundary_active_pairs'])
    C = [x[2 * i:2 * i + 2] for i in range(nc)]
    t = x[system['t_index']]
    if not (t.lo > 0 and C[system['gauge_center']][0].lo * C[system['gauge_center']][0].hi > 0):
        raise ValueError('Certificate acceptance condition failed at line 128')
    distance_alpha = [x[nbase + 3 * k + l] for k in range(nt) for l in range(3)]
    distance_alpha += [x[nbase + 3 * nt + 3 * k + l] for k in range(nb) for l in range(2)]
    if not all((a.lo > 0 for a in distance_alpha)):
        raise ValueError('Certificate acceptance condition failed at line 131')
    lambdas = [sum((x[nbase + 3 * k + l] for l in range(3)), I(0)) for k in range(nt)]
    lambdas += [sum((x[nbase + 3 * nt + 3 * k + l] for l in range(2)), I(0)) for k in range(nb)]
    if not all((a.lo > 0 for a in lambdas)):
        raise ValueError('Certificate acceptance condition failed at line 134')
    if '--generate' in sys.argv:
        generate(system, mid)
    prop = json.loads((ROOT / 'algebraic_cover_topology.json').read_text())
    hull = prop['hull']
    tri = prop['triangles']
    if not (len(set(hull)) == len(hull) and len(tri) == 2 * nc - 2 - len(hull)):
        raise ValueError('Certificate acceptance condition failed at line 137')
    if not all((type(ci) is int and 0 <= ci < nc for ci in hull)):
        raise ValueError('Certificate acceptance condition failed at line 138')
    if not all((len(row) == 3 and len(set(row)) == 3 and all((type(ci) is int and 0 <= ci < nc for ci in row)) for row in tri)):
        raise ValueError('Certificate acceptance condition failed at line 139')
    min_separation = None
    for i in range(nc):
        for j in range(i):
            d = dist(C[i], C[j])
            if not d.lo > 0:
                raise ValueError('Certificate acceptance condition failed at line 143')
            min_separation = d.lo if min_separation is None else min(min_separation, d.lo)
    hullmargin = None
    for a, b in zip(hull, hull[1:] + hull[:1]):
        for j in range(nc):
            if j in (a, b):
                continue
            v = cross(sub(C[b], C[a]), sub(C[j], C[a]))
            if not v.lo > 0:
                raise ValueError('Certificate acceptance condition failed at line 150')
            hullmargin = v.lo if hullmargin is None else min(hullmargin, v.lo)
    if not len({tuple(sorted(v)) for v in tri}) == len(tri):
        raise ValueError('Certificate acceptance condition failed at line 152')
    active = {tuple(sorted(v)) for v in system['internal_active_triangles']}
    if not len(active) == nt:
        raise ValueError('Certificate acceptance condition failed at line 154')
    empty_margin = None
    inside_margin = None
    outside_margin = None
    slack = None
    area = None
    cnt = {'active_interior': 0, 'inactive_interior': 0, 'outside_vertices': 0}
    for a, b, c in tri:
        p, r2, det = circumcenter(C[a], C[b], C[c])
        if not det.lo > 0:
            raise ValueError('Certificate acceptance condition failed at line 160')
        area = det.lo if area is None else min(area, det.lo)
        for j in range(nc):
            if j in (a, b, c):
                continue
            v = dist(C[j], p) - r2
            if not v.lo > 0:
                raise ValueError('Certificate acceptance condition failed at line 164')
            empty_margin = v.lo if empty_margin is None else min(empty_margin, v.lo)
        pn = norm(p)
        if pn.hi < 1:
            inside_margin = 1 - pn.hi if inside_margin is None else min(inside_margin, 1 - pn.hi)
            if tuple(sorted((a, b, c))) in active:
                cnt['active_interior'] += 1
            else:
                gap = t.lo - r2.hi
                if not gap > 0:
                    raise ValueError('Certificate acceptance condition failed at line 174')
                slack = gap if slack is None else min(slack, gap)
                cnt['inactive_interior'] += 1
        else:
            if not pn.lo > 1:
                raise ValueError('Certificate acceptance condition failed at line 177')
            outside_margin = pn.lo - 1 if outside_margin is None else min(outside_margin, pn.lo - 1)
            if not tuple(sorted((a, b, c))) not in active:
                raise ValueError('Certificate acceptance condition failed at line 179')
            cnt['outside_vertices'] += 1
    if not cnt['active_interior'] == len(active):
        raise ValueError('Certificate acceptance condition failed at line 180')
    pairs = system['boundary_active_pairs']
    nint = len(active)
    nb = len(pairs)
    Qb = [x[ncoord + 2 * (nint + k):ncoord + 2 * (nint + k) + 2] for k in range(nb)]
    order = prop['boundary_order']
    if not sorted(order) == list(range(nb)):
        raise ValueError('Certificate acceptance condition failed at line 186')
    quads = [quadrant(Qb[k]) for k in order]
    if not (quads == sorted(quads) and set(quads) == {0, 1, 2, 3}):
        raise ValueError('Certificate acceptance condition failed at line 187')
    arcmargin = None
    centerdot = None
    near_margin = None
    for ka, kb in zip(order, order[1:] + order[:1]):
        a, b = (Qb[ka], Qb[kb])
        cr = cross(a, b)
        dd = dot(a, b)
        if not (cr.lo > 0 and dd.lo > 0):
            raise ValueError('Certificate acceptance condition failed at line 190')
        arcmargin = cr.lo if arcmargin is None else min(arcmargin, cr.lo)
        common = set(pairs[ka]) & set(pairs[kb])
        if not len(common) == 1:
            raise ValueError('Certificate acceptance condition failed at line 192')
        ci = common.pop()
        va = dot(C[ci], a)
        vb = dot(C[ci], b)
        if not (va.lo > 0 and vb.lo > 0):
            raise ValueError('Certificate acceptance condition failed at line 193')
        cd = min(va.lo, vb.lo)
        centerdot = cd if centerdot is None else min(centerdot, cd)
    for k, (a, b) in enumerate(pairs):
        for j in range(nc):
            if j in (a, b):
                continue
            v = dist(C[j], Qb[k]) - t
            if not v.lo > 0:
                raise ValueError('Certificate acceptance condition failed at line 198')
            near_margin = v.lo if near_margin is None else min(near_margin, v.lo)
    radius_den = 10 ** 65
    rl = math.isqrt(t.lo.numerator * radius_den ** 2 // t.lo.denominator)
    rh = math.isqrt(t.hi.numerator * radius_den ** 2 // t.hi.denominator) + 1
    report = {'verified': True, 'statement': 'The unique certified algebraic KKT-root centers cover the entire closed unit disk with radius sqrt(t). No global optimum or rigorous second-order minimality claim.', 'triangles': len(tri), 'hull_vertices': len(hull), 'boundary_arcs': nb, 'counts': cnt, 'n_variables': system['n_variables'], 'gauge_center': system['gauge_center'], 'anchored_unused_centers': sorted(set((a['center'] for a in system['unused_center_anchors']))), 'certified_radius_interval': [[rl, radius_den], [rh, radius_den]], 'radius_interval_decimal': [decimal(Q(rl, radius_den)), decimal(Q(rh, radius_den))], 'certified_kkt_multiplier_lower_bounds': {'distance_alpha': decimal(min((a.lo for a in distance_alpha))), 'collapsed_contact_lambda': decimal(min((a.lo for a in lambdas)))}, 'strict_margin_lower_bounds': {'center_squared_separation': decimal(min_separation), 'hull_orientation': decimal(hullmargin), 'triangle_orientation': decimal(area), 'empty_circumcircle_squared_distance': decimal(empty_margin), 'inside_vertex_unit_margin': decimal(inside_margin), 'outside_vertex_unit_margin': decimal(outside_margin), 'inactive_vertex_squared_radius_slack': decimal(slack), 'successive_boundary_cross_product': decimal(arcmargin), 'arc_center_dot_endpoint': decimal(centerdot), 'boundary_other_center_squared_distance': decimal(near_margin)}, 'elapsed_seconds': time.time() - started}
    (ROOT / 'algebraic_cover_verification.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)
    return report
if __name__ == '__main__':
    main()
