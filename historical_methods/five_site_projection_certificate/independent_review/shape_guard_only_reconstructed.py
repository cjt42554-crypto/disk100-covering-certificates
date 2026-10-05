# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Certify an exact algebraic five-site tie near a rational configuration.

The two added centers are projected radially onto the primary triangle's
exact circle. Full-disk coverage is checked at the original rational centers
and transported using an exact rational upper bound on the displacement.
"""
import argparse
from fractions import Fraction as F
from functools import cmp_to_key
import hashlib
import importlib.util
import json
from math import isqrt
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
VERIFIER = ROOT/'global_lower_sync_1132/source/disk100_lower_bound_0p1132_20261004/upper_configuration/verify_cover.py'
EXPECTED = '06806070a7c185bc0c710629561d4c2d713eaf231efaadb9f7b9924a05ecb514'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def pair(x):
    return [x.numerator, x.denominator]


def cross(a, b):
    return a[0]*b[1]-a[1]*b[0]


def norm(a):
    return a[0]*a[0]+a[1]*a[1]


def sqrt_interval(x, bits=100):
    require(x >= 0, 'negative square-root argument')
    Q = 1 << bits
    k = isqrt((x.numerator*Q*Q)//x.denominator)
    lo = F(k, Q)
    hi = lo if k*k*x.denominator == x.numerator*Q*Q else F(k+1, Q)
    require(lo*lo <= x <= hi*hi, 'square-root isolation failed')
    return lo, hi


def circle(points, triangle):
    a, b, c = [points[i] for i in triangle]
    dx, dy = b[0]-a[0], b[1]-a[1]
    ex, ey = c[0]-a[0], c[1]-a[1]
    det = 2*(dx*ey-dy*ex)
    require(det != 0, 'primary triangle is collinear')
    bn, cn = norm(b)-norm(a), norm(c)-norm(a)
    q = ((bn*ey-dy*cn)/det, (dx*cn-bn*ex)/det)
    r2 = norm((q[0]-a[0], q[1]-a[1]))
    require(r2 > 0, 'primary radius is zero')
    require(all(norm((q[0]-points[i][0], q[1]-points[i][1])) == r2
                for i in triangle), 'exact primary circle identity failed')
    return q, r2


def ray_order(vectors):
    def half(v):
        return 0 if v[1] > 0 or (v[1] == 0 and v[0] >= 0) else 1
    def compare(i, j):
        a, b = vectors[i], vectors[j]
        if half(a) != half(b):
            return -1 if half(a) < half(b) else 1
        x = cross(a, b)
        if x:
            return -1 if x > 0 else 1
        require(i == j, 'two selected sites have the same circle ray')
        return 0
    order = sorted(vectors, key=cmp_to_key(compare))
    require(all(cross(vectors[a],vectors[b]) > 0
                for a,b in zip(order,order[1:]+order[:1])),
            'circle center is not strictly inside the five-site polygon')
    return order


def certify(csv, primary, added, radius):
    require(__debug__ and sys.flags.no_site, 'Use unoptimized python3 -S -B')
    require(hashlib.sha256(VERIFIER.read_bytes()).hexdigest() == EXPECTED,
            'full-cover acceptor source changed')
    spec = importlib.util.spec_from_file_location('cover_acceptor', VERIFIER)
    acceptor = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(acceptor)
    rows, source_sha = acceptor.read_centers(str(csv))
    require(all(len(row) == 2 for row in rows),
            'each center must have exactly two coordinates')
    points = [tuple(F(x) for x in p) for p in rows]
    require(len(points) == 100, 'exactly100 input centers required')
    sites = primary + added
    require(len(primary)==3 and len(added)==2 and len(set(sites))==5,
            'three primary and two distinct added sites required')
    require(all(0 <= i < 100 for i in sites), 'site index outside range')
    q, r2 = circle(points, primary)
    require(norm(q) < 1, 'common circle center is outside the open unit disk')
    vectors = {i:(points[i][0]-q[0],points[i][1]-q[1]) for i in sites}
    require(all(norm(v)>0 for v in vectors.values()), 'zero circle ray')
    order = ray_order(vectors)
    other = [i for i in range(100) if i not in sites]
    gaps = [(norm((points[i][0]-q[0],points[i][1]-q[1]))-r2,i) for i in other]
    gap, gap_site = min(gaps)
    require(gap > 0, 'another center blocks the five-site nearest circle')
    Rlo, Rhi = sqrt_interval(r2)
    centers = []
    D = F(0)
    for i in added:
        d2 = norm(vectors[i])
        lo, hi = sqrt_interval(d2)
        displacement_upper = max(Rhi-lo,hi-Rlo,F(0))
        D = max(D,displacement_upper)
        k = r2/d2
        require(k > 0 and k*d2 == r2, 'projection identity failed')
        centers.append({'site':i,'rational_base_q':[pair(x) for x in q],
                        'rational_ray_vector':[pair(x) for x in vectors[i]],
                        'scale_squared':pair(k),'positive_scale_branch':True,
                        'formula':'q+sqrt(scale_squared)*rational_ray_vector',
                        'displacement_upper':pair(displacement_upper)})
    original_cover = acceptor.verify(rows,radius)
    require(original_cover['status']=='PASS', 'source full-disk cover failed')
    U = F(radius)
    cap = U+D
    return {'status':'PASS','source_csv_sha256':source_sha,
            'full_cover_acceptor_sha256':EXPECTED,'source_closed_disk_cover':original_cover,
            'primary_sites':primary,'added_sites':added,
            'common_circle_center':[pair(x) for x in q],
            'common_circle_radius_squared':pair(r2),
            'exact_algebraic_added_centers':centers,
            'remaining95_centers_unchanged_rational':True,
            'distinct_convex_five_sites':True,'circle_center_inside_five_polygon':True,
            'cyclic_ray_order':order,'other_strictly_farther_centers':len(other),
            'minimum_other_site_squared_gap':pair(gap),'minimum_gap_site':gap_site,
            'maximum_center_displacement_upper':pair(D),'cover_radius':pair(cap),
            'cover_radius_integer_polynomial_ascending':[-cap.numerator,cap.denominator],
            'cover_radius_approx':float(cap),'projected_closed_unit_disk_cover':True,
            'argument':'Each projected center moves by at most D; the nearest-distance function therefore rises by at most D everywhere on the closed unit disk.',
            'scope':'Specified projected five-site algebraic event configuration only; no path, stationary-point, record improvement or global optimality claim.'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('coordinates',type=Path)
    p.add_argument('--primary',nargs=3,type=int,required=True)
    p.add_argument('--added',nargs=2,type=int,required=True)
    p.add_argument('--radius',required=True)
    p.add_argument('--output',type=Path,required=True)
    a = p.parse_args()
    require(not a.output.exists(), 'refuse overwriting an existing certificate')
    result = certify(a.coordinates,a.primary,a.added,a.radius)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','cover_radius_approx':result['cover_radius_approx'],
                      'other_strictly_farther_centers':95}))


if __name__ == '__main__':
    main()
