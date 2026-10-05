# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Exact generic four-site event identity, not a covering/optimality proof.

For A=orient(a,b,c), B=orient(a,b,d), L=|b-a|^2, and
Nc=(c-a).(c-b), Nd=(d-a).(d-b), standard incircle D satisfies
D=Nc B-Nd A.  Clearing denominators in circumradii gives
4 A^2 B^2 (g_abc-g_abd)=L D (Nc B+Nd A).
Consequently, on D=0 and A B !=0, their gradient difference is
beta grad(D), beta=L(Nc B+Nd A)/(4 A^2 B^2).
This explains a redundant epigraph row at an event.  It does not
certify a nearest empty circle, a regular KKT point or a path.
All generic polynomial coefficients are recomputed over integers.
"""
from fractions import Fraction
from itertools import permutations
from pathlib import Path
import hashlib
import json

N = 8
ZERO = (0,) * N


def scalar(v):
    return {} if v == 0 else {ZERO: v}


def variable(i):
    e = list(ZERO)
    e[i] = 1
    return {tuple(e): 1}


def add(a, b):
    r = a.copy()
    for k, v in b.items():
        r[k] = r.get(k, 0) + v
        if not r[k]:
            del r[k]
    return r


def neg(a):
    return {k: -v for k, v in a.items()}


def sub(a, b):
    return add(a, neg(b))


def mul(a, b):
    r = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(x + y for x, y in zip(ka, kb))
            r[k] = r.get(k, 0) + va * vb
    return {k: v for k, v in r.items() if v}


def sq(a):
    return mul(a, a)


def dot(a, b):
    return add(mul(a[0], b[0]), mul(a[1], b[1]))


def vecsub(a, b):
    return [sub(x, y) for x, y in zip(a, b)]


def cross(a, b):
    return sub(mul(a[0], b[1]), mul(a[1], b[0]))


def determinant(rows):
    ans = {}
    for p in permutations(range(len(rows))):
        inv = sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))
        term = scalar((-1) ** inv)
        for i, j in enumerate(p):
            term = mul(term, rows[i][j])
        ans = add(ans, term)
    return ans


def evaluate(p, x):
    total = Fraction(0)
    for e, coefficient in p.items():
        term = Fraction(coefficient)
        for v, power in zip(x, e):
            term *= v ** power
        total += term
    return total


def derivative(p, i):
    r = {}
    for e, coefficient in p.items():
        if e[i]:
            k = list(e)
            k[i] -= 1
            r[tuple(k)] = coefficient * e[i]
    return r


def check_equal(a, b, label):
    difference = sub(a, b)
    if difference:
        raise ValueError(label + ': nonzero generic integer polynomial')


def gradient_radius(L, area, numerator, x):
    l, o, n = (evaluate(p, x) for p in (L, area, numerator))
    if not o:
        raise ValueError('degenerate triangle')
    ans = []
    for i in range(N):
        dl, do, dn = (evaluate(derivative(p, i), x) for p in (L, area, numerator))
        ans.append(dl / 4 + dl * n * n / (4 * o * o)
                   + l * n * dn / (2 * o * o)
                   - l * n * n * do / (2 * o * o * o))
    return ans


def main():
    a, b, c, d = [[variable(2 * i), variable(2 * i + 1)] for i in range(4)]
    e = vecsub(b, a)
    L = dot(e, e)
    A = cross(e, vecsub(c, a))
    B = cross(e, vecsub(d, a))
    Nc = dot(vecsub(c, a), vecsub(c, b))
    Nd = dot(vecsub(d, a), vecsub(d, b))
    D = determinant([[p[0], p[1], dot(p, p), scalar(1)] for p in (a, b, c, d)])
    check_equal(D, sub(mul(Nc, B), mul(Nd, A)), 'incircle sign identity')
    for point, area, num in ((c, A, Nc), (d, B, Nd)):
        left = mul(dot(vecsub(point, a), vecsub(point, a)),
                   dot(vecsub(point, b), vecsub(point, b)))
        check_equal(left, add(sq(area), sq(num)), 'circumradius numerator identity')
    S = add(mul(Nc, B), mul(Nd, A))
    cleared_difference = mul(L, sub(mul(sq(Nc), sq(B)), mul(sq(Nd), sq(A))))
    factorization = mul(mul(L, D), S)
    check_equal(cleared_difference, factorization, 'event factorization')
    # These exact cyclic examples cover beta=0 and beta!=0.  They are
    # formula checks, not proposed 100-center coverings or physical events.
    examples = [
        [Fraction(v) for v in (1, 0, -1, 0, 0, 1, 0, -1)],
        [Fraction(1), Fraction(0), Fraction(-3, 5), Fraction(4, 5),
         Fraction(0), Fraction(1), Fraction(0), Fraction(-1)],
    ]
    checks = []
    for x in examples:
        if evaluate(D, x) != 0:
            raise ValueError('example not exactly cyclic')
        aa, bb = evaluate(A, x), evaluate(B, x)
        beta = evaluate(L, x) * evaluate(S, x) / (4 * aa * aa * bb * bb)
        ga = gradient_radius(L, A, Nc, x)
        gb = gradient_radius(L, B, Nd, x)
        gd = [evaluate(derivative(D, i), x) for i in range(N)]
        if any(u - v != beta * w for u, v, w in zip(ga, gb, gd)):
            raise ValueError('exact gradient relation failed')
        checks.append({'beta': str(beta), 'gradient_components_checked': N,
                       'triangles_noncollinear': aa != 0 and bb != 0,
                       'opposite_sides_of_shared_edge': aa * bb < 0})
    coefficients = [[list(k), v] for k, v in sorted(factorization.items())]
    result = {'status': 'PASS', 'integer_polynomial_identities': 4,
              'generic_variables': ['ax', 'ay', 'bx', 'by', 'cx', 'cy', 'dx', 'dy'],
              'standard_incircle_monomials': len(D),
              'cleared_factorization_monomials': len(factorization),
              'coefficient_sha256': hashlib.sha256(json.dumps(coefficients, separators=(',', ':')).encode()).hexdigest(),
              'exact_cyclic_gradient_checks': checks,
              'domain': 'orient(a,b,c)*orient(a,b,d) != 0',
              'normalized_event_row_coefficient': 'beta*D0, with E=D/D0; divide by t0 if holes divided by t0',
              'scope': 'Generic exact algebraic redundancy only. No nearest-circle qualification, KKT regularity, continuous covering or global optimum.'}
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
