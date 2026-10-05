"""Exact integer verifier for a local quadratic-system root enclosure.

Every midpoint has denominator D=10^60, every preconditioner entry denominator
P=2^40, and the uniform root box radius is rho=10^-58. The acceptance path uses
only Python's standard library and integers/Fraction. No floating number
participates in the acceptance inequalities. --generate additionally uses
NumPy/SciPy to propose a dyadic approximate inverse, then fully checks it.
"""
import os
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import json, time, sys, hashlib, gzip
from fractions import Fraction
from decimal import Decimal, localcontext, ROUND_HALF_EVEN
from pathlib import Path
ROOT = Path(__file__).resolve().parent
D = 10 ** 60
P = 2 ** 40
RHODEN = 10 ** 58

def validate_system(system):
    n = system['n_variables']
    if not (type(n) is int and n > 0):
        raise ValueError('Certificate acceptance condition failed at line 21')
    if not len(system['polynomials']) == n:
        raise ValueError('Certificate acceptance condition failed at line 22')
    for terms in system['polynomials']:
        if not type(terms) is list:
            raise ValueError('Certificate acceptance condition failed at line 24')
        for term in terms:
            if not (type(term) is list and len(term) == 2):
                raise ValueError('Certificate acceptance condition failed at line 26')
            coeff, vs = term
            if not (type(vs) is list and len(vs) <= 2):
                raise ValueError('Certificate acceptance condition failed at line 28')
            if not all((type(v) is int and 0 <= v < n for v in vs)):
                raise ValueError('Certificate acceptance condition failed at line 29')
            if type(coeff) is list:
                if not (len(coeff) == 2 and all((type(v) is int for v in coeff))):
                    raise ValueError('Certificate acceptance condition failed at line 31')
                if not (coeff[1] > 0 and (not vs) and (D * D % coeff[1] == 0)):
                    raise ValueError('Certificate acceptance condition failed at line 32')
            elif not type(coeff) is int:
                raise ValueError('Certificate acceptance condition failed at line 34')

def frac_text(f):
    with localcontext() as c:
        c.prec = 35
        return str(Decimal(f.numerator) / Decimal(f.denominator))

def evaluate_integer(system, x):
    n = system['n_variables']
    f = [0] * n
    J = [{} for _ in range(n)]
    lip = [0] * n
    for row, terms in enumerate(system['polynomials']):
        for coeff, vs in terms:
            if isinstance(coeff, list):
                if not (not vs and D * D % coeff[1] == 0):
                    raise ValueError('Certificate acceptance condition failed at line 46')
                f[row] += coeff[0] * (D * D // coeff[1])
                continue
            if not vs:
                f[row] += coeff * D * D
            elif len(vs) == 1:
                f[row] += coeff * x[vs[0]] * D
            else:
                f[row] += coeff * x[vs[0]] * x[vs[1]]
            for k, col in enumerate(vs):
                others = vs[:k] + vs[k + 1:]
                value = coeff * (D if not others else x[others[0]])
                J[row][col] = J[row].get(col, 0) + value
                if others:
                    lip[row] += abs(coeff)
    return (f, J, lip)

def generate(system, x, J):
    import numpy as np
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import splu
    n = system['n_variables']
    row = []
    col = []
    data = []
    for i, entries in enumerate(J):
        for j, v in entries.items():
            row.append(i)
            col.append(j)
            data.append(float(Fraction(v, D)))
    lu = splu(coo_matrix((data, (row, col)), shape=(n, n)).tocsc())
    approx = lu.solve(np.eye(n))
    scaled = np.rint(approx * P)
    if not np.max(np.abs(scaled)) < 2 ** 62:
        raise ValueError('Certificate acceptance condition failed at line 67')
    R = scaled.astype(np.int64)
    np.savez_compressed(ROOT / 'krawczyk_preconditioner.npz', R=R)
    payload = json.dumps({'denominator': P, 'matrix': R.tolist()}, separators=(',', ':')).encode()
    (ROOT / 'krawczyk_preconditioner.json.gz').write_bytes(gzip.compress(payload, mtime=0))

def verify(generate_requested=None, emit=True):
    started = time.time()
    system = json.loads((ROOT / 'polynomial_system.json').read_text())
    validate_system(system)
    n = system['n_variables']
    high = json.loads((ROOT / 'root_high_precision.json').read_text())
    if not (type(high['variables']) is list and len(high['variables']) == n):
        raise ValueError('Certificate acceptance condition failed at line 76')
    if not all((type(v) is str for v in high['variables'])):
        raise ValueError('Certificate acceptance condition failed at line 77')
    with localcontext() as c:
        c.prec = 130
        x = [int((Decimal(v) * D).to_integral_value(rounding=ROUND_HALF_EVEN)) for v in high['variables']]
    (ROOT / 'rational_midpoint.json').write_text(json.dumps({'denominator': D, 'numerators': [str(v) for v in x]}, separators=(',', ':')))
    f, J, lip = evaluate_integer(system, x)
    if generate_requested is None:
        generate_requested = '--generate' in sys.argv
    if generate_requested:
        generate(system, x, J)
    witness = json.loads(gzip.decompress((ROOT / 'krawczyk_preconditioner.json.gz').read_bytes()))
    if not witness['denominator'] == P:
        raise ValueError('Certificate acceptance condition failed at line 86')
    R = witness['matrix']
    if not (len(R) == n and all((len(row) == n for row in R))):
        raise ValueError('Certificate acceptance condition failed at line 87')
    if not all((type(v) is int for row in R for v in row)):
        raise ValueError('Certificate acceptance condition failed at line 88')
    PD = P * D
    error_rows = [0] * n
    columns = [[] for _ in range(n)]
    for i, entries in enumerate(J):
        for j, v in entries.items():
            if v:
                columns[j].append((i, v))
    for j, entries in enumerate(columns):
        for row in range(n):
            value = sum((R[row][i] * v for i, v in entries)) - (PD if row == j else 0)
            error_rows[row] += abs(value)
    eps = Fraction(int(max(error_rows)), PD)
    rf = [sum((v * ff for v, ff in zip(row, f))) for row in R]
    eta = Fraction(max(map(abs, rf)), P * D * D)
    beta = Fraction(max((sum((abs(v) * l for v, l in zip(row, lip))) for row in R)), P)
    rho = Fraction(1, RHODEN)
    q = eps + rho * beta
    image = eta + q * rho
    if not q < 1:
        raise ValueError('Contraction norm must be less than 1')
    if not image < rho:
        raise ValueError('Krawczyk image must lie strictly inside the root box')
    report = {'verified': True, 'statement': 'Unique real root of the explicit local KKT polynomial system in the certified box; no global optimum claim.', 'n_variables': n, 'midpoint_denominator': D, 'inverse_denominator': P, 'box_radius': [1, RHODEN], 'exact_inequalities': {'epsilon': [eps.numerator, eps.denominator], 'eta': [eta.numerator, eta.denominator], 'jacobian_lipschitz_preconditioned': [beta.numerator, beta.denominator], 'contraction': [q.numerator, q.denominator], 'image_radius': [image.numerator, image.denominator]}, 'decimal_diagnostics': {'epsilon': frac_text(eps), 'eta': frac_text(eta), 'beta': frac_text(beta), 'contraction': frac_text(q), 'image_radius': frac_text(image), 'box_radius': frac_text(rho)}, 'elapsed_seconds': time.time() - started, 'system_sha256': hashlib.sha256((ROOT / 'polynomial_system.json').read_bytes()).hexdigest(), 'preconditioner_sha256': hashlib.sha256((ROOT / 'krawczyk_preconditioner.json.gz').read_bytes()).hexdigest()}
    (ROOT / 'krawczyk_verification.json').write_text(json.dumps(report, indent=2))
    if emit:
        print(json.dumps(report, indent=2), flush=True)
    return report
if __name__ == '__main__':
    verify()
