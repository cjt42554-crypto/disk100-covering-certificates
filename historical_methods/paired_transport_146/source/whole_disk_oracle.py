# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Independent numerical disk-cover scorer including all boundary extrema.

Enumerates every pair bisector/unit-circle intersection and every center's
circle antipode, plus interior Delaunay vertices. Scores each candidate by the
actual nearest of ALL centers. This is numerical verification, not an exact
arithmetic certificate. No soft radius or released offsets enter the scorer.
"""
import os
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
             'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
from scipy.spatial import Delaunay


def true_cover(C):
    C = np.asarray(C, dtype=float)
    if C.ndim != 2 or C.shape[1] != 2 or not np.all(np.isfinite(C)):
        raise ValueError('Invalid coordinate array')
    n = len(C)
    tri = Delaunay(C).simplices
    P = C[tri]
    U = P[:, 1:] - P[:, 0, None, :]
    rhs = np.sum(U * U, axis=2)
    det = np.linalg.det(U)
    usable = np.abs(det) > 1e-14
    rel = np.linalg.solve(2 * U[usable], rhs[usable])
    vertices = P[usable, 0] + rel
    vertices = vertices[np.sum(vertices * vertices, axis=1) <= 1]

    i, j = np.triu_indices(n, 1)
    u = C[j] - C[i]
    ss = np.sum(u * u, axis=1)
    h = (np.sum(C[j] * C[j], axis=1) - np.sum(C[i] * C[i], axis=1)) / 2
    discriminant = ss - h * h
    valid = (ss > 1e-20) & (discriminant >= 0)
    u, ss, h, discriminant = u[valid], ss[valid], h[valid], discriminant[valid]
    perpendicular = np.column_stack((-u[:, 1], u[:, 0]))
    circle_base = h[:, None] * u / ss[:, None]
    circle_offset = np.sqrt(discriminant)[:, None] * perpendicular / ss[:, None]
    boundary = np.concatenate((circle_base + circle_offset, circle_base - circle_offset))
    if len(boundary) and np.max(np.abs(np.sum(boundary * boundary, axis=1) - 1)) > 1e-7:
        raise ValueError('Unstable pair/circle intersection')

    norm = np.linalg.norm(C, axis=1)
    nonzero = norm > 1e-14
    directions = C[nonzero] / norm[nonzero, None]
    antipodes = -directions
    axes = np.array([[1., 0.], [0., 1.], [-1., 0.], [0., -1.]])
    # Both directions are harmless disk points; antipodes are the potential
    # interior-of-arc maxima. Axes also cover the zero-center special case.
    all_points = np.concatenate((vertices, boundary, antipodes, directions, axes))
    best_squared = -1.
    best_point = None
    best_nearest = None
    best_index = None
    for begin in range(0, len(all_points), 512):
        chunk = all_points[begin:begin + 512]
        distances = np.sum((chunk[:, None, :] - C[None, :, :]) ** 2, axis=2)
        nearest = np.argmin(distances, axis=1)
        values = distances[np.arange(len(chunk)), nearest]
        at = int(np.argmax(values))
        if values[at] > best_squared:
            best_squared = float(values[at])
            best_point = chunk[at].tolist()
            best_nearest = int(nearest[at])
            best_index = begin + at
    count_interior, count_boundary = len(vertices), len(boundary)
    if best_index < count_interior:
        kind = 'interior_delaunay_vertex'
    elif best_index < count_interior + count_boundary:
        kind = 'all_pair_unit_circle_intersection'
    else:
        kind = 'antipode_or_auxiliary_circle_point'
    return float(np.sqrt(best_squared)), {
        'method': 'all-pair-circle+all-antipodes+interior-Delaunay-nearest-distance',
        'all_pairs_examined': len(i), 'interior_candidates': count_interior,
        'boundary_intersection_candidates': count_boundary,
        'antipodes': len(antipodes), 'auxiliary_circle_points': len(directions) + 4,
        'maximum_candidate_kind': kind, 'maximum_point': best_point,
        'nearest_center_at_maximum': best_nearest, 'maximum_squared_distance': best_squared,
        'numerical_only': True}
