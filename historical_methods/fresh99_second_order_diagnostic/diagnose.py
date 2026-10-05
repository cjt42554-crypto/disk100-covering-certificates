# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Three-source numerical 99-core KKT and physical second-order diagnosis.

No full100 optimizer. Free-center coordinates are exact unconstrained flat
coordinates of the local active-branch model, never anchors adding curvature.
Only this tree is written. All Hessians use t=R^2 objective multipliers.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import hashlib
import importlib.util
import json
from pathlib import Path
import time
import numpy as np
from scipy.linalg import null_space, eigh
from scipy.optimize import nnls

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent
FREES = [24, 73, 80]
SHAS = ['e053e7012512613df84d2329e549f6a9722d263eb7154116a6bd5771b79017ce',
        '2119c57582b6142b7650285b8329cbcf7f132b752cc79b18150a3096f9f13db5',
        'e482be7b4d0f6eb8d94316eee003c371d1017bea70efa9e3ba74548e7e591367']


def mod(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


geom = mod('fresh99_condensed_geometry', BASE / 'topology_patch_hybrid_146/source/geometry.py')
oracle = mod('fresh99_whole_domain', BASE / 'topology_patch_hybrid_146/source/whole_disk_oracle.py')


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write(p, x):
    Path(p).write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


class ContactBranches:
    def __init__(self, C, tol=1e-12):
        T, E, sg, F = geom.topology(C)
        assert not len(F), 'source has active singleton boundary branch requiring new implementation'
        tg, p, _ = geom.triangles(C, T)
        bg, q, _ = geom.boundary(C, E, sg)
        g = np.r_[tg, bg]
        active = np.flatnonzero(g.max() - g < tol)
        ti = active[active < len(T)]
        bi = active[active >= len(T)] - len(T)
        self.T, self.E, self.sg = T[ti], E[bi], sg[bi]
        self.keys = ([['t', *sorted(map(int, t))] for t in self.T]
                     + [['b', int(e[0]), int(e[1]), int(s)] for e, s in zip(self.E, self.sg)])
        self.m = len(self.keys)
        self.inventory_top = (T, E, sg, F)
        self.initial_g, self.initial_active = g, active

    def evaluate(self, C):
        tg, p, tj = geom.triangles(C, self.T)
        bg, q, bj = geom.boundary(C, self.E, self.sg)
        J = np.zeros((self.m, 200))
        for k, (ids, gradient) in enumerate(zip(self.T, tj)):
            cols = (2 * ids[:, None] + np.array([0, 1])).ravel()
            J[k, cols] = gradient.ravel()
        for k, (ids, gradient) in enumerate(zip(self.E, bj), len(self.T)):
            cols = (2 * ids[:, None] + np.array([0, 1])).ravel()
            J[k, cols] = gradient.ravel()
        return np.r_[tg, bg], np.vstack((p, q)), J

    def physical_hessian(self, C, weights):
        _, points, _ = self.evaluate(C)
        H = np.zeros((200, 200))
        minimum_bary = 1.
        minimum_IFT = np.inf
        betas = []
        local_records = []
        for k, (p, lam) in enumerate(zip(points, weights)):
            inside = k < len(self.T)
            ids = self.T[k] if inside else self.E[k - len(self.T)]
            P = C[ids]
            d = len(ids)
            RHS = np.zeros((2, 2 * d))
            if inside:
                B = P[1:] - P[0]
                for j in [1, 2]:
                    RHS[j - 1, :2] = -(P[0] - p)
                    RHS[j - 1, 2 * j:2 * j + 2] = P[j] - p
                w = np.linalg.solve(np.vstack((P.T, np.ones(3))), np.r_[p, 1.])
                beta = 0.
            else:
                B = np.vstack((P[1] - P[0], p))
                RHS[0, :2] = -(P[0] - p)
                RHS[0, 2:4] = P[1] - p
                tangent = np.array([-p[1], p[0]])
                w0 = P[1] @ tangent / ((P[1] - P[0]) @ tangent)
                w = np.array([w0, 1 - w0])
                beta = -lam * (1 - (w @ P) @ p)
                betas.append(beta)
            Q = np.linalg.solve(B, RHS)
            Hk = 2 * beta * Q.T @ Q
            for j in range(d):
                D = Q.copy()
                D[:, 2 * j:2 * j + 2] -= np.eye(2)
                Hk += 2 * lam * w[j] * D.T @ D
            cols = (2 * ids[:, None] + np.array([0, 1])).ravel()
            H[np.ix_(cols, cols)] += Hk
            minimum_bary = min(minimum_bary, w.min())
            minimum_IFT = min(minimum_IFT, abs(np.linalg.det(B)))
            local_records.append(dict(key=self.keys[k], lambda_contact=float(lam),
                                      collapsed_distance_weights=(lam * w).tolist(),
                                      barycentric_weights=w.tolist(), beta_unit=float(beta),
                                      hole_IFT_determinant=float(np.linalg.det(B)),
                                      hole=points[k].tolist()))
        return H, dict(minimum_barycentric_weight=float(minimum_bary),
                       minimum_hole_IFT_abs_determinant=float(minimum_IFT),
                       boundary_unit_multiplier_range=[float(min(betas)), float(max(betas))],
                       contacts=local_records)


def root_correct(C, branches, free):
    labels = np.delete(np.arange(100), free)
    corecols = (2 * labels[:, None] + np.array([0, 1])).ravel()
    reference = C.ravel()[corecols].copy()
    w = np.column_stack((-C[labels, 1], C[labels, 0])).ravel()
    w /= np.linalg.norm(w)
    g, _, J = branches.evaluate(C)
    A = J[:, corecols]
    weights, nnres = nnls(np.vstack((A.T, np.ones(branches.m))), np.r_[np.zeros(198), 1.], maxiter=10000)
    weights /= weights.sum()
    t, eta = float(g.mean()), 0.
    initial = dict(stationarity_norm=float(np.linalg.norm(A.T @ weights)),
                   contact_equality_range=float(g.max() - g.min()),
                   smallest_contact_singular_values=np.linalg.svd(A, compute_uv=False)[-5:].tolist(),
                   minimum_multiplier=float(weights.min()), multiplier_sum=float(weights.sum()))
    records = []
    corrected = C.copy()
    for outer in range(8):
        g, _, J = branches.evaluate(corrected)
        A = J[:, corecols]
        H = branches.physical_hessian(corrected, weights)[0][np.ix_(corecols, corecols)]
        residual = np.r_[A.T @ weights + eta * w, 1 - weights.sum(), g - t, w @ (corrected.ravel()[corecols] - reference)]
        records.append(dict(iteration=outer, residual_max=float(np.max(np.abs(residual))),
                            stationarity_norm=float(np.linalg.norm(A.T @ weights)),
                            radius=float(np.sqrt(t)), minimum_multiplier=float(weights.min())))
        if np.max(np.abs(residual)) < 2e-14:
            break
        K = np.block([[H, np.zeros((198, 1)), A.T, w[:, None]],
                      [np.zeros((1, 199)), -np.ones((1, branches.m)), np.zeros((1, 1))],
                      [A, -np.ones((branches.m, 1)), np.zeros((branches.m, branches.m + 1))],
                      [w[None], np.zeros((1, branches.m + 2))]])
        step = np.linalg.solve(K, -residual)
        records[-1]['bordered_KKT_condition2'] = float(np.linalg.cond(K))
        records[-1]['core_Newton_max_center_motion'] = float(np.linalg.norm(step[:198].reshape(99, 2), axis=1).max())
        records[-1]['predicted_t_step'] = float(step[198])
        scale = 1.
        for damping in range(12):
            trial = corrected.copy()
            trial[labels] += scale * step[:198].reshape(99, 2)
            tw = weights + scale * step[199:199 + branches.m]
            tt = t + scale * step[198]
            te = eta + scale * step[-1]
            gg, _, JJ = branches.evaluate(trial)
            rr = np.r_[JJ[:, corecols].T @ tw + te * w, 1 - tw.sum(), gg - tt,
                       w @ (trial.ravel()[corecols] - reference)]
            if np.linalg.norm(rr) < np.linalg.norm(residual) and tw.min() > 0 and tt > 0:
                corrected, weights, t, eta = trial, tw, tt, te
                records[-1]['damping'] = scale
                break
            scale *= .5
        else:
            records[-1]['failure'] = 'Newton line search failed; original and last iterate retained'
            break
    return corrected, weights, t, w, corecols, initial, records


def main():
    (ROOT / 'inputs').mkdir(exist_ok=True)
    (ROOT / 'matrices').mkdir(exist_ok=True)
    report = dict(schema='fresh99-core-second-order-numerical-v1', numerical_only=True,
                  no_strict_local_or_global_claim=True, no_full100_optimizer=True,
                  objective='t=R^2', source_reports=[])
    for i, free in enumerate(FREES):
        source = BASE / 'fresh_rattler_geometry_bridge/inputs' / f'basin{i:02d}.npy'
        assert sha(source) == SHAS[i]
        (ROOT / 'inputs' / source.name).write_bytes(source.read_bytes())
        C0 = np.load(source)
        branches = ContactBranches(C0)
        initial_g, _, initial_J = branches.evaluate(C0)
        assert np.max(np.abs(initial_J[:, 2 * free:2 * free + 2])) == 0
        C, weights, t, gauge, cols, initial, newton = root_correct(C0, branches, free)
        corrected_file = ROOT / 'inputs' / f'corrected_source{i:02d}.npy'
        np.save(corrected_file, C)
        g, points, Aall = branches.evaluate(C)
        H, qual = branches.physical_hessian(C, weights)
        A = Aall[:, cols]
        singular = np.linalg.svd(A, compute_uv=False)
        N = null_space(np.vstack((A, gauge)), rcond=1e-9)
        core_H = H[np.ix_(cols, cols)]
        reduced = N.T @ core_H @ N
        values, vectors = eigh(reduced)
        full_N = np.zeros((200, N.shape[1] + 2))
        full_N[cols, :N.shape[1]] = N
        full_N[2 * free, -2] = 1.
        full_N[2 * free + 1, -1] = 1.
        full_values = np.linalg.eigvalsh(full_N.T @ H @ full_N)
        fd_results = []
        fd_matrices = []
        for h in [1e-5, 1e-6, 1e-7]:
            Hfd = np.zeros((200, 200))
            for k in cols:
                Cp, Cm = C.copy(), C.copy()
                Cp.ravel()[k] += h
                Cm.ravel()[k] -= h
                gp = branches.evaluate(Cp)[2].T @ weights
                gm = branches.evaluate(Cm)[2].T @ weights
                Hfd[:, k] = (gp - gm) / (2 * h)
            asym = np.max(np.abs(Hfd - Hfd.T))
            Hfd = (Hfd + Hfd.T) / 2
            fd_values = np.linalg.eigvalsh(N.T @ Hfd[np.ix_(cols, cols)] @ N)
            error = np.linalg.norm(N.T @ (Hfd - H)[np.ix_(cols, cols)] @ N, 2)
            fd_results.append(dict(h=h, maximum_unsymmetrized_asymmetry=float(asym),
                                   full_H_entry_max_error=float(np.max(np.abs(Hfd - H))),
                                   reduced_H_operator2_error=float(error),
                                   minimum_eigenvalue=float(fd_values[0]), maximum_eigenvalue=float(fd_values[-1]),
                                   negative_count=int(np.sum(fd_values < -1e-8)),
                                   zero_count=int(np.sum(np.abs(fd_values) < 1e-8))))
            fd_matrices.append(Hfd)
        actual_r, actual_info = oracle.true_cover(C)
        corrected_top = geom.topology(C)
        same_topology = geom.top_key(corrected_top) == geom.top_key(branches.inventory_top)
        allT, allE, allsg, allF = corrected_top
        allg = np.r_[geom.triangles(C, allT)[0], geom.boundary(C, allE, allsg)[0]]
        gaps = allg.max() - allg
        inactive = gaps[gaps > 1e-12]
        support_clearances = []
        for key, p in zip(branches.keys, points):
            support = key[1:] if key[0] == 't' else key[1:3]
            nonsupport = np.delete(C, support, axis=0)
            support_clearances.append(float(np.min(np.sum((nonsupport - p) ** 2, axis=1)) - np.sum((C[support[0]] - p) ** 2)))
        modes = np.zeros((200, min(6, len(values))))
        modes[cols] = N @ vectors[:, :modes.shape[1]]
        matrix_file = ROOT / 'matrices' / f'source{i:02d}.npz'
        np.savez(matrix_file, C_original=C0, C=C, contact_values=g, contact_points=points,
                 contact_J=Aall, contact_multipliers=weights, rotation_gauge=gauge,
                 H_analytic_physical=H, H_fd_1e5=fd_matrices[0], H_fd_1e6=fd_matrices[1], H_fd_1e7=fd_matrices[2],
                 core_columns=cols, core_N=N, core_reduced_H=reduced, core_eigenvalues=values,
                 complete_N=full_N, complete_eigenvalues=full_values, modes=modes)
        summary = dict(source=i, input_sha256=SHAS[i], free_label=free, active_count=branches.m,
                       active_internal=len(branches.T), active_boundary=len(branches.E),
                       initial_numerical_KKT=initial, newton_records=newton,
                       corrected_file=str(corrected_file.relative_to(ROOT)), corrected_file_sha256=sha(corrected_file),
                       maximum_core_correction=float(np.linalg.norm(C-C0, axis=1).max()),
                       free_center_bitwise_unchanged=bool(np.array_equal(C[free], C0[free])),
                       stationarity_norm=float(np.linalg.norm(A.T @ weights)),
                       multiplier_sum=float(weights.sum()), minimum_multiplier=float(weights.min()),
                       maximum_active_equality_residual=float(np.max(np.abs(g-t))),
                       smallest_contact_singular_values=singular[-5:].tolist(),
                       contact_rank_rcond1e9=int(np.sum(singular > singular[0] * 1e-9)),
                       all_active_core_tangent_dimension=N.shape[1], free_exact_flat_dimensions=2,
                       maximum_tangent_null_residual=float(np.max(np.abs(A @ N))),
                       rotation_slice_null_residual=float(np.max(np.abs(gauge @ N))),
                       free_J_column_max=float(np.max(np.abs(Aall[:, 2 * free:2 * free + 2]))),
                       free_H_rowcol_max=float(np.max(np.abs(H[2 * free:2 * free + 2]))),
                       core_eigenvalues=values.tolist(), core_minimum_eigenvalue=float(values[0]),
                       core_maximum_eigenvalue=float(values[-1]),
                       core_negative_count=int(np.sum(values < -1e-8)),
                       complete_tangent_eigenvalues=full_values.tolist(),
                       complete_tangent_zero_count=int(np.sum(np.abs(full_values) < 1e-8)),
                       finite_difference_checks=fd_results, physical_geometry_qualification=qual,
                       active_non_support_squared_clearance_min=float(min(support_clearances)),
                       minimum_inactive_squared_radius_slack=float(min(inactive)),
                       corrected_topology_same_as_raw=bool(same_topology),
                       actual_whole_disk_radius=actual_r, model_radius=float(np.sqrt(t)),
                       actual_minus_model_radius=float(actual_r - np.sqrt(t)),
                       actual_whole_disk_oracle=actual_info,
                       matrices_file=str(matrix_file.relative_to(ROOT)), matrices_file_sha256=sha(matrix_file),
                       conclusion='Numerical strict-positive-stress core diagnosis; full100 has two genuine rattler flats. Numerical positive all-active curvature is not an interval SOSC certificate.',
                       no_negative_direction_trials_because_core_positive=bool(values[0] > 1e-8))
        report['source_reports'].append(summary)
        write(ROOT / 'CHECKPOINT.json', report)
        print('SOURCE', i, 'coremin', values[0], 'res', summary['stationarity_norm'], 'coremove', summary['maximum_core_correction'], flush=True)
    write(ROOT / 'RESULT.json', report)


if __name__ == '__main__':
    main()
