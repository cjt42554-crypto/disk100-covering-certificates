# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Independent matrix-level sanity and error audit; not an SOSC certificate."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.linalg import null_space

ROOT = Path(__file__).resolve().parent


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    result = json.loads((ROOT / 'RESULT.json').read_text())
    rows = []
    for report in result['source_reports']:
        i, free = report['source'], report['free_label']
        z = np.load(ROOT / report['matrices_file'])
        C, A, H, N, cols = z['C'], z['contact_J'], z['H_analytic_physical'], z['core_N'], z['core_columns']
        weights = z['contact_multipliers']
        gauge = z['rotation_gauge']
        B = np.vstack((np.column_stack((A[:, cols], -np.ones(len(A)))), np.r_[gauge, 0.]))
        bs = np.linalg.svd(B, compute_uv=False)
        ungauged = null_space(A, rcond=1e-9)
        ungauged_values = np.linalg.eigvalsh(ungauged.T @ H @ ungauged)
        rotation = np.column_stack((-C[:, 1], C[:, 0])).ravel()
        rotation /= np.linalg.norm(rotation)
        flat = np.eye(200)[:, [2 * free, 2 * free + 1]]
        norms = dict(H_rotation_norm=float(np.linalg.norm(H @ rotation)),
                     A_rotation_norm=float(np.linalg.norm(A @ rotation)),
                     free_A_max=float(np.max(np.abs(A @ flat))),
                     free_H_max=float(np.max(np.abs(H @ flat))))
        numerical_error = max(float(np.linalg.norm(N.T @ (z[f'H_fd_1e{h}'] - H)[np.ix_(cols, cols)] @ N, 2))
                              for h in [5, 6, 7])
        mineig = float(np.linalg.eigvalsh(N.T @ H[np.ix_(cols, cols)] @ N)[0])
        assert norms['free_A_max'] == 0 and norms['free_H_max'] == 0
        assert max(norms['H_rotation_norm'], norms['A_rotation_norm']) < 1e-12
        assert len(bs) == len(A) + 1 and bs.min() > .001
        assert np.count_nonzero(np.abs(ungauged_values) < 1e-8) == 3
        assert weights.min() > 1e-4 and abs(weights.sum() - 1) < 1e-12
        assert mineig > 10000 * numerical_error and mineig > .001
        assert report['corrected_topology_same_as_raw']
        assert abs(report['actual_minus_model_radius']) < 2e-13
        assert sha(ROOT / report['corrected_file']) == report['corrected_file_sha256']
        assert sha(ROOT / report['matrices_file']) == report['matrices_file_sha256']
        assert sha(ROOT / 'inputs' / f'basin{i:02d}.npy') == report['input_sha256']
        raw = ROOT.parent / 'fresh_rattler_geometry_bridge/inputs' / f'basin{i:02d}.npy'
        assert sha(raw) == report['input_sha256']
        rows.append(dict(source=i, pass_numerical_audit=True,
                         epigraph_plus_gauge_full_row_rank_min_singular=float(bs.min()),
                         ungauged_tangent_dimension=ungauged.shape[1],
                         ungauged_zero_count=3,
                         numerical_core_minimum_eigenvalue=mineig,
                         maximum_reduced_finite_difference_error=numerical_error,
                         positive_to_error_ratio=float(mineig / numerical_error), **norms))
    out = dict(schema='fresh99-second-order-audit-v1', passed=True, numerical_only=True,
               all_three_sources_no_robust_negative_direction=True,
               all_three_sources_positive_core_and_two_free_flats=True,
               no_negative_direction_proposals_generated=True,
               no_strict_local_or_global_claim=True, no_interval_root_or_Hessian_enclosure=True,
               limitations=['Rank, multipliers, curvature and geometry eligibility are binary64 numerical checks.',
                            'Positive all-active tangent alone is insufficient when a contact multiplier is zero.',
                            'Here all contact multipliers are numerically positive; cone=tangent would require a validated root and positive multiplier enclosure.',
                            'Core positivity cannot make the original100 configuration strict: its free center has two genuine flat directions.'],
               source_audits=rows)
    (ROOT / 'AUDIT.json').write_text(json.dumps(out, indent=2) + '\n')


if __name__ == '__main__':
    main()
