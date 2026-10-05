#!/usr/bin/env python3
"""Recompute both proofs, then show this fixed configuration is nonglobal.

The companion is a bundled dependency. No stored report or success flag
is used as acceptance of its root existence or whole closed-disk coverage.
"""
from contextlib import redirect_stdout
from fractions import Fraction as F
import hashlib
import io
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
SNAPSHOT=HERE/'companion_snapshot'
EXPECTED_COMPANION_DEPENDENCIES={'stdioverify.py': '2e1006536cbb5f8a2ab0624f67f711647666bd3f56de90a7edcc53cbaa57f597', 'exact_algebraic_cover.py': '0f12003fd5620bfd21f3a0a8ff3a11237b1da726c20e63798e6651a805e61232', 'exact_krawczyk.py': '84720038ce84f3a9d6839f25735771328b41e2f309898529a4e4efdce8939225', 'source_configuration.json': '84fb2a5dd9869e641f18250b8839ceb44bfe1263fa54bf588cfdb4867fb28a36', 'polynomial_system.json': 'c2da24e7a5109c776c81d49ae3a82c36ba0d891ad0de553a5d0d4e9223205b44', 'root_high_precision.json': '618ab635ddae5689aee23ca829c614f0927c0c797326a8bdd574c2cc4a822ec3', 'krawczyk_preconditioner.json.gz': '131ddfbcb2d29d2e0e74ee4a561b3a8ba883f33f89955081d574506c284e97a3', 'algebraic_cover_topology.json': 'd8c40cf0541aed6a1d793e02e652a570094acc4d2c9069d1b34069a5578906a4', 'LOCAL_PROOF.md': 'd703c32429871485e8907445b70ff57e809f6293c17c5b626892c177b76eeae1'}
sys.path.insert(0,str(HERE))
from verify_exact import verify_files
from critical_candidates import require

def main():
    fixed=verify_files(HERE/'configuration.json',HERE/'exact_certificate.json',HERE/'exact_result.json')
    source=json.loads((SNAPSHOT/'SNAPSHOT_SHA256.json').read_text())
    require(source['files']==EXPECTED_COMPANION_DEPENDENCIES,'companion source inventory differs from pinned dependencies')
    for name,expected in source['files'].items():
        require(hashlib.sha256((SNAPSHOT/name).read_bytes()).hexdigest()==expected,'companion dependency identity mismatch: '+name)
    require(source['frozen_companion_manifest_sha256']=='4bc4d4647c5b5e8e6e7bb01d1ce10f50fae1823139ce0333efffa15f329472e2','wrong formal companion provenance')
    require(hashlib.sha256((SNAPSHOT/'polynomial_system.json').read_bytes()).hexdigest()=='c2da24e7a5109c776c81d49ae3a82c36ba0d891ad0de553a5d0d4e9223205b44','wrong companion polynomial system')
    system=json.loads((SNAPSHOT/'polynomial_system.json').read_text())
    require(system['n_variables']==932 and system['n_centers']==100,'wrong companion configuration size')
    # main() rechecks system semantics, proves an exact Krawczyk contraction
    # and root inclusion, and proves full closed-disk coverage by rational
    # interval geometry. Its freshly returned interval is the only one used.
    sys.path.insert(0,str(SNAPSHOT))
    from exact_algebraic_cover import main as accept_companion
    with redirect_stdout(io.StringIO()):
        companion=accept_companion()
    detail=json.loads((HERE/'exact_result.json').read_text())
    actual=companion['certified_radius_interval']
    require(actual==detail['external_companion_comparison']['given_companion_radius_interval'],'fresh companion interval differs from comparison data')
    fixed_lo=F(*fixed['isolating_interval'][0])
    smaller_hi=F(*actual[1])
    require(smaller_hi<fixed_lo,'companion radius is not strictly smaller')
    gap=fixed_lo-smaller_hi
    require([gap.numerator,gap.denominator]==detail['external_companion_comparison']['strict_radius_gap_lower_bound'],'incorrect strict comparison gap')
    result={'fixed_radius_exact_verified':True,'companion_root_and_whole_closed_disk_cover_freshly_recomputed':True,
            'companion_n_variables':932,'companion_n_centers':100,
            'fixed_worst_candidate':fixed['unique_worst_candidate'],
            'companion_radius_interval':actual,'fixed_radius_interval':fixed['isolating_interval'],
            'strict_radius_gap_lower_bound':[gap.numerator,gap.denominator],
            'statement':'This fixed rational configuration is strictly nonglobal: the independently verified companion is a 100-center full closed-unit-disk cover with strictly smaller radius. Its fixed-radius polynomial is not a global optimum polynomial.'}
    print(json.dumps(result,indent=2))
    print('PASS: both exact proofs freshly recomputed; fixed configuration is strictly nonglobal.')
    return result

if __name__=='__main__':main()
