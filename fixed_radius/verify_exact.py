#!/usr/bin/env python3
"""Exact covering radius of the pinned rational configuration; stdlib only."""
from fractions import Fraction as F
import hashlib
import json
from math import gcd,isqrt
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from critical_candidates import candidates,rational_dist2,require,verify

EXPECTED_CONFIGURATION_SHA256='84fb2a5dd9869e641f18250b8839ceb44bfe1263fa54bf588cfdb4867fb28a36'
EXPECTED_CHECKER_SHA256='93b202d1893f5f3aaa2ab013c2d6a623ffc182a471a6ae255d54a3cb1d4d9c3e'

def verify_files(config_path,certificate_path,detail_path):
    cfgpath=Path(config_path);certpath=Path(certificate_path);detailpath=Path(detail_path)
    cfg=json.loads(cfgpath.read_text());cert=json.loads(certpath.read_text());detail=json.loads(detailpath.read_text())
    sha=hashlib.sha256(cfgpath.read_bytes()).hexdigest()
    require(sha==EXPECTED_CONFIGURATION_SHA256==detail['configuration_sha256'],'fixed configuration identity mismatch')
    require(hashlib.sha256((HERE/'critical_candidates.py').read_bytes()).hexdigest()==EXPECTED_CHECKER_SHA256,'source math checker changed')
    require(set(cert)=={'bound_squared','covering_indices','exception','exact_radius_squared','root_interval'},'unexpected certificate fields')
    for pair in (cert['bound_squared'],cert['exact_radius_squared'],*cert['root_interval']):
        require(type(pair) is list and len(pair)==2 and all(type(z) is int for z in pair) and pair[1]>0,'invalid exact rational encoding')
    # Regenerates ALL candidates. The exception is independently proved to
    # have exactly the listed nearest centers; every other candidate <= gap.
    result=verify(cfgpath,certpath)
    N,D=cert['exact_radius_squared']
    require(0<N<D and gcd(N,D)==1,'invalid reduced radius squared')
    require(detail['radius_squared']==[N,D],'detail radius mismatch')
    require(detail['polynomial_coefficients_descending']==[D,0,-N],'polynomial coefficients do not encode exact radius')
    require(detail['polynomial_degree']==2,'incorrect polynomial degree')
    require(detail['polynomial_primitive'] is True and gcd(D,N)==1,'polynomial is not primitive')
    require(not(isqrt(N)**2==N and isqrt(D)**2==D),'claimed quadratic is reducible over rationals')
    require(detail['polynomial_irreducible_over_Q'] is True,'missing exact irrationality claim')
    require(detail['maximum_key']==cert['exception'],'maximum key mismatch')
    require(detail['nearest_centers_zero_based']==cert['exception'][1:],'nearest centers mismatch')
    require(cert['exception'][0]=='T','this quadratic certificate requires an interior triple maximum')
    require(detail['isolating_interval']==cert['root_interval'],'interval metadata mismatch')
    lo=F(*cert['root_interval'][0]);hi=F(*cert['root_interval'][1])
    grid=detail['isolating_interval_common_denominator']
    nums=detail['isolating_interval_common_numerators']
    require(type(grid) is int and grid==10**60,'wrong isolation grid')
    require(type(nums) is list and len(nums)==2 and all(type(z) is int for z in nums) and nums[1]==nums[0]+1,'wrong adjacent grid endpoints')
    require([lo,hi]==[F(n,grid) for n in nums],'grid/interval inconsistency')
    require(0<lo<hi and D*lo*lo-N<0<D*hi*hi-N,'incorrect polynomial root isolation signs')
    # Since P'(r)=2Dr>0 for r>0, the signs isolate exactly one positive root.
    require(F(*detail['all_other_candidate_squared_threshold'])==F(*cert['bound_squared']),'gap threshold mismatch')
    max_point=None
    for key,P in candidates(cfg['centers_numerators'],cfg['coordinate_denominator']):
        if list(key)==cert['exception']:
            require(P[0]=='R','maximum is not rational')
            require(rational_dist2(P,cfg['centers_numerators'][key[1]],cfg['coordinate_denominator'])==F(N,D),'wrong exact radius at the proposed maximum')
            require(P[1]**2+P[2]**2<P[3]**2,'maximum is not strictly inside the disk')
            point=[F(P[1],P[3]),F(P[2],P[3])]
            require([[z.numerator,z.denominator] for z in point]==detail['point_rational'],'wrong maximum point metadata')
            max_point=point
            break
    require(max_point is not None,'maximum candidate absent')
    comparison=detail.get('external_companion_comparison')
    if comparison:
        companion_hi=F(*comparison['given_companion_radius_interval'][1])
        companion_lo=F(*comparison['given_companion_radius_interval'][0])
        require(0<companion_lo<companion_hi<lo,'given companion interval is not strictly smaller')
        require(comparison['companion_coverage_acceptance_external'] is True,'companion acceptance must be explicitly external')
        require(F(*comparison['strict_radius_gap_lower_bound'])==lo-companion_hi,'incorrect conditional radius gap')
    result.update({'fixed_radius_squared':[N,D],
                   'polynomial_coefficients_descending':[D,0,-N],
                   'isolating_interval':cert['root_interval'],
                   'scope':'Exact radius and unique worst point of this fixed rational configuration only.',
                   'external_companion_interval_strictly_smaller':bool(comparison),
                   'external_companion_coverage_not_replayed_here':bool(comparison)})
    return result

def main():
    result=verify_files(HERE/'configuration.json',HERE/'exact_certificate.json',HERE/'exact_result.json')
    print(json.dumps(result,indent=2))
    print('PASS: exact fixed-configuration radius, primitive irreducible quadratic, positive root isolation, full closed disk coverage, and unique worst point.')
    if result['external_companion_interval_strictly_smaller']:
        print('Exact interval comparison PASS; strict nonglobal conclusion additionally uses the separately accepted companion full-disk covering proof.')

if __name__=='__main__':main()
