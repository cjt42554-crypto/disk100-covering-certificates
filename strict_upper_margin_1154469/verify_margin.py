#!/usr/bin/env python3
"""Pure exact covering + positive rational-margin acceptance; stdlib only."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from critical_candidates import require,verify

EXPECTED_CONFIG_SHA='84fb2a5dd9869e641f18250b8839ceb44bfe1263fa54bf588cfdb4867fb28a36'
EXPECTED_CHECKER_SHA='93b202d1893f5f3aaa2ab013c2d6a623ffc182a471a6ae255d54a3cb1d4d9c3e'
EXPECTED_V=[1154469548855829,10**16]
EXPECTED_U=[115446954885583,10**15]

def rational_pair(value):
    require(type(value) is list and len(value)==2 and all(type(x) is int for x in value) and value[1]>0,'invalid integer rational encoding')
    return F(*value)

def verify_files(config_path,certificate_path,ledger_path):
    config_path=Path(config_path);certificate_path=Path(certificate_path);ledger_path=Path(ledger_path)
    require(hashlib.sha256(config_path.read_bytes()).hexdigest()==EXPECTED_CONFIG_SHA,'frozen configuration identity mismatch')
    require(hashlib.sha256((HERE/'critical_candidates.py').read_bytes()).hexdigest()==EXPECTED_CHECKER_SHA,'source math checker identity mismatch')
    ledger=json.loads(ledger_path.read_text());cert=json.loads(certificate_path.read_text())
    require(set(cert)=={'bound_squared','covering_indices'},'unexpected coverage certificate fields')
    V=rational_pair(ledger['local_cover_radius_bound'])
    U=rational_pair(ledger['reference_radius_bound'])
    squared=rational_pair(cert['bound_squared'])
    require(V>0 and U>0,'nonpositive radius')
    require(V==F(*EXPECTED_V),'unexpected local covering radius')
    require(V*V==squared,'explicit local radius squared binding mismatch')
    require(ledger['source_configuration_sha256']==EXPECTED_CONFIG_SHA,'ledger source configuration mismatch')
    # configuration.json is copied byte-for-byte. Its old radius_bound is
    # provenance only; the explicit V^2 above is the actual acceptance bound.
    covered=verify(config_path,certificate_path)
    vn,vd=ledger['local_cover_radius_bound'];un,ud=ledger['reference_radius_bound']
    require(vn*ud<un*vd,'strict rational margin inequality failed')
    require(U==F(*EXPECTED_U),'wrong reference radius')
    margin=rational_pair(ledger['strict_radius_margin'])
    require(margin==U-V and margin>0,'incorrect positive rational margin')
    require(ledger['configuration_radius_bound_is_metadata_only'] is True,'old configuration radius must be explicitly excluded from acceptance')
    return {'covering_exact_verified':True,'strict_margin_exact_verified':True,
            'local_cover_radius_bound':[V.numerator,V.denominator],
            'reference_radius_bound':[U.numerator,U.denominator],
            'strict_radius_margin':[margin.numerator,margin.denominator],
            'configuration_sha256':EXPECTED_CONFIG_SHA,
            'candidate_counts':covered['candidate_counts'],'total':covered['total'],
            'actual_covering_squared_bound':covered['upper_squared'],
            'statement':'The explicit frozen 100-center rational configuration covers the whole closed unit disk at V, and V<U. Its covering radius, and consequently the arbitrary-configuration optimum, is strictly below U.',
            'cloud_statement_not_changed':'The separate cloud result remains its reported non-strict <= bound.'}

def main():
    print(json.dumps(verify_files(HERE/'configuration.json',HERE/'cover_certificate.json',HERE/'margin_ledger.json'),indent=2))
    print('PASS: exact whole closed disk cover at V and pure integer strict comparison V<U.')

if __name__=='__main__':main()
