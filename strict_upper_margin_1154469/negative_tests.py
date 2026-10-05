#!/usr/bin/env python3
"""Reject truncated coverage, an insufficient radius, and a non-strict gap."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from critical_candidates import verify
from verify_margin import verify_files

def main():
    cert=json.loads((HERE/'cover_certificate.json').read_text())
    ledger=json.loads((HERE/'margin_ledger.json').read_text())
    shortened=copy.deepcopy(cert);shortened['covering_indices'].pop()
    smaller=copy.deepcopy(cert)
    too_small=F(*ledger['local_cover_radius_bound'])-F(2,10**16)
    squared=too_small*too_small
    smaller['bound_squared']=[squared.numerator,squared.denominator]
    equality=copy.deepcopy(ledger)
    equality['reference_radius_bound']=equality['local_cover_radius_bound'][:]
    equality['strict_radius_margin']=[0,1]
    cases=[('truncated_covering_witness_list',shortened,ledger,'cover','truncated certificate'),
           ('radius_reduced_by_2e-16',smaller,ledger,'cover','uncovered candidate'),
           ('reference_bound_equal_to_V',cert,equality,'margin','strict rational margin inequality failed')]
    result=[]
    with tempfile.TemporaryDirectory(prefix='negative_',dir=HERE) as directory:
        cp=Path(directory)/'cover_certificate.json';lp=Path(directory)/'margin_ledger.json'
        for name,c,l,mode,expected in cases:
            cp.write_text(json.dumps(c));lp.write_text(json.dumps(l))
            try:
                if mode=='cover':verify(HERE/'configuration.json',cp)
                else:verify_files(HERE/'configuration.json',cp,lp)
            except ValueError as error:
                if expected not in str(error):raise RuntimeError(f'{name}: wrong rejection: {error}')
                result.append({'case':name,'rejected':True,'reason':str(error)})
            else:raise RuntimeError(f'{name}: corrupted mathematical claim accepted')
    print(json.dumps({'all_rejected':True,'cases':result},indent=2))

if __name__=='__main__':main()
