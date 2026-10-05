#!/usr/bin/env python3
"""Build a small cover certificate without KKT, floats or numerical ranking."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from verify_margin import verify_files,EXPECTED_CONFIG_SHA,EXPECTED_CHECKER_SHA

def main():
    V=F(1154469548855829,10**16);U=F(115446954885583,10**15)
    proposals=json.loads((HERE/'savedsource'/'witness_proposal.json').read_text())
    squared=V*V
    certificate={'bound_squared':[squared.numerator,squared.denominator],
                 'covering_indices':proposals['covering_indices']}
    ledger={'source_configuration_sha256':EXPECTED_CONFIG_SHA,
            'source_checker_sha256':EXPECTED_CHECKER_SHA,
            'local_cover_radius_bound':[V.numerator,V.denominator],
            'reference_radius_bound':[U.numerator,U.denominator],
            'strict_radius_margin':[(U-V).numerator,(U-V).denominator],
            'configuration_radius_bound_is_metadata_only':True,
            'acceptance_bound':'V is explicitly bound to cover_certificate.json bound_squared; no old configuration radius is used.',
            'proof_scope':'Independent local whole-closed-disk cover and strict margin. Geometry is unchanged; no optimum polynomial or material new geometry is claimed.',
            'witness_proposal_sha256':hashlib.sha256((HERE/'savedsource'/'witness_proposal.json').read_bytes()).hexdigest(),
            'witness_proposal_is_not_acceptance':True,
            'cloud_statement_not_changed':'The cloud result remains its own non-strict <= statement; this is a separate locally recomputed strict proof.'}
    (HERE/'cover_certificate.json').write_text(json.dumps(certificate,separators=(',',':'))+'\n')
    (HERE/'margin_ledger.json').write_text(json.dumps(ledger,indent=2)+'\n')
    result=verify_files(HERE/'configuration.json',HERE/'cover_certificate.json',HERE/'margin_ledger.json')
    (HERE/'verification_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
