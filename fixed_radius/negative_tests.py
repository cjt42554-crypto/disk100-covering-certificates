#!/usr/bin/env python3
"""Reject a false polynomial, a false exact radius, and a false root box."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from verify_exact import verify_files

def main():
    cert=json.loads((HERE/'exact_certificate.json').read_text())
    detail=json.loads((HERE/'exact_result.json').read_text())
    cases=[]
    false_poly=copy.deepcopy(detail)
    false_poly['polynomial_coefficients_descending'][2]+=1
    cases.append(('constant_polynomial_coefficient_plus_one',cert,false_poly,'polynomial coefficients'))
    false_radius=copy.deepcopy(cert);false_detail=copy.deepcopy(detail)
    false_radius['exact_radius_squared'][0]+=1
    false_detail['radius_squared']=false_radius['exact_radius_squared'][:]
    false_detail['polynomial_coefficients_descending'][2]-=1
    cases.append(('exact_radius_squared_numerator_plus_one',false_radius,false_detail,'wrong exact radius'))
    wrong_box=copy.deepcopy(cert);wrong_detail=copy.deepcopy(detail)
    grid=wrong_detail['isolating_interval_common_denominator']
    wrong_detail['isolating_interval_common_numerators']=[n+2 for n in wrong_detail['isolating_interval_common_numerators']]
    box=[F(n,grid) for n in wrong_detail['isolating_interval_common_numerators']]
    wrong_box['root_interval']=[[x.numerator,x.denominator] for x in box]
    wrong_detail['isolating_interval']=wrong_box['root_interval']
    cases.append(('root_interval_shifted_two_grid_units',wrong_box,wrong_detail,'invalid positive root isolation'))
    results=[]
    with tempfile.TemporaryDirectory(prefix='negative_',dir=HERE) as directory:
        cp=Path(directory)/'exact_certificate.json';dp=Path(directory)/'exact_result.json'
        for name,c,d,expected in cases:
            cp.write_text(json.dumps(c));dp.write_text(json.dumps(d))
            try:verify_files(HERE/'configuration.json',cp,dp)
            except ValueError as error:
                if expected not in str(error):raise RuntimeError(f'{name}: unexpected rejection: {error}')
                results.append({'case':name,'rejected':True,'reason':str(error)})
            else:raise RuntimeError(f'{name}: corrupted exact claim accepted')
    print(json.dumps({'all_rejected':True,'cases':results},indent=2))

if __name__=='__main__':main()
