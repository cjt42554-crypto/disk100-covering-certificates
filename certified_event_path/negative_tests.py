#!/usr/bin/env python3
"""Tamper tests exercise mathematical checks, explicitly bypassing SHA checks."""
from pathlib import Path
from tempfile import TemporaryDirectory
from fractions import Fraction as F
import copy,json,shutil
from verify_path import HERE,verify_path

def read(p):return json.loads(p.read_text())
def write(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def pair(x):return [x.numerator,x.denominator]

def short_copy(destination):
    path=read(HERE/'path.json')
    path['nodes']=[n for n in path['nodes'] if n['index'] in [18,19]]
    path['segments']=[read(HERE/'short_segment.json')]
    path['path_coverage_cap']=path['segments'][0]['coverage_cap']
    path['build_mode']='short event segment'
    write(destination/'path.json',path)
    for filename in ['event_geometry.json','short_segment.json']:
        shutil.copyfile(HERE/filename,destination/filename)
    for i in [18,19]:
        (destination/'nodes').mkdir(exist_ok=True)
        shutil.copytree(HERE/'nodes'/str(i),destination/'nodes'/str(i))

def alter_polynomial(d):
    p=d/'event_geometry.json';x=read(p);x['incircle_primitive_coefficients_low_to_high'][0]+=1;write(p,x)
def alter_bracket(d):
    p=d/'event_geometry.json';x=read(p);x['root_interval']=[[1,10],[1,5]];write(p,x)
def alter_cap(d):
    p=d/'path.json';x=read(p);x['segments'][0]['coverage_cap']=pair(F(*x['segments'][0]['coverage_cap'])-F(1,10**15));write(p,x)
def alter_empty(d):
    p=d/'event_geometry.json';x=read(p);x['interval_evidence']['strict_empty_gap_bounds'].pop();write(p,x)
def alter_endpoint_bound(d):
    U=F(11,100)
    p=d/'nodes/18/configuration.json';x=read(p);x['radius_bound']=pair(U);write(p,x)
    p=d/'nodes/18/upper_certificate.json';x=read(p);x['bound_squared']=pair(U*U);write(p,x)
    p=d/'path.json';x=read(p);x['nodes'][0]['radius_bound']=pair(U);write(p,x)

def main():
    result=[]
    for name,mutate in [('event_polynomial_tamper',alter_polynomial),('nonisolating_root_bracket',alter_bracket),
                        ('lowered_continuous_coverage_cap',alter_cap),('missing_other_site_empty_check',alter_empty),
                        ('false_endpoint_cover_radius',alter_endpoint_bound)]:
        with TemporaryDirectory(prefix='.negative-',dir=HERE) as td:
            folder=Path(td);short_copy(folder);mutate(folder)
            try:verify_path(folder,check_manifest=False,emit=False)
            except ValueError as e:result.append({'test':name,'status':'REJECTED','reason':str(e)})
            else:raise RuntimeError('tampered mathematical certificate unexpectedly passed: '+name)
    write(HERE/'negative_results.json',result)
    print(json.dumps(result,indent=2));print('PASS: all five mathematical tamper tests rejected.')

if __name__=='__main__':main()
