#!/usr/bin/env python3
"""One corrected lambda14 scalar replay; full JSONMATCH and all-source pins."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import importlib.util
import json
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent
TARGET=RESEARCH/'quadratic_lambda14_conjugate'
PRIMITIVE=RESEARCH/'coupled_edges/verify_support.py'
PRIMITIVE_PIN='90fd6e3426d8c113bacedd55f9a13f5e47981e9a48d21bfef48183ec4602efec'
IDENTITY={
 'PROOF.md':'aaaad3b1629ba98fd44907735b75ae111fecaa45c8d77f078f4e28acf4713ea3',
 'verify_conjugate.py':'e95db0a6ead5ace03884ab54743b6238f8fe4cd242b1563617b789fce9cb57f4',
 'certificate.json':'f325d42bfe651599881b7dc652567729a1586389ee385c38060840bc41696bc2',
 'replay.log':'e1c3535b2023a350366024fb734dd5f936956a8eb21ad9e604d6571af3cafb3a',
 'SHA256.json':'4606ec6e225fa4c2e00da75c1eb1835e721d666125d14a2ed79520b361335b9e',
 'history_factor2_correction/PROOF.md':'266dc36cf80d739dd2f99bcf1703fb77589a88644228f3adde64128b780913e4',
 'history_factor2_correction/verify_conjugate.py':'a5573868e6af596d48eefdd3fc0a15fc9f53285e4a30de42b9614622f8cae960',
 'history_factor2_correction/certificate.json':'f325d42bfe651599881b7dc652567729a1586389ee385c38060840bc41696bc2',
 'history_factor2_correction/replay.log':'12b42809fcf4b3aa142d28d9e5bf9c8182956f670540ab7f19f2e7c4054a862d',
 'history_factor2_correction/SHA256.json':'a8c001b5655bd4eeb42728b9210e184895093f808367eaddb4be2c1dfd57184e',
 'history_factor2_correction/STATUS.json':'f496cc9823ce38f76ba05ee659e50aa9b87ad4c84abe174a82d9303afe70cdaf',
}
def need(ok,message):
 if not ok:raise RuntimeError(message)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def interval(enc):return Q(*enc[0]),Q(*enc[1])
def verify():
 need(__debug__,'optimized Python is not an accepting execution')
 originals={n:(TARGET/n).read_bytes() for n in IDENTITY}
 for n,pin in IDENTITY.items():need(sha(originals[n])==pin,'reviewed source identity '+n)
 manifest=json.loads(originals['SHA256.json'])
 for n,pin in manifest.items():need(IDENTITY[n]==pin,'target manifest consistency '+n)
 status=json.loads(originals['history_factor2_correction/STATUS.json'])
 need(status['status']=='SUPERSEDED_ANALYTIC_GUARD_FACTOR2' and
      status['independent_fresh_calls']==0,'old formula is preserved as superseded, not accepted')
 primitive=PRIMITIVE.read_bytes();need(sha(primitive)==PRIMITIVE_PIN,'frozen arithmetic primitive identity')
 saved=json.loads(originals['certificate.json'])
 spec=importlib.util.spec_from_file_location('new_corrected_lambda14_replay',TARGET/'verify_conjugate.py')
 module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 result=module.verify()
 need(result==saved,'complete corrected fresh dictionary matches saved certificate')
 need(result['lambda_fraction']==[7,50] and result['tilt_fraction']==[1071883,20000000],
      'fixed lambda and positive tilt')
 C=interval(result['conjugate']);Co=interval(result['obtuse_branch']['value'])
 C0=interval(result['zero_tilt_conjugate']);Co0=interval(result['zero_tilt_obtuse'])
 need(C[0]>Co[1] and C0[0]>Co0[1],'both complete branch dominances')
 need(Co0[0]<=Q(25,56)<=Co0[1],'exact zero-tilt obtuse maximum')
 for branch in ['acute_branch','obtuse_branch','zero_tilt_acute']:
  a,b=interval(result[branch]['peak_angle']);lo,hi=interval(result[branch]['value'])
  need(0<a<b and b-a<Q(1,10**24) and 0<hi-lo<Q(1,10**23),'certified unique-peak enclosure '+branch)
 need(result['primitive_sha256']==PRIMITIVE_PIN and result['old_batches_run']==0 and result['optimizers_run']==0,
      'primitive and no historical/search acceptance')
 for n,raw in originals.items():need((TARGET/n).read_bytes()==raw,'target bytes unchanged '+n)
 need(PRIMITIVE.read_bytes()==primitive,'primitive bytes unchanged')
 return dict(status='PASS_INDEPENDENT_CORRECTED_LAMBDA14_CONJUGATE_AUDIT',
  analytic_review='PASS_CORRECTED_SOURCE',initial_factor2_issue_preserved_and_repaired=True,
  independent_fresh_verify_calls=1,complete_saved_JSON_match=True,child_checker_calls=0,
  old_checker_batches_run=0,old_primitive_verify_calls=0,optimizers_run=0,
  target_source_sha256=IDENTITY,primitive_sha256=PRIMITIVE_PIN,source_bytes_unchanged=True,
  lambda_fraction=[7,50],tilt_fraction=[1071883,20000000],
  conjugate=result['conjugate'],zero_tilt_conjugate=result['zero_tilt_conjugate'],
  exact_zero_tilt_obtuse_maximum=[25,56],acute_obtuse_margin=result['acute_obtuse_margin'],
  own_source_sha256={n:sha((HERE/n).read_bytes()) for n in ['REVIEW.md','verify_audit.py']},
  scope='Exact unrestricted cap-anchor full-apex single-face conjugate at lambda=7/50 and fixed positive tilt, plus zero tilt. No actual shared-anchor joint epsilon, arc bootstrap, global covering-radius bound or optimum.')
if __name__=='__main__':
 result=verify();(HERE/'audit_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
 print(result['status'])
 print('One corrected new scalar verify: complete saved JSONMATCH; all current/history source and primitive bytes unchanged.')
 print('Full-apex single-face support only; historical verify calls and optimizers zero.')
