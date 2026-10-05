# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Integer target-face audit for sealed seeds; leaves the handoff untouched."""
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[n]='1'
import json,time,hashlib,sys
from pathlib import Path
import numpy as np
from exact_current_dt import orient,dot
ROOT=Path(__file__).resolve().parent
FREEZE=ROOT/'frozen_handoff'
if not FREEZE.exists():FREEZE=ROOT.parent/'frozen_handoff'


def main():
    tic=time.monotonic();ledger=json.loads((FREEZE/'seed_ledger.json').read_text());cat=json.loads((FREEZE/'catalog.json').read_text())
    graphs={g['independent_index']:g for g in cat['graphs']};records=[]
    for row in ledger['records']:
        path=FREEZE/row['coordinate_relative_file'];assert hashlib.sha256(path.read_bytes()).hexdigest()==row['file_sha256']
        C=np.load(path);rr=[[float(v).as_integer_ratio() for v in c] for c in C];Q=max(d for c in rr for n,d in c)
        P=[tuple(n*(Q//d) for n,d in c) for c in rr]
        negative=[];ties=[];inverted=[];outside=[];checks=0
        graph=graphs[row['independent_index']]
        for fi,f in enumerate(graph['faces']):
            a,b,c=(P[j] for j in f);D=2*orient(a,b,c)
            if D<=0:inverted.append(fi);continue
            u=(b[0]-a[0],b[1]-a[1]);v=(c[0]-a[0],c[1]-a[1]);uu=dot(u,u);vv=dot(v,v)
            pr=(uu*v[1]-vv*u[1],u[0]*vv-v[0]*uu);N=(a[0]*D+pr[0],a[1]*D+pr[1]);r2=dot(pr,pr)
            if dot(N,N)>=(Q*D)**2:outside.append(fi)
            for j,p in enumerate(P):
                if j in f:continue
                dp=(p[0]*D-N[0],p[1]*D-N[1]);sg=dot(dp,dp)-r2;checks+=1
                if sg<0:negative.append([fi,j])
                elif sg==0:ties.append([fi,j])
        exact_pass=not negative and not ties and not inverted and not outside and all(dot(p,p)<Q*Q for p in P)
        records.append(dict(independent_index=row['independent_index'],coordinate_file=str(path),file_sha256=row['file_sha256'],
                            numerical_target_match=row['actual_target_fingerprint_equal'],integer_target_faces_strict_pass=exact_pass,
                            negative_emptycircle_signs=negative,exact_extra_site_ties=ties,nonpositive_oriented_faces=inverted,
                            nonstrict_inside_circumcenter_faces=outside,integer_emptycircle_checks=checks))
    out=dict(status='COMPLETE',scope='Exact target-face emptiness, orientation, inside and center-unitdisk checks on binary64 harmonic proposals. No coverage-radius certificate.',
             handoff_manifest_sha256=hashlib.sha256((FREEZE/'MANIFEST_SHA256.json').read_bytes()).hexdigest(),
             exact_strict_pass=sum(r['integer_target_faces_strict_pass'] for r in records),count=len(records),seconds=time.monotonic()-tic,records=records)
    target=Path(sys.argv[1]) if len(sys.argv)>1 else ROOT/'exact_seed_face_audit.json'
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='records'},indent=2))
    print('NUMERICAL_MATCH_BUT_EXACT_FAILURE',[(r['independent_index'],Path(r['coordinate_file']).name,len(r['negative_emptycircle_signs']),len(r['exact_extra_site_ties'])) for r in records if r['numerical_target_match'] and not r['integer_target_faces_strict_pass']])

if __name__=='__main__':main()
