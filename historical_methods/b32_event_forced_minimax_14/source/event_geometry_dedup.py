# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Finite old-six and inter-pilot physical geometry comparison, no optimizer."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,hashlib,itertools
import numpy as np
ROOT=Path(__file__).resolve().parent.parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def match(P,Q):
    best=None
    for perm in itertools.permutations(range(4)):
        for reflect in (False,True):
            A=P.copy();A[:,1]*=-1 if reflect else 1
            U,s,V=np.linalg.svd(A.T@Q[list(perm)]);S=np.eye(2);S[-1,-1]=np.linalg.det(U@V);O=U@S@V
            err=np.linalg.norm(A@O-Q[list(perm)],axis=1)
            row=dict(maximum_support_distance=float(err.max()),rms_support_distance=float(np.sqrt(np.mean(err**2))),permutation=list(perm),reflection=reflect,origin_fixed_rotation=O.tolist())
            if best is None or row['maximum_support_distance']<best['maximum_support_distance']:best=row
    return best
def main():
    h=json.loads((ROOT/'PHYSICAL_EVENT_HANDOFF.json').read_text());oldpath=ROOT.parent/'event_geometry_146/FINAL_BASIN_HANDOFF.json';old=json.loads(oldpath.read_text());records=[];new=[]
    for row in h['events']:
        tie=next(r for r in row['frames'] if r['fraction']==1);C=np.load(ROOT/tie['file']);assert sha(ROOT/tie['file'])==tie['file_sha256'];P=C[row['labels']];new.append((row,C,P))
        comparisons=[]
        for r in old:
            p=Path(r['event_npz']);Z=np.load(p);O=Z['C'];Q=O[r['event_quad']];local=match(P,Q);radial=float(np.max(np.abs(np.sort(np.linalg.norm(C,axis=1))-np.sort(np.linalg.norm(O,axis=1)))))
            comparisons.append(dict(old_route_id=r['event_id'],old_quad=r['event_quad'],old_event_npz=str(p),old_event_npz_sha256=sha(p),local_quad=local,all100_sorted_radius_assignment_lower_bound=radial,no_all100_origin_rotation_reflection_match_at2e6=radial>2e-6,local_quad_match_at2e6=local['maximum_support_distance']<2e-6))
        records.append(dict(target_id=row['target_id'],tie_file=tie['file'],tie_file_sha256=tie['file_sha256'],comparisons=comparisons,no_old6_local_quad_geometry_match=not any(r['local_quad_match_at2e6'] for r in comparisons)))
    pair=match(new[0][2],new[1][2]);pair['no_interpilot_quad_geometry_match_at2e6']=pair['maximum_support_distance']>2e-6
    report=dict(status='PASS' if all(r['no_old6_local_quad_geometry_match'] for r in records) and pair['no_interpilot_quad_geometry_match_at2e6'] else 'FAIL',scope='Actual floating event-frame local4 support compared all24 assignments, origin-fixed rotation/reflection; all100 radial multiset lower bound. Finite old6 only, no global historical absence claim.',prior_handoff_sha256=sha(oldpath),records=records,interpilot_comparison=pair)
    p=ROOT/'audit/PHYSICAL_GEOMETRY_DEDUP.json';p.write_text(json.dumps(report,indent=2)+'\n');h['physical_geometry_dedup_status']=report['status'];h['physical_geometry_dedup_file']='audit/PHYSICAL_GEOMETRY_DEDUP.json';h['physical_geometry_dedup_sha256']=sha(p);(ROOT/'PHYSICAL_EVENT_HANDOFF.json').write_text(json.dumps(h,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],nearest_old6_localquad_residual=min(x['local_quad']['maximum_support_distance'] for r in records for x in r['comparisons']),interpilot=pair['maximum_support_distance']),indent=2))
if __name__=='__main__':main()
