# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Read-only old metric hash inventory, including all local B32 profile arrays."""
import engine as run
from pathlib import Path
import json
import numpy as np
ROOT=run.ROOT;BASE=ROOT.parent;PREVIOUS=BASE/'b32_metric_seed_search'

def unordered_hash(C):
    C=np.asarray(C,dtype='<f8');return run.csha(C[np.lexsort((C[:,1],C[:,0]))])

def load_record(r):
    p=(BASE/r['path']).resolve()
    if not p.is_relative_to(BASE):raise ValueError('Other-project coordinate read prohibited: '+r['path'])
    field=r['field'];schema=r['schema']
    if p.suffix=='.npy':
        C=np.load(p,allow_pickle=False)
        if'slice_index'in r:C=C[r['slice_index']]
    elif p.suffix=='.npz':
        with np.load(p,allow_pickle=False)as z:
            C=z[field[:-6]][:200].reshape(100,2)if schema=='state_first_200_center_variables'else z[field]
    else:
        d=json.loads(p.read_text())
        if schema=='rational_numerators/denominator':C=np.array([[int(x)/int(d[r['denominator_field']])for x in row]for row in d[field[2:]]])
        elif schema=='quadratic_system_first_200_center_variables':
            C=np.array([float(x)for x in d['variables'][:200]]).reshape(100,2)if'variables'in field else np.array([int(x)/int(d['denominator'])for x in d['numerators'][:200]]).reshape(100,2)
        elif field=='$':C=np.array(d,dtype=float)
        else:raise ValueError('Unsupported saved coordinate field: '+field)
    C=np.asarray(C,dtype=float)
    if C.shape==(200,):C=C.reshape(100,2)
    elif C.shape==(201,):C=C[:200].reshape(100,2)
    if C.shape!=(100,2)or run.csha(C)!=r['coordinate_sha256']:raise ValueError('Old coordinate inventory SHA mismatch')
    return C

def inventory():
    ordered=json.loads((PREVIOUS/'prior_ordered_hashes.json').read_text());unordered=json.loads((PREVIOUS/'prior_unordered_hashes.json').read_text());newcount=0;readchecks=[];errors=[]
    # Permanent541-file metric archive contains all87 new starts,134 stages,
    # 96 raw harmonic proposals, old reference copies, and pinned poses.
    manifest=json.loads((PREVIOUS/'MANIFEST_SHA256.json').read_text())
    for relative,h in manifest.items():
        p=PREVIOUS/relative
        if p.suffix not in('.npy','.npz'):continue
        assert run.sha(p)==h
        try:
            if p.suffix=='.npy':values=[('',np.load(p,allow_pickle=False))]
            else:
                with np.load(p,allow_pickle=False)as z:values=[(k,z[k])for k in z.files]
            for k,C in values:
                if C.shape==(200,):C=C.reshape(100,2)
                elif C.shape==(201,):C=C[:200].reshape(100,2)
                if C.shape!=(100,2):continue
                label='b32_metric_seed_search/'+relative+(':'+k if k else'');ordered.setdefault(run.csha(C),[]).append(label);unordered.setdefault(unordered_hash(C),[]).append(label);newcount+=1
        except Exception as e:errors.append(dict(path=str(p),error=str(e)))
    history=PREVIOUS/'inputs/history_snapshot';records=json.loads((history/'retained_coordinate_records.json').read_text());matches=json.loads((history/'historical_b32_profile_matches.json').read_text())
    # All34868 ordered historical hashes are retained. Read only the local484
    # matching profile geometries for a complete unordered target-graph check;
    # other projects are never opened. Non-B32 arrays cannot be exact unordered
    # copies of an accepted strict B32 target realization.
    for r in records:ordered.setdefault(r['coordinate_sha256'],[]).append('history:'+r['path']+':'+r['field'])
    for m in matches:
        r=records[m['record_index']]
        try:
            C=load_record(r);label='history:'+r['path']+':'+r['field'];unordered.setdefault(unordered_hash(C),[]).append(label);readchecks.append(dict(path=r['path'],field=r['field'],coordinate_sha256=r['coordinate_sha256'],unordered_sha256=unordered_hash(C),matched_index=m['independent_catalog_index']))
        except Exception as e:errors.append(dict(path=r['path'],field=r['field'],error=str(e)))
    return ordered,unordered,dict(previous_new_array_count=newcount,local_historical_profile_reads=len(readchecks),old_ordered_unique=len(ordered),old_unordered_unique=len(unordered),readchecks=readchecks,errors=errors,scope='Old actual traces plus frozen harmonic87; local B32-profile sorted geometry hashes supplement complete old ordered inventory. No other-project array read and no optimizer.')
