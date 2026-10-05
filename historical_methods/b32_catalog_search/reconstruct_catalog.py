# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
"""Independent combinatorial reconstruction from supplied strip regions.

The source regions come from the older N26/F106 catalog.  This result makes no
completeness claim; expected 16 count/family names must be checked against the
parent's exact catalog when it arrives.
"""
import sys,json,hashlib
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from graph_tools import canonical_disk
ROOT=Path(__file__).resolve().parent
SRC=ROOT.parent/'github_archive_staging/repository/conditional/metric_family_11316/lower_bounds_extension_20261004'
if not SRC.exists():SRC=ROOT/'source_lattice'
sys.dont_write_bytecode=True
sys.path.insert(0,str(SRC))
from lattice_mesh_utils import region_faces,canonical,mesh


def deleted_states(faces,neighbors,bcounts,count):
    allmask=(1<<len(faces))-1
    boundary=sum(1<<j for j,c in enumerate(bcounts) if c)
    states={0}
    sizes=[1]
    for level in range(count):
        nxt=set()
        for removed in states:
            choices=boundary
            rm=removed
            while rm:
                bit=rm&-rm;rm-=bit;j=bit.bit_length()-1
                choices|=neighbors[j]
            choices&=allmask^removed
            while choices:
                bit=choices&-choices;choices-=bit
                nxt.add(removed|bit)
        states=nxt;sizes.append(len(states))
    return states,sizes


def full_graph(core):
    ring=core['boundary_ids'];p=[4-q for q in core['boundary_face_incidences']]
    assert len(ring)==26 and sum(p)==32 and max(p)<=2 and min(p)>=0
    connectors=[];q1=[];offset=0
    for j,v in enumerate(ring):
        chain=[68+(offset+k)%32 for k in range(p[j]+1)]
        q1.extend([[a,b,v] for a,b in zip(chain,chain[1:])])
        connectors.append([v,chain[-1],ring[(j+1)%26]])
        offset+=p[j]
    faces=core['triangles']+connectors+q1
    ca=canonical_disk(faces)
    return dict(core_q0_triangles=core['triangles'],ring_q0_triangles=connectors,q1_triangles=q1,
                faces=faces,outer_boundary_ids=list(range(68,100)),interior_ids=list(range(68)),ring_ids=ring,
                q1_fan_sizes_in_ring_order=p,fingerprint=ca['fingerprint'],vertex_to_canonical=ca['vertex_to_canonical'])


def main():
    data=json.loads((SRC/'flat_N26_core_catalog.json').read_text())
    cores={};stats=[]
    for ri,row in enumerate(data['containing_regions']):
        count=row['triangle_count']-108
        if count<0:continue
        faces,neighbors,bcounts,boundary=region_faces(row)
        states,sizes=deleted_states(faces,neighbors,bcounts,count)
        valid=0
        for removed in states:
            # Count the remaining boundary exactly before the stronger mesh checks.
            ep=defaultdict(int)
            remaining=[]
            for j,f in enumerate(faces):
                if (removed>>j)&1:continue
                remaining.append(f)
                for k in range(3):ep[tuple(sorted((f[k],f[(k+1)%3])))]+=1
            if sum(v==1 for v in ep.values())!=26:continue
            key,oriented=canonical(remaining)
            if key in cores:continue
            core=mesh(oriented,None,68,108,26,{2,3,4})
            if core is not None:
                core.update(source_region=ri,removed_face_indices=[j for j in range(len(faces)) if (removed>>j)&1])
                cores[key]=core;valid+=1
        stats.append(dict(region=ri,triangle_count=row['triangle_count'],delete=count,state_counts=sizes,new_valid_classes=valid))
    graphs=[]
    for j,(key,core) in enumerate(sorted(cores.items())):
        g=full_graph(core)
        g.update(independent_index=j,core=core)
        graphs.append(g)
    assert len({g['fingerprint'] for g in graphs})==len(graphs),'nonlattice_isomorphic_duplicates'
    out=dict(status='PASS',scope='Independent N68/F108/P26 reconstruction from older supplied containing regions; no completeness claim',
             source_regions_path=str(SRC/'flat_N26_core_catalog.json'),source_regions_sha256=hashlib.sha256((SRC/'flat_N26_core_catalog.json').read_bytes()).hexdigest(),
             graph_count=len(graphs),statistics=stats,graphs=graphs)
    (ROOT/'independent_catalog_numerical_index.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='graphs'},indent=2))
    current=json.loads((ROOT/'current_graph_numerical.json').read_text())
    print('CURRENT_MATCHES',[(g['independent_index'],g['fingerprint']) for g in graphs if g['fingerprint']==current['fingerprint']])

if __name__=='__main__':main()
