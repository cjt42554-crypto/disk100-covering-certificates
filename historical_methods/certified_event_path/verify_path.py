# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Self-contained standard-library exact replay. Works with -S and -O."""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json
from critical_candidates import verify as verify_coverage
from exact_path import (require, coords, primitive, incircle, root_count, peval,
                        circumcenter, orient, interval, sq, add, sub, scale, mul,
                        distance_difference, divrem, displacement, cap, derivative, triangle_witness)
HERE=Path(__file__).resolve().parent

def read(path):return json.loads(path.read_text())
def rational(value):
    require(isinstance(value,list) and len(value)==2 and all(type(x) is int for x in value),'bad rational pair')
    require(value[1]>0,'nonpositive denominator')
    return F(*value)

def box_pair(box):return [rational(v) for v in box]

def verify_geometry(G,A,B,Q):
    sites=G['sites_in_incircle_order'];cyclic=G['convex_cyclic_order']
    require(sites==[65,30,64,27] and cyclic==[65,64,30,27],'unexpected event sites/order')
    CP=coords(A,B);P=primitive(incircle([CP[i] for i in sites]))
    require(len(P)<=5,'event polynomial degree exceeds four')
    require(P==G['incircle_primitive_coefficients_low_to_high'],'event coefficient mismatch')
    lo,hi=[rational(v) for v in G['root_interval']]
    require(0<lo<hi<1,'root box must lie strictly inside the segment')
    require(peval(P,F(0))*peval(P,F(1))<0,'event endpoint signs do not differ')
    require(peval(P,lo)*peval(P,hi)<0,'root box does not cross the event')
    require(root_count(P,F(0),F(1))==1 and G['root_count_0_1']==1,'wrong segment root count')
    require(root_count(P,lo,hi)==1,'root box does not isolate one root')
    nx,ny,den=circumcenter(*(CP[i] for i in sites[:3]),Q)
    function={'x_numerator_coefficients':nx,'y_numerator_coefficients':ny,'denominator_coefficients':den}
    require(function==G['circumcenter'],'circumcenter rational function mismatch')
    evidence=G['interval_evidence']
    derivative_box=interval(derivative(P),lo,hi)
    require(list(derivative_box)==box_pair(evidence['incircle_derivative']) and derivative_box[0]>0,'event root is not simple/transverse')
    actual_den=interval(den,lo,hi)
    require(list(actual_den)==box_pair(evidence['denominator']) and actual_den[0]>0,'degenerate original triangle/denominator')
    actual_inside=interval(sub(sq(den),add(sq(nx),sq(ny))),lo,hi)
    require(list(actual_inside)==box_pair(evidence['inside_disk_numerator']) and actual_inside[0]>0,'event circumcenter not proved inside disk')
    for j,a in enumerate(cyclic):
        b=cyclic[(j+1)%4];c=cyclic[(j+2)%4]
        obox=interval(orient(CP[a],CP[b],CP[c]),lo,hi)
        require(list(obox)==box_pair(evidence['cyclic_orientation_bounds'][j]) and obox[1]<0,'quadrilateral not strictly convex')
        dx,dy=[sub(x,y) for x,y in zip(CP[b],CP[a])]
        qdx=sub(scale(nx,Q),mul(den,CP[a][0]));qdy=sub(scale(ny,Q),mul(den,CP[a][1]))
        qbox=interval(sub(mul(dx,qdy),mul(dy,qdx)),lo,hi)
        require(list(qbox)==box_pair(evidence['q_inside_convex_quad_edge_bounds'][j]) and qbox[1]<0,'hole point not strictly inside quadrilateral')
    for i in sites[1:3]:
        require(distance_difference(CP[sites[0]],CP[i],nx,ny,den,Q)==[0],'three-site equal-distance identity failed')
    H=distance_difference(CP[sites[0]],CP[sites[3]],nx,ny,den,Q)
    ratio,rem=divrem(H,P)
    require(rem==[0] and ratio==[rational(v) for v in G['fourth_distance_difference_multiple_of_P']],'fourth site is not tied at the event root')
    require(len(ratio)==1 and ratio[0]!=0,'invalid fourth-site crossing multiplicity')
    radius_num=add(sq(sub(scale(nx,Q),mul(CP[sites[0]][0],den))),
                   sq(sub(scale(ny,Q),mul(CP[sites[0]][1],den))))
    require(G['common_radius_squared']=={'numerator_coefficients':radius_num,'denominator_coefficients':scale(sq(den),Q*Q)},'common radius expression mismatch')
    gap_entries=evidence['strict_empty_gap_bounds']
    others=[i for i in range(100) if i not in sites]
    require([x['site'] for x in gap_entries]==others,'missing/repeated other-site emptiness checks')
    for entry in gap_entries:
        gap=interval(distance_difference(CP[sites[0]],CP[entry['site']],nx,ny,den,Q),lo,hi)
        require(list(gap)==[rational(entry['lower']),rational(entry['upper'])] and gap[0]>0,'common circumcircle not strictly empty')
    bracket={'before':[triangle_witness(CP,lo,t,Q) for t in [(65,30,64),(65,30,27)]],
             'after':[triangle_witness(CP,hi,t,Q) for t in [(65,64,27),(30,64,27)]]}
    require(G['bracket_delaunay_triangles']==bracket,'incorrect bracket Delaunay triangles')
    return {'degree':len(P)-1,'roots_on_open_segment':1,'strict_empty_other_sites':len(others),
            'inside_closed_disk':True,'inside_convex_quadrilateral':True,'root_interval':G['root_interval']}

def verify_path(directory=HERE,check_manifest=True,emit=True):
    directory=Path(directory)
    if check_manifest:
        manifest=read(directory/'MANIFEST_SHA256.json')
        require(type(manifest) is dict and manifest,'empty manifest')
        for filename,digest in manifest.items():
            require(not Path(filename).is_absolute() and '..' not in Path(filename).parts,'unsafe manifest path')
            require(hashlib.sha256((directory/filename).read_bytes()).hexdigest()==digest,f'file SHA mismatch: {filename}')
    path=read(directory/'path.json');Q=path['coordinate_denominator']
    require(type(Q) is int and Q==10**15,'unexpected coordinate denominator')
    nodes=path['nodes'];require(2<=len(nodes)<=30,'node budget violated')
    indices=[n['index'] for n in nodes]
    require(len(set(indices))==len(indices),'repeated node')
    if path['build_mode']=='all accepted historical nodes':require(indices==[-1]+list(range(22)),'full path must include the previously certified baseline and all accepted nodes in source order')
    else:require(indices==[18,19] and path['build_mode']=='short event segment','unexpected short path')
    C={};U={};counts=[]
    for n in nodes:
        i=n['index'];require(n['directory']==f'nodes/{i:02d}','bad node directory')
        folder=directory/n['directory'];cfg=read(folder/'configuration.json');cert=read(folder/'upper_certificate.json')
        require(cfg['coordinate_denominator']==Q,'node denominator mismatch')
        C[i]=cfg['centers_numerators'];U[i]=rational(cfg['radius_bound'])
        require(U[i]==rational(n['radius_bound']) and U[i]*U[i]==rational(cert['bound_squared']),'radius binding mismatch')
        require(all(x*x+y*y<Q*Q for x,y in C[i]),'center outside unit disk')
        source=cfg['source'];require(source['snapshot']==n['snapshot'],'source snapshot binding mismatch')
        if i==-1:
            require(source['kind']=='previously certified rational configuration' and
                    source['snapshot']=='certified_candidates/U1155018241946/configuration.json','bad baseline provenance')
            require(U[i]==F(577509120973,5000000000000),'wrong previously certified baseline bound')
            old=dict(cfg);old.pop('source')
            old_digest=hashlib.sha256((json.dumps(old,indent=2)+'\n').encode()).hexdigest()
            require(old_digest==source['sha256']=='024276d929e0325c6ad7447bb55779bcb17da8aa33e8463cd8fa9e4869415461','old baseline configuration identity mismatch')
            require(hashlib.sha256((folder/'upper_certificate.json').read_bytes()).hexdigest()==source['certificate_sha256']=='014144d0458126fd5dddd0145c39dc7dc7f6458d9db1d525795cd7f74a2cf146','old baseline certificate identity mismatch')
        else:
            require(source['rounding']=='nearest integer, exact rational ties to even','wrong rounding rule')
            raw=source['binary64_rationals']
            require(len(raw)==100 and all(len(c)==2 for c in raw),'missing source coordinates')
            for center,original in zip(C[i],raw):
                require([round(rational(x)*Q) for x in original]==center,'exact source rounding mismatch')
        result=verify_coverage(folder/'configuration.json',folder/'upper_certificate.json')
        counts.append(result['total'])
        if emit:print(json.dumps({'node':i,'U':str(U[i]),'closed_disk_candidates':result['total'],'exact_cover':'PASS'}),flush=True)
    segments=path['segments'];require(len(segments)==len(nodes)-1,'wrong segment count')
    caps=[]
    for segment,left,right in zip(segments,indices,indices[1:]):
        require((segment['left'],segment['right'])==(left,right),'wrong segment adjacency')
        D,sqmax=displacement(C[left],C[right],Q)
        require(segment['maximum_squared_integer_displacement']==sqmax,'wrong maximum displacement square')
        require(rational(segment['displacement_bound'])==D,'wrong upward integer square-root displacement bound')
        upper=cap(U[left],U[right],D)
        require(rational(segment['coverage_cap'])==upper,'wrong continuous segment coverage cap')
        caps.append(upper)
    require(rational(path['path_coverage_cap'])==max(caps),'wrong path coverage cap')
    require(path['event_segment']==[18,19],'wrong event segment')
    require(read(directory/'short_segment.json')==next(s for s in segments if s['left']==18 and s['right']==19),'short-segment receipt mismatch')
    geometry=verify_geometry(read(directory/'event_geometry.json'),C[18],C[19],Q)
    result={'nodes':len(nodes),'segments':len(segments),'total_exact_candidates':sum(counts),
            'path_radius_cap':str(max(caps)),'event':geometry,
            'claim':'Explicit rational PL path covers the whole closed unit disk at the given cap and crosses a certified empty convex inside-D cocircular event. Global optimality and the original numerical curve are not certified.'}
    if emit:print(json.dumps(result,indent=2));print('PASS: exact closed-disk path and algebraic flip certificate.',flush=True)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--directory',type=Path,default=HERE);args=parser.parse_args()
    verify_path(args.directory)
