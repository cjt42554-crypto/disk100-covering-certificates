"""Rejection checks for deliberately corrupted local proof witnesses."""
import sys,json,gzip,tempfile,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import exact_krawczyk as k
import exact_algebraic_cover as c

FILES=['polynomial_system.json','root_high_precision.json',
       'krawczyk_preconditioner.json.gz','algebraic_cover_topology.json','source_configuration.json']

def run(name,mutate,cover):
    with tempfile.TemporaryDirectory(prefix='disk100_negative_') as td:
        tmp=Path(td)
        for f in FILES:shutil.copyfile(ROOT/f,tmp/f)
        k.ROOT=tmp;c.ROOT=tmp;mutate(tmp)
        try:
            if cover:c.main()
            else:k.verify(generate_requested=False,emit=False)
        except ValueError:
            print('PASS rejected:',name,flush=True)
        else:raise RuntimeError('Corrupt witness accepted: '+name)

def alter_json(tmp,name,change):
    p=tmp/name;obj=json.loads(p.read_text());change(obj);p.write_text(json.dumps(obj))

def bad_inverse(tmp):
    p=tmp/'krawczyk_preconditioner.json.gz';obj=json.loads(gzip.decompress(p.read_bytes()))
    obj['matrix'][0][0]+=obj['denominator']*10**10
    p.write_bytes(gzip.compress(json.dumps(obj).encode(),mtime=0))

def bad_root(tmp):
    alter_json(tmp,'root_high_precision.json',lambda o:o['variables'].__setitem__(0,'0.08410828486037'))
def bad_system(tmp):
    alter_json(tmp,'polynomial_system.json',lambda o:o['polynomials'][0][0].__setitem__(0,2))
def bad_degree(tmp):
    alter_json(tmp,'polynomial_system.json',lambda o:o['polynomials'][0][0][1].append(0))
def missing_anchor(tmp):
    alter_json(tmp,'polynomial_system.json',lambda o:(o['unused_center_anchors'].pop() if o['unused_center_anchors'] else o['unused_center_anchors'].append({'center':0,'dimension':0,'value':[0,1]})))
def bad_dimension(tmp):
    alter_json(tmp,'root_high_precision.json',lambda o:o['variables'].pop())

def missing_face(tmp):
    alter_json(tmp,'algebraic_cover_topology.json',lambda o:o['triangles'].pop())
def bad_arc(tmp):
    def change(o):o['boundary_order'][0],o['boundary_order'][10]=o['boundary_order'][10],o['boundary_order'][0]
    alter_json(tmp,'algebraic_cover_topology.json',change)

if __name__=='__main__':
    for name,fn,cover in [('inverse perturbation',bad_inverse,False),('midpoint perturbation',bad_root,False),
        ('altered geometric equation',bad_system,True),('missing Delaunay face',missing_face,True),
        ('broken boundary arc order',bad_arc,True),('cubic malformed system',bad_degree,False),
        ('corrupted unused-center anchors',missing_anchor,True),('wrong root dimension',bad_dimension,False)]:run(name,fn,cover)
