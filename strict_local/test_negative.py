"""Meaningful corruption tests for the exact certificate acceptance path."""
import copy,gzip,json,tempfile,shutil,subprocess,sys
from pathlib import Path
from fractions import Fraction as Q
from verify_sosc import verify_second_order
ROOT=Path(__file__).resolve().parent

def main():
    system=json.loads((ROOT/'source/polynomial_system.json').read_text())
    mid=json.loads((ROOT/'source/rational_midpoint.json').read_text())
    witness=json.loads(gzip.decompress((ROOT/'sosc_witness.json.gz').read_bytes()))
    rho=Q(1,10**58);results=[]
    def rejection(name,w= None,mm=None):
        try:verify_second_order(system,mid if mm is None else mm,witness if w is None else w,rho)
        except (ValueError,KeyError,TypeError) as exc:
            results.append({'test':name,'rejected':True,'reason':str(exc)});return
        raise RuntimeError('Corruption was accepted: '+name)
    w=copy.deepcopy(witness);w['inverse_matrix'][0]=[0]*system['n_constraints'];rejection('zero inverse row',w)
    w=copy.deepcopy(witness);w['cholesky_transform']=[[0]*54 for _ in range(54)];rejection('zero Cholesky transform',w)
    w=copy.deepcopy(witness);w['pivot_columns'][1]=w['pivot_columns'][0];rejection('duplicate pivot column',w)
    w=copy.deepcopy(witness);w['inverse_denominator']+=1;rejection('non-dyadic inverse denominator',w)
    w=copy.deepcopy(witness);w['n_base']+=1;rejection('wrong witness dimension',w)
    mm=copy.deepcopy(mid)
    for i in range(system['n_base'],system['n_variables']):mm['numerators'][i]=str(-int(mm['numerators'][i]))
    rejection('negated Lagrangian Hessian',mm=mm)
    with tempfile.TemporaryDirectory(prefix='sosc-negative-') as temp:
        dest=Path(temp)/'certificate';shutil.copytree(ROOT,dest,ignore=shutil.ignore_patterns('__pycache__','portable_certificate.tar.gz'))
        w=copy.deepcopy(witness);w['source_system_sha256']='0'*64
        (dest/'sosc_witness.json.gz').write_bytes(gzip.compress(json.dumps(w).encode(),mtime=0))
        result=subprocess.run([sys.executable,'-S','-O',str(dest/'verify_sosc.py')],capture_output=True,text=True)
        require_rejection(result,'witness source SHA mismatch',results)
        shutil.copyfile(ROOT/'sosc_witness.json.gz',dest/'sosc_witness.json.gz')
        source=dest/'source/polynomial_system.json';data=json.loads(source.read_text());data['polynomials'][0][0][0]+=1;source.write_text(json.dumps(data))
        result=subprocess.run([sys.executable,'-S','-O',str(dest/'verify_sosc.py')],capture_output=True,text=True)
        require_rejection(result,'altered geometric polynomial source',results)
    report={'all_rejected':True,'number_tests':len(results),'results':results}
    (ROOT/'negative_tests.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))

def require_rejection(result,name,results):
    if result.returncode==0:raise RuntimeError('Corruption accepted: '+name)
    payload=json.loads(result.stdout)
    if payload.get('verified') is not False:raise RuntimeError('Expected explicit rejection: '+name)
    results.append({'test':name,'rejected':True,'reason':payload['error']})

if __name__=='__main__':main()
