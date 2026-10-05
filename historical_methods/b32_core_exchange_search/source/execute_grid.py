# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""P4/BLAS1 full100 original-problem execution, only after prepared READY grid."""
import engine as run
import os,json,time,traceback,multiprocessing as mp,queue as qm
from pathlib import Path
import numpy as np
from metric_validation import validate
ROOT=run.ROOT

def pidinfo():
    return dict(pid=os.getpid(),namespace_pid_status=[s for s in Path('/proc/self/status').read_text().splitlines()if s.startswith('NSpid:')],affinity=sorted(os.sched_getaffinity(0)))

def worker_main(worker,rows,queue,cpus):
    os.sched_setaffinity(0,{cpus[worker%len(cpus)]});folder=ROOT/'workers'/f'worker{worker:02d}';records=[];tic=time.monotonic();run.write(folder/'runtime.json',dict(worker=worker,BLAS_threads=1,**pidinfo()))
    if(folder/'checkpoint.json').exists():records=json.loads((folder/'checkpoint.json').read_text())['records']
    done={r['grid_id']:r for r in records}
    with open(folder/'records.jsonl','a',buffering=1)as fp:
        for row in rows:
            if row['grid_id']in done:queue.put(dict(type='record',record=done[row['grid_id']]));continue
            tt=time.monotonic();i=row['grid_id']
            try:
                C=np.load(folder/'starts'/f'{i:04d}.npy');R,C,initial,rounds=run.optimize(C,folder,row);p=folder/'finals'/f'{i:04d}.npy';np.save(p,C);R,d=run.true_cover(C);v=validate(C,row['target_graph']);v.pop('actual_graph',None)
                result=dict(row,initial_radius=float(np.sqrt(initial['maximum_squared_distance'])),initial_whole_disk_score=initial,final_radius=R,final_whole_disk_score=d,final_diagnostic=run.diagnose(C),final_target_and_actual_metric=v,final_coordinate_file=str(p.relative_to(ROOT)),final_file_sha256=run.sha(p),final_coordinate_sha256=run.csha(C),rounds=rounds,seconds=time.monotonic()-tt,numerical_only=True)
            except Exception as e:result=dict(row,error=str(e),traceback=traceback.format_exc(),seconds=time.monotonic()-tt)
            records.append(result);fp.write(json.dumps(run.clean(result),allow_nan=False)+'\n');run.write(folder/'checkpoint.json',dict(worker=worker,completed=len(records),requested=len(rows),records=records,best_radius=min((r.get('final_radius',1)for r in records),default=1),seconds=time.monotonic()-tic,BLAS_threads=1,**pidinfo()));queue.put(dict(type='record',record=result))
    queue.put(dict(type='done',worker=worker))

def main():
    checks=json.loads((ROOT/'input_checks.json').read_text());assert checks['status']=='READY_VERIFIED';assert checks['executable']<=96
    grid=json.loads((ROOT/'grid.json').read_text());rows=[r for r in grid if not r['skip']];assert len(rows)==checks['executable'];assert rows
    cpus=sorted(os.sched_getaffinity(0));runtime=dict(source_sha256={str(p.relative_to(ROOT)):run.sha(p)for p in sorted((ROOT/'source').glob('*.py'))},worker_cpu_assignments=[cpus[w%len(cpus)]for w in range(4)],workers=4,BLAS_threads_per_worker=1,**pidinfo());run.write(ROOT/'runtime_provenance.json',runtime)
    ctx=mp.get_context('spawn');queue=ctx.Queue();procs=[ctx.Process(target=worker_main,args=(w,[r for r in rows if r['worker']==w],queue,cpus))for w in range(4)]
    for p in procs:p.start()
    runtime['worker_pids']=[p.pid for p in procs];run.write(ROOT/'runtime_provenance.json',runtime);print('EXCHANGE_WORKERS',json.dumps(runtime),flush=True)
    records=[];done=set();best=float('inf');tic=time.monotonic()
    while len(done)<4:
        try:item=queue.get(timeout=20)
        except qm.Empty:
            dead=[p for p in procs if p.exitcode not in(None,0)]
            if dead:raise RuntimeError('worker failed: '+str([(p.pid,p.exitcode)for p in dead]))
            continue
        if item['type']=='done':done.add(item['worker']);continue
        r=item['record'];records.append(r);R=r.get('final_radius',1)
        if R<run.THRESHOLD and run.REFERENCE-R>run.MATERIAL:
            p=ROOT/'frozen_proposals'/f'grid{r["grid_id"]:04d}_R{R:.17f}.npy';np.save(p,np.load(ROOT/r['final_coordinate_file']));run.write(p.with_suffix('.json'),dict(r,frozen_file_sha256=run.sha(p),frozen_coordinate_sha256=run.csha(np.load(p))));print('FROZEN_PROPOSAL',r['grid_id'],R,str(p),run.sha(p),flush=True)
        if R<best:
            best=R;np.save(ROOT/'best.npy',np.load(ROOT/r['final_coordinate_file']));run.write(ROOT/'best.json',r);print('EXCHANGE_BEST',r['grid_id'],best,flush=True)
        run.write(ROOT/'checkpoint.json',dict(completed=len(records),requested=len(rows),budget=96,skipped=len(grid)-len(rows),workers_completed=sorted(done),best_radius=best,records=records,seconds=time.monotonic()-tic))
        if len(records)%8==0:print('EXCHANGE_PROGRESS',len(records),len(rows),'best',best,'seconds',time.monotonic()-tic,flush=True)
    for p in procs:p.join()
    run.write(ROOT/'completion.json',dict(completed=len(records),requested=len(rows),best_radius=best,errors=sum('error'in r for r in records),worker_exitcodes=[p.exitcode for p in procs],seconds=time.monotonic()-tic));print('EXCHANGE_DONE',len(records),best,'seconds',time.monotonic()-tic,flush=True)

if __name__=='__main__':main()
