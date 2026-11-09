"""Run a clean, committed B2 protocol once into a new external directory."""
import argparse,hashlib,json,os,platform,subprocess,sys,time,traceback
from datetime import datetime,timezone
from pathlib import Path

sys.dont_write_bytecode=True
os.environ['PYTHONDONTWRITEBYTECODE']='1'
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
import scipy
from model import HERE,OLD
from precheck import precheck
from analyze import analyze,massive
from point import write

def run(out):
    root=HERE.parents[2]
    def git(*args):return subprocess.check_output(['git','-C',str(root),*args],text=True).strip()
    if git('status','--porcelain'):raise RuntimeError('Production requires a clean committed checkout')
    if root==out or root in out.parents:raise ValueError('Output must be outside the source checkout')
    card=json.loads((HERE/'protocol_b2_01.json').read_text());out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic();limit=card['resources'];events=[]
    paths=list(HERE.glob('*.py'))+[HERE/'protocol_b2_01.json',OLD/'continuum.py',OLD/'lattice.py']
    # The inherited lattice module dynamically imports the original Schwinger implementation.
    paths+=[HERE.parent/'grut_qed_detector_histories_20261107/ed/schwinger.py']
    hashes={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    provenance={'freeze_commit':git('rev-parse','HEAD'),'utc_started':datetime.now(timezone.utc).isoformat(),
        'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
        'platform':platform.platform(),'source_sha256':hashes,'thread_environment':{k:os.environ[k] for k in
            ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')},'limits':limit}
    write(out/'provenance.json',provenance)
    def point(label,n,refined=False,mass=0.):
        remaining=limit['total_wall_seconds']-(time.monotonic()-start)
        if remaining<=0:raise TimeoutError('Total resource cap reached')
        cmd=[sys.executable,str(HERE/'point.py'),'--N',str(n),'--mass',str(mass),'--output',str(out/(label+'.json'))]
        if refined:cmd.append('--refined')
        print(json.dumps({'event':'starting','point':label}),flush=True)
        tick=time.monotonic()
        with (out/(label+'.log')).open('x',encoding='utf-8') as log:
            result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=min(limit['point_wall_seconds'],remaining))
        event={'point':label,'returncode':result.returncode,'wall_seconds':time.monotonic()-tick};events.append(event)
        print(json.dumps(event),flush=True)
        if result.returncode:raise RuntimeError('Point failed: '+label)
    try:
        checks=precheck(card);write(out/'precheck.json',checks)
        if not checks['passed']:raise RuntimeError('Precheck failed')
        for n in card['N']:point(f'N{n}',n)
        point('N14_refined',card['refined_N'],True)
        assessment=analyze(out);write(out/'assessment.json',assessment)
        if assessment['massless_pass']:point('massive_N20',20,mass=.5)
        transfer=massive(out);write(out/'massive_status.json',transfer)
        status={'execution_complete':True,'massless_pass':assessment['massless_pass'],'massive_status':transfer['status']}
    except Exception as exc:
        status={'execution_complete':False,'error':str(exc),'traceback':traceback.format_exc()}
        print(json.dumps(status),flush=True)
    status['events']=events;status['seconds']=time.monotonic()-start
    write(out/'execution.json',status);print(json.dumps(status),flush=True)
    return 0 if status['execution_complete'] else 1

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();sys.exit(run(args.output.resolve()))
