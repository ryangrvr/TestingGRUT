"""Single frozen N22 experiment, including refinement at N22; no massive path."""
import argparse,hashlib,json,os,platform,subprocess,sys,time,traceback
from datetime import datetime,timezone
from pathlib import Path
sys.dont_write_bytecode=True;os.environ['PYTHONDONTWRITEBYTECODE']='1'
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
import numpy as np
import scipy
from shared import HERE,B2,card,precheck,write
from b3_analyze import analyze
from resource_limit import bounded_run

def run(out):
    root=HERE.parents[2]
    def git(*a):return subprocess.check_output(['git','-C',str(root),*a],text=True).strip()
    if git('status','--porcelain'):raise RuntimeError('Clean committed source required')
    if out==root or root in out.parents:raise ValueError('Output must be external')
    cfg=card();limits=cfg['resources'];out.mkdir(parents=True,exist_ok=False);start=time.monotonic();events=[]
    paths=list(HERE.glob('*.py'))+[HERE/'protocol.json']+list(B2.glob('*.py'))+[B2/'protocol_b2_01.json']
    paths+=[B2.parent/'grut_detector_matching_20261107'/p for p in ('continuum.py','lattice.py')]
    paths+=[B2.parent/'grut_qed_detector_histories_20261107/ed/schwinger.py']+[B2/'raw'/f'N{n}.json' for n in (10,14,18,20)]
    hashes={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    write(out/'provenance.json',{'freeze_commit':git('rev-parse','HEAD'),'utc_started':datetime.now(timezone.utc).isoformat(),
        'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),
        'source_sha256':hashes,'limits':limits,'thread_environment':{k:os.environ[k] for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')}})
    try:
        checks=precheck(cfg);write(out/'precheck.json',checks)
        if not checks['passed']:raise RuntimeError('Precheck failed')
        for refined in (False,True):
            label='N22_refined' if refined else 'N22'
            remaining=limits['total_wall_seconds']-(time.monotonic()-start)
            if remaining<=0:raise TimeoutError('Total resource cap')
            command=[sys.executable,str(HERE/'point_n22.py'),'--output',str(out/(label+'.json'))]
            if refined:command.append('--refined')
            print(json.dumps({'event':'starting','point':label}),flush=True)
            event=bounded_run(command,out/(label+'.log'),label,min(remaining,limits['point_wall_seconds']),limits,out/(label+'_resources.json'))
            events.append(event)
            if event['returncode']!=0 or event['error']:raise RuntimeError('Point failed: '+label)
        assessment=analyze(out);write(out/'assessment.json',assessment)
        status={'execution_complete':True,'verdict':assessment['verdict'],'nominal_pass':assessment['nominal_pass'],'margin_pass':assessment['margin_pass']}
    except Exception as exc:
        status={'execution_complete':False,'error':str(exc),'traceback':traceback.format_exc()}
        if not (out/'assessment.json').exists():write(out/'assessment.json',analyze(out))
    status.update(seconds=time.monotonic()-start,events=events,massive_status='withheld_no_transfer_authorized')
    write(out/'execution.json',status);print(json.dumps(status),flush=True)
    return 0 if status['execution_complete'] else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    sys.exit(run(args.output.resolve()))
