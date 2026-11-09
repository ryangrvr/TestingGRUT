"""One bounded, immutable R1 source-response diagnostic suite."""
import argparse,hashlib,importlib.util,json,os,platform,subprocess,sys,time,traceback
from datetime import datetime,timezone
from pathlib import Path
sys.dont_write_bytecode=True;os.environ['PYTHONDONTWRITEBYTECODE']='1'
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
import numpy as np,scipy
from shared import HERE,B2,D1,N22,card,write,bounded_run
from r1_reference import synthetic_precheck
from analyze import analyze

def run(inputs,out):
    root=HERE.parents[2]
    def git(*args):return subprocess.check_output(['git','-C',str(root),*args],text=True).strip()
    if git('status','--porcelain'):raise RuntimeError('Clean committed source required')
    if out==root or root in out.parents:raise ValueError('Use external outputs')
    cfg=card();out.mkdir(parents=True,exist_ok=False);start=time.monotonic();events=[]
    schwinger=B2.parent/'grut_qed_detector_histories_20261107/ed/schwinger.py'
    sources=list(HERE.glob('*.py'))+[HERE/'protocol.json',HERE/'INPUTS.json',HERE/'PROTOCOL.md',HERE/'DERIVATION.md',HERE/'EXPERIMENT_PROTOCOL.md',D1/'reference.py',D1/'source_precheck.py',B2/'model.py',B2.parent/'grut_detector_matching_20261107/lattice.py',B2.parent/'grut_detector_matching_20261107/continuum.py',N22/'resource_limit.py',schwinger]
    write(out/'provenance.json',{'freeze_commit':git('rev-parse','HEAD'),'utc_started':datetime.now(timezone.utc).isoformat(),
      'source_sha256':{p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
      'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),
      'resources':cfg['resources'],'input_manifest':json.loads((HERE/'INPUTS.json').read_text())})
    try:
        spec=importlib.util.spec_from_file_location('d1_source_precheck',D1/'source_precheck.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
        source=mod.check(schwinger);synthetic=synthetic_precheck()
        pre={'source':source,'synthetic':synthetic,'passed':source['passed'] and synthetic['passed']};write(out/'prechecks.json',pre)
        if not pre['passed']:raise RuntimeError('Precheck failed')
        for n,refined in [(n,False) for n in cfg['N']]+[(22,True)]:
            label=f'N{n}'+('_refined' if refined else '')
            remaining=cfg['resources']['total_wall_seconds']-(time.monotonic()-start)
            if remaining<=0:raise TimeoutError('Total resource cap')
            command=[sys.executable,str(HERE/'point.py'),'--N',str(n),'--inputs',str(inputs),'--output',str(out/(label+'.json'))]
            if refined:command.append('--refined')
            print(json.dumps({'event':'starting','point':label}),flush=True)
            event=bounded_run(command,out/(label+'.log'),label,min(remaining,cfg['resources']['point_wall_seconds']),cfg['resources'],out/(label+'_resources.json'));events.append(event)
            if event['returncode']!=0 or event['error']:raise RuntimeError('Point failed: '+label)
        assessment=analyze(out);write(out/'assessment.json',assessment)
        status={'execution_complete':True,'verdict':assessment['verdict'],'numerical_pass':assessment['pass']}
    except Exception as exc:
        status={'execution_complete':False,'error':str(exc),'traceback':traceback.format_exc()}
        if not (out/'assessment.json').exists():
            try:assessment=analyze(out)
            except Exception as assessment_error:assessment={'complete':False,'pass':False,'verdict':'analysis_failed','error':str(assessment_error),'massive_status':'withheld'}
            write(out/'assessment.json',assessment)
    status.update(seconds=time.monotonic()-start,events=events,massive_status='withheld')
    write(out/'execution.json',status);print(json.dumps(status),flush=True)
    return 0 if status['execution_complete'] else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();sys.exit(run(a.inputs.resolve(),a.output.resolve()))
