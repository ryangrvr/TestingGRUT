"""One bounded W1 validation suite, requiring clean frozen Git source."""
import argparse,hashlib,json,os,platform,subprocess,sys,time,traceback
from datetime import datetime,timezone
from pathlib import Path
sys.dont_write_bytecode=True;os.environ['PYTHONDONTWRITEBYTECODE']='1'
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
import numpy as np,scipy
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];N22=HERE.parent/'grut_n22_refinement_20261109'
sys.path.insert(1,str(N22))
from resource_limit import bounded_run

def write(path,data):
    with path.open('x',encoding='utf-8') as f:json.dump(data,f,indent=2,allow_nan=False);f.write('\n')
def main(out):
    def git(*args):return subprocess.check_output(['git','-C',str(ROOT),*args],text=True).strip()
    if git('status','--porcelain'):raise RuntimeError('Clean frozen source required')
    if out==ROOT or ROOT in out.parents:raise ValueError('Use external output directory')
    out.mkdir(parents=True,exist_ok=False);cfg=json.loads((HERE/'protocol.json').read_text());start=time.monotonic();events=[]
    sources=[p for p in HERE.iterdir() if p.is_file()]+[N22/'resource_limit.py',HERE.parent/'grut_response_separation_20261109/r1_reference.py']
    write(out/'provenance.json',{'freeze_commit':git('rev-parse','HEAD'),'utc_started':datetime.now(timezone.utc).isoformat(),
      'source_sha256':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
      'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'resources':cfg['resources'],
      'actual_experimental_data':False,'new_large_lattice_evolution':False})
    try:
        for label,source in [('finite_model','finite_model.py'),('autonomous_parent','autonomous_parent.py'),('spectral','spectral_run.py')]:
            remaining=cfg['resources']['total_wall_seconds']-(time.monotonic()-start)
            if remaining<=0:raise TimeoutError('Total resource cap')
            command=[sys.executable,str(HERE/source),'--output',str(out/(label+'.json'))]
            print(json.dumps({'starting':label}),flush=True)
            event=bounded_run(command,out/(label+'.log'),label,min(remaining,cfg['resources']['point_wall_seconds']),cfg['resources'],out/(label+'_resources.json'));events.append(event)
            if event['returncode']!=0 or event['error']:raise RuntimeError('Registered calculation failed: '+label)
        finite=json.loads((out/'finite_model.json').read_text());parent=json.loads((out/'autonomous_parent.json').read_text());spectral=json.loads((out/'spectral.json').read_text())
        assessment={'execution_complete':True,'finite_parametric_model':finite['status']=='theoretical_internal_checks_passed',
          'autonomous_parent':parent['passed'],'spectral_implementation':spectral['numerically_valid'],
          'exact_mathematical_truncation_bound':spectral['primary_mathematical_uniform_certificate']['passed'],
          'primary_mathematical_bound':spectral['primary_mathematical_uniform_certificate']['bound'],
          'primary_measured_new_time_error':spectral['comparisons']['1024']['max_canonical_primary_error'],
          'physical_dictionary_status':'unresolved beyond the historical qualified scope',
          'empirical_status':'not_run; device adapter and calibration data missing',
          'distinctive_fundamental_law_status':'not supplied',
          'scope':'Working finite operational model and adopted-continuum response predictor, with separate proof and numerical status. No aggregate empirical theory pass.',
          'massive_status':'withheld'}
    except Exception as exc:
        assessment={'execution_complete':False,'error':str(exc),'traceback':traceback.format_exc(),'empirical_status':'not_run','massive_status':'withheld'}
    write(out/'assessment.json',assessment);write(out/'execution.json',{'seconds':time.monotonic()-start,'events':events,'execution_complete':assessment['execution_complete']})
    print(json.dumps(assessment),flush=True)
    return 0 if assessment['execution_complete'] else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args();raise SystemExit(main(args.output.resolve()))
