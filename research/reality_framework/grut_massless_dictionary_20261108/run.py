"""Execute only frozen B1 tasks and preserve the conditional-transfer decision."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
import scipy
from analyze import analyze,analyze_massive
from dictionary import verify,OLD
from point import write

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();out=args.output_dir.resolve()
    root=Path(subprocess.check_output(['git','rev-parse','--show-toplevel'],cwd=HERE,text=True).strip())
    if out.is_relative_to(root):raise ValueError('Outputs must be outside frozen checkout')
    if subprocess.check_output(['git','status','--porcelain'],cwd=HERE,text=True):raise ValueError('Clean source required')
    out.mkdir(parents=True,exist_ok=False)
    card=json.loads((HERE/'protocol_b1_01.json').read_text())
    inputs=list(HERE.glob('*.py'))+list(HERE.glob('*.md'))+[HERE/'protocol_b1_01.json']
    inputs+=[OLD/f for f in ('continuum.py','lattice.py')]+[OLD.parent/'grut_qed_detector_histories_20261107/ed/schwinger.py']
    inputs+=[OLD/'raw_m0_02'/f'point_N{n}_primary.json' for n in card['N']]
    inputs+=[HERE.parent/'grut_parent_dictionary_20261108/raw/E14_q128.json']
    inputs+=[OLD/'raw_m1_01/lattice_N20.json']
    write(out/'provenance.json',{'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=HERE,text=True).strip(),
        'source_clean':True,'card':card,'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
        'source_and_input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}})
    began=time.monotonic();checks=verify(card);write(out/'verification.json',checks)
    if max(checks.values())>card['gates']['algebraic_check']:raise ValueError('Algebraic verification failed')
    events=[];env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    def task(label,n,mass=0,refined=False):
        start=time.monotonic();left=card['resources']['total_wall_seconds']-(start-began);event={'label':label}
        if left<=0:event['status']='total_wall_cap'
        else:
            cmd=[sys.executable,str(HERE/'point.py'),'--N',str(n),'--mass',str(mass),'--output',str(out/(label+'.json'))]
            if refined:cmd.append('--refined')
            try:
                r=subprocess.run(cmd,env=env,capture_output=True,text=True,timeout=min(left,card['resources']['point_wall_seconds']))
                event.update(status='complete' if r.returncode==0 else 'exception',returncode=r.returncode,stdout=r.stdout,stderr=r.stderr)
            except subprocess.TimeoutExpired as e:event.update(status='wall_cap',stdout=str(e.stdout),stderr=str(e.stderr))
        event['seconds']=time.monotonic()-start;events.append(event);write(out/('event_'+label+'.json'),event)
        print(json.dumps(event),flush=True)
    for n in card['N']:task(f'N{n}',n)
    task('N14_refined',14,refined=True)
    result=analyze(out);write(out/'analysis.json',result)
    print(json.dumps({'verdict':result['verdict'],'gates':result.get('gates')}),flush=True)
    if result['massless_pass']:task('massive_N20',20,mass=.5)
    write(out/'massive_stage.json',analyze_massive(out));write(out/'events.json',events)


if __name__=='__main__':main()
