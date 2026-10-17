"""Sequential remaining frozen G1 workers and immutable S7 comparison.

Requires the first G1 worker to have terminated. Never overwrites an outcome.
"""
from pathlib import Path
import json
import os
import subprocess
import sys
import campaign

BASE=Path(__file__).resolve().parent
S7=BASE.parent/'wm2_positive_residual_20261017'


def worker(kind,counts,filename):
    output=BASE/'results'/filename
    assert not output.exists(),'Refuse to overwrite '+str(output)
    if kind=='s7':
        args=[sys.executable,str(S7/'worker.py'),'--method','s7','--counts',counts,'--output',str(output)]
        cwd=S7
    else:
        args=[sys.executable,str(BASE/'worker_global.py'),'--counts',counts,'--output',str(output)]
        cwd=BASE
    subprocess.run(args,cwd=cwd,env=dict(os.environ,OPENBLAS_NUM_THREADS='1'),check=True)


def main():
    assert (BASE/'results/g1_N5_root0.json').exists()
    worker('g1','3,1','g1_N5_root1.json')
    worker('s7','2,2','s7_N5_root0.json')
    worker('s7','3,1','s7_N5_root1.json')
    try:campaign.gate_six()
    except AssertionError:
        print(json.dumps(dict(six_event='NOT_RUN',reason='frozen G1 recovery gate failed')),flush=True)
    else:worker('g1','2,3','g1_N6_root0.json')
    print(json.dumps(campaign.summary(),indent=2),flush=True)


if __name__=='__main__':main()
