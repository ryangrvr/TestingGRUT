"""Exactly two predeclared fresh workers; no fallback, N6 or overwrites."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys

BASE=Path(__file__).resolve().parent


def main():
    assert json.loads((BASE/'SWEEP_CONTROLS.json').read_text())['status']=='PASS'
    for name,value in json.loads((BASE/'EXECUTED_SOURCES.json').read_text())['sha256'].items():
        assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==value,name
    for root,counts in enumerate(['2,2','3,1']):
        output=BASE/f'results/S9_N5_root{root}.json'
        assert not output.exists(),'Refuse to overwrite '+str(output)
        subprocess.run([sys.executable,str(BASE/'worker.py'),'--counts',counts,
            '--output',str(output),'--seconds','1200'],cwd=BASE,
            env=dict(os.environ,OPENBLAS_NUM_THREADS='1'),check=True)
    print('Both frozen components recorded. No further trajectory is authorized in S9.',flush=True)


if __name__=='__main__':main()
