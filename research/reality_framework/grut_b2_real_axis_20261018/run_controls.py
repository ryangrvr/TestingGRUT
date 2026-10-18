"""One bounded primary or replay; never changes the frozen source."""
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

p=Path(__file__).resolve().parent
mode=sys.argv[1]
assert mode in ('primary','replay')
digest=lambda data:hashlib.sha256(data).hexdigest()
frozen=json.loads((p/'EXECUTED_SOURCES.json').read_text())
for name,h in frozen['files'].items():
    assert digest((p/name).read_bytes())==h, name
meta=p/('PRIMARY_WORKER.json' if mode=='primary' else 'REPLAY.json')
assert not meta.exists(), 'This bounded worker has already executed.'
if mode=='replay':
    primary=json.loads((p/'PRIMARY_WORKER.json').read_text())
    assert primary['exit_code']==0, 'No rescue replay after a failed primary.'

def limits():
    resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
    resource.setrlimit(resource.RLIMIT_CPU,(25,25))
    resource.setrlimit(resource.RLIMIT_FSIZE,(25*1024**2,25*1024**2))

env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
t0=time.monotonic()
try:
    r=subprocess.run([sys.executable,str(p/'controls.py')],capture_output=True,
                     timeout=25,env=env,preexec_fn=limits)
    output,error,code=r.stdout,r.stderr,r.returncode
except subprocess.TimeoutExpired as e:
    output,error,code=e.stdout or b'',e.stderr or b'',124
wall=time.monotonic()-t0
if mode=='primary':
    (p/'CONTROLS.json').write_bytes(output)
elif output!=(p/'CONTROLS.json').read_bytes():
    (p/'REPLAY_DIFFERENCE.json').write_bytes(output)
    code=code or 3
if error:
    (p/(mode.upper()+'_STDERR.txt')).write_bytes(error)
record={'mode':mode,'exit_code':code,'wall_seconds':wall,
        'max_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
        'output_sha256':digest(output),'source_freeze_commit':frozen['commit'],
        'timeout_seconds':25,'memory_limit_bytes':1024**3,
        'replay_identical':None if mode=='primary' else output==(p/'CONTROLS.json').read_bytes()}
meta.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record))
raise SystemExit(code)
