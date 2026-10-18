"""One frozen primary plus identical replay. Never launches a replacement."""
import hashlib,json,os,resource,subprocess,sys,time
from pathlib import Path
p=Path(__file__).resolve().parent
mode=sys.argv[1]
assert mode in ('primary','replay')
sha=lambda b:hashlib.sha256(b).hexdigest()
frozen=json.loads((p/'EXECUTED_SOURCES.json').read_text())
for name,digest in frozen['files'].items():assert sha((p/name).read_bytes())==digest,name
meta=p/('PRIMARY_WORKER.json' if mode=='primary' else 'REPLAY.json')
assert not meta.exists(),'Frozen worker already ran.'
if mode=='replay':assert json.loads((p/'PRIMARY_WORKER.json').read_text())['exit_code']==0
def limits():
    resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
    resource.setrlimit(resource.RLIMIT_CPU,(25,25))
    resource.setrlimit(resource.RLIMIT_FSIZE,(25*1024**2,25*1024**2))
t=time.monotonic()
try:
    r=subprocess.run([sys.executable,str(p/'diagnostic.py')],capture_output=True,
                     timeout=25,preexec_fn=limits,env=dict(os.environ,PYTHONHASHSEED='0'))
    output,error,code=r.stdout,r.stderr,r.returncode
except subprocess.TimeoutExpired as e:
    output,error,code=e.stdout or b'',e.stderr or b'',124
wall=time.monotonic()-t
if mode=='primary':(p/'DIAGNOSTIC.json').write_bytes(output)
elif output!=(p/'DIAGNOSTIC.json').read_bytes():
    (p/'REPLAY_DIFFERENCE.json').write_bytes(output);code=code or 3
if error:(p/(mode.upper()+'_STDERR.txt')).write_bytes(error)
record={'mode':mode,'exit_code':code,'wall_seconds':wall,
        'max_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
        'output_sha256':sha(output),'source_freeze_commit':frozen['commit'],
        'timeout_seconds':25,'memory_limit_bytes':1024**3,
        'replay_identical':None if mode=='primary' else output==(p/'DIAGNOSTIC.json').read_bytes()}
meta.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record))
raise SystemExit(code)
