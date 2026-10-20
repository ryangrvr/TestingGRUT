"""Run the fixed Gaussian-cloud suite with two numerical workers."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import subprocess
import time

ROOT = Path(__file__).resolve().parent
manifest = json.loads((ROOT / 'FREEZE_MANIFEST.json').read_text())
for name, expected in manifest.items():
    if hashlib.sha256((ROOT.parent / name).read_bytes()).hexdigest() != expected:
        raise RuntimeError('Changed cloud input: ' + name)
protocol = json.loads((ROOT / 'protocol.json').read_text())
logs = ROOT / 'logs'
logs.mkdir(exist_ok=True)


def execute(spec):
    start = time.monotonic()
    with (logs / (spec['id'] + '.log')).open('w') as log:
        result = subprocess.run(['node', 'run.cjs', spec['id']], cwd=ROOT,
                                stdout=log, stderr=subprocess.STDOUT)
    return dict(id=spec['id'], returncode=result.returncode,
                wall_seconds=time.monotonic() - start)


outcomes = []
with ThreadPoolExecutor(max_workers=2) as pool:
    for future in as_completed([pool.submit(execute, spec) for spec in protocol['runs']]):
        result = future.result()
        outcomes.append(result)
        print(json.dumps(result), flush=True)
        (ROOT / 'execution_status.json').write_text(json.dumps(outcomes, indent=2) + '\n')
if any(x['returncode'] for x in outcomes):
    raise SystemExit(1)
