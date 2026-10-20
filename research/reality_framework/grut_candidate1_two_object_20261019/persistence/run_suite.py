"""Run the fixed suite with three independent numerical workers (standard library)."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import subprocess
import time

ROOT = Path(__file__).resolve().parent


def main():
    manifest = json.loads((ROOT / 'FREEZE_MANIFEST.json').read_text())
    for name, expected in manifest.items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
            raise RuntimeError('Frozen input changed: ' + name)
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
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = [pool.submit(execute, spec) for spec in protocol['runs']]
        for future in as_completed(futures):
            result = future.result()
            outcomes.append(result)
            print(json.dumps(result), flush=True)
            (ROOT / 'execution_status.json').write_text(json.dumps(outcomes, indent=2) + '\n')
    if any(result['returncode'] for result in outcomes):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
