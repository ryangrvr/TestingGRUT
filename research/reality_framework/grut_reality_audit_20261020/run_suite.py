"""Run ten specified tests with two local numerical workers."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import hashlib
import json
from full_field import execute

ROOT=Path(__file__).resolve().parent

if __name__=='__main__':
    manifest=json.loads((ROOT/'LOCAL_FREEZE.json').read_text())
    for name,digest in manifest.items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    protocol=json.loads((ROOT/'protocol.json').read_text())
    outcomes=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(execute,spec) for spec in protocol['runs']]):
            outcomes.append(future.result())
            print(json.dumps(outcomes[-1]),flush=True)
            (ROOT/'execution_status.json').write_text(json.dumps(outcomes,indent=2)+'\n')
