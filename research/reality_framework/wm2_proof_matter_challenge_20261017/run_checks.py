"""Run only the bounded challenge controls; never launches a formation solver."""
from pathlib import Path
import hashlib
import json
import platform
import time
import numpy as np
import scipy
from charge_controls import run as charge
from dyson_controls import run as dyson
from matter_controls import run as matter
from resource_controls import run as resource

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def main():
    pins = json.loads((HERE / 'SOURCE_PINS.json').read_text())
    for relative, digest in pins['sha256'].items():
        assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == digest, relative
    names = ['run_checks.py', 'charge_controls.py', 'dyson_controls.py', 'matter_controls.py', 'resource_controls.py']
    executed = {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in names}
    (HERE / 'EXECUTED_SOURCES.json').write_text(json.dumps(executed, indent=2, sort_keys=True) + '\n')
    result = {'base_commit': pins['base_commit'], 'charter_commit': 'a0c53ef6df3667036e160202b8da76d80e5e49c1',
              'provenance': 'same Builder; fresh code, not independent external review',
              'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__}}
    for name, work in [('charge_blocks', charge), ('interaction_words', dyson), ('distinct_matter', matter), ('receiver_dependence', resource)]:
        start = time.monotonic()
        result[name] = work()
        result[name]['elapsed_seconds'] = time.monotonic() - start
        print(name, result[name]['status'], f"{result[name]['elapsed_seconds']:.3f}s", flush=True)
    result['status'] = 'PASS_DIAGNOSTICS_NOT_EXTERNAL_PROOF_REVIEW'
    (HERE / 'CHECKS.json').write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], flush=True)


if __name__ == '__main__':
    main()
