"""Run the single frozen finite-diagnostic suite. No formation simulations."""
from pathlib import Path
import hashlib
import json
import resource
import time
from diagnostics import (qubits, full_kernel_controls, prefix_control,
                         configuration_count)


def run():
    started = time.perf_counter()
    here = Path(__file__).resolve().parent
    repo = here.parents[2]
    pins = json.loads((here / 'SOURCE_PINS.json').read_text())
    assert hashlib.sha256((here / 'CHARTER.md').read_bytes()).hexdigest() == pins['charter_sha256']
    for path, digest in pins['files'].items():
        assert hashlib.sha256((repo / path).read_bytes()).hexdigest() == digest, path
    result = {'status': 'RUNNING',
              'precision': 'ordinary complex128, no outward rounding',
              'receiver_scope': 'two-mode diagnostic models only; no physical-L formation computation',
              'qubits': qubits(),
              'physical_kernel': full_kernel_controls(),
              'independent_prefix': [prefix_control(7, 1), prefix_control(8, 2)],
              'configuration_counts': [configuration_count(7, 6, 4, 0, 1),
                                       configuration_count(9, 8, 5, 0, 2),
                                       configuration_count(11, 10, 7, 1, 2)]}
    elapsed = time.perf_counter() - started
    peak_kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    assert elapsed < 60, 'Frozen worker wall-time stop exceeded'
    assert peak_kib < 1024 ** 2, 'Frozen 1 GiB resident-memory stop exceeded'
    result.update(status='PASS_FINITE_DIAGNOSTICS', elapsed_seconds=elapsed,
                  peak_resident_MiB=peak_kib / 1024,
                  mathematical_review='AUTHOR_DERIVED_INDEPENDENT_REVIEW_PENDING',
                  derived_universal_locality='NOT_ESTABLISHED',
                  inherited_formation='UNCHANGED_AUTHOR_DERIVED_S9',
                  new_formation_trajectories=0)
    for path, digest in pins['files'].items():
        assert hashlib.sha256((repo / path).read_bytes()).hexdigest() == digest, path
    return result


if __name__ == '__main__':
    print(json.dumps(run(), indent=2, allow_nan=False))
