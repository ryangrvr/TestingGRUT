"""Hostile freshness checks of the lane reporter; no manifest mutation."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent


def main():
    cases = ('unchanged', 'missing_results', 'wrong_result', 'wrong_table', 'changed_code', 'changed_proof')
    outcomes = []
    for case in cases:
        with tempfile.TemporaryDirectory(prefix='grut-physical-fixture-') as directory:
            copy = Path(directory)/'lab'
            shutil.copytree(ROOT, copy, ignore=shutil.ignore_patterns('evidence', '__pycache__'))
            if case == 'missing_results':
                (copy/'RESULTS.json').unlink()
            elif case == 'wrong_result':
                data = json.loads((copy/'RESULTS.json').read_text())
                data['controls']['C02_quantum_memory_bound']['perfect_quantum_dimension_lower_bound'] = 2
                (copy/'RESULTS.json').write_text(json.dumps(data))
            elif case == 'wrong_table':
                with (copy/'TABLES.md').open('a') as handle:
                    handle.write('unexpected table\n')
            elif case == 'changed_code':
                with (copy/'controls.py').open('a') as handle:
                    handle.write('\n# unhashed amendment\n')
            elif case == 'changed_proof':
                with (copy/'PROOFS.md').open('a') as handle:
                    handle.write('\nUNREVIEWED AMENDMENT\n')
            process = subprocess.run([sys.executable, str(copy/'run.py'), '--check'],
                                     capture_output=True, text=True, timeout=30)
            wanted = 0 if case == 'unchanged' else 1
            correct = process.returncode == wanted
            if wanted == 1:
                correct = correct and 'FAIL: missing/stale' in process.stderr
            outcomes.append({'fixture': case, 'returncode': process.returncode,
                             'expected_returncode': wanted, 'pass': correct})
    print(json.dumps(outcomes, indent=2))
    return 0 if all(item['pass'] for item in outcomes) else 1


if __name__ == '__main__':
    sys.exit(main())
