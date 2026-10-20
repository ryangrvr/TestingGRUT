"""Freeze before evaluation; then anchor the public source commit. Never overwrite."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

P = Path(__file__).resolve().parent
ROOT = P.parents[2]
BASE = P.parent / 'grut_3d_empirical_20261020'
BASE_COMMIT = '479f111117ea7162bd7888ff8f3792b2f6c45fdb'
FROZEN = ['README.md', 'PROTOCOL.json', 'DERIVATION.md', 'ASSUMPTIONS.md',
          'PARAMETER_PROVENANCE.md', 'FEASIBILITY.md', 'maser_table1.json',
          'freeze.py', 'audit.py']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def dump_new(path, value):
    with path.open('x') as f:
        f.write(json.dumps(value, indent=2, sort_keys=True) + '\n')


def main():
    if len(sys.argv) > 1 and sys.argv[1] == 'anchor':
        commit = sys.argv[2]
        if not re.fullmatch('[0-9a-f]{40}', commit):
            raise ValueError('Expected an exact commit SHA')
        rel = (P / 'FREEZE.json').relative_to(ROOT).as_posix()
        published = subprocess.check_output(['git', 'show', commit + ':' + rel], cwd=ROOT)
        assert hashlib.sha256(published).hexdigest() == sha(P / 'FREEZE.json')
        dump_new(P / 'PUBLIC_FREEZE.json', {
            'commit': commit,
            'url': 'https://github.com/ryangrvr/TestingGRUT/commit/' + commit,
            'freeze_sha256': sha(P / 'FREEZE.json'),
            'role': 'Public source freeze before new numerical evaluation; data already known.',
            'recorded_utc': datetime.now(timezone.utc).isoformat(),
        })
        print('Public source anchor recorded:', commit)
        return
    subprocess.run(['git', 'diff', '--quiet', BASE_COMMIT, '--',
                    BASE.relative_to(ROOT).as_posix(),
                    'research/reality_framework/grut_candidate_1_20261018'],
                   cwd=ROOT, check=True)
    protocol = json.loads((P / 'PROTOCOL.json').read_text())
    refs = [BASE / 'SOURCE_MANIFEST.json', BASE / 'RAW_MANIFEST.json',
            P / protocol['action'], P / protocol['calibration_file'],
            P / protocol['parameter_file'],
            P.parent / 'grut_candidate_1_20261018/code/results.json']
    dump_new(P / 'FREEZE.json', {
        'id': protocol['id'], 'base_commit': BASE_COMMIT,
        'local_freeze_utc': datetime.now(timezone.utc).isoformat(),
        'client_timezone': 'America/Chicago',
        'public_preregistration': False,
        'disclosure': 'Observations and prior rejection already known. Public source commit to follow before new numerical evaluation.',
        'source_sha256': {name: sha(P / name) for name in FROZEN},
        'baseline_input_sha256': {
            p.resolve().relative_to(ROOT).as_posix(): sha(p) for p in refs},
    })
    print('Local source freeze recorded; publish before running audit.py.')


if __name__ == '__main__':
    main()
