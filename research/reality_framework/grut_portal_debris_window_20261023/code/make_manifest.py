#!/usr/bin/env python3
"""Pin this additive packet plus the exact consumed parent inputs."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PARENTS=[
    '../grut_unified_working_theory_20261022/RESULTS.json',
    '../grut_unified_working_theory_20261022/data/profile_Q39944.csv',
    '../grut_unified_working_theory_20261022/data/profile_Q79988.csv',
    '../grut_portal_empirical_test_20261022/SOURCE_INPUTS.json',
    '../grut_portal_empirical_test_20261022/data/public_inputs.npz',
    '../grut_portal_empirical_test_20261022/code/model.py',
    '../grut_portal_empirical_test_20261022/code/detector.py',
    '../grut_fragmentation_transport_20261022/code/transport.py',
]


def fingerprint(p):
    return dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def payload():
    own={str(p.relative_to(ROOT)):fingerprint(p) for p in sorted(ROOT.rglob('*'))
         if p.is_file() and '__pycache__' not in p.parts and p.name!='MANIFEST.json'}
    return dict(parent_commit='ac96f9054b06a3d3c8cff1c1724b01e67fca1a4e',
                packet_files=own,consumed_parent_files={n:fingerprint(ROOT/n) for n in PARENTS})


def main(check):
    p=ROOT/'MANIFEST.json'
    actual=payload()
    if check:
        assert json.loads(p.read_text())==actual,'Manifest mismatch'
    else:
        p.write_text(json.dumps(actual,indent=2)+'\n')
    print('Manifest matches' if check else 'Manifest written',len(actual['packet_files']),
          'packet files plus',len(PARENTS),'consumed parent inputs; manifest is additional.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true')
    main(parser.parse_args().check)
