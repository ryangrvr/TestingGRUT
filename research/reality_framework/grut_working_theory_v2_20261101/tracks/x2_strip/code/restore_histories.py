#!/usr/bin/env python3
"""Restore the exact verified scan histories from the durable packet."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('work',type=Path);a=ap.parse_args()
    pack=HERE.parent/'results'/'restored_histories'
    meta=json.loads((pack/'provenance.json').read_text())
    for node,row in meta['histories'].items():
        root=a.work.resolve()/'histories'/node
        dst=root/'hist'/f'bench_{node}.json'
        data=gzip.decompress((pack/(dst.name+'.gz')).read_bytes())
        assert hashlib.sha256(data).hexdigest()==row['bench_sha256']
        assert row['rebuild_check']['passed']
        dst.parent.mkdir(parents=True,exist_ok=True)
        if dst.exists() and dst.read_bytes()!=data:
            raise FileExistsError(f'{dst} differs from the verified history; use a new work directory')
        dst.write_bytes(data)
        (root/'rebuild_check.json').write_text(json.dumps(row['rebuild_check'],indent=2)+'\n')
        print('Restored',node, row['bench_sha256'])


if __name__=='__main__':main()
