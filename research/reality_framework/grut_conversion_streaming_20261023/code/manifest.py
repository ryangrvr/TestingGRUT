#!/usr/bin/env python3
import argparse,hashlib,json
from pathlib import Path
from streaming import ROOT,PARENT,POINT,FIXED


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
    files=[f for f in sorted(ROOT.rglob('*')) if f.is_file() and f.name!='MANIFEST.json' and '__pycache__' not in f.parts]
    inputs=[PARENT/'RATE_RESULTS.json',PARENT/'MATHEMATICS.md',POINT,FIXED]
    d=dict(parent_commit='bb211d1680d8ee8540eb7a8996ddb9b2b8c1546e',
           files=[dict(path=str(f.relative_to(ROOT)),bytes=f.stat().st_size,sha256=sha(f)) for f in files],
           inputs=[dict(path=str(f.relative_to(ROOT.parent)),sha256=sha(f)) for f in inputs])
    dest=ROOT/'MANIFEST.json'
    if args.check:assert json.loads(dest.read_text())==d;print('MANIFEST_PASS',len(files),'files;',len(inputs),'inputs')
    else:dest.write_text(json.dumps(d,indent=2)+'\n')


if __name__=='__main__':main()
