#!/usr/bin/env python3
import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def main():
    files=[]
    for p in sorted(ROOT.rglob('*')):
        if p.is_file() and 'data' not in p.relative_to(ROOT).parts and '__pycache__' not in p.parts and p.name!='MANIFEST.json':
            b=p.read_bytes();files.append(dict(path=str(p.relative_to(ROOT)),bytes=len(b),sha256=hashlib.sha256(b).hexdigest()))
    inputs=[ROOT.parent/'grut_soft_receiver_20261023/code/trace_audit.py',ROOT.parent/'grut_soft_receiver_20261023/BOUNDARY_RESULTS.json',ROOT.parent/'grut_core_conversion_20261023/code/radial_nodal.cpp',ROOT.parent/'grut_core_conversion_20261023/code/refinement.py',ROOT.parent/'grut_core_conversion_20261023/N00625_RESULTS.json']
    pins=[dict(path=str(p.relative_to(ROOT.parent)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in inputs]
    out=dict(parent_commit='9636950b4b56b582663c2805737629907454d95d',files=files,inputs=pins,
        raw_policy='Lossless raw arrays remain in the builder workspace and are reproducible; their numeric hashes are included in every run result. No bulk traces are duplicated in this compact publication.')
    if '--check' in sys.argv:
        assert out==json.loads((ROOT/'MANIFEST.json').read_text());print('MANIFEST_PASS',len(files),'files;',len(pins),'inputs')
    else:(ROOT/'MANIFEST.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
