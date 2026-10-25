#!/usr/bin/env python3
import hashlib,json,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    files=[]
    for p in sorted(ROOT.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts and p.name!='MANIFEST.json':
            raw=p.read_bytes();files.append(dict(path=str(p.relative_to(ROOT)),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
    relative=['grut_cstar_audit_20261024/CHARGE_LEAKAGE_RESULTS.json','grut_cstar_audit_20261024/MATHEMATICS.md',
              'grut_cstar_audit_20261024/GLOBAL_CHARGE.md','grut_portal_debris_window_20261023/FIXED_YIELD_RESULTS.json']
    inputs=[dict(path=r,sha256=hashlib.sha256((ROOT.parent/r).read_bytes()).hexdigest()) for r in relative]
    out=dict(parent_commit='8986d348f59a654bf329dfab914fea180fdd30e7',read_only_WIP_review='fdc7204225e984258231b5b9290b804938a3b699',files=files,inputs=inputs)
    if '--check' in sys.argv:
        assert out==json.loads((ROOT/'MANIFEST.json').read_text());print('MANIFEST_PASS',len(files),'files;',len(inputs),'inputs')
    else:(ROOT/'MANIFEST.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
