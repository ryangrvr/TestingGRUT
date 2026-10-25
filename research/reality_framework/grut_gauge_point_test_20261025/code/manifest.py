#!/usr/bin/env python3
import hashlib,json,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    files=[]
    for p in sorted(ROOT.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts and p.name!='MANIFEST.json':
            data=p.read_bytes();files.append(dict(path=str(p.relative_to(ROOT)),bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
    rel=['grut_z6_gauge_bridge_20261025/MATCHING_RESULTS.json','grut_cstar_audit_20261024/CHARGE_LEAKAGE_RESULTS.json',
         'grut_core_conversion_20261023/code/scattering.py','grut_portal_debris_window_20261023/FIXED_YIELD_RESULTS.json',
         'grut_portal_debris_window_20261023/data/profile_Q200000_lambda0.03_refined.csv']
    inputs=[dict(path=r,sha256=hashlib.sha256((ROOT.parent/r).read_bytes()).hexdigest()) for r in rel]
    out=dict(parent_commit='7737487a890c4261e544d067700f9a01cdba9904',review_commit='2e9f2968e30b2d6d011dd35289ab7aa56129465a',files=files,inputs=inputs)
    if '--check' in sys.argv:
        assert out==json.loads((ROOT/'MANIFEST.json').read_text());print('MANIFEST_PASS',len(files),'files;',len(inputs),'inputs')
    else:(ROOT/'MANIFEST.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
