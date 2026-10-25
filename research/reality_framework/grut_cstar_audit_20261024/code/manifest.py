#!/usr/bin/env python3
import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def main():
    files=[]
    for path in sorted(ROOT.rglob('*')):
        if path.is_file() and '__pycache__' not in path.parts and path.name!='MANIFEST.json':
            data=path.read_bytes();files.append(dict(path=str(path.relative_to(ROOT)),bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
    base=ROOT.parent
    paths=[base/'grut_core_conversion_20261023/code/scattering.py',base/'grut_core_conversion_20261023/MATHEMATICS.md',
        base/'grut_portal_debris_window_20261023/data/profile_Q200000_lambda0.03_refined.csv',base/'grut_portal_debris_window_20261023/data/profile_Q200000_lambda0.03.csv',
        base/'grut_portal_debris_window_20261023/FIXED_YIELD_RESULTS.json',base/'grut_candidate_1_20261018/audit/CONSTRUCTOR_NOTES.md',
        base/'grut_q_stage2_build_20261021/REPORT.md',base/'grut_unified_working_theory_20261022/WORKING_THEORY.md']
    inputs=[dict(path=str(p.relative_to(base)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths]
    out=dict(parent_commit='ff49e5d0b4b40be1f7beb472c8bfa9d7070cdef2',files=files,inputs=inputs)
    if '--check' in sys.argv:
        assert out==json.loads((ROOT/'MANIFEST.json').read_text());print('MANIFEST_PASS',len(files),'files;',len(inputs),'inputs')
    else:(ROOT/'MANIFEST.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
