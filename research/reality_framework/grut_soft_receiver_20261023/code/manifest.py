#!/usr/bin/env python3
import hashlib,json,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    files=[]
    for p in sorted(ROOT.rglob('*')):
        if p.is_file() and p.name!='MANIFEST.json' and '__pycache__' not in p.parts:
            b=p.read_bytes();files.append(dict(path=str(p.relative_to(ROOT)),bytes=len(b),sha256=hashlib.sha256(b).hexdigest()))
    parent=ROOT.parent/'grut_core_conversion_20261023'
    inputs=[parent/'code/radial_nodal.cpp']+[parent/f'{name}_RESULTS.json' for name in ['N025','N0125','N00625','N0125_far','N0125_box','N0125_abs']]
    pins=[dict(path=str(p.relative_to(ROOT.parent)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in inputs]
    out=dict(parent_commit='59a1a712495c5205cc5df042154468225602d3d6',files=files,inputs=pins,
        raw_input_note='Raw numeric trace/final/history hashes and performed checks are in TRACE_AUDIT_RESULTS and VERIFICATION_RESULTS; raw arrays are not duplicated here.')
    if '--check' in sys.argv:
        assert out==json.loads((ROOT/'MANIFEST.json').read_text());print('MANIFEST_PASS',len(files),'files;',len(pins),'inputs')
    else:(ROOT/'MANIFEST.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
