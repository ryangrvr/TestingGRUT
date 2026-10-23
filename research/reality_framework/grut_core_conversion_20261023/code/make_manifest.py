#!/usr/bin/env python3
import argparse
import hashlib
import json
from pathlib import Path
from scattering import ROOT

PARENT=ROOT.parent
INPUTS=[
 'grut_unified_working_theory_20261022/WORKING_THEORY.md',
 'grut_portal_debris_window_20261023/FIXED_YIELD_RESULTS.json',
 'grut_portal_debris_window_20261023/data/profile_Q200000_lambda0.03.csv',
 'grut_portal_debris_window_20261023/data/profile_Q200000_lambda0.03_refined.csv',
 'grut_daughter_core_20261023/PARAMETER_POINT.json',
 'grut_daughter_core_20261023/CORE_RESULTS.json',
 'grut_daughter_core_20261023/SPECTRAL_RESULTS.json',
 'grut_daughter_core_20261023/code/run_core.py',
 'grut_daughter_core_20261023/code/radial.cpp',
]


def record(path,base):
    b=path.read_bytes()
    return dict(path=str(path.relative_to(base)),bytes=len(b),sha256=hashlib.sha256(b).hexdigest())


def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
    files=[x for x in sorted(ROOT.rglob('*')) if x.is_file() and x.name!='MANIFEST.json'
           and '__pycache__' not in x.parts and 'data' not in x.relative_to(ROOT).parts]
    data=dict(publication_scope='Compact report, code and recorded results; raw arrays separately inventoried in RAW_DATA_MANIFEST.json.',parent_commit='c4b99f3db8c6394b29f952198f8fd9bb7640d071',
              files=[record(x,ROOT) for x in files],consumed_parent_inputs=[record(PARENT/x,PARENT) for x in INPUTS])
    path=ROOT/'MANIFEST.json'
    if args.check:
        assert json.loads(path.read_text())==data,'Manifest differs'
        print('MANIFEST_PASS',len(files),'packet files and',len(INPUTS),'consumed parent inputs')
    else:path.write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':main()
