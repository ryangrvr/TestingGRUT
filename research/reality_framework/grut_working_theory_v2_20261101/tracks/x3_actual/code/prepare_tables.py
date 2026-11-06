#!/usr/bin/env python3
"""Rebuild exactly the four archived tables, without importing other solvers."""
import concurrent.futures
import json
import os
import sys
from pathlib import Path

work=Path(sys.argv[1]).resolve()
os.environ['EK_SCRATCH']=str(work/'original_ptables')
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'e2e_kinetic/code'))
import ek_common as C
import ek_ptable as P

C.compile_engine()
reference=json.loads((C.RES/'ptable_summary.json').read_text())['tables']
masses=[.41667,1.6,.8,.48]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    values=list(pool.map(P.build,masses))
checks={}
for mass,rows in zip(masses,values):
    archived=reference[f'{mass:.5f}']
    assert len(rows)==len(archived)
    checks[str(mass)]=max(abs(a[k]-b[k])/max(abs(b[k]),1e-30)
                         for a,b in zip(rows,archived) for k in a)
result={'passed':max(checks.values())<2e-6,'relative_errors':checks}
(C.SP/'check.json').write_text(json.dumps(result,indent=2)+'\n')
print(result)
if not result['passed']:raise SystemExit(1)
