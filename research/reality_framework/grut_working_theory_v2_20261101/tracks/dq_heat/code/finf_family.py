#!/usr/bin/env python3
"""dq_heat (Claude, own code): sensitivity of the inflationary Higgs norm F_inf along the prepared family.

The GMR modes near k_* left the horizon N_* = ln(1/k_*) = 2.7-3.8 e-folds before the slice (e = 1e-5 ... 1e-6). If the
inflaton was larger then (s = phi_I/phi_I(0) > 1, fixed H as in fx_holder's s-scan), the conditional minimum and hence F_inf
differ; the superhorizon Goldstone angle is H/(2 pi F_exit), so N_A scales as (F_inf(s=1)/F_inf(s))^2.  [C on the static
family; which s applies at horizon exit is a property of the supplied inflaton history, not computed here.]
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import dq_model as fc

p = fc.params(0.5)
Q = np.array([1.0, 3.0, 6.0])
rows = []
guess = None
for s in (1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0):
    u = p['mu'] + (fc.U_PREP - p['mu']) * s
    st = fc.stationary(u, p, guess=guess)
    guess = st['z']
    F = float(np.sqrt(2 * np.sum(Q ** 2 * np.abs(st['z']) ** 2)))
    rows.append(dict(s=s, amps=np.abs(st['z']).tolist(), F_inf=F, lightest_m2=float(sorted(st['eig'])[1]), grad=st['grad']))
F1 = rows[0]['F_inf']
for r in rows:
    r['N_A_factor_vs_s1'] = (F1 / r['F_inf']) ** 2
    print(r)
json.dump(dict(method=__doc__, rows=rows), open(Path(__file__).resolve().parent.parent / 'results' / 'finf_family.json', 'w'), indent=1)
