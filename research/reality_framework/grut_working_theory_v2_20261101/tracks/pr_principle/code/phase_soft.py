#!/usr/bin/env python3
"""pr_principle (Claude): candidate B (maximal CP torque: held X-Sigma relative phase = -90 deg) and the post-hoc
'softest held mode' point. Same solver and holder line as edge.py. [C] on the model."""
import sys, json, math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'e2e_chain' / 'skeptic' / 'code'))
import sk_common as S  # noqa: E402
mc = S.E.mc; fc, P, MU = mc.fc, mc.P05, mc.MU
G = 10j - MU
def st(c, d):
    return fc.stationary(MU + c * G * np.exp(1j * math.radians(d)), P, R=12.0)
out = {'torque': {}, 'soft': {}}
# B: relative phase (stationary 'delta' = arg z_Sigma - 6 arg z_X) as a function of delta at fixed c0
for c in (1.0, 1.40, 1.45, 1.50, 1.55, 1.60):
    f = lambda d: math.degrees(st(c, d)['delta']) + 90.0
    ds = np.arange(-40, 40.1, 5.0); vals = [f(d) for d in ds]
    root = None
    for i in range(len(ds) - 1):
        if vals[i] * vals[i + 1] < 0:
            root = brentq(f, ds[i], ds[i + 1], xtol=1e-4); break
    out['torque'][f'{c:.2f}'] = dict(rel_phase_at_delta0=math.degrees(st(c, 0.0)['delta']), delta_max_torque=root)
    print('c0 %.2f: rel phase at delta=0: %.2f deg; max-torque delta = %s' % (c, math.degrees(st(c, 0.0)['delta']), None if root is None else '%.2f' % root))
# post-hoc: minimum of the lightest non-gauge held eigenvalue along c0 (delta = 0, +-4)
for d in (-4.0, 0.0, 4.0):
    cs = np.arange(0.9, 1.5001, 0.01); lam = []
    for c in cs:
        ev = np.sort(st(c, d)['eig']); lam.append(ev[1])
    i = int(np.argmin(lam))
    out['soft'][f'{d:+.0f}'] = dict(c0_min=float(cs[i]), m2_min=float(lam[i]))
    print('delta %+.0f: softest held mode m2 = %.3f H^2 at c0 = %.2f' % (d, lam[i], cs[i]))
json.dump(out, open(HERE.parent / 'results' / 'phase_soft.json', 'w'), indent=1)
