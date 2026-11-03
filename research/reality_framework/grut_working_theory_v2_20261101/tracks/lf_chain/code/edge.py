#!/usr/bin/env python3
"""pr_principle (Claude): candidate E, the edge of stability of the held hidden state along the holder line.
Held state = all-direction stationary point at R = 12 (constant-H inflation), from fx_holder's fc.stationary via the
e2e_chain skeptic's import path (read-only). u = mu + c0 * G * e^{i delta}, G = 10i - mu (the M2 definition).
We track the stationary branch from c0 = 1 upward and record the Hessian spectrum; the edge is where the lightest
non-gauge eigenvalue crosses zero (or the held branch is lost). Labels: [C] on the model."""
import sys, json, math
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'e2e_chain' / 'skeptic' / 'code'))
import sk_common as S  # noqa: E402  (read-only; scratch redirected by sk_common)
mc = S.E.mc
fc, P, MU = mc.fc, mc.P05, mc.MU
G = 10j - MU

def scan(delta_deg, cs):
    out, prev = [], None
    for c in cs:
        u = MU + c * G * np.exp(1j * math.radians(delta_deg))
        st = fc.stationary(u, P, R=12.0, guess=prev)
        if st['grad'] > 1e-6:
            st = fc.stationary(u, P, R=12.0)
        ev = np.sort(st['eig'])
        out.append(dict(c0=float(c), abs_u=abs(u), arg_u=math.degrees(np.angle(u)), grad=st['grad'], V=st['V'],
                        eig=[float(x) for x in ev], amps=[float(a) for a in st['amps']]))
        prev = st['z']
    return out

res = {}
cs = np.round(np.arange(1.0, 2.6001, 0.02), 3)
for d in (-4.0, 0.0, 4.0):
    rows = scan(d, cs)
    res[f'{d:+.0f}'] = rows
    print(f'delta {d:+.0f}:')
    for r in rows[::5]:
        print('  c0 %.2f |u| %6.2f arg %6.1f grad %.1e eig %s amps %s' % (r['c0'], r['abs_u'], r['arg_u'], r['grad'],
              ' '.join('%.3g' % e for e in r['eig'][:3]), ' '.join('%.2f' % a for a in r['amps'])))
json.dump(res, open(HERE.parent / 'results' / 'edge.json', 'w'), indent=1)
