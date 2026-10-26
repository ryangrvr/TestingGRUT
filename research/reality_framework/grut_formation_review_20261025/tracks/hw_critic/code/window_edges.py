#!/usr/bin/env python3
"""hw_critic (Claude, own code): lower edge m_edge of the heavy window (R_central(m) = R_L) for each threshold variant,
from the hw_verify ladder R_central(m) [Q, own track] interpolated as a power law between rungs, for the mediator
debris alone and with the inflationary vector piece R_A(e) added (R_A = 48.47 (1e-7/e), m-independent at fixed yield
normalisation to within the 2 % spread of E_X along the ladder [E]). Stability edge m = 42.8 (kappa = kappa_max) [Q].
Usage: python3 -B window_edges.py ../results/window_edges.json
"""
import json, math, sys
from pathlib import Path
LADDER = [(10, 144.52), (15, 34.234), (20, 8.9396), (30, 1.1440), (40, 0.25109)]
THR = json.load(open(Path(__file__).resolve().parents[1] / 'results' / 'thresholds.json'))['thresholds']


def m_edge(target):
    for (m1, r1), (m2, r2) in zip(LADDER, LADDER[1:]):
        if r1 >= target > r2:
            s = math.log(r2 / r1) / math.log(m2 / m1)
            return m1 * (target / r1) ** (1 / s)
    return None if target < LADDER[-1][1] else float('nan')


out = {}
for k, RL in THR.items():
    row = dict(R_L=RL, mediator_only=m_edge(RL))
    for e in (1e-7, 1e-6, 3e-6, 1e-5):
        RA = 48.47 * (1e-7 / e)
        row['with_vector_e=%.0e' % e] = m_edge(RL - RA) if RL > RA else 'closed (R_A = %.1f > R_L)' % RA
    out[k] = row
Path(sys.argv[1]).write_text(json.dumps(out, indent=1))
for k, v in out.items():
    print('%-26s R_L %6.2f  m_edge(med) %5s  ' % (k, v['R_L'], ('%.1f' % v['mediator_only']) if v['mediator_only'] else 'none'),
          '  '.join('%s: %s' % (kk.split('=')[1], (('%.1f' % vv) if isinstance(vv, float) else vv)) for kk, vv in v.items() if kk.startswith('with')))
