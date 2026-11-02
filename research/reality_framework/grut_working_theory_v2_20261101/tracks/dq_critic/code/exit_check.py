#!/usr/bin/env python3
"""dq_critic quick check (Claude, own driver, 2 Nov 2026) of the one contradiction between dq_heat and its skeptic:
does a Hubble rate that grows with phi_I cancel the horizon-exit sensitivity of the vector heat?

dq_heat [E]: factor s^2 (F(1)/F(s))^2 = 1.06-1.15 (fixed-H condensate).  SKEPTIC [C]: re-solve the held minimum at
R = 12 h_r^2 -> f = [h_r F(1,1)/F(s, h_r)]^2 = 0.77-0.46 along H ~ phi_I.
Here: dq_heat's own model (code/dq_model.py, imported read-only, unchanged) at lambda = 0.5, stationary point at
u = mu + s (10i - mu), R = 12 h_r^2, F = (2 sum q^2 |z|^2)^{1/2}; plus the physical reason: |X|^2 at the held minimum
for xi_X R |X|^2 with xi_X = -5/6, i.e. |X| ~ H when X dominates its own stabilisation.
Output: ../results/exit_check.json
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', 'dq_heat', 'code')))
import dq_model as fc

Q = np.array([1.0, 3.0, 6.0])


def F_at(p, s, hr, guess=None):
    u = p['mu'] + (fc.U_PREP - p['mu']) * s
    st = fc.stationary(u, p, R=12.0 * hr * hr, guess=guess)
    return float(np.sqrt(2 * np.sum(Q ** 2 * np.abs(st['z']) ** 2))), st


def main():
    p = fc.params(0.5)
    F11, st0 = F_at(p, 1.0, 1.0)
    out = dict(method=__doc__, F11=F11, z11=np.abs(st0['z']).tolist(), rows=[])
    g = {}
    for s, hr in ((1.0, 1.5), (1.0, 2.0), (1.5, 1.0), (2.0, 1.0), (1.5, 1.5), (2.0, 2.0)):
        F, st = F_at(p, s, hr, g.get(s))
        g[s] = st['z']
        Ffix, _ = F_at(p, s, 1.0, None) if hr != 1.0 else (F, st)
        row = dict(s=s, h_r=hr, absz=np.abs(st['z']).tolist(), F=F, F_over_Hexit=F / hr, f=(hr * F11 / F) ** 2,
                   dq_heat_estimate=(s * s * (F11 / Ffix) ** 2) if s == hr else None, grad=float(st['grad']),
                   lightest_m2_over_Hexit2=float(sorted(st['eig'])[1] / hr ** 2))
        out['rows'].append(row)
        print('s=%.2f h_r=%.2f |X|,|Y|,|S| = %s  F/H_exit = %.2f  f = %.3f  dq_heat est = %s' % (
            s, hr, ' '.join('%.2f' % x for x in row['absz']), row['F_over_Hexit'], row['f'],
            '%.3f' % row['dq_heat_estimate'] if row['dq_heat_estimate'] else '-'))
    json.dump(out, open(os.path.join(HERE, '..', 'results', 'exit_check.json'), 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
