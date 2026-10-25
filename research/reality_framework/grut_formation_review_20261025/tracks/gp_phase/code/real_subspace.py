#!/usr/bin/env python3
"""gp_phase (a): numerical checks of the real-subspace statement with our own integrator
(release.c, conformal variables, Dormand-Prince 5(4)), independent of Codex's full_release.py.

1. Exactly real initial data (own-potential alignment, zero velocity): all charges vanish
   identically, with and without an extra operator whose coefficient is real (delta = 0, pi).
2. Codex's misaligned alternatives (X phase -/+ pi/24) recomputed at a = 14.177 (Codex t=100):
   comparison with RELEASE_RESULTS.json numbers quoted below.
3. Stability of the real subspace: a CP-odd seed e0 (X phase) grows; a3 n_X / e0 versus tau
   for several e0 measures the transverse (CP-odd) growth rate in the linear regime.
"""
import json
from pathlib import Path
import numpy as np
import gauge_point as g

OUT = Path(__file__).resolve().parents[1] / 'results'
Z = np.array([23.824169767032206 + 0j, -1.891386373567112 + 0j, -0.26697595666823387 + 0j])
CODEX_PI24 = [-16.261387466138004, 6.395904448609672, -0.4877209799485097]  # RELEASE_RESULTS.json rows[1]


def main():
    out = dict(scope='Homogeneous classical X,Y,Sigma; Codex release model (Hubble mass off, radiation era).')
    # 1. exact zero
    zero = []
    for op, delta in [(0, 0.0), (4, 0.0), (4, np.pi), (2, 0.0), (1, np.pi)]:
        if op:
            Cabs = g.coefficient_for_ratio(op, 0.1 if op == 4 else 1.0, Z)
            C = g.SIGN[op] * Cabs * np.exp(1j * delta)
            C = complex(np.real(C), 0.0)  # real coefficient: no CP phase
            z, v, gn, nd = g.inflation_minimum(op, C)
            z = np.array([complex(x.real, 0.0) for x in z])  # real minimum (imag parts are roundoff)
        else:
            C, z = 0j, Z.copy()
        rows, end = g.run_release(z, op, C, tau_end=1000.0, rtol=1e-11, atol=1e-13, nsamp=30)
        zero.append(dict(op=op, delta=delta, max_abs_a3n=max(max(abs(r['a3nX']), abs(r['a3nY']), abs(r['a3nS'])) for r in rows),
                         end=end))
    out['exact_zero_runs'] = zero
    # 2. Codex comparison
    comp = []
    for ang in (np.pi / 24, -np.pi / 24):
        zz = Z.copy(); zz[0] = Z[0] * np.exp(-1j * ang)
        rows, end = g.run_release(zz, tau_end=14.177446878757825, nsamp=10, rtol=1e-12, atol=1e-14)
        mine = [rows[-1]['a3nX'], rows[-1]['a3nY'], rows[-1]['a3nS']]
        sign = 1 if ang > 0 else -1
        comp.append(dict(angle=ang, mine=mine, codex=[sign * c for c in CODEX_PI24],
                         max_rel_diff=max(abs(m - sign * c) / abs(c) for m, c in zip(mine, CODEX_PI24)), end=end))
    out['codex_pi24_comparison'] = comp
    # 3. transverse growth
    growth = {}
    for e0 in (1e-15, 1e-12, 1e-9, 1e-6):
        zz = Z.copy(); zz[0] = Z[0] * np.exp(-1j * e0)
        rows, end = g.run_release(zz, tau_end=300.0, nsamp=60, rtol=1e-12, atol=1e-14)
        growth[repr(e0)] = dict(tau=[r['tau'] for r in rows], a3nX_over_e0=[r['a3nX'] / e0 for r in rows],
                                epsX=[r['epsX'] for r in rows])
    out['cp_odd_growth'] = growth
    # growth-rate fit on the e0=1e-15 envelope, linear regime |a3 n_X| < 1
    t = np.array(growth['1e-15']['tau']); y = np.abs(np.array(growth['1e-15']['a3nX_over_e0'])) * 1e-15
    lin = (t > 5) & (y < 1.0) & (y > 0)
    env_t, env_y = [], []
    for k in range(len(t)):
        if lin[k] and (not env_y or np.log(y[k]) > max(env_y)):
            env_t.append(t[k]); env_y.append(np.log(y[k]))
    if len(env_t) > 3:
        slope = np.polyfit(env_t, env_y, 1)[0]
        out['cp_odd_envelope_growth_rate_per_unit_conformal_time'] = float(slope)
        out['envelope_points'] = dict(tau=env_t, ln_abs_a3nX=env_y)
    (OUT / 'real_subspace.json').write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(dict(zero=[(z['op'], z['delta'], z['max_abs_a3n']) for z in zero],
                          codex=[c['max_rel_diff'] for c in comp],
                          growth_rate=out.get('cp_odd_envelope_growth_rate_per_unit_conformal_time')), indent=1))
    for k, v in growth.items():
        print(k, ' '.join('%.3g:%.2e' % (a, b) for a, b in list(zip(v['tau'], v['a3nX_over_e0']))[::5]))


if __name__ == '__main__':
    main()
