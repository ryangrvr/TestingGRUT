#!/usr/bin/env python3
"""Item 4, strongly excited cores (own code). Compare the frozen split's 'free' charge at time t,
Q(t) - Q_bound(t), with the charge that is actually free in the sense of escaping by the end of
the run: Q(t) - Q_r(T_end) (charge on the grid now minus what is still inside radius r at T_end).
Cases: soliton amplitude x1.3, x0.75 (in phase), and a prepared Gaussian that collapses.
"""
import json, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from soliton import solve_fixed_norm
from evolve import Evolver, frozen_bound_split
from check_bound import embed

OUT = Path(__file__).resolve().parents[1] / 'results'


def run(kind, T=4000., dt=.025, L=300., R=240., h=.25, every=5.0, sol=None, scale=1.0):
    n = int(round(L / h))
    ev = Evolver(n, h, R)
    r = ev.r
    if kind == 'soliton_scaled':
        U, V = embed(sol, n, h)
        U[0] *= scale; V[1] *= scale
    else:   # prepared Gaussian, V = -i U (as in the packet's preparation, smaller and denser)
        f = .6 * np.exp(-r ** 2 / (2 * 12. ** 2))
        U = np.zeros((3, n + 1)); V = np.zeros((3, n + 1))
        U[0] = r * f; V[1] = -r * f
        U[:, 0] = U[:, n] = 0; V[:, 0] = V[:, n] = 0
    Y = [U, V, np.zeros((3, n)), np.zeros((3, n + 1))]
    rows = []
    steps = int(round(T / dt)); rec = int(round(every / dt))
    for j in range(steps + 1):
        if j % rec == 0:
            U, V = Y[0], Y[1]
            d = frozen_bound_split(ev, U, V)
            rows.append(dict(t=j * dt, Q=ev.Q(U, V), Q30=ev.Q_within(U, V, 30.), Q60=ev.Q_within(U, V, 60.),
                             Qbound=d['Qbound'], nb=len(d['omegas']), omegas=d['omegas'], mode_charges=d['mode_charges']))
        if j < steps:
            Y = ev.step(Y, dt)
    t = np.array([x['t'] for x in rows]); Q = np.array([x['Q'] for x in rows]); Qb = np.array([x['Qbound'] for x in rows])
    Q30 = np.array([x['Q30'] for x in rows])
    split_free = Q - Qb
    escaped_by_end = Q - Q30[-1]               # charge on grid now that is not inside r<30 at T_end
    win = (t >= 1000) & (t <= T - 1000)
    err = split_free[win] - escaped_by_end[win]
    # late drift of the core charge (to show it is quasi-stationary)
    late = t >= T - 1000
    slope = np.polyfit(t[late], Q30[late], 1)[0]
    return dict(kind=kind, scale=scale, T=T, Q0=float(Q[0]), Q30_end=float(Q30[-1]),
                window='1000 <= t <= T-1000',
                split_free_window_min=float(split_free[win].min()), split_free_window_max=float(split_free[win].max()),
                escaped_by_end_window_min=float(escaped_by_end[win].min()), escaped_by_end_window_max=float(escaped_by_end[win].max()),
                split_minus_escaped_min=float(err.min()), split_minus_escaped_max=float(err.max()),
                split_minus_escaped_rms=float(np.sqrt(np.mean(err ** 2))),
                Qbound_over_Q30_window_min=float((Qb[win] / Q30[win]).min()), Qbound_over_Q30_window_max=float((Qb[win] / Q30[win]).max()),
                n_bound_modes_seen=sorted(set(int(x['nb']) for x in rows if x['t'] >= 1000)),
                late_dQ30_dt=float(slope),
                series_every_50=[dict(t=x['t'], Q=x['Q'], Q30=x['Q30'], Qbound=x['Qbound'], nb=x['nb'], omegas=x['omegas']) for x in rows[::10]])


def main():
    sol = solve_fixed_norm(1029.)
    out = {}
    t0 = time.time()
    for kind, scale in [('soliton_scaled', 1.3), ('soliton_scaled', .75), ('gaussian_collapse', 1.0)]:
        key = f'{kind}_{scale}'
        out[key] = run(kind, sol=sol, scale=scale)
        print(key, {k: v for k, v in out[key].items() if not k.startswith('series')}, time.time() - t0, flush=True)
        (OUT / 'bound_split_strong.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
