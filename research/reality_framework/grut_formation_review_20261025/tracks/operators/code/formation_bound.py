#!/usr/bin/env python3
"""Item 4 (production context). Own solver: rerun the prepared radial formation at the
P025_time settings (h=0.25, dt=0.015625, L=450, R=360, sigma=20, alpha=0.1) and apply my
re-implementation of the frozen-snapshot bound split as a TIME SERIES, next to the
charge inside r=60. Independent of Codex's code (initial data rebuilt from REPORT/run.py
description: x = 0.18 exp(-r^2/(2*60^2)), V = -i U, S static equilibrium given x).
"""
import json, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from evolve import Evolver, frozen_bound_split
from lattice import L_matrix

OUT = Path(__file__).resolve().parents[1] / 'results'
SCR = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/operators')


def static_S(r, f, h, n):
    """-L Z + r (3 (S^2-1) + f^2) S = 0 with S = 1 + Z/r (discrete equilibrium), Newton."""
    import scipy.sparse as sp, scipy.sparse.linalg as spl
    L = sp.csr_matrix(L_matrix(n, h))
    Z = r * (np.sqrt(1 - f * f / 3) - 1)
    for _ in range(30):
        S = 1 + Z / r
        F = -L @ Z + r * (3 * (S * S - 1) + f * f) * S
        if np.max(abs(F)) < 1e-13:
            break
        J = -L + sp.diags(9 * S * S - 3 + f * f)
        Z = Z - spl.spsolve(J.tocsc(), F)
    return Z


def main(dt=.015625, T=12000., h=.25, L=450., R=360., tag='P025_time_own'):
    n = int(round(L / h))
    ev = Evolver(n, h, R)
    r = ev.r
    f = .18 * np.exp(-r ** 2 / (2 * 60. ** 2))
    U = np.zeros((3, n + 1)); V = np.zeros((3, n + 1))
    U[0] = r * f; V[1] = -r * f
    Z = static_S(r[1:n], f[1:n], h, n)
    U[2, 1:n] = Z
    U[:, 0] = U[:, n] = 0; V[:, 0] = V[:, n] = 0
    Y = [U, V, np.zeros((3, n)), np.zeros((3, n + 1))]
    Q0 = ev.Q(U, V); E0 = ev.energy(U, V)
    steps = int(round(T / dt))
    rec_every = int(round(10 / dt))
    rows = []
    t0 = time.time()
    for j in range(steps + 1):
        t = j * dt
        if j % rec_every == 0:
            U, V = Y[0], Y[1]
            row = dict(t=t, Q=ev.Q(U, V), Q60=ev.Q_within(U, V, 60.), Q30=ev.Q_within(U, V, 30.))
            if t >= 8000 or j % (20 * rec_every) == 0:
                d = frozen_bound_split(ev, U, V)
                row.update(Qbound=d['Qbound'], omegas=d['omegas'], mode_charges=d['mode_charges'], inside60=d['inside60'])
            rows.append(row)
            if j % (100 * rec_every) == 0:
                print(f"t={t:.0f} Q60/Q0={row['Q60']/Q0:.6f} Qb={row.get('Qbound')} om={row.get('omegas')} wall={time.time()-t0:.0f}s", flush=True)
        if j < steps:
            Y = ev.step(Y, dt)
    U, V = Y[0], Y[1]
    late = [x for x in rows if x['t'] >= 9600]
    out = dict(tag=tag, h=h, dt=dt, L=L, R=R, T=T, Q0=Q0, E0=E0,
               late_mean_Q60_fraction=float(np.mean([x['Q60'] for x in late]) / Q0),
               final_Q60_fraction=float(rows[-1]['Q60'] / Q0),
               final_Q_physical_field=float(rows[-1]['Q']), rows=rows)
    qb = np.array([x['Qbound'] for x in late]); q60 = np.array([x['Q60'] for x in late])
    nb = np.array([len(x['omegas']) for x in late])
    out['late_window_summary'] = dict(
        Qbound_min=float(qb.min()), Qbound_max=float(qb.max()), Qbound_mean=float(qb.mean()),
        Q60_min=float(q60.min()), Q60_max=float(q60.max()),
        Qbound_over_Q60_min=float((qb / q60).min()), Qbound_over_Q60_max=float((qb / q60).max()),
        fraction_snapshots_with_1_bound_mode=float(np.mean(nb == 1)),
        fraction_snapshots_with_2_bound_modes=float(np.mean(nb == 2)),
        fraction_snapshots_with_3plus=float(np.mean(nb >= 3)),
        final_snapshot=dict(Qbound=rows[-1]['Qbound'], omegas=rows[-1]['omegas'], mode_charges=rows[-1]['mode_charges']))
    (OUT / f'formation_bound_{tag}.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out['late_window_summary'], indent=1), out['late_mean_Q60_fraction'], out['final_Q60_fraction'])


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'quick':
        main(dt=.03125, T=400., tag='quick')
    else:
        main()
