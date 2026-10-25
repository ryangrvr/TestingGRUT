#!/usr/bin/env python3
"""Driver for radial_spec.c (fast port of spectral_radial.py).  Prepares the initial data with the
Python solver (prepared Gaussian + spectral static neutral solve), runs the C binary, then post-processes
the history and the snapshots (exterior particle/antiparticle projection with Radial.mode_numbers).

usage: drive_c.py run  --h .5 --L 6400 --dt .0625 --tmax 12000 --out PREFIX [--scheme 4] [--sponge Rs,D,g,p]
       drive_c.py post --h .5 --L 6400 --out PREFIX
"""
import argparse, json, subprocess, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from spectral_radial import Radial, T

RAD = [30, 60, 120, 240, 480]
COLS = ['tau', 'E', 'Q', 'Erec', 'Qrec', 'E_local_sum', 'S0']
for R in RAD:
    COLS += [f'Q{R}', f'Q{R}_trap', f'E{R}', f'Ech_out{R}', f'Ene_out{R}', f'pad{R}']
BIN = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/independent_solver/radial_spec')


def run(a):
    sim = Radial(a.h, a.L, a.dt, 'yoshida4' if a.scheme == 4 else 'strang', True, None)
    sim.prepared_gaussian(A=a.amp)
    q, p = sim.physical()
    init = Path(a.out + '_init.bin')
    np.concatenate([q.ravel(), p.ravel()]).astype('<f8').tofile(init)
    sp = [float(x) for x in a.sponge.split(',')] if a.sponge else [-1, 1, 0, 1]
    cmd = [str(BIN), str(init), a.out, repr(a.h), repr(a.L), repr(a.dt), repr(a.tmax), str(a.scheme),
           repr(a.diag)] + [repr(x) for x in sp] + [a.snaps or '-', repr(a.ckpt)]
    t0 = time.time()
    with open(a.out + '.log', 'a') as log:
        log.write(' '.join(cmd) + '\n'); log.flush()
        subprocess.run(cmd, stderr=log, check=True)
    Path(a.out + '_meta.json').write_text(json.dumps(dict(args=vars(a), seconds=time.time() - t0,
                                                          static_iterations=sim.static_iterations)))
    post(a)


def history(prefix):
    H = np.fromfile(prefix + '_history.bin', dtype='<f8')
    return H.reshape(-1, len(COLS))


def post(a):
    H = history(a.out)
    c = {k: i for i, k in enumerate(COLS)}
    tau = H[:, 0]; E0, Q0 = H[0, c['E']], H[0, c['Q']]
    res = dict(args=vars(a), E0=E0, Q0=Q0, tau_end=float(tau[-1]),
               max_rel_energy_error=float(np.max(np.abs((H[:, c['E']] + H[:, c['Erec']]) / E0 - 1))),
               max_rel_charge_error=float(np.max(np.abs((H[:, c['Q']] + H[:, c['Qrec']]) / Q0 - 1))),
               max_rel_local_energy_mismatch=float(np.max(np.abs(H[:, c['E_local_sum']] / H[:, c['E']] - 1))))
    for lo, hi in [(9600, 12000), (19200, 24000), (28800, 36000), (38400, 48000), (11000, 12000)]:
        m = (tau > lo) & (tau <= hi + 1e-9)
        if m.sum() > 2 and tau[-1] >= hi - 1e-6:
            y = H[m, c['Q60']] / Q0
            res[f'Q60_mean_{lo}_{hi}'] = float(y.mean())
            res[f'Q60trap_mean_{lo}_{hi}'] = float(np.mean(H[m, c['Q60_trap']]) / Q0)
            res[f'Q60_std_{lo}_{hi}'] = float(y.std())
            res[f'Q60_slope_per_1000_{lo}_{hi}'] = float(np.polyfit(tau[m], y, 1)[0] * 1000)
    # samples of the escaped quantities (region r > 240) and the receivers
    samples = {}
    for t in [2000, 4000, 6000, 8000, 10000, 12000, 16000, 20000, 24000, 30000, 36000, 42000, 48000]:
        i = np.argmin(abs(tau - t))
        if abs(tau[i] - t) < 1e-6:
            row = H[i]
            samples[t] = dict(Q60=row[c['Q60']] / Q0, Q30=row[c['Q30']] / Q0, Q120=row[c['Q120']] / Q0,
                              Q_out240=(row[c['Q']] - row[c['Q240']]), Q_out240_plus_rec=(row[c['Q']] - row[c['Q240']] + row[c['Qrec']]),
                              E_out240=(row[c['E']] - row[c['E240']]), E_out240_plus_rec=(row[c['E']] - row[c['E240']] + row[c['Erec']]),
                              Ech_out240=row[c['Ech_out240']], Ene_out240=row[c['Ene_out240']],
                              Ech_out480=row[c['Ech_out480']], Ene_out480=row[c['Ene_out480']],
                              Q_out480=(row[c['Q']] - row[c['Q480']]), Erec=row[c['Erec']], Qrec=row[c['Qrec']])
    res['samples'] = samples
    # snapshots: exterior projection
    modes = []
    sim = Radial(a.h, a.L, a.dt, 'strang', True, None)
    for f in sorted(Path(a.out).parent.glob(Path(a.out).name + '_snap*.bin')):
        d = np.fromfile(f, dtype='<f8')
        n = sim.n
        q = d[:3 * n].reshape(3, n); p = d[3 * n:6 * n].reshape(3, n)
        sim.Qh, sim.Ph = T(q), T(p); sim.tau = d[6 * n]
        for rmin in (240.0, 480.0, 1000.0):
            if rmin < a.L - 100:
                m = sim.mode_numbers(rmin); m['tau'] = float(d[6 * n]); modes.append(m)
        dg = sim.diagnostics()
        # spectral tail: free-energy fraction carried by sine modes above fractions of k_max
        e_k = 0.5 * (sim.Ph ** 2 + sim.omega ** 2 * sim.Qh ** 2)
        kmax = sim.k[-1]
        tail = {f'frac_E_k_gt_{f}kmax': float(e_k[:, sim.k > f * kmax].sum() / e_k.sum()) for f in (0.25, 0.5, 0.8)}
        modes.append(dict(tau=float(d[6 * n]), check_Q60_python=dg['Q60'] / Q0, check_E_python=dg['E'],
                          kmax=float(kmax), **tail))
    res['modes'] = modes
    Path(a.out + '_RESULT.json').write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k not in ('samples', 'modes', 'args')}, indent=1))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['run', 'post'])
    ap.add_argument('--h', type=float, default=0.5)
    ap.add_argument('--L', type=float, default=6400.0)
    ap.add_argument('--dt', type=float, default=0.0625)
    ap.add_argument('--tmax', type=float, default=12000.0)
    ap.add_argument('--scheme', type=int, default=4)
    ap.add_argument('--diag', type=float, default=4.0)
    ap.add_argument('--sponge', default='')
    ap.add_argument('--snaps', default='')
    ap.add_argument('--ckpt', type=float, default=0.0)
    ap.add_argument('--amp', type=float, default=0.18)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    run(a) if a.mode == 'run' else post(a)
