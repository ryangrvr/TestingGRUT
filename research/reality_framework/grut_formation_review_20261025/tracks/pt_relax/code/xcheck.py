"""Track pt_relax (Claude): integrator cross-check of the late, UV-filled regime.

From the saved s = 894 state of the long IR12_s1_w1 run, advance one 10-s window three ways:
  (a) Codex's FD4/Yoshida engine, hmax 0.004 (as in the long run; reproduction tool);
  (b) the same engine, hmax 0.002 (time-step control);
  (c) our own independent spectral-FD4 Blanes-Moan SRKN6b integrator (tracks/en_indep/code/indep3d.py), h 0.004.
Compare the windowed negative-frequency action (own estimator, long_run.Obs.window), energy and final fields.
"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'en_indep', 'code'))
import long_run as LR          # own estimator; imports Codex engine from the private copy
import pt_common as P
import indep3d as I            # own independent integrator (en_indep)

RUN = os.environ['RUNDIR']
R = np.load(os.path.join(RUN, 'IR12_s1_w1_restart_s894.npz'))
z0 = np.ascontiguousarray(R['z']); nQ0 = float(R['nQ0']); s0 = float(R['s'])
N, L = 32, 12.
ob = LR.Obs(N, L, nQ0)
out = {}
finals = {}
for name in ['codex_h004', 'codex_h002', 'bm6_h004']:
    t = time.monotonic()
    z = z0.copy()
    buf = np.zeros((LR.M, N, N, N), complex)
    row, psi = ob.equal_time(z); buf[0] = psi; rho0 = row['rho']
    if name.startswith('codex'):
        P.set_box(L); lib = P.library(); h = 0.004 if name.endswith('004') else 0.002
        for j in range(1, LR.M):
            P.step(lib, z, LR.DS, h, threads=2)
            row, psi = ob.equal_time(z); buf[j] = psi
    else:
        lat = I.Lattice(N, L); q = z[:3].copy(); p = z[3:].copy(); F = None
        for j in range(1, LR.M):
            q, p, F = I.step_bm6(lat, q, p, LR.DS / 10, 10, F)
            z = np.ascontiguousarray(np.r_[q, p])
            row, psi = ob.equal_time(z); buf[j] = psi
    w = ob.window(buf); w.pop('shell_minus')
    w.update(energy_rel_change=row['rho'] / rho0 - 1, Q_rel_change=row['Q'] / ob.nQ0 - 1, seconds=time.monotonic() - t)
    out[name] = w; finals[name] = z.copy()
    print(name, {k: (round(v, 7) if isinstance(v, float) else v) for k, v in w.items() if k in
                 ['coherent_minus', 'band_minus', 'uv_minus', 'total_minus_all', 'energy_rel_change', 'seconds']}, flush=True)
for a, b in [('codex_h004', 'codex_h002'), ('codex_h004', 'bm6_h004')]:
    d = finals[a] - finals[b]
    rms = np.sqrt(np.mean(finals[a] ** 2, axis=(1, 2, 3)))
    out[f'field_maxdiff_over_rms_{a}_vs_{b}'] = (np.max(np.abs(d), axis=(1, 2, 3)) / rms).tolist()
    out[f'total_minus_all_reldiff_{a}_vs_{b}'] = out[a]['total_minus_all'] / out[b]['total_minus_all'] - 1
    print(a, b, out[f'field_maxdiff_over_rms_{a}_vs_{b}'], out[f'total_minus_all_reldiff_{a}_vs_{b}'])
out['window_s'] = [s0, s0 + 10]
out['long_run_value_same_window'] = dict(coherent_minus=0.18551, band_minus=0.54563, total_minus_all=0.88569)
json.dump(out, open(os.path.join(HERE, '..', 'results', 'xcheck_integrators.json'), 'w'), indent=1)
