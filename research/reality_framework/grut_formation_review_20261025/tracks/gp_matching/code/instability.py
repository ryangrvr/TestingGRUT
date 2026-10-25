#!/usr/bin/env python3
"""Angular instability of the homogeneous release and its consequence for the phase -> charge map.

1. Linear regime: delta0 = 1e-14 (the real subspace is exactly invariant, so the imaginary components
   evolve linearly without round-off seeding). Record N_X/delta0 and |Im u|/delta0 versus a, two tolerances.
2. Sensitivity: pairs of runs at O(1) delta0 differing by 1e-9 rad, and the same delta0 at two tolerances;
   record when they decorrelate.
3. Ensemble: 40 delta0 in [0.6, 1.0] -> distribution of the late X charge fraction eps_X (window a=600-1000).
Outputs results/instability.json.
"""
import json, sys
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from release import params, rhs_factory, initial, run, summarize, OUT

POINT = json.loads((OUT / 'point_full.json').read_text())['baseline']
P = params()


def linear(rtol):
    d0 = 1e-14
    s0 = initial(POINT, d0)
    a_eval = np.concatenate([np.linspace(1, 300, 3000), np.linspace(300.5, 1000, 1400)])
    sol = solve_ivp(rhs_factory(P), (1, 1000), s0, method='DOP853', rtol=rtol, atol=1e-30, t_eval=a_eval)
    im = np.sqrt(sol.y[1]**2 + sol.y[3]**2 + sol.y[5]**2 + sol.y[7]**2 + sol.y[9]**2 + sol.y[11]**2) / d0
    NX = -2 * (sol.y[0] * sol.y[7] - sol.y[1] * sol.y[6]) / d0
    return sol.t, im, NX


def pair(args):
    d0, kind = args
    if kind == 'tol':
        h1, _ = run(POINT, P, d0, 1000., rtol=1e-10, n_eval=2000)
        h2, _ = run(POINT, P, d0, 1000., rtol=1e-12, n_eval=2000)
    else:
        h1, _ = run(POINT, P, d0, 1000., rtol=1e-11, n_eval=2000)
        h2, _ = run(POINT, P, d0 + 1e-9, 1000., rtol=1e-11, n_eval=2000)
    rel = np.abs(h1[:, 1] - h2[:, 1]) / (np.abs(h1[:, 1]) + np.abs(h2[:, 1]) + 1e-300)
    j = np.argmax(rel > 0.1) if np.any(rel > 0.1) else -1
    w1 = summarize(h1)['windows']['600-1000']; w2 = summarize(h2)['windows']['600-1000']
    return dict(delta0=d0, kind=kind, a_decorrelate=float(h1[j, 0]) if j >= 0 else None,
                eps_X_1=w1['eps_X_mean'], eps_X_2=w2['eps_X_mean'])


def ens(d0):
    h, _ = run(POINT, P, d0, 1000., rtol=1e-10, n_eval=1500)
    w = summarize(h)['windows']['600-1000']
    return dict(delta0=d0, eps_X=w['eps_X_mean'], N_X=w['N_X_mean'], N_eff=w['N_eff_mean'], EY_over_EX=w['EY_over_EX_mean'])


def main():
    out = {}
    lin = {}
    for rtol in [1e-10, 1e-12]:
        t, im, NX = linear(rtol)
        sel = {a: int(np.argmin(abs(t - a))) for a in [10, 20, 35, 50, 70, 100, 150, 200, 300, 500, 1000]}
        lin[str(rtol)] = dict(Im_norm_over_delta0={str(a): float(im[j]) for a, j in sel.items()},
                              N_X_over_delta0={str(a): float(NX[j]) for a, j in sel.items()})
        k1, k2 = sel[35], sel[150]
        lin[str(rtol)]['mean_growth_rate_per_conformal_time_35_150'] = float(np.log(im[k2] / im[k1]) / (t[k2] - t[k1]))
        print('linear', rtol, lin[str(rtol)], flush=True)
    out['linear'] = lin
    with Pool(4) as pool:
        pr = pool.map(pair, [(np.pi / 4, 'tol'), (np.pi / 2, 'tol'), (0.3, 'tol'), (2.0, 'tol'),
                             (np.pi / 4, 'shift'), (np.pi / 2, 'shift'), (0.3, 'shift'), (2.0, 'shift')])
        for r in pr:
            print('pair', r, flush=True)
        out['sensitivity_pairs'] = pr
        grid = np.linspace(0.6, 1.0, 40)
        er = pool.map(ens, grid)
    eps = np.array([r['eps_X'] for r in er])
    out['ensemble'] = dict(rows=er, frac_positive=float(np.mean(eps > 0)), mean_eps=float(eps.mean()),
                           rms_eps=float(np.sqrt(np.mean(eps**2))), median_abs_eps=float(np.median(abs(eps))),
                           min_abs_eps=float(abs(eps).min()), max_abs_eps=float(abs(eps).max()),
                           sign_changes=int(np.sum(np.diff(np.sign(eps)) != 0)))
    print('ensemble', {k: v for k, v in out['ensemble'].items() if k != 'rows'})
    (OUT / 'instability.json').write_text(json.dumps(out, indent=1, default=float) + '\n')


if __name__ == '__main__':
    main()
