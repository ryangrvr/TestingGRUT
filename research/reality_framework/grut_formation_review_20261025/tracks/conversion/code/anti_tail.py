#!/usr/bin/env python3
"""Late net anti-X abundance after inverse conversion (Claude review, track `conversion`).

First-order estimate, own method.  A positive particle on its final approach to
the inverse threshold passes through delta = E_+ - (1+2 Omega) -> 0.  While it
does, an inverse event X -> anti-X occurs at rate lam_+(a, delta) and makes an
anti-X with p_- = sqrt(delta (delta+2)).  That anti-X survives to today with
probability S_-(a, p_-) = exp(-int lam_- dln a) along its own redshift.  Then
    f_anti(today) ~= < int lam_+ S_- dln a >   over escaping particles,
taking the whole band 0 < delta < delta_0 as an upper envelope (particles
created at smaller delta traverse less of it).  Second-order chains (anti-X ->
X -> anti-X again) are neglected; they are suppressed by further factors of
lam_+ S_-.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import transport_mc as T


def survival_anti(rate, a, p):
    q = a * p
    f = lambda y: T.lam_minus(rate, np.exp(y), q / np.exp(y))
    return np.exp(-quad(f, np.log(a), 0.0, epsrel=1e-8, limit=400)[0])


def tail_for_escape(rate, a_esc, delta0=np.sqrt(1.5) - 1):
    """Integrate backwards from the threshold crossing at a_esc to delta0."""
    # positive particle: E_+ = EB + delta, comoving q = a_esc * UESC
    q = a_esc * T.UESC

    def integrand(y):   # y = ln a
        a = np.exp(y)
        pp = q / a
        Ep = np.sqrt(1 + pp * pp)
        d = Ep - T.EB
        if d <= 0:
            return 0.0
        pm = np.sqrt(d * (d + 2))
        return T.lam_plus(rate, a, np.array([pp]))[0] * survival_anti(rate, a, pm)
    # a at which delta = delta0: p = sqrt((EB+delta0)^2-1)
    p0 = np.sqrt((T.EB + delta0) ** 2 - 1)
    a0 = q / p0
    ys = np.log(np.geomspace(a0, a_esc, 200))
    # split integration at points: dense near threshold
    tot = 0.0; parts = []
    edges = np.r_[np.log(a0), np.log(a_esc) - np.geomspace(np.log(a_esc / a0), 1e-9, 40)[1:], np.log(a_esc)]
    edges = np.unique(edges)
    for x0, x1 in zip(edges[:-1], edges[1:]):
        v = quad(integrand, x0, x1, epsrel=1e-6, limit=100)[0]
        tot += v
    return tot


def main():
    table = sys.argv[1] if len(sys.argv) > 1 else 'codex'
    rate = T.codex_table() if table == 'codex' else T.my_table()
    out = {}
    for k in ['z1e6', 'T1keV', 'T20MeV']:
        af = T.EPOCHS[k]
        # representative escape: median from MC if available, else equilibrium
        f = sorted(T.SCR.glob(f'mc_aesc_{k}_{table}_*.npy'))
        if f:
            ae = np.load(f[0])
            sample = np.quantile(ae, [0.05, 0.25, 0.5, 0.75, 0.95])
        else:
            sample = np.array([T.equilibrium_drift(rate, af)['a_escape']])
        vals = [tail_for_escape(rate, a) for a in sample]
        out[k] = dict(a_escape_sample=sample.tolist(), anti_today_per_escaping_particle=vals,
                      mean_estimate=float(np.mean(vals)))
        print(k, out[k], flush=True)
    (T.RES / f'anti_tail_{table}.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
