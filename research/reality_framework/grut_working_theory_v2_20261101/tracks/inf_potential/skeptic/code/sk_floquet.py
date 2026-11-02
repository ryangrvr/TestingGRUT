#!/usr/bin/env python3
"""inf_potential SKEPTIC (Claude, own code): Minkowski Floquet exponents of inflaton fluctuations on the periodic orbit of
the benchmark potentials at the track's first turning-point amplitude (and fractions of it), via the monodromy matrix
over one exact period (track: modes.py floquet).  Units m_I = 1.  Output ../results/sk_floquet.json.  Label [C]."""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sk_cmb as cm
OUT = Path(__file__).resolve().parent.parent / 'results'

def floquet(name, frac_list=(1.0, 0.5, 0.2, 0.1)):
    d = json.load(open(cm.TRK / 'results' / ('bench_%s.json' % name)))['by_Hi'][0]
    P = cm.Pot(name, d['params']); pmin = P.phimin()
    m2 = P.d2V(pmin)
    V = lambda f: (P.V(f + pmin)) / m2; dV = lambda f: P.dV(f + pmin) / m2; d2V = lambda f: P.d2V(f + pmin) / m2
    A1 = d['osc_gamma0']['turning'][0]['phiI_turn']           # first turning point of the real background (M_P)
    res = []
    for fr in frac_list:
        A = A1 * fr; E = V(A)
        # other turning point at the same energy
        sgn = -np.sign(A)
        B = brentq(lambda f: V(f) - E, sgn * 1e-9, sgn * 50 * abs(A) if name != 'HT6t' else sgn * min(50 * abs(A), 0.99 * (pmin - (P.vv * (-(P.b / P.p) ** 0.2)))))
        ev = lambda t, y: y[1]; ev.direction = 0
        def rhs(t, y):
            return [y[1], -dV(y[0])]
        # period: integrate from A at rest until the velocity returns to zero twice
        s = solve_ivp(rhs, (0, 1e3), [A, 0.0], method='DOP853', rtol=1e-12, atol=1e-15, events=ev, dense_output=True)
        tz = [t for t in s.t_events[0] if t > 1e-6]
        T = tz[1]
        mus = []
        for k in np.linspace(0.0, 2.0, 201):
            def r2(t, y):
                f = s.sol(t)[0]; w2 = k * k + d2V(f)
                return [y[1], -w2 * y[0], y[3], -w2 * y[2]]
            q = solve_ivp(r2, (0, T), [1, 0, 0, 1], method='DOP853', rtol=1e-11, atol=1e-13)
            Mo = np.array([[q.y[0, -1], q.y[2, -1]], [q.y[1, -1], q.y[3, -1]]])
            lam = np.linalg.eigvals(Mo)
            mus.append(max(np.log(np.abs(lam))) / T)
        mus = np.array(mus); i = int(np.argmax(mus))
        res.append(dict(frac=fr, A_MP=A, other_turn_MP=B, period_over_2pi_inv_m=T / (2 * math.pi), omega_over_m=2 * math.pi / T,
                        mu_max_over_m=float(mus[i]), k_at_max_over_m=float(np.linspace(0, 2, 201)[i]),
                        mu_max_over_Hi=float(mus[i]) * d['m_I_over_Hi']))
        print('%s A=%.4f (frac %.1f) other=%.4f omega/m=%.3f mu_max/m=%.4f at k/m=%.2f -> mu/H_i=%.2f' % (
            name, A, fr, B, 2 * math.pi / T, mus[i], np.linspace(0, 2, 201)[i], mus[i] * d['m_I_over_Hi']), flush=True)
    return res

if __name__ == '__main__':
    allr = {}
    for n in sys.argv[1:] or ['T1', 'HT6t', 'E1']:
        allr[n] = floquet(n)
    json.dump(dict(method=__doc__, runs=allr), open(OUT / 'sk_floquet.json', 'w'), indent=1)
