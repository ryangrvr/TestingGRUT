#!/usr/bin/env python3
"""inf_reheat (Claude, own code): numerical check of non-perturbative fermion production at inflaton zero crossings.

Fermion with mass m_N(t) = y phi(t), phi = Phi cos(m t) (Minkowski, back-reaction off): for each helicity and momentum k,
    i d/dt (u+, u-) = [[m_N, k], [k, -m_N]] (u+, u-),
started in the instantaneous positive-frequency state at the turning point t = 0; n_k = |<v_-(t)|u(t)>|^2 with v_- the
instantaneous negative-frequency eigenvector.  Compared with Greene-Kofman (hep-ph/9807339, hep-ph/0003018):
n_k = exp(-pi k^2 / (y |phidot_0|)) after one crossing, |phidot_0| = Phi m.  The second crossing shows the Pauli-bounded
(interference) evolution n_k in [0, 1].  The energy at the next turning point is compared with n (y Phi), the estimate used in
channels.py.  Units: m (inflaton mass) and y Phi in H_i.  Labels: [C].
Usage: fermion_modes.py -> ../results/fermion_modes.json
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rh_common as rc


def evolve(k, yPhi, m, t_end):
    def rhs(t, z):
        u = z[:2] + 1j * z[2:]
        mn = yPhi * math.cos(m * t)
        du = -1j * np.array([mn * u[0] + k * u[1], k * u[0] - mn * u[1]])
        return np.concatenate([du.real, du.imag])
    mn0 = yPhi
    w = math.hypot(k, mn0)
    vp = np.array([w + mn0, k]) / math.hypot(w + mn0, k)          # +omega eigenvector (mn0 > 0)
    z0 = np.concatenate([vp, np.zeros(2)])
    sol = solve_ivp(rhs, (0, t_end), z0, method='DOP853', rtol=1e-10, atol=1e-12, max_step=0.05 / max(yPhi, 1.0) * 20,
                    dense_output=True)
    return sol


def occupation(sol, t, k, yPhi, m):
    z = sol.sol(t); u = z[:2] + 1j * z[2:]
    mn = yPhi * math.cos(m * t); w = math.hypot(k, mn)
    vm = np.array([-k, w + mn]) if w + mn >= w - mn else np.array([mn - w, k])   # both span the -omega eigenvector
    vm = vm / np.linalg.norm(vm)
    return float(abs(np.vdot(vm, u)) ** 2), w


def case(yPhi, m, nk=40):
    kgk = math.sqrt(yPhi * m / math.pi)
    ks = np.linspace(0.0, 3.0 * kgk, nk)
    t1, t2 = math.pi / m, 2 * math.pi / m          # next turning points (after 1 and 2 crossings)
    rows = []
    for k in ks:
        sol = evolve(k, yPhi, m, t2)
        n1, w1 = occupation(sol, t1, k, yPhi, m)
        n2, w2 = occupation(sol, t2, k, yPhi, m)
        rows.append(dict(k=float(k), n_after_1=n1, n_GK=math.exp(-math.pi * k * k / (yPhi * m)), n_after_2=n2, omega_turn=w1))
    k = np.array([r['k'] for r in rows]); n1 = np.array([r['n_after_1'] for r in rows]); ng = np.array([r['n_GK'] for r in rows])
    w = np.array([r['omega_turn'] for r in rows]); n2 = np.array([r['n_after_2'] for r in rows])
    # Majorana: 2 helicities; d^3k/(2 pi)^3 = k^2 dk / (2 pi^2)
    N1 = 2 * np.trapezoid(k * k * n1, k) / (2 * math.pi ** 2); N2 = 2 * np.trapezoid(k * k * n2, k) / (2 * math.pi ** 2)
    E1 = 2 * np.trapezoid(k * k * n1 * w, k) / (2 * math.pi ** 2)
    Nan = (yPhi * m) ** 1.5 / (4 * math.pi ** 3)
    return dict(yPhi=yPhi, m=m, k_GK=kgk, max_abs_dev_n1_vs_GK=float(np.max(np.abs(n1 - ng))), N_after_1=N1, N_analytic=Nan,
                N_ratio=N1 / Nan, N_after_2=N2, E_after_1=E1, E_estimate_n_yPhi=Nan * yPhi, E_ratio=E1 / (Nan * yPhi),
                max_n_after_2=float(n2.max()), rows=rows)


def main():
    res = dict(method=__doc__, cases=[])
    for yPhi, m in ((60.0, 3.0), (300.0, 3.0), (1500.0, 3.0)):
        c = case(yPhi, m)
        print('yPhi=%g m=%g: max|n1 - GK| = %.3e  N1/N_GK = %.4f  E1/(N yPhi) = %.4f  N2/N1 = %.3f  max n2 = %.3f' % (
            yPhi, m, c['max_abs_dev_n1_vs_GK'], c['N_ratio'], c['E_ratio'], c['N_after_2'] / c['N_after_1'], c['max_n_after_2']))
        res['cases'].append(c)
    json.dump(res, open(rc.RES / 'fermion_modes.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
