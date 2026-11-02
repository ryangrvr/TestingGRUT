"""mf_gemini: the end of the mild waterfall and the chi oscillation that follows [C].

Once chi >> chi_trap, sigma is pinned near 0 (m_sigma = g chi >> H) and carries < 1e-6 of the energy
(waterfall_dn.py end_state.rho_sig_frac).  The rest is single-field hilltop rolling of chi on
lam/4 (chi^2 - v^2)^2 with v^2 = 8/3, lam = 27/16 (code units M_P = H0 = 1).  This script integrates that
roll from chi = 1e-6 through the end of inflation (eps_H = 1) and 40 H0^-1 beyond, and records:
  - H_end/H0, m_chi/H_end, the slice values (a = 1 at eps_H = 1 here);
  - the first swings of chi about v (overshoot chi_max - v, return v - chi_min), the decay of the
    oscillation amplitude, the number of oscillations per Hubble time;
  - w(t) averaged over oscillations, and K(t) = a^2 H/H_end that a sudden decay at time t would give
    (radiation thereafter): K = a^2 H at the decay time, in H_end units;
  - the chi(t) table (thinned) for preheat_lin.py.
Also a check of the classical end state against waterfall_dn's two-stage runs.

Run: python3 -B waterfall_bg.py -> results/waterfall_bg.json
"""
import json
import os
import sys
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
v2 = 8.0 / 3.0
v = np.sqrt(v2)
lam = 27.0 / 16.0


def V(c):
    return 0.25 * lam * (c**2 - v2)**2


def dV(c):
    return lam * (c**2 - v2) * c


def rhs(t, y):
    N, c, dc = y
    h = np.sqrt((0.5 * dc**2 + V(c)) / 3.0)
    return [h, dc, -3 * h * dc - dV(c)]


def main():
    c0 = 1e-6
    s = (-3 + np.sqrt(9 + 4 * 0.5 * lam * v2)) / 2      # tachyonic growth rate at chi = 0 (H = 1)
    y0 = [0.0, c0, s * c0]

    def ev_end(t, y):
        h = np.sqrt((0.5 * y[2]**2 + V(y[1])) / 3.0)
        return 0.5 * y[2]**2 / h**2 - 1.0
    ev_end.terminal = True
    ev_end.direction = 1
    s1 = solve_ivp(rhs, (0, 500), y0, method='DOP853', rtol=1e-12, atol=1e-16, events=ev_end)
    te = s1.t_events[0][0]
    ye = s1.y_events[0][0]
    he = np.sqrt((0.5 * ye[2]**2 + V(ye[1])) / 3.0)
    Ne = ye[0]
    s2 = solve_ivp(rhs, (te, te + 40.0 / he * he), ye, method='DOP853', rtol=1e-12, atol=1e-16,
                   dense_output=True, max_step=0.01)
    t = np.linspace(te, s2.t[-1], 40001)
    N, c, dc = s2.sol(t)
    h = np.sqrt((0.5 * dc**2 + V(c)) / 3.0)
    rho = 3 * h**2
    p = 0.5 * dc**2 - V(c)
    # turning points of chi (dc changes sign)
    sg = np.sign(dc)
    idx = np.where(sg[1:] != sg[:-1])[0]
    turns = [dict(t_after_end_Hend=float((t[i] - te) * he), chi_over_v=float(c[i] / v),
                  dev_over_v=float((c[i] - v) / v)) for i in idx[:12]]
    # oscillation-averaged w between successive maxima
    wavg = []
    for i0, i1 in zip(idx[:-2:2], idx[2::2]):
        sl = slice(i0, i1)
        wavg.append(dict(t_mid_Hend=float((t[(i0 + i1) // 2] - te) * he),
                         w=float(np.trapezoid(p[sl], t[sl]) / np.trapezoid(rho[sl], t[sl]))))
    # K(t) if the field decayed suddenly at t into radiation: K = a^2 H in H_end units (a = 1 at end)
    a = np.exp(N - Ne)
    Kt = a**2 * h / he
    Ksamp = {f'{x:g}': float(np.interp(x, (t - te) * he, Kt)) for x in [0.5, 1, 2, 4, 8, 16]}
    # crossings of chi = v per Hubble time
    cr = np.where(np.sign(c[1:] - v) != np.sign(c[:-1] - v))[0]
    mchi_over_Hend = np.sqrt(2 * lam * v2) / he
    # amplitude envelope |chi - v| at turning points vs a^{-3/2}
    devs = np.array([abs(x['dev_over_v']) for x in turns])
    tt_turn = np.array([x['t_after_end_Hend'] for x in turns])
    out = dict(
        units='code: M_P = H0 = 1, H0 = sqrt(V0/3); "Hend" = H at eps_H = 1',
        N_from_chi_1em6_to_end=float(Ne), H_end_over_H0=float(he), m_chi_over_H_end=float(mchi_over_Hend),
        chi_end_over_v=float(ye[1] / v), dchi_end=float(ye[2]),
        representative=dict(osc_amp_over_v_first=float(abs(turns[0]['dev_over_v'])),
                            first_return_dev_over_v=float(abs(turns[1]['dev_over_v'])) if len(turns) > 1 else None,
                            note='first swing after the end of inflation: overshoot of chi above v, then return'),
        turning_points=turns, w_oscillation_averaged=wavg[:10],
        K_if_sudden_decay_at_t_Hend=Ksamp,
        chi_v_crossings_per_Hend_time=float(len(cr) / ((t[-1] - te) * he)),
        chi_ever_returns_below_zero=bool((c < 0).any()), chi_min_after_end_over_v=float(c.min() / v),
        table=dict(t_after_end_H0=[float(x) for x in (t[::40] - te)], N_after_end=[float(x) for x in (N[::40] - Ne)],
                   chi=[float(x) for x in c[::40]], dchi=[float(x) for x in dc[::40]], h=[float(x) for x in h[::40]]))
    with open(os.path.join(RES, 'waterfall_bg.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: out[k] for k in out if k != 'table'}, indent=1))


if __name__ == '__main__':
    main()
