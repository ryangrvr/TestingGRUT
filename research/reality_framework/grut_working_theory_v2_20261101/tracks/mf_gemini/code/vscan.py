"""mf_gemini: m_chi/H_end and the post-inflation amplitude versus v for a chi-dominated waterfall [C].

V = lam/4 (chi^2 - v^2)^2 with lam v^4/4 = 3 (code units M_P = H0 = 1).  m_chi/H0 = sqrt(2 lam) v = sqrt(24)/v.
The proposal's m_chi = 3 H0 is v^2 = 8/3.  M1r's H_i is the Hubble rate at the end-of-inflation slice, so the
quiet condition reads m_chi/H_end; this scan gives both, the e-folds from chi = H0/2pi (hub = 5.34e-8) to the
end at the tachyonic rate, and the first overshoot amplitude about v.
Run: python3 -B vscan.py -> results/vscan.json
"""
import json, os
import numpy as np
from scipy.integrate import solve_ivp
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')


def one(v2, hub=1.3e11 / 2.435e18):
    lam = 12.0 / v2**2
    v = np.sqrt(v2)
    V = lambda c: 0.25 * lam * (c**2 - v2)**2
    dV = lambda c: lam * (c**2 - v2) * c
    def rhs(t, y):
        h = np.sqrt((0.5 * y[2]**2 + V(y[1])) / 3)
        return [h, y[2], -3 * h * y[2] - dV(y[1])]
    s = (-3 + np.sqrt(9 + 4 * lam * v2)) / 2
    c0 = hub / (2 * np.pi)
    def ev(t, y):
        h = np.sqrt((0.5 * y[2]**2 + V(y[1])) / 3)
        return 0.5 * y[2]**2 / h**2 - 1
    ev.terminal = True; ev.direction = 1
    s1 = solve_ivp(rhs, (0, 1e4), [0, c0, s * c0], method='DOP853', rtol=1e-11, atol=1e-20, events=ev)
    ye = s1.y_events[0][0]
    he = np.sqrt((0.5 * ye[2]**2 + V(ye[1])) / 3)
    s2 = solve_ivp(rhs, (s1.t_events[0][0], s1.t_events[0][0] + 30), ye, method='DOP853', rtol=1e-11, atol=1e-16,
                   max_step=0.01, dense_output=True)
    tt = np.linspace(s2.t[0], s2.t[-1], 30001)
    c = s2.sol(tt)[1]
    dc = s2.sol(tt)[2]
    i = np.where(np.sign(dc[1:]) != np.sign(dc[:-1]))[0]
    amp = abs(c[i[0]] - v) / v if i.size else None
    return dict(v_over_MP=float(v), m_over_H0=float(np.sqrt(24) / v), m_over_Hend=float(np.sqrt(24) / v / he),
                Hend_over_H0=float(he), eta_top=float(-4 / v2), growth_rate=float(s),
                N_quantum_to_end=float(ye[0]), first_overshoot_over_v=float(amp) if amp else None,
                first_overshoot_over_MP=float(amp * v) if amp else None, chi_end_over_v=float(ye[1] / v))


out = [one(v2) for v2 in [8.0 / 3.0, 3.5, 4.5, 6.0, 8.0, 12.0, 24.0]]
with open(os.path.join(RES, 'vscan.json'), 'w') as f:
    json.dump(out, f, indent=1)
for r in out:
    print(r)
