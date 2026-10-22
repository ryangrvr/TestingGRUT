"""v_dissipation, script 2: independent reproduction of KV 1601.00601 (KV16a) numbers quoted by the track (NOTES sec 3).

Two methods, neither taken from ../../dissipation/code:
  (i)  ODE in tau for the state (h, r_V) only (f is reconstructed afterwards by inverting r_V(f)), radiation only:
           hdot = 2 r_V - 2 h^2 ,   r_V-dot = -gamma |h| (hdot + 2 h^2)^2      (E_* units; v1 check B)
       integrated with LSODA and Radau at two tolerances.
  (ii) the exact App. B solution in xi = ln a: r_V = 1/(4 gamma (xi - xi0)), h^2 = (C1 e^{-4xi} + e^{-y} Ei(y))/gamma,
       y = 4(xi - xi0), with xi0, C1 fixed by the data at tau = 1 (a(1) = 1), and tau(xi) = 1 + int_0^xi dxi'/h(xi')
       done with mpmath; the root tau(xi) = 1e5 is found by bracketing.  No ODE solver is involved.
Paper values [Q] (1601.00601v6 Figs 1-3 captions, sec 7): {r_V(1), r_M(1)} = {8.97333, 27.0267} (lambda = 1),
{6.97333, 29.0267} (lambda = -1), r_V(1e5) = 3e-3 for both; gamma = 1/100; {h, f}(1) = {6, 3/5}; corners {6 +- 1, 3/5 +- 1/5}.
Track values [track NOTES sec 3]: r_V(1e5) = 3.0285e-3 (lambda = 1), 3.0284e-3 (lambda = -1, 0), corners 3.0284e-3..3.0285e-3,
App. B agreement 3.0e-12, (7.1) ratio 1 - 2.8e-5 at tau = 1e5, f(1e5) = 1.51685 (fbar 1.51749), 1.14091 (1.14139).
Output: ../results/v2_kv16a.json.  Run: python3 -I v2_kv16a.py
"""
import json, os
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

mp.mp.dps = 30


def rV_of_f(f, lam):
    return lam - f**2 + 3 / f**2


def f_of_rV(r, lam):           # r_V decreasing in x = f^2: x^2 + (r - lam) x - 3 = 0
    b = r - lam
    return np.sqrt((-b + np.sqrt(b * b + 12)) / 2)


def rhs(tau, s, gam):
    h, r = s
    hd = 2 * r - 2 * h * h
    k = hd + 2 * h * h
    return [hd, -gam * abs(h) * k * k]


def ode_run(h1, r1, gam, taus, method, rtol):
    sol = solve_ivp(rhs, (1.0, taus[-1]), [h1, r1], args=(gam,), method=method, rtol=rtol, atol=1e-16,
                    t_eval=taus, dense_output=False)
    return sol.y


def appB(h1, r1, gam, tau_target):
    """exact solution; returns r_V, h, xi at tau_target."""
    gam = mp.mpf(gam); r1 = mp.mpf(r1); h1 = mp.mpf(h1)
    xi0 = -1 / (4 * gam * r1)                        # xi(1) = 0
    C1 = gam * h1**2 - mp.e**(4 * xi0) * mp.ei(-4 * xi0)

    def h_of(xi):
        y = 4 * (xi - xi0)
        return mp.sqrt((C1 * mp.e**(-4 * xi) + mp.e**(-y) * mp.ei(y)) / gam)

    def tau_of(xi):
        return 1 + mp.quad(lambda s: 1 / h_of(s), [0, min(xi, 5), min(xi, 50), min(xi, 500), xi])

    # bracket: on the attractor xi ~ (3 tau/(4 sqrt gamma))^(2/3)
    xg = (3 * tau_target / (4 * mp.sqrt(gam)))**(mp.mpf(2) / 3)
    lo, hi = xg * 0.9, xg * 1.1
    xs = mp.findroot(lambda s: tau_of(s) - tau_target, (lo, hi), solver='anderson')
    return float(1 / (4 * gam * (xs - xi0))), float(h_of(xs)), float(xs)


out = {'paper_[Q]': {'Fig1': {'lambda': 1, 'rV1': 8.97333, 'rM1': 27.0267, 'rV(1e5)': 3e-3},
                     'Fig2': {'lambda': -1, 'rV1': 6.97333, 'rM1': 29.0267, 'rV(1e5)': 3e-3}},
       'track': {'Fig1 rV(1e5)': 3.0285e-3, 'Fig2 rV(1e5)': 3.0284e-3, 'lambda0 rV(1e5)': 3.0284e-3,
                 'f(1e5) lambda=1': 1.51685, 'f(1e5) lambda=-1': 1.14091}}
gam = 0.01
taus = np.array([10.0, 100.0, 1e3, 1e4, 1e5])
res = {}
for lam in (1.0, -1.0, 0.0):
    h1, f1 = 6.0, 0.6
    r1 = rV_of_f(f1, lam)
    rM1 = h1**2 - r1
    entry = {'rV1': r1, 'rM1': rM1}
    for method, rtol in (('LSODA', 1e-11), ('Radau', 1e-12), ('DOP853', 1e-13)):
        Y = ode_run(h1, r1, gam, taus, method, rtol)
        entry['ode_%s_rV' % method] = list(map(float, Y[1]))
    rB, hB, xB = appB(h1, r1, gam, 1e5)
    entry['appB_rV(1e5)'] = rB
    entry['appB_h(1e5)'] = hB
    entry['appB_xi(1e5)'] = xB
    entry['ode_vs_appB_rel'] = abs(entry['ode_DOP853_rV'][-1] - rB) / rB
    entry['asym_(6 gamma tau)^-2/3'] = (6 * gam * 1e5)**(-2 / 3)
    entry['rV/asym - 1'] = rB / (6 * gam * 1e5)**(-2 / 3) - 1
    entry['f(1e5)'] = float(f_of_rV(rB, lam))
    entry['fbar'] = float(f_of_rV(0.0, lam))
    res['lambda=%g' % lam] = entry
out['KV16a_figs'] = res

# Fig 3 corners
corners = []
for lam in (1.0, -1.0):
    for h1 in (5.0, 7.0):
        for f1 in (0.4, 0.8):
            r1 = rV_of_f(f1, lam)
            rB, hB, xB = appB(h1, r1, gam, 1e5)
            corners.append({'lambda': lam, 'h1': h1, 'f1': f1, 'rV1': r1, 'rV(1e5)_appB': rB, 'f(1e5)': float(f_of_rV(rB, lam))})
out['Fig3_corners'] = corners
out['Fig3_range_rV(1e5)'] = [min(c['rV(1e5)_appB'] for c in corners), max(c['rV(1e5)_appB'] for c in corners)]

# (7.1) ratio for several gamma, from the exact solution
ratios = {}
for g_ in (1.0, 1e-3, 1 / (32 * np.pi**2), 0.01):
    rB, hB, xB = appB(6.0, rV_of_f(0.6, 1.0), g_, 1e5)
    ratios['%.6g' % g_] = {'rV/(6 g tau)^-2/3 - 1': rB / (6 * g_ * 1e5)**(-2 / 3) - 1,
                           'h/(6 g tau)^-1/3 - 1': hB / (6 * g_ * 1e5)**(-1 / 3) - 1}
out['eq7.1_ratios_at_tau=1e5'] = ratios

with open(os.path.join(os.path.dirname(__file__), '..', 'results', 'v2_kv16a.json'), 'w') as fh:
    json.dump(out, fh, indent=1)
for k, v in res.items():
    print(k, 'rV1=%.6g rM1=%.6g  appB rV(1e5)=%.6e  ode=%.6e  rel=%.1e  f=%.6f fbar=%.6f  asym-1=%.2e' % (
        v['rV1'], v['rM1'], v['appB_rV(1e5)'], v['ode_DOP853_rV'][-1], v['ode_vs_appB_rel'], v['f(1e5)'], v['fbar'], v['rV/asym - 1']))
    print('   LSODA', v['ode_LSODA_rV'][-1], 'Radau', v['ode_Radau_rV'][-1])
print('corners range', out['Fig3_range_rV(1e5)'])
print(json.dumps(ratios, indent=1))
