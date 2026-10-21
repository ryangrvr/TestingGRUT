"""Analytic / numerical checks of the independent code (v_charged_exit).
Writes results/checks.json."""
import json, os, time
import numpy as np
from scipy.integrate import solve_ivp
import vlat
from vlat import Lattice, A_G, B_X, NBAR_C

HERE = os.path.dirname(os.path.abspath(__file__))
res = {}

# --- thresholds (closed form) ---
Yc2 = 1.0 / (2 * A_G)                         # |Y_c|^2 where phi mass^2 at phi=0 vanishes: -1 + 2A|Y|^2 = 0
res['nbar_c'] = NBAR_C                          # sqrt(8B)|Y_c|^3 in sqrt(lambda_phi) v^3
res['n_c_over_v3'] = NBAR_C * 0.1              # sqrt(lambda_phi) = 0.1
F2s1, S0s1, oms1 = vlat.branch1_state(0.9 * NBAR_C)
res['branch1_at_0.9nc'] = dict(F2=F2s1, S0=S0s1, omega=oms1)
# branch-1 spinodal n_s1 = max_F2 omega(F2) F2 on 0<F2<1
F2g = np.linspace(1e-6, 1 - 1e-9, 2000001)
ng = np.sqrt(A_G * (1 - A_G * F2g) + B_X * F2g) * F2g
res['nbar_s1'] = float(ng.max()); res['F2_at_s1'] = float(F2g[ng.argmax()])

# --- 1. exact homogeneous phi=0 orbit in radiation era, kappa = 2 and 4 ---
for kappa in (2.0, 4.0):
    lat = Lattice(4, 1.0, 40.0)
    info = vlat.init_charged_exit(lat, kappa, 0, phi_std=0.0, x_std=0.0)
    lat.S[...] = 0.0
    lat.q1[...] = np.sqrt(2) * info['Y0']; lat.q2[...] = 0.0
    lat.p1[...] = -1 / 40.0 * lat.q1; lat.p2[...] = -info['Omega'] * lat.q1
    Q0 = lat.charge()
    devs = []; ndev = []
    for a_t in np.linspace(1.1, 3.0, 20):
        lat.evolve(40.0 * a_t, 0.01)
        Y2 = float(lat.F2().mean()) / 2
        n = float(lat.nphys().mean())
        pred_n = kappa * NBAR_C / a_t ** 3
        devs.append(abs(Y2 / (n ** 2 / (8 * B_X)) ** (1 / 3) - 1))
        ndev.append(abs(n / pred_n - 1))
    res['orbit_kappa%g' % kappa] = dict(max_rel_dev_Y2_vs_n_relation=max(devs), max_rel_dev_n_vs_ni_a3=max(ndev),
                                        max_abs_S=float(np.abs(lat.S).max()), Q_rel=abs(lat.charge() / Q0 - 1),
                                        a_c_pred=kappa ** (1 / 3))

# physical-time check that |X| ~ 1/a is exact: |X|'' + 3H|X|' - thetadot^2 |X| + 2 B |X|^3 = 0 with
# H = 1/(2t): (1/a)'' + 3H (1/a)' = (3H^2 - 3H^2)/a = 0 and thetadot^2 = 2B|Y|^2  (closed form, see VERIFY.md)

# --- 2. linear tachyonic growth of the k=0 phi mode through n_c: lattice vs ODE ---
kappa, ti = 2.0, 40.0
lat = Lattice(4, 1.0, ti)
info = vlat.init_charged_exit(lat, kappa, 0, phi_std=0.0, x_std=0.0)
Y0, Om = info['Y0'], info['Omega']
lat.q1[...] = np.sqrt(2) * Y0; lat.q2[...] = 0.0
lat.p1[...] = -1 / ti * lat.q1; lat.p2[...] = -Om * lat.q1
eps = 1e-10
lat.S[...] = eps; lat.pS[...] = 0.0
# ODE for delta S in conformal time: (a^2 dS')' = -a^4 (-1 + F2) dS, F2 = 2 Y0^2/a^2
def rhs(t, y):
    a = t / ti
    F2 = 2 * Y0 ** 2 / a ** 2
    return [y[1] / a ** 2, -a ** 4 * (-1 + A_G * F2) * y[0]]
taus = np.linspace(ti, 2.0 * ti, 41)
sol = solve_ivp(rhs, (ti, 2.0 * ti), [eps, 0.0], t_eval=taus, rtol=1e-11, atol=1e-30, method='DOP853')
lg = []
for t in taus[1:]:
    lat.evolve(t, 0.01)
    lg.append(float(lat.S.mean()))
lg = np.array(lg); od = sol.y[0][1:]
res['tachyon_growth'] = dict(kappa=kappa, tau_i=ti, a_c=kappa ** (1 / 3),
                            growth_lattice=float(lg[-1] / eps), growth_ode=float(od[-1] / eps),
                            max_rel_dev=float(np.max(np.abs(lg / od - 1))),
                            a_min_of_S=float(taus[1:][np.argmin(lg)] / ti))
# sign of phi mass^2 at the closed-form crossing
for a_t in (kappa ** (1 / 3) * 0.99, kappa ** (1 / 3) * 1.01):
    res.setdefault('mass2_phi_at_phi0', {})['a=%.4f' % a_t] = -1 + 2 * Y0 ** 2 / a_t ** 2

# --- 3. charge conservation + dt convergence on a strongly nonlinear small run ---
def small(dt, N=96, dx=0.25, kappa=2.0, ti=10.0, a_end=3.0):
    lat = Lattice(N, dx, ti)
    vlat.init_charged_exit(lat, kappa, 7, phi_std=5e-2, x_std=5e-2)
    Q0 = lat.charge(); qd = 0.0
    for a_t in np.linspace(1.2, a_end, 10):
        lat.evolve(ti * a_t, dt)
        qd = max(qd, abs(lat.charge() / Q0 - 1))
    return lat, qd
t0 = time.time()
l1, qd1 = small(0.05)
l2, qd2 = small(0.025)
l3, qd3 = small(0.0125)
d12 = float(np.abs(l1.S - l2.S).max()); d23 = float(np.abs(l2.S - l3.S).max())
res['nonlinear_small'] = dict(max_charge_rel_dev=[qd1, qd2, qd3], maxdiff_S_dt05_vs_025=d12,
                              maxdiff_S_dt025_vs_0125=d23, ratio=d12 / d23, S_range=[float(l3.S.min()), float(l3.S.max())],
                              F2_max=float(l3.F2().max()), wall_s=time.time() - t0)

# --- 4. independent friction-form RK (solve_ivp) on a 24^2 grid vs the canonical KDK ---
N, dx, ti = 24, 0.5, 10.0
lat = Lattice(N, dx, ti)
vlat.init_charged_exit(lat, 2.0, 3, phi_std=0.3, x_std=0.1)
f0 = [lat.S.copy(), lat.q1.copy(), lat.q2.copy()]
v0 = [lat.pS.copy() / lat.a ** 2, lat.p1.copy() / lat.a ** 2, lat.p2.copy() / lat.a ** 2]
def lapn(f):
    return (np.roll(f, 1, 0) + np.roll(f, -1, 0) + np.roll(f, 1, 1) + np.roll(f, -1, 1) - 4 * f) / dx ** 2
def rhs2(t, y):
    a = t / ti; ap = 1 / ti
    S, q1, q2, vS, v1, v2 = y.reshape(6, N, N)
    F2 = q1 ** 2 + q2 ** 2
    aS = lapn(S) - a ** 2 * (S * (S * S - 1) + A_G * S * F2) - 2 * ap / a * vS
    m2 = A_G * S * S + B_X * F2
    a1 = lapn(q1) - a ** 2 * m2 * q1 - 2 * ap / a * v1
    a2 = lapn(q2) - a ** 2 * m2 * q2 - 2 * ap / a * v2
    return np.concatenate([vS.ravel(), v1.ravel(), v2.ravel(), aS.ravel(), a1.ravel(), a2.ravel()])
y0 = np.concatenate([x.ravel() for x in f0 + v0])
sol = solve_ivp(rhs2, (ti, 2.5 * ti), y0, rtol=1e-10, atol=1e-12, method='DOP853')
yS = sol.y[:, -1].reshape(6, N, N)
lat.evolve(2.5 * ti, 0.005)
res['friction_form_vs_canonical'] = dict(max_abs_diff_S=float(np.abs(yS[0] - lat.S).max()),
                                         max_abs_diff_q1=float(np.abs(yS[1] - lat.q1).max()),
                                         S_range=[float(lat.S.min()), float(lat.S.max())],
                                         Q_friction=float(np.sum(yS[2] * yS[3 + 1] - yS[1] * yS[3 + 2]) * (2.5) ** 2 * dx * dx),
                                         Q_canonical=lat.charge())

# --- 5. diagnostics on synthetic configurations ---
N, dx = 256, 0.25
x = (np.arange(N) + 0.5) * dx; X, Yg = np.meshgrid(x, x, indexing='ij')
lat = Lattice(N, dx, 10.0, with_X=False)
Lb = N * dx
lat.S[...] = np.tanh((X - Lb / 4) / np.sqrt(2)) * -np.tanh((X - 3 * Lb / 4) / np.sqrt(2))
d = vlat.diagnostics(lat)
res['synthetic_two_walls'] = dict(L_basin=d['L_basin'] * Lb * Lb, exact=2 * Lb,
                                  note='Crofton factor pi/4 makes axis-aligned lines read 4/pi too low; isotropic networks are unbiased')
r = np.sqrt((X - Lb / 2) ** 2 + (Yg - Lb / 2) ** 2)
lat.S[...] = np.tanh((r - 16) / np.sqrt(2))
d = vlat.diagnostics(lat)
res['synthetic_circle_R16'] = dict(L_basin=d['L_basin'] * Lb * Lb, exact=2 * np.pi * 16)

with open(os.path.join(HERE, 'results', 'checks.json'), 'w') as fh:
    json.dump(res, fh, indent=1)
print(json.dumps(res, indent=1))
