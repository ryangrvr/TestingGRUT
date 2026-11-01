"""fx_formation item 1a (Claude, 1 Nov 2026): the M1 object family E(Q), Omega(Q), radii, Q_min .. 1e8.

Solver: tracks/m1_vacuum/code/qball_m1.py, copied unchanged to code/qball_m1_copy.py (Numerov + Newton, fixed-Q mode
with Omega as the extra unknown; fixed-Omega mode for the branch near Q_min).  Units as there: rho = m_X r,
Qbar = g^2 Q, E/(Q m_X) = E_b/Qbar, n/m_X^3 = 100 Omega F^2.
Outputs results/family.json: rows (Q, Omega, E/Qm, S0, F0, r50, r90, r99), the stable branch near Q_min, the fit
E/Q = Omega_b + A Q^-1/3 + B Q^-2/3 + C Q^-1 on Q >= 3e5, and first-law checks dE/dQ vs Omega from the table.
"""
import json, math, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import qball_m1_copy as qb

OM_BULK = 0.7758647074595397


def main(out):
    t0 = time.time()
    rows = []
    # primary seed at Q = 2e5
    sv, u, w, Om, ob = qb.solve_fixed_Q(200000, 60.0, 0.025)
    seed0 = (sv, u, w, Om)
    rows.append(ob)
    print('Q=2e5 Om=%.10f E/Qm=%.10f' % (Om, ob['E_over_Qm']), flush=True)
    # downward
    cur = seed0
    for Q in (150000, 120000, 100000, 90000, 80000, 75000, 71262.76, 70000, 65000, 62000, 60000):
        svq, uq, wq, Omq, ob = qb.solve_fixed_Q(Q, 60.0, 0.025, seed=cur)
        cur = (svq, uq, wq, Omq); rows.append(ob)
        print('Q=%.6g Om=%.8f E/Qm=%.8f r90=%.3f t=%.0fs' % (Q, Omq, ob['E_over_Qm'], ob['r90'], time.time() - t0), flush=True)
    # upward
    cur = seed0
    for Q, L in ((3e5, 60.), (5e5, 60.), (1e6, 60.), (2e6, 70.), (5e6, 80.), (1e7, 90.), (2e7, 100.), (3e7, 110.),
                 (5e7, 120.), (7e7, 130.), (1e8, 140.)):
        svq, uq, wq, Omq, ob = qb.solve_fixed_Q(Q, L, 0.025, seed=cur)
        cur = (svq, uq, wq, Omq); rows.append(ob)
        print('Q=%.6g Om=%.8f E/Qm=%.8f r90=%.3f t=%.0fs' % (Q, Omq, ob['E_over_Qm'], ob['r90'], time.time() - t0), flush=True)
    rows.sort(key=lambda d: d['Q'])
    # stable branch near Q_min by fixed Omega (L = 100, h = 0.05), Omega from 0.955 to 0.9728
    svb = qb.Solver(100.0, 0.05)
    s60, u60, w60, Om60, _ = qb.solve_fixed_Q(60000, 60.0, 0.025, seed=seed0)
    ub, wb = svb.interp_from(s60, u60, w60)
    ub, wb, _, _, _ = svb.solve(ub, wb, Om60)
    branch = []
    for Omv in np.arange(round(Om60, 3), 0.99, 0.001):
        ub, wb, _, res, it = svb.solve(ub, wb, float(Omv))
        o = svb.observables(ub, wb, float(Omv))
        branch.append(dict(Omega=float(Omv), Q=o['Q'], E_over_Qm=o['E_over_Qm'], r90=o['r90']))
    Qs = np.array([b['Q'] for b in branch]); Os = np.array([b['Omega'] for b in branch])
    i = int(np.argmin(Qs)); p = np.polyfit(Os[i - 1:i + 2], Qs[i - 1:i + 2], 2)
    Om_min = -p[1] / (2 * p[0]); Q_min = float(np.polyval(p, Om_min))
    pe = np.polyfit(Os[i - 1:i + 2], [branch[j]['E_over_Qm'] for j in range(i - 1, i + 2)], 2)
    EQ_min = float(np.polyval(pe, Om_min))
    # fit on Q >= 3e5
    Q = np.array([r['Q'] for r in rows]); EQ = np.array([r['E_over_Qm'] for r in rows]); OmA = np.array([r['Omega'] for r in rows])
    m = Q >= 3e5
    X = np.vstack([Q[m] ** (-1 / 3), Q[m] ** (-2 / 3), Q[m] ** (-1.0)]).T
    coef, *_ = np.linalg.lstsq(X, EQ[m] - OM_BULK, rcond=None)
    fit_res = EQ[m] - OM_BULK - X @ coef
    X2 = X[:, :2]
    coef2, *_ = np.linalg.lstsq(X2, EQ[m] - OM_BULK, rcond=None)
    fit2_res = EQ[m] - OM_BULK - X2 @ coef2
    # first-law check from the table: d(E)/dQ by spline vs Omega
    from scipy.interpolate import CubicSpline
    lq = np.log(Q); Etot = EQ * Q
    cs = CubicSpline(lq, Etot)
    dEdQ = cs(lq, 1) / Q
    # Omega implied by the 3-term fit: d(E)/dQ = Om_b + (2/3)A Q^-1/3 + (1/3) B Q^-2/3
    om_fit = OM_BULK + (2 / 3) * coef[0] * Q ** (-1 / 3) + (1 / 3) * coef[1] * Q ** (-2 / 3)
    res = dict(rows=[{k: r[k] for k in ('Q', 'Omega', 'E_over_Qm', 'S0', 'F0', 'r50', 'r90', 'r99', 'L', 'h', 'virial_residual')} for r in rows],
               branch_near_Qmin=branch, Q_min=Q_min, Omega_at_Qmin=float(Om_min), E_over_Qm_at_Qmin=EQ_min,
               fit3=dict(Omega_b=OM_BULK, A=coef[0], B=coef[1], C=coef[2], max_abs_residual=float(np.max(np.abs(fit_res))), range='Q>=3e5'),
               fit2=dict(Omega_b=OM_BULK, A=coef2[0], B=coef2[1], max_abs_residual=float(np.max(np.abs(fit2_res)))),
               first_law=[dict(Q=float(q), dEdQ_spline=float(d), Omega=float(o), rel=float(d / o - 1), Omega_fit=float(of))
                          for q, d, o, of in zip(Q, dEdQ, OmA, om_fit)],
               runtime_s=time.time() - t0)
    print('Q_min=%.1f Om=%.5f E/Qm=%.6f  fit3 A=%.4f B=%.3f C=%.1f maxres=%.1e; fit2 A=%.4f B=%.3f maxres=%.1e' % (
        Q_min, Om_min, EQ_min, coef[0], coef[1], coef[2], res['fit3']['max_abs_residual'], coef2[0], coef2[1], res['fit2']['max_abs_residual']))
    for r in res['first_law']:
        print('  Q=%.4g dE/dQ=%.6f Omega=%.6f rel=%.1e Omega_fit=%.6f' % (r['Q'], r['dEdQ_spline'], r['Omega'], r['rel'], r['Omega_fit']))
    json.dump(res, open(out, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main(sys.argv[1])
