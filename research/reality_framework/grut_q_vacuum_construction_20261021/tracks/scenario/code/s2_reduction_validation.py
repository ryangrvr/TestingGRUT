"""Part 2: validated reduction. Full QV1 equations (q dynamics included) at compressed hierarchies, against the
adiabatic reduction "GR + V_m(t) - c*(mu) + adiabatic q-correction".

Units: K0 = 1, time unit 1/H_c where H_c ~ H at the transition centre N = 0; then m = m_q/H_c = r.
Matter: radiation + dust + a matter-vacuum step DeltaV(N) = Delta (1 - tanh(N/w))/2 (supplied history [P]); under
P_q, c*(mu) = final matter vacuum, absorbed into eps_inf, so the gravitating matter-vacuum term is DeltaV(N).
q-sector shape: KV 0806.2805 (b3 = 4/3, b4 = 1/3, beta = 1).

Reductions compared at equal N = ln a (matter densities are identical functions of N in all versions):
  R0: 3 H^2 = rho_m(N)                                      (GR + V_m(t) - c*(mu))
  R1: 3 H^2 (1 + (S + S_N)/(3 m^2)) = rho_m + S^2/(12 m^2)  (+ adiabatic q-correction, O(H^2/m^2))
Diagnostics: N-space window averages of H over k periods (k chosen so that the window is ~0.002-0.01 e-folds),
the Einstein-frame oscillation energy about the quasi-static track (times a^3), <rho_V> against S^2/(12 m^2),
the Friedmann constraint C and its Bianchi propagation C a^4 = const.

Run: python3 -I s2_reduction_validation.py [outdir]   (about 15-25 min on 4 cores)
"""
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qv1_core import KV_SHAPE, QV1, ToyMatter, integrate_full, osc_energy, quasi_static_toy, reduced_toy  # noqa: E402

FD_DEFAULT, FDUST_DEFAULT = 0.05, 0.02


def setup(r, w, fD, fdust, Ni, Delta_on=True):
    q = QV1(r, **KV_SHAPE)
    Delta = 3 * fD if Delta_on else 0.0
    mat = ToyMatter(Delta, 0.0, w)
    rd_i = 3 * fdust * math.exp(-3 * Ni)
    rr_i = 3 * (1 - fdust - fD) * math.exp(-4 * Ni)
    z_i = (rr_i, rd_i)
    H = math.sqrt(mat.rho(z_i, Ni) / 3)
    for _ in range(4):
        Yq, Vq, _, _ = quasi_static_toy(q, mat, Ni, rr_i, rd_i, H)
        H = q.H_from_constraint(Yq, Vq, mat.rho(z_i, Ni))
    return q, mat, z_i, Yq, Vq, H


def run_toy(r, w, fD=FD_DEFAULT, fdust=FDUST_DEFAULT, Ni=-1.0, Nf=1.5, rtol=1e-12, n_windows=160, n_sub=16,
            Delta_on=True, keep_series=False, k_periods=None):
    t0w = time.time()
    q, mat, z_i, Y_i, V_i, H_i = setup(r, w, fD, fdust, Ni, Delta_on)
    red0, H0f = reduced_toy(mat, q, z_i, Ni, Nf + 0.05, order=0)
    red1, H1f = reduced_toy(mat, q, z_i, Ni, Nf + 0.05, order=1)
    t_f = red1.sol(Nf)[2]
    period = 2 * math.pi / q.m
    if k_periods is None:
        k_periods = max(1, int(round(0.002 * r)))
    win = k_periods * period
    lo, hi = 3 * win, t_f - win
    centres = np.exp(np.linspace(math.log(lo), math.log(hi), n_windows))
    pts = [c + win * (np.arange(n_sub + 1) / n_sub - 0.5) for c in centres]
    keep = [0]
    for i in range(1, len(pts)):
        if pts[i][0] > pts[keep[-1]][-1]:
            keep.append(i)
    pts = np.array([pts[i] for i in keep])
    tt = pts.ravel()
    y0 = [Y_i, V_i, H_i, Ni, 0.0, z_i[0], z_i[1]]
    scale = [abs(Y_i), abs(V_i) + abs(Y_i) * H_i, H_i, 1.0, 1.0, z_i[0], z_i[1]]
    nf = [0]
    out = integrate_full(q, mat, y0, 0.0, tt, rtol=rtol, atol=rtol * 1e-4, scale=scale, counter=nf)
    out = out.reshape(pts.shape[0], pts.shape[1], 7)
    rows = []
    C0 = None
    for k in range(pts.shape[0]):
        o = out[k]
        Y, V, H, N, J, rr, rd = o.T
        N1, N2, J1, J2 = N[0], N[-1], J[0], J[-1]
        Hbar = (J2 - J1) / (N2 - N1)
        u1, u2 = red0.sol(N1), red0.sol(N2)
        Hbar0 = (u2[3] - u1[3]) / (N2 - N1)
        v1, v2 = red1.sol(N1), red1.sol(N2)
        Hbar1 = (v2[3] - v1[3]) / (N2 - N1)
        jc = n_sub // 2
        Nc_, Hc_ = N[jc], H[jc]
        H0p = H0f(Nc_, *red0.sol(Nc_)[:2])
        H1p = H1f(Nc_, *red1.sol(Nc_)[:2])
        # quasi-static subtraction (O(m^-4)) and oscillation energy, averaged over the window
        Eo, rV, Sq, Kk = [], [], [], []
        for j in range(len(N)):
            Yq, Vq, _, _ = quasi_static_toy(q, mat, N[j], rr[j], rd[j], H1f(N[j], *red1.sol(N[j])[:2]))
            Eo.append(osc_energy(q, Y[j], V[j], Yq, Vq) * math.exp(3 * N[j]))
            rV.append(q.rhoV(Y[j]))
            Sq.append(mat.S((rr[j], rd[j]), N[j]) ** 2 / (12 * q.m2))
            # KV kick formula (9/8) kappa^2 H^4 K'^2/eps'' = (3/4) kappa^2 H^4 K0/m^2, kappa = S/rho_m
            kap = mat.S((rr[j], rd[j]), N[j]) / mat.rho((rr[j], rd[j]), N[j])
            Kk.append(0.75 * kap ** 2 * H[j] ** 4 / q.m2)
        Eo, rV, Sq, Kk = np.array(Eo[:-1]), np.array(rV[:-1]), np.array(Sq[:-1]), np.array(Kk[:-1])
        C, Crel = q.constraint(Y[jc], V[jc], H[jc], mat.rho((rr[jc], rd[jc]), N[jc]))
        Ca4 = C * math.exp(4 * N[jc])
        if C0 is None:
            C0 = (C, Ca4)
        rows.append(dict(N=float(0.5 * (N1 + N2)), t=float(pts[k][jc]), H=float(Hbar),
                         eps0_win=float(Hbar / Hbar0 - 1), eps1_win=float(Hbar / Hbar1 - 1),
                         pred_corr_win=float(Hbar1 / Hbar0 - 1),
                         eps2_win=float(Hbar / math.sqrt(Hbar1 ** 2 + Eo.mean() * math.exp(-1.5 * (N1 + N2)) / 3) - 1),
                         eps0_pt=float(Hc_ / H0p - 1), eps1_pt=float(Hc_ / H1p - 1), pred_corr_pt=float(H1p / H0p - 1),
                         Eosc_a3=float(Eo.mean()), Eosc_a3_spread=float(Eo.std()),
                         rhoV_mean=float(rV.mean()), S2_12m2_mean=float(Sq.mean()), kick_mean=float(Kk.mean()),
                         C_rel=float(Crel), C_a4=float(Ca4), m_over_H=float(q.m / Hbar)))
    # transition-centre H and summary
    Hc = float(H1f(0.0, *red1.sol(0.0)[:2]))
    res = dict(r=r, w=w, fD=fD, fdust=fdust, Ni=Ni, Nf=Nf, rtol=rtol, Delta=mat.Delta, m=q.m, H_c=Hc,
               m_over_Hc=q.m / Hc, k_periods=k_periods, n_rhs=nf[0], wall_s=time.time() - t0w,
               x=math.pi * q.m * (w / Hc) / 2)
    res['max_abs_C_rel'] = float(max(abs(rw['C_rel']) for rw in rows))
    res['C_a4_rel_drift'] = float(max(abs(rw['C_a4'] - rows[0]['C_a4']) for rw in rows) /
                                  (3 * H_i ** 2 * math.exp(4 * Ni)))
    if keep_series:
        res['series'] = rows
    res['rows'] = rows
    return res


def summarize_reduction(res, Nmin=-0.25):
    rows = [rw for rw in res['rows'] if rw['N'] >= Nmin]
    pick = lambda key: float(max(abs(rw[key]) for rw in rows))
    at = lambda Nq, key: float(np.interp(Nq, [rw['N'] for rw in res['rows']], [rw[key] for rw in res['rows']]))
    s = dict(r=res['r'], m_over_Hc=res['m_over_Hc'], w=res['w'], x=res['x'], rtol=res['rtol'],
             max_eps0_win=pick('eps0_win'), max_eps1_win=pick('eps1_win'), max_eps2_win=pick('eps2_win'),
             max_pred_corr=pick('pred_corr_win'),
             max_eps0_pt=pick('eps0_pt'), max_eps1_pt=pick('eps1_pt'),
             eps0_win_at={str(n): at(n, 'eps0_win') for n in (0.0, 0.5, 1.0)},
             eps1_win_at={str(n): at(n, 'eps1_win') for n in (0.0, 0.5, 1.0)},
             eps2_win_at={str(n): at(n, 'eps2_win') for n in (0.0, 0.5, 1.0)},
             eps0_pt_at={str(n): at(n, 'eps0_pt') for n in (0.0, 0.5, 1.0)},
             eps1_pt_at={str(n): at(n, 'eps1_pt') for n in (0.0, 0.5, 1.0)},
             pred_corr_at={str(n): at(n, 'pred_corr_win') for n in (0.0, 0.5, 1.0)},
             max_abs_C_rel=res['max_abs_C_rel'], C_a4_rel_drift=res['C_a4_rel_drift'],
             n_rhs=res['n_rhs'], wall_s=res['wall_s'])
    late = [rw for rw in res['rows'] if rw['N'] > 0.6]
    ratio = [rw['rhoV_mean'] / rw['S2_12m2_mean'] for rw in late]
    s['rhoV_over_S2_12m2_late'] = [float(min(ratio)), float(max(ratio))]
    near = [rw for rw in res['rows'] if abs(rw['N']) < 0.3]
    ratio2 = [rw['rhoV_mean'] / rw['S2_12m2_mean'] for rw in near]
    s['rhoV_over_S2_12m2_transition'] = [float(min(ratio2)), float(max(ratio2))] if ratio2 else None
    ratio3 = [rw['rhoV_mean'] / rw['kick_mean'] for rw in res['rows'] if rw['N'] > -0.3]
    s['rhoV_over_KV_kick_formula'] = [float(min(ratio3)), float(max(ratio3))]
    return s


def deposit_plateau(res, w):
    Nmin = max(0.4, 6 * w + 0.2)
    late = [rw for rw in res['rows'] if rw['N'] > Nmin]
    E = np.array([rw['Eosc_a3'] for rw in late])
    return float(E.mean()), float(E.std() / E.mean()) if E.mean() else None, len(E)


def job(args):
    kind, kw = args
    res = run_toy(**kw)
    if kind == 'red':
        return kind, kw, summarize_reduction(res), res
    if kind == 'dep':
        return kind, kw, deposit_plateau(res, kw['w']), res
    return kind, kw, None, res


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                                                                  'results')
    T0 = time.time()
    jobs = []
    # (A) reduction-error scan: adiabatic transition (w = 1: x = pi r/2 >= 15)
    # r >= 3e4 is not used: in these variables dV/dt is a difference S/3 - m^2 F(Y) of O(1) terms that cancel to
    # O(H^2/m^2), so round-off (1e-16) exceeds the requested 1e-12 accuracy and DOP853 stalls (limit, NOTES 2.5)
    R_LIST = [10, 30, 100, 300, 1000, 3000, 10000]
    NI = {10: -0.25, 30: -0.6}
    for r in R_LIST:
        Ni = NI.get(r, -1.0)
        # rtol 1e-12: at 1e-13 DOP853 reaches the round-off limit for r >= 3000 (step count x8, same results)
        jobs.append(('red', dict(r=r, w=1.0, Ni=Ni, Nf=1.5, rtol=1e-12)))
    # (A') the same at w = 0.3 (deposit visible at r = 10)
    for r in (10, 30, 100):
        jobs.append(('red', dict(r=r, w=0.3, Ni=NI.get(r, -1.0), Nf=1.5)))
    # (B) deposit scan: x = pi m tau/2 with tau = w/H_c
    X_LIST = [0.1, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0]
    for r in (30, 100, 300, 1000):
        for x in X_LIST:
            w = 2 * x / (math.pi * r)
            jobs.append(('dep', dict(r=r, w=w, Ni=-0.5 if r <= 30 else -0.7, Nf=max(1.0, 6 * w + 0.9), n_windows=80,
                                     n_sub=16)))
        jobs.append(('ctl', dict(r=r, w=2.0 / (math.pi * r), Ni=-0.5 if r <= 30 else -0.7, Nf=1.0, n_windows=80,
                                 Delta_on=False)))
    # Delta^2 scaling check
    jobs.append(('dep', dict(r=100, w=2 * 1.0 / (math.pi * 100), fD=0.2, Ni=-0.7, Nf=1.0, n_windows=80)))
    jobs.append(('dep', dict(r=100, w=2 * 1.0 / (math.pi * 100), fD=0.0125, Ni=-0.7, Nf=1.0, n_windows=80)))
    # (C) step-size / tolerance study
    for rtol in (1e-9, 1e-10, 1e-11, 1e-12, 1e-13):
        jobs.append(('conv_red', dict(r=100, w=1.0, Ni=-1.0, Nf=1.5, rtol=rtol)))
        jobs.append(('conv_dep', dict(r=100, w=2 * 2.0 / (math.pi * 100), Ni=-0.7, Nf=1.0, rtol=rtol, n_windows=80)))
    # (D) window-length check (k periods) at r = 300
    for k in (1, 3):
        jobs.append(('conv_win', dict(r=300, w=1.0, Ni=-1.0, Nf=1.5, k_periods=k)))
    jobs.sort(key=lambda j: -(j[1]['r'] if j[0] == 'red' else 0))
    results = {'reduction': [], 'reduction_w03': [], 'deposit': [], 'control': [], 'delta_scaling': [],
               'convergence_reduction': [], 'convergence_deposit': [], 'window_check': []}
    series = {}
    with ProcessPoolExecutor(max_workers=4) as ex:
        for kind, kw, summ, res in ex.map(job, jobs):
            print(kind, kw.get('r'), kw.get('w'), kw.get('rtol'), 'wall %.1f s' % res['wall_s'], flush=True)
            if kind == 'red':
                (results['reduction'] if kw['w'] == 1.0 else results['reduction_w03']).append(summ)
                if kw['r'] in (10, 100, 1000) and kw['w'] == 1.0:
                    series['reduction_r%d' % kw['r']] = [{k: rw[k] for k in ('N', 'H', 'eps0_win', 'eps1_win', 'eps2_win',
                                                                             'pred_corr_win', 'C_rel', 'C_a4',
                                                                             'rhoV_mean', 'S2_12m2_mean')}
                                                         for rw in res['rows']]
            elif kind == 'dep':
                E, spread, n = summ
                Hc, m = res['H_c'], res['m']
                x = res['x']
                FF = (x / math.sinh(x)) ** 2
                Esudden = (4.0 / 3.0) * res['Delta'] ** 2 / m ** 2
                d = dict(r=kw['r'], m_over_Hc=res['m_over_Hc'], w=kw['w'], fD=res['fD'], x=x,
                         tau_Hc=kw['w'], E_a3=E, E_rel_spread=spread, n_windows_used=n, E_sudden_formula=Esudden,
                         FF_measured=E / Esudden, FF_formula=FF, ratio_meas_formula=E / (Esudden * FF),
                         max_abs_C_rel=res['max_abs_C_rel'], C_a4_rel_drift=res['C_a4_rel_drift'])
                if kw.get('fD', FD_DEFAULT) != FD_DEFAULT:
                    results['delta_scaling'].append(d)
                else:
                    results['deposit'].append(d)
                if kw['r'] == 100 and abs(x - 2.0) < 0.05:
                    series['deposit_r100_x2'] = [{k: rw[k] for k in ('N', 'Eosc_a3', 'Eosc_a3_spread', 'H')}
                                                 for rw in res['rows']]
            elif kind == 'ctl':
                E, spread, n = deposit_plateau(res, kw['w'])
                results['control'].append(dict(r=kw['r'], E_a3_floor=E,
                                               E_sudden_formula_ref=(4.0 / 3.0) * (3 * FD_DEFAULT) ** 2 / res['m'] ** 2,
                                               floor_over_sudden=E / ((4.0 / 3.0) * (3 * FD_DEFAULT) ** 2 /
                                                                      res['m'] ** 2)))
            elif kind == 'conv_red':
                results['convergence_reduction'].append(summarize_reduction(res))
            elif kind == 'conv_dep':
                E, spread, n = deposit_plateau(res, kw['w'])
                results['convergence_deposit'].append(dict(rtol=kw['rtol'], E_a3=E, E_rel_spread=spread,
                                                           n_rhs=res['n_rhs'], x=res['x'],
                                                           max_abs_C_rel=res['max_abs_C_rel']))
            elif kind == 'conv_win':
                s = summarize_reduction(res)
                s['k_periods'] = kw['k_periods']
                results['window_check'].append(s)
    # fits
    red = sorted(results['reduction'], key=lambda d: d['r'])
    lr = np.log([d['m_over_Hc'] for d in red])
    fits = {}
    for key in ('max_eps0_win', 'max_eps1_win', 'max_eps2_win', 'max_eps0_pt', 'max_eps1_pt'):
        vals = np.array([d[key] for d in red])
        mH = np.array([d['m_over_Hc'] for d in red])
        floor = 3e-13
        sel = (vals > floor) & (mH >= 90)
        if sel.sum() >= 2:
            sl, ic = np.polyfit(lr[sel], np.log(vals[sel]), 1)
            fits[key] = dict(slope=float(sl), prefactor=float(math.exp(ic)), n_points=int(sel.sum()),
                             fit_range_m_over_H=[float(mH[sel].min()), float(mH[sel].max())],
                             note='log-log fit for m/H_c >= 90 and values above the %.0e integration floor' % floor)
        sel2 = (vals > floor) & (mH < 90)
        if sel2.sum() >= 2:
            sl, ic = np.polyfit(lr[sel2], np.log(vals[sel2]), 1)
            fits[key + '_small_r'] = dict(slope=float(sl), n_points=int(sel2.sum()))
    dep = results['deposit']
    for r in sorted(set(d['r'] for d in dep)):
        dd = sorted([d for d in dep if d['r'] == r], key=lambda d: d['x'])
        xs = np.array([d['x'] for d in dd if d['x'] >= 2.0])
        Es = np.array([d['FF_measured'] for d in dd if d['x'] >= 2.0])
        if len(xs) >= 3:
            # (x/sinh x)^2 -> 4 x^2 e^{-2x}: fit ln(FF/x^2) = a - kappa x ; formula kappa = 2 (i.e. E ~ exp(-pi m tau))
            sl, ic = np.polyfit(xs, np.log(Es / xs ** 2), 1)
            FFf = (xs / np.sinh(xs)) ** 2
            slf, icf = np.polyfit(xs, np.log(FFf / xs ** 2), 1)
            fits['deposit_r%d' % r] = dict(kappa_fit=float(-sl), kappa_formula=2.0, prefactor_fit=float(math.exp(ic)),
                                           prefactor_formula=4.0, kappa_fit_of_formula_values=float(-slf),
                                           prefactor_fit_of_formula_values=float(math.exp(icf)),
                                           note='fit ln(FF/x^2) = ln A - kappa x for x >= 2; E ~ exp(-kappa x)')
    results['fits'] = fits
    results['meta'] = dict(units='K0 = 1, time unit 1/H_c (H at transition centre); m = m_q/H_c',
                           shape='KV 0806.2805: b3 = 4/3, b4 = 1/3, beta = 1',
                           matter='radiation + dust (3*0.02 at N=0) + step Delta = 3*0.05 (tanh in N, width w)',
                           wall_total_s=time.time() - T0)
    results['series'] = series
    with open(os.path.join(outdir, 'reduction_validation.json'), 'w') as fh:
        json.dump(results, fh, indent=1)
    print(json.dumps(results['fits'], indent=1))
    print('total wall', time.time() - T0)


if __name__ == '__main__':
    main()
