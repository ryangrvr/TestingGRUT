#!/usr/bin/env python3
"""dq_heat (Claude, own code), step 0: the M1r background on the fixed engine copy rptri_dq.c.

M1r: lambda_Y = lambda_Sigma = 0.5, inflaton trilinear holder (g_I e^{i alpha_I} phi_I Sigma* Y^2), (m_I, Gamma_I) = (3, 0.5) H_i,
prepared at u = 10 i (the conditional minimum, fx_holder section 4).  Two runs of the same Lagrangian and initial state:
  coarse: tau = H_i t in [0, 2000], samples every 0.2 (validation against fx_holder: N_X(400-1000), E_X, N_YL, N_YH, R_r);
  dense : tau in [0, 200], samples every 0.01 (the Higgs norm F(t) and a(t) for the vector-mode solver).
Records |X|, |Y|, |Sigma|(t), F^2 = 2 (|X|^2 + 9|Y|^2 + 36|Sigma|^2), a, H, the conformal time eta = int dt/a and the late
radiation constant K = a^2 H (a' = da/deta = K in radiation).  Output: ../results/background.json, scratch npz.
"""
import json, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import dq_model as fc

RES = Path(__file__).resolve().parent.parent / 'results'
B, GAM, LAM = 3.0, 0.5, 0.5
Q = np.array([1.0, 3.0, 6.0])


def series(D):
    """D columns: t, y[0..18], H, R, Hdot.  y: 0 ln a, 1 I, 2 J, 3..8 fields (re, im), 9..14 velocities, 15 rho_r, 16 Haux, 17 W, 18 rhoI."""
    t = D[:, 0]; y = D[:, 1:20]; H, R, Hd = D[:, 20], D[:, 21], D[:, 22]
    z = np.stack([y[:, 3] + 1j * y[:, 4], y[:, 5] + 1j * y[:, 6], y[:, 7] + 1j * y[:, 8]], axis=1)
    v = np.stack([y[:, 9] + 1j * y[:, 10], y[:, 11] + 1j * y[:, 12], y[:, 13] + 1j * y[:, 14]], axis=1)
    a = np.exp(y[:, 0])
    F2 = 2 * np.sum(Q ** 2 * np.abs(z) ** 2, axis=1)
    # conformal time by the trapezoid on the dense samples (a smooth; error << 1e-6 at dt = 0.01)
    eta = np.concatenate([[0.0], np.cumsum(0.5 * (1 / a[1:] + 1 / a[:-1]) * np.diff(t))])
    return dict(t=t, a=a, H=H, R=R, Hdot=Hd, z=z, v=v, F2=F2, eta=eta, I=y[:, 1], J=y[:, 2], rho_r=y[:, 15])


def main():
    t0 = time.time()
    p = fc.params(LAM)
    st = fc.stationary(fc.U_PREP, p)
    zh = st['z']
    G, I0, ph0 = fc.trilinear_G(zh, p, B)
    out = dict(prep=dict(z=[[c.real, c.imag] for c in zh], amps=np.abs(zh).tolist(), hessian_eig=st['eig'].tolist(), V=st['V'],
                         G=[G.real, G.imag], g_I=abs(G), alpha_I=float(np.angle(G)), I0=I0, phihat0=ph0,
                         F2_inf=float(2 * np.sum(Q ** 2 * np.abs(zh) ** 2)), F_inf=float(np.sqrt(2 * np.sum(Q ** 2 * np.abs(zh) ** 2))),
                         F_vac=6 * p['h'], X_part_of_F2=float(2 * abs(zh[0]) ** 2 / (2 * np.sum(Q ** 2 * np.abs(zh) ** 2)))))
    rc = fc.run('m1r_coarse', p, B, GAM, z0=zh, G=G, end=2000.0, nsamp=10000, binw=20.0)
    rd = fc.run('m1r_dense', p, B, GAM, z0=zh, G=G, end=200.0, nsamp=20000, binw=20.0)
    # validation: fx_holder conventions (400-1000 window, 80-wide running windows)
    nx = fc.window_stats(rc, 'NX', 400.0, 1000.0, 80.0)
    val = dict(NX_400_1000=nx['mean'], NX_run80=[nx['run_min'], nx['run_max']], NX_1000_2000=fc.window_mean(rc, 'NX', 1000.0, 2000.0),
               checks=rc['END'])
    for k in ('NYL', 'NYH', 'Nsig', 'EX', 'eps', 'NS', 'NY'):
        val[k] = fc.window_mean(rc, k, 400.0, 1000.0)
    val['EX_1000_2000'] = fc.window_mean(rc, 'EX', 1000.0, 2000.0)
    val['NYL_1000_2000'] = fc.window_mean(rc, 'NYL', 1000.0, 2000.0)
    Y, ratio, Rr = fc.yield_over_target(rc, val['NX_400_1000'])
    val.update(yield_over_YQ=ratio, Rr_end=Rr, K_late=float(np.sqrt(Rr / 3.0)))
    fxh = json.load(open(Path(__file__).resolve().parents[2] / 'fx_holder' / 'results' / 'tri_runs.json'))['runs']['lam0.5_b3_g0.5']['trilinear']
    val['fx_holder_reference'] = dict(NX_400_1000=fxh['NX_400_1000'], EX=fxh['EX'], NYL=fxh['NYL'], NYH=fxh['NYH'], Nsig=fxh['Nsig'],
                                      eps=fxh['eps'], Rr_end=fxh['Rr_end'], yield_over_YQ=fxh['yield_over_YQ_400_1000'])
    val['rel_diff_vs_fx_holder'] = {k: float(val[k] / fxh[k] - 1) for k in ('NX_400_1000', 'EX', 'NYL', 'NYH', 'eps')}
    out['validation'] = val
    # early N_YL / N_YH with the fixed diagnostic (the old one was wrong while Im u != 0): bins of 20
    S = rc['S']; C = fc.COLS
    early = {}
    for tt in (0.0, 0.5, 1, 2, 3, 5, 10, 20, 50, 100, 200, 400):
        i = int(np.argmin(abs(S[:, 0] - tt)))
        early[str(tt)] = dict(NYL=S[i, C.index('NYL')], NYH=S[i, C.index('NYH')], Nsig=S[i, C.index('Nsig')], NX=S[i, C.index('NX')],
                              ure=S[i, C.index('ure')], uim=S[i, C.index('uim')], mL2_forced=S[i, C.index('mL2')], a=S[i, C.index('a')])
    out['early_debris_fixed_diag'] = early
    B20 = rc['bins']
    out['bins20'] = dict(t0=B20[:, 0].tolist(), NYL=fc.bin_col(rc, 'NYL').tolist(), NYH=fc.bin_col(rc, 'NYH').tolist(),
                         Nsig=fc.bin_col(rc, 'Nsig').tolist(), EX=fc.bin_col(rc, 'EX').tolist(), NX=fc.bin_col(rc, 'NX').tolist())
    # series for the mode solvers
    sd, sc = series(rd['D']), series(rc['D'])
    np.savez(fc.SP / 'runs' / 'm1r_background.npz', **{'d_' + k: v for k, v in sd.items()}, **{'c_' + k: v for k, v in sc.items()})
    # summary of F(t), |X|, |Y|, |Sigma|
    tab = []
    for tt in (0, 0.25, 0.5, 1, 1.5, 2, 3, 5, 7, 10, 15, 20, 30, 50, 100, 200):
        i = int(np.argmin(abs(sd['t'] - tt)))
        tab.append(dict(tau=float(sd['t'][i]), a=float(sd['a'][i]), H=float(sd['H'][i]), eta=float(sd['eta'][i]),
                        absX=float(abs(sd['z'][i, 0])), absY=float(abs(sd['z'][i, 1])), absS=float(abs(sd['z'][i, 2])),
                        F=float(np.sqrt(sd['F2'][i])), aF=float(sd['a'][i] * np.sqrt(sd['F2'][i])), a2H=float(sd['a'][i] ** 2 * sd['H'][i])))
    out['F_table'] = tab
    # conformal X amplitude chi = a|X| late (for the smooth X part of F^2 after the run): mean of chi^2 over the last 400
    sel = sc['t'] > 1600
    chi2 = np.mean((sc['a'][sel] * np.abs(sc['z'][sel, 0])) ** 2)
    out['late'] = dict(a_end=float(sc['a'][-1]), H_end=float(sc['H'][-1]), eta_end=float(sc['eta'][-1]), K_end=float(sc['a'][-1] ** 2 * sc['H'][-1]),
                       chi2_mean_late=float(chi2), absS_late=float(np.mean(np.abs(sc['z'][sel, 2]))), h_over_sqrt2=p['h'] / np.sqrt(2),
                       dense_vs_coarse_a_at_200=float(sd['a'][-1] / sc['a'][np.argmin(abs(sc['t'] - 200))] - 1))
    out['seconds'] = time.time() - t0
    json.dump(out, open(RES / 'background.json', 'w'), indent=1, default=float)
    print(json.dumps(dict(val={k: v for k, v in val.items() if k != 'fx_holder_reference'}, late=out['late'], prep=out['prep']), indent=1, default=float))
    for r in tab:
        print(r)
    print('seconds', out['seconds'])


if __name__ == '__main__':
    main()
