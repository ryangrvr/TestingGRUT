#!/usr/bin/env python3
"""dq_entropy item 1 (Claude, own code): X <-> Sigma exchange to freeze-out from the long rpdq.c runs.

For each run: conformal-time windows (width DETA in eta) of N_X, N_Sigma, N_Y, E_X; least-squares extrapolation of the
window means to a -> infinity with N(a) = N_inf + c2 a^-2 (the exchange operator X^6 Sigma* enters the conformal X equation
with a factor a^-2), plus the variants N_inf + c1/a + c2/a^2 and a constant, over several late ranges; the frozen
yield Y/Y_Q at the M1 H_i (fx_holder's convention: Y = |N_X| H^3/(a^3 s), s from the radiation at the end of the run,
Y_Q = 5.07941423e-9); E_X late (conformal) and eps = N_X/N_c(E_X) with N_c = 2 sqrt(2 lX)(E/(3 lX))^{3/4}; the decay ledger
of fx_holder (fx_common.ledger, read-only) on the late N_YL, N_YH, N_sigma; comparison with the 400-1000 window.
Usage: freeze_analysis.py OUT.json TAG [TAG ...]
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import dq_common as dq
fc = dq.fc
C = dq.COLS
DETA = 20.0


def series(tag):
    d = np.load(dq.SP / 'runs' / (tag + '.npz'))
    meta = json.load(open(dq.SP / 'runs' / (tag + '_meta.json')))
    return d, meta


def windows(Cb, names, deta=DETA, start=20.0):
    eta0, eta1 = Cb[:, 0], Cb[:, 1]; de = eta1 - eta0
    edges = np.arange(start, eta1[-1] + 1e-9, deta)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (eta0 >= lo - 1e-9) & (eta1 <= hi + 1e-9)
        if s.sum() < 4:
            continue
        w = de[s]; r = dict(eta_lo=float(lo), eta_hi=float(hi))
        for nm in names + ['a']:
            v = Cb[s, 4 + C.index(nm)] if nm in C else None
            r[nm] = float(np.sum(v * w) / np.sum(w))
            r[nm + '_binstd'] = float(np.std(v))
        rows.append(r)
    return rows


def fit(rows, nm, amin, model):
    a = np.array([r['a'] for r in rows]); y = np.array([r[nm] for r in rows])
    s = a >= amin
    a, y = a[s], y[s]
    if model == 'a-2':
        X = np.c_[np.ones_like(a), a ** -2.0]
    elif model == 'a-1,a-2':
        X = np.c_[np.ones_like(a), a ** -1.0, a ** -2.0]
    else:
        X = np.c_[np.ones_like(a)]
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ c
    dof = max(len(y) - X.shape[1], 1)
    cov = np.linalg.inv(X.T @ X) * float(res @ res) / dof
    return dict(model=model, a_min=float(amin), n=int(len(y)), N_inf=float(c[0]), N_inf_err=float(np.sqrt(cov[0, 0])), coef=c.tolist(),
                rms_resid=float(np.sqrt(np.mean(res ** 2))))


def main():
    out = {}
    for tag in sys.argv[2:]:
        d, meta = series(tag)
        Cb, B, S = d['cbins'], d['bins'], d['S']
        rows = windows(Cb, ['NX', 'NS', 'NY', 'EX', 'NYL', 'NYH', 'Nsig', 'eps'])
        fits = {}
        a_end_run = float(Cb[-1, 4 + C.index('a')])
        A1, A2, A3 = (round(a_end_run / k, 1) for k in (6.0, 3.0, 2.0))     # 105, 210, 315 for the 1e5 runs
        for nm in ('NX', 'NS', 'EX'):
            fits[nm] = [fit(rows, nm, amin, mdl) for amin in (A1, A2, A3) for mdl in ('a-2', 'a-1,a-2', 'const')]
        # adopted: a^-2 fit on a >= a_end/3; uncertainty = max deviation among {a-2 fits a >= a_end/6, a_end/2;
        # a-1,a-2 on a >= a_end/6; const on a >= a_end/2} combined with its own standard error
        frozen = {}
        for nm in ('NX', 'NS', 'EX'):
            F = {(f['model'], f['a_min']): f for f in fits[nm]}
            ad = F[('a-2', A2)]
            alts = [F[('a-2', A1)], F[('a-2', A3)], F[('a-1,a-2', A1)], F[('const', A3)]]
            spread = max(abs(f['N_inf'] - ad['N_inf']) for f in alts)
            frozen[nm] = dict(value=ad['N_inf'], stat_err=ad['N_inf_err'], model_spread=spread,
                              err=float(np.hypot(ad['N_inf_err'], spread)), c2=ad['coef'][1])
        late = [r for r in rows if r['a'] >= A3]
        NY_late = float(np.mean([r['NY'] for r in late]))
        gauss = frozen['NX']['value'] + 6 * frozen['NS']['value'] + 3 * NY_late
        # yield at the M1 H_i (fx_holder convention) with the end-of-run radiation
        last = dict(zip(C, S[-1]))
        a_end, rr_end = last['a'], last['rho_r']
        rho_rad = rr_end * fc.MPL ** 2 * fc.H_GEV ** 2
        T = (30 * rho_rad / (np.pi ** 2 * fc.GSTAR)) ** 0.25
        s = 2 * np.pi ** 2 / 45 * fc.GSTAR * T ** 3
        Yfac = fc.H_GEV ** 3 / (a_end ** 3 * s)
        NXf, dNX = frozen['NX']['value'], frozen['NX']['err']
        y_over = abs(NXf) * Yfac / fc.YQ_TARGET
        # 400-1000 physical window (fx_holder's late convention) for comparison
        nx_400 = dq.window_mean(B, 'NX', 400.0, 1000.0)
        ex_400 = dq.window_mean(B, 'EX', 400.0, 1000.0)
        EXf = frozen['EX']['value']
        lx = 1e-3
        Nc = 2 * np.sqrt(2 * lx) * (EXf / (3 * lx)) ** 0.75
        eps = abs(NXf) / Nc
        c_const = abs(NXf) / (eps * EXf ** 0.75)
        late_deb = {k: float(np.mean([r[k] for r in late])) for k in ('NYL', 'NYH', 'Nsig')}
        p = fc.params(0.5)
        Rr = rr_end * a_end ** 4
        led = {('%g' % e): fc.ledger(late_deb['NYH'], late_deb['NYL'], late_deb['Nsig'], EXf, Rr, p, e=e) for e in (1e-7, 1e-6, 2e-6, 3e-6, 5e-6)}
        o = dict(meta=dict(cport_units=meta['cport_units'], end=meta['end'], eta_end=meta['eta_end'], a_end=a_end, END=meta['END'],
                           prep_amps=meta['prep_amps'], prep_eig=meta['prep_eig'][:4], du_prep=meta['du_prep'], seconds=meta['seconds']),
                 windows=rows, fits=fits, frozen=frozen, NY_late=NY_late, gauss_frozen=gauss,
                 Rr_end=float(Rr), yield_factor=float(Yfac),
                 yield_over_YQ_frozen=float(y_over), yield_over_YQ_err=float(dNX * Yfac / fc.YQ_TARGET),
                 NX_400_1000=nx_400, yield_over_YQ_400_1000=float(abs(nx_400) * Yfac / fc.YQ_TARGET), EX_400_1000=ex_400,
                 frozen_over_400_1000=NXf / nx_400,
                 H_new_Q2e5_relic_GeV=float(fc.H_GEV * y_over ** (-2 / 3)),
                 E_X_late=EXf, eps_late=float(eps), E_over_Ec=float(eps ** (-4 / 3)), c_NcE34=float(c_const),
                 energy_per_charge_conf=float(EXf / abs(NXf)), debris_late=late_deb, ledger=led)
        out[tag] = o
        print('%s: N_X frozen = %.4f +- %.4f (stat %.4f, model %.4f); N_S = %.4f; N_Y late = %.2e; gauss %.2e' % (
            tag, NXf, dNX, frozen['NX']['stat_err'], frozen['NX']['model_spread'], frozen['NS']['value'], NY_late, gauss))
        print('   400-1000: %.4f  ratio frozen/400-1000 = %.6f;  Y/Y_Q frozen = %.6f +- %.6f (400-1000: %.6f); H(Q=2e5 relic) = %.5e' % (
            nx_400, NXf / nx_400, y_over, o['yield_over_YQ_err'], o['yield_over_YQ_400_1000'], o['H_new_Q2e5_relic_GeV']))
        print('   E_X late = %.3f +- %.3f (400-1000: %.3f); eps = %.6f; E/E_c = %.4f; c = %.6f; E_X/N_X = %.4f; debris %s; R(1e-6) = %.4f (Y %.4f casc %.4f)' % (
            EXf, frozen['EX']['err'], ex_400, eps, eps ** (-4 / 3), c_const, EXf / abs(NXf), late_deb, led['1e-06']['R_central'], led['1e-06']['R_Y'], led['1e-06']['R_cascade']))
    if 'base_1e5' in out and 'base_2e4_rtol12' in out:
        b, r12 = out['base_1e5'], out['base_2e4_rtol12']
        wb = {(w['eta_lo'], w['eta_hi']): w for w in b['windows']}
        diffs = [abs(w['NX'] - wb[(w['eta_lo'], w['eta_hi'])]['NX']) for w in r12['windows'] if (w['eta_lo'], w['eta_hi']) in wb]
        out['rtol_control'] = dict(max_abs_window_NX_diff=float(max(diffs)), n_windows=len(diffs),
                                   frozen_diff=r12['frozen']['NX']['value'] - b['frozen']['NX']['value'])
        print('rtol 1e-12 control: max |dNX| over common eta windows = %.3e; frozen(2e4, rtol12) - frozen(1e5) = %.4f' % (
            max(diffs), out['rtol_control']['frozen_diff']))
    if 'cp_plus1_2e4' in out and 'cp_minus1_2e4' in out:
        pp, mm = out['cp_plus1_2e4'], out['cp_minus1_2e4']
        out['portal_derivative_per_unit'] = dict(
            dlnNX=(np.log(abs(pp['frozen']['NX']['value'])) - np.log(abs(mm['frozen']['NX']['value']))) / 2,
            dE_over_Ec=(pp['E_over_Ec'] - mm['E_over_Ec']) / 2,
            dR_med_1e6=(pp['ledger']['1e-06']['R_central'] - mm['ledger']['1e-06']['R_central']) / 2)
        print('per unit ln Q:', out['portal_derivative_per_unit'])
    if 'base_1e5' in out and 'planck_1e5' in out:
        b, q = out['base_1e5'], out['planck_1e5']
        out['planck_over_base'] = dict(NX=q['frozen']['NX']['value'] / b['frozen']['NX']['value'],
                                       yield_ratio=q['yield_over_YQ_frozen'] / b['yield_over_YQ_frozen'],
                                       H_ratio_same_relic=(q['yield_over_YQ_frozen'] / b['yield_over_YQ_frozen']) ** (-2 / 3))
        print('Planck/base:', out['planck_over_base'])
    json.dump(out, open(sys.argv[1], 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
