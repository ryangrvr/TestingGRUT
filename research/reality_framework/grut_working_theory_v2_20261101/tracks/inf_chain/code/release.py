#!/usr/bin/env python3
"""inf_chain item 1 (Claude, own code): the M1r release on each benchmark inflaton, run to freeze-out.

For a benchmark (ch_common registry):
  * the held preparation is re-solved at the slice (u_slice, R = 6(2 - eps_H) H_i^2), g_I re-fixed so that mu_eff = 10i at
    the calibration point (S1 unless stated);
  * DENSE run: tau <= T_dense with fine S rows (background for the vector modes, release diagnostics);
  * LONG run: tau <= T_long with conformal bins (step cap relaxed after the inflaton has decayed, CH1) for the X <-> Sigma
    exchange: conformal windows (width 20 in eta) fitted with N_inf + c2 a^-2 (dq_entropy's method), alternatives
    a^-1 + a^-2 and a constant; uncertainty = largest deviation among them;
  * late values: E_X (conformal), eps = N_X/N_c(E_X), E/E_c = eps^{-4/3}, K = a^2 H (radiation; Gamma_I = 0.5 runs),
    rho_I a^3 (matter era; Gamma_I ~ 0 runs), N_YL, N_YH, N_sigma (late means and spread), Gauss/H/work checks.
Usage: python3 -B release.py KEY [T_dense T_long]     Output: ../results/release_<KEY>.json (series in scratch)
Labels: [C] given the model choices.
"""
import json, math, sys
import numpy as np
import ch_common as ch
from benchmarks_reg import get_bench

C = ch.COLS


def windows(Cb, names, deta=4.0, start=8.0):
    eta0, eta1 = Cb[:, 0], Cb[:, 1]; de = eta1 - eta0
    edges = np.arange(start, eta1[-1] + 1e-9, deta)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (eta0 >= lo - 1e-9) & (eta1 <= hi + 1e-9)
        if s.sum() < 4:
            continue
        w = de[s]; r = dict(eta_lo=float(lo), eta_hi=float(hi))
        for nm in names + ['a', 'tau']:
            v = Cb[s, 4 + C.index(nm)]
            r[nm] = float(np.sum(v * w) / np.sum(w)); r[nm + '_binstd'] = float(np.std(v))
        rows.append(r)
    return rows


def fit(rows, nm, amin, model):
    a = np.array([r['a'] for r in rows]); y = np.array([r[nm] for r in rows])
    s = a >= amin; a, y = a[s], y[s]
    X = {'a-2': np.c_[np.ones_like(a), a ** -2.0], 'a-1,a-2': np.c_[np.ones_like(a), a ** -1.0, a ** -2.0],
         'const': np.c_[np.ones_like(a)]}[model]
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ c; dof = max(len(y) - X.shape[1], 1)
    cov = np.linalg.pinv(X.T @ X) * float(res @ res) / dof
    return dict(model=model, a_min=float(amin), n=int(len(y)), N_inf=float(c[0]), N_inf_err=float(np.sqrt(abs(cov[0, 0]))))


def frozen(rows, nm, a_end):
    A1, A2, A3 = a_end / 6, a_end / 3, a_end / 2
    F = {(m, A): fit(rows, nm, A, m) for A in (A1, A2, A3) for m in ('a-2', 'a-1,a-2', 'const')}
    ad = F[('a-2', A2)]
    alts = [F[('a-2', A1)], F[('a-2', A3)], F[('a-1,a-2', A1)], F[('const', A3)]]
    dev = max(abs(x['N_inf'] - ad['N_inf']) for x in alts)
    return dict(value=ad['N_inf'], err=float(math.hypot(dev, ad['N_inf_err'])), fits=[ad] + alts)


def slope_std(rows, nm):
    a = np.array([r['a'] for r in rows]); s = np.array([r[nm + '_binstd'] for r in rows])
    sel = (a > 60) & (s > 0)
    if sel.sum() < 4:
        return None
    return float(np.polyfit(np.log(a[sel]), np.log(s[sel]), 1)[0])


def main(key, Td=100.0, Tl=3.0e4):
    b = get_bench(key)
    st = b.prep()
    z = st['z']
    info = dict(bench=b.summary(), prep=dict(amps=[float(x) for x in np.abs(z)], z=[[c.real, c.imag] for c in z], F=ch.F_of(z),
                                             eig=[float(x) for x in st['eig']], grad=st['grad']))
    gam = b.gam
    # dense run
    rd = reuse(key + '_dense', Td) or ch.run_engine(b, key + '_dense', gam, Td, int(round(Td / 0.005)), z=z)
    S = rd['S']
    lateD = S[S[:, 0] >= 0.75 * Td]
    info['dense'] = dict(T=Td, END=rd['END'], seconds=rd['seconds'], I0=rd['I0'],
                         late_mean={nm: float(np.mean(lateD[:, C.index(nm)])) for nm in ('NX', 'EX', 'eps', 'NYL', 'NYH', 'Nsig', 'a', 'H')},
                         checks=dict(max_gauss=rd['END']['max_gauss'], max_H=rd['END']['max_H_consistency'], max_work=rd['END']['max_hidden_work']))
    tl = [dict(tau=float(r[0]), a=float(r[1]), H=float(r[2]), NX=float(r[4]), EX=float(r[8]), NYL=float(r[11]), NYH=float(r[10]),
               Nsig=float(r[12]), absX=float(r[13]), absY=float(r[14]), Sv=float(r[15]), rhoI=float(r[17]), radfrac=float(r[18]),
               ure=float(r[26]), uim=float(r[27])) for r in S[::max(1, len(S) // 80)]]
    info['dense']['timeline'] = tl
    print('%s dense tau<=%g: %.1fs, steps %d; late NX=%.5g EX=%.5g eps=%.4f NYL=%.4g NYH=%.4g Nsig=%.3g | gauss %.1e H %.1e work %.1e' % (
        key, Td, rd['seconds'], rd['END']['steps'], *[info['dense']['late_mean'][k] for k in ('NX', 'EX', 'eps', 'NYL', 'NYH', 'Nsig')],
        rd['END']['max_gauss'], rd['END']['max_H_consistency'], rd['END']['max_hidden_work']), flush=True)
    # long run
    capT = 60.0 if gam >= 0.2 else None
    rl = reuse(key + '_long', Tl) or ch.run_engine(b, key + '_long', gam, Tl, 4000, cbin=0.5, capT=capT, z=z)
    S, Cb = rl['S'], rl['C']
    rows = windows(Cb, ['NX', 'NS', 'NY', 'EX', 'NYL', 'NYH', 'Nsig', 'eps', 'rho_r', 'rho_I', 'H'])
    a_end = float(Cb[-1, 4 + C.index('a')])
    fr = {nm: frozen(rows, nm, a_end) for nm in ('NX', 'NS', 'EX', 'NYL', 'NYH', 'Nsig')}
    last = dict(zip(C, S[-1]))
    late = S[S[:, 0] >= 0.5 * Tl]
    lm = {nm: float(np.mean(late[:, C.index(nm)])) for nm in ('NX', 'EX', 'eps', 'NYL', 'NYH', 'Nsig')}
    lsd = {nm: float(np.std(late[:, C.index(nm)]) / max(1e-300, abs(np.mean(late[:, C.index(nm)])))) for nm in ('NYL', 'NYH', 'Nsig')}
    NX = fr['NX']['value']; EX = fr['EX']['value']
    Nc = 2 * math.sqrt(2 * P_lx()) * (EX / (3 * P_lx())) ** 0.75
    eps = abs(NX) / Nc
    K = last['a'] ** 2 * last['H']
    out = dict(info)
    out['long'] = dict(T=Tl, END=rl['END'], seconds=rl['seconds'], a_end=a_end, windows=rows, frozen=fr, late_mean=lm, late_spread=lsd,
                       std_slope_NX=slope_std(rows, 'NX'),
                       last=dict(a=last['a'], H=last['H'], rho_r=last['rho_r'], rho_I=last['rho_I'], radfrac=last['radfrac']))
    res = dict(N_X=NX, N_X_err=fr['NX']['err'], N_Sigma=fr['NS']['value'], E_X=EX, E_X_err=fr['EX']['err'], eps=eps, E_over_Ec=eps ** (-4 / 3),
               N_YL=fr['NYL']['value'], N_YH=fr['NYH']['value'], N_sigma=fr['Nsig']['value'],
               N_X_over_window_400_1000=None, gauss=rl['END']['max_gauss'], H_consistency=rl['END']['max_H_consistency'],
               hidden_work=rl['END']['max_hidden_work'], normalisation='a = 1 at the slice (%s)' % ('S1' if b.Nslice == 1 else b.Nslice))
    if gam >= 0.2:
        res['K'] = K; res['Rr'] = last['rho_r'] * last['a'] ** 4
        res['radiation_fraction_end'] = last['radfrac']
        # yield at the M1 H_i in fx_holder's convention: Y = |N_X| H^3/(a^3 s), s from rho_r
        rho_rad = last['rho_r'] * ch.MP ** 2 * ch.H_M1 ** 2
        T = (30 * rho_rad / (math.pi ** 2 * ch.GSTAR)) ** 0.25; s = 2 * math.pi ** 2 / 45 * ch.GSTAR * T ** 3
        res['Y_over_YQ_at_HM1_lastrow'] = abs(last['NX']) * ch.H_M1 ** 3 / (last['a'] ** 3 * s) / ch.YQ_TARGET
        res['Y_over_YQ_at_HM1_frozen'] = res['Y_over_YQ_at_HM1_lastrow'] * abs(NX) / abs(last['NX'])
    else:
        # matter era: x0 = <rho_I a^3> over the last windows (Gamma ~ 0 in the engine; the decay is added in the fluid)
        xw = late[:, C.index('rho_I')] * late[:, C.index('a')] ** 3     # per sample (constant once w averages to 0)
        res['x0'] = float(np.mean(xw)); res['x0_spread'] = float(np.std(xw) / np.mean(xw))
        res['w_note'] = 'Gamma_I = %g in the engine; K and the yield follow from the fluid continuation (vectors.py / chain.py)' % gam
    out['result'] = res
    json.dump(out, open(ch.RES / ('release_%s.json' % key), 'w'), indent=1, default=float)
    print('%s long tau<=%g: %.1fs steps %d | N_X=%.6g+-%.2g (std slope %.2f) E_X=%.6g eps=%.5f E/Ec=%.4f N_YL=%.5g N_YH=%.5g N_sig=%.4g %s | spreads %s | gauss %.1e H %.1e work %.1e' % (
        key, Tl, rl['seconds'], rl['END']['steps'], NX, fr['NX']['err'], out['long']['std_slope_NX'] or float('nan'), EX, eps, eps ** (-4 / 3),
        res['N_YL'], res['N_YH'], res['N_sigma'], ('K=%.5f Y/YQ(HM1)=%.5f' % (K, res['Y_over_YQ_at_HM1_frozen'])) if gam >= 0.2 else ('x0=%.5g' % res['x0']),
        {k: '%.1e' % v for k, v in lsd.items()}, rl['END']['max_gauss'], rl['END']['max_H_consistency'], rl['END']['max_hidden_work']), flush=True)


def reuse(tag, T):
    """reuse a saved run of the same length (set CH_RERUN=1 to force a new run)"""
    import os
    f = ch.RUNS / (tag + '.npz'); mf = ch.RUNS / (tag + '_meta.json')
    if os.environ.get('CH_RERUN') or not f.exists() or not mf.exists():
        return None
    meta = json.load(open(mf))
    if abs(meta['end'] - T) > 1e-9:
        return None
    d = np.load(f)
    return dict(S=d['S'], C=d['C'], END=meta['END'], I0=meta['I0'], W=meta['W'], seconds=meta['seconds'])


def P_lx():
    return ch.P05['lx']


if __name__ == '__main__':
    k = sys.argv[1]
    Td = float(sys.argv[2]) if len(sys.argv) > 2 else 100.0
    Tl = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0e4
    main(k, Td, Tl)
