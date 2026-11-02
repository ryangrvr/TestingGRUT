#!/usr/bin/env python3
"""mf_ideal item 1 (Claude, own code): the M1r release on an idealised fast history with the holder on a spectator s.

A history (model choice) is named by a tag and specified by
  bg  : 'fluid'  -- constant-H inflation to the slice (a = 1, H = H_i = 1, R = 12 just before), then dust rho_m = 3
                    (R = rho_m) decaying into radiation at the rate gam(K) tuned to K = a^2 H (K = 1: pure radiation,
                    R = 0 from the slice).  The step of R at the slice is the heavy-waterfall (m_chi >> H) limit;
        'table:<file>' -- any MODE 5 table (bench histories, oscillator backgrounds, smooth ends; see README);
  s   : m_s (H_i), t_on (H_i^-1, smootherstep turn-on of m_s^2 from 0), Gamma_s (H_i), s_* (GeV; back-reaction
        (H_i/s_*)^2), c0 = s(slice)/s_* (1 = the calibrated holder, mu_eff = 10i at the slice; c0 != 1 is the
        isocurvature variation at fixed g_s).
The preparation is re-solved at u = mu + (10i - mu) c0 and R = 12 (constant-H inflation, s frozen).
Runs (rpmfi.c MODE 5): DENSE tau <= 100 (dt = 0.005; background and F(t) for the vectors) and LONG tau <= T_long with
conformal bins (0.5 in eta) for the X <-> Sigma exchange, fitted as N_inf + c2 a^-2 over conformal windows (inf_chain /
dq_entropy method, alternatives a^-1 + a^-2 and a constant; the uncertainty is the largest deviation).
Usage: python3 -B release_mfi.py TAG K MS TON [GS SSTAR C0 BG T_LONG]
Output: ../results/release/release_<TAG>.json (+ series in scratch).  Labels: [C] given the model choices.
"""
import json, math, sys
import numpy as np
import mfi_common as mc

C = mc.COLS
OUT = mc.RES / 'release'


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


def spec_env(K, ms, ton, gs, sstar, c0, bg, xis=0.0):
    if bg == 'fluid':
        env, gam = mc.history_env(K, ms, ton, gs=gs, br=mc.br_of(sstar), c0=c0)
    else:
        kind, path = bg.split(':', 1)            # 'table:' (s evolves) or 'tables:' (s from the table's c, cdot columns)
        env = dict(MFI_TABLE=path, MFI_MS2='%.17g' % (ms * ms), MFI_TON='%.17g' % ton, MFI_GS='%.17g' % gs,
                   MFI_BR='%.17g' % mc.br_of(sstar), MFI_C0='%.17g' % c0)
        if kind == 'tables':
            env['MFI_STAB'] = '1'
        gam = 0.0
    if xis:
        env['MFI_XIS'] = '%.17g' % xis
    return env, gam


def release(tag, K, ms, ton, gs=0.0, sstar=1e16, c0=1.0, bg='fluid', Td=100.0, Tl=1.0e4, bcap=None, save_dense=True, xis=0.0):
    u0 = mc.MU + (mc.U_PREP - mc.MU) * c0
    st = mc.prep(u0, 12.0)
    z = st['z']
    G = mc.U_PREP - mc.MU
    env, gam = spec_env(K, ms, ton, gs, sstar, c0, bg, xis)
    b = max(ms, 1.0) if bcap is None else bcap
    info = dict(tag=tag, spec=dict(K_target=K, m_s=ms, t_on=ton, Gamma_s=gs, s_star_GeV=sstar, c0=c0, bg=bg, gam_bg=gam, xi_s=xis,
                                   ramp='smootherstep m_s^2(t) = m_s^2 S(t/t_on)', step_cap_b=b),
                prep=dict(u=[u0.real, u0.imag], amps=[float(x) for x in np.abs(z)], z=[[c.real, c.imag] for c in z], F=mc.F_of(z),
                          eig=[float(x) for x in st['eig']], grad=st['grad']),
                env=env)
    rd = mc.run(tag + '_dense', 5, b, gam, z, G, Td, int(round(Td / 0.005)), env=env, save=save_dense)
    S = rd['S']
    lateD = S[S[:, 0] >= 0.75 * Td]
    info['dense'] = dict(T=Td, END=rd['END'], seconds=rd['seconds'],
                         late_mean={nm: float(np.mean(lateD[:, C.index(nm)])) for nm in ('NX', 'EX', 'eps', 'NYL', 'NYH', 'Nsig', 'a', 'H')},
                         timeline=[dict(tau=float(r[0]), a=float(r[1]), H=float(r[2]), R=float(r[3]), NX=float(r[4]), EX=float(r[8]),
                                        NYL=float(r[11]), NYH=float(r[10]), Nsig=float(r[12]), absX=float(r[13]), absY=float(r[14]),
                                        c=float(r[28]), ure=float(r[26]), uim=float(r[27]), mL2=float(r[29]))
                                   for r in S[np.unique(np.r_[np.arange(0, min(len(S), 400), 20), np.geomspace(1, len(S) - 1, 60).astype(int)])]])
    rl = mc.run(tag + '_long', 5, b, gam, z, G, Tl, 2000, env=env, cbin=0.5, capT=60.0, save=True)
    S, Cb = rl['S'], rl['C']
    rows = windows(Cb, ['NX', 'NS', 'NY', 'EX', 'NYL', 'NYH', 'Nsig', 'eps'])
    a_end = float(Cb[-1, 4 + C.index('a')])
    fr = {nm: frozen(rows, nm, a_end) for nm in ('NX', 'NS', 'EX', 'NYL', 'NYH', 'Nsig')}
    last = dict(zip(C, S[-1]))
    NX = fr['NX']['value']; EX = fr['EX']['value']
    lx = mc.P05['lx']
    Nc = 2 * math.sqrt(2 * lx) * (EX / (3 * lx)) ** 0.75
    eps = abs(NX) / Nc
    K_end = last['a'] ** 2 * last['H']
    res = dict(N_X=NX, N_X_err=fr['NX']['err'], N_Sigma=fr['NS']['value'], E_X=EX, E_X_err=fr['EX']['err'], eps=eps, E_over_Ec=eps ** (-4 / 3),
               N_YL=fr['NYL']['value'], N_YH=fr['NYH']['value'], N_sigma=fr['Nsig']['value'], K=K_end,
               charge_balance=float(abs(NX + 6 * fr['NS']['value']) / abs(NX)),
               gauss=rl['END']['max_gauss'], H_consistency=rl['END']['max_H_consistency'], hidden_work=rl['END']['max_hidden_work'],
               min_mL2=min(rl['END']['min_mL2'], rd['END']['min_mL2']), a_end=a_end, T_long=Tl,
               normalisation='a = 1 at the slice (end of constant-H inflation), H_i = 1')
    info['long'] = dict(END=rl['END'], seconds=rl['seconds'], frozen=fr, windows=rows[::3])
    info['result'] = res
    OUT.mkdir(parents=True, exist_ok=True)
    info["long"].pop("windows", None)
    info['dense']['timeline'] = [{k: float('%.7g' % v) for k, v in row.items()} for row in info['dense']['timeline']]
    json.dump(info, open(OUT / ("release_%s.json" % tag), "w"), indent=0, default=float)
    print('%s: K=%.4f N_X=%.6g+-%.2g E_X=%.6g eps=%.5f E/Ec=%.3f N_YL=%.4g N_YH=%.4g N_sig=%.3g | gauss %.1e work %.1e Hc %.1e | %.0f+%.0f s' % (
        tag, K_end, NX, fr['NX']['err'], EX, eps, eps ** (-4 / 3), res['N_YL'], res['N_YH'], res['N_sigma'], res['gauss'], res['hidden_work'],
        res['H_consistency'], rd['seconds'], rl['seconds']), flush=True)
    return info


if __name__ == '__main__':
    a = sys.argv[1:]
    tag, K, ms, ton = a[0], float(a[1]), float(a[2]), float(a[3])
    gs = float(a[4]) if len(a) > 4 else 0.0
    sstar = float(a[5]) if len(a) > 5 else 1e16
    c0 = float(a[6]) if len(a) > 6 else 1.0
    bg = a[7] if len(a) > 7 else 'fluid'
    Tl = float(a[8]) if len(a) > 8 else 1.0e4
    release(tag, K, ms, ton, gs, sstar, c0, bg, Tl=Tl)
