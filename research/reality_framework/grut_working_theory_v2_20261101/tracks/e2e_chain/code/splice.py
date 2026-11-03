#!/usr/bin/env python3
"""e2e_chain item 1 (Claude): the final spliced histories of HFa / HFb in the bench_*.json format of mf_ideal's
bench_run.py.

PIECES (all read-only):
  * t < 0: the mf_hybrid bench rows (constant-H inflation; exit-factor check only), tracks/mf_hybrid/results/bench_<name>.json.
  * 0 <= t <= t_sw: e2e_lattice's 3D coherent stage (box means), results/hist_<name>.json, closure key (P0 = local Pauli
    blocking, P0GM = the mf_hybrid skeptic's saturation on the box-mean condensate; H84_* at 8.4e10 GeV).  Columns used:
    a, H, R_cyc (cycle-mean Ricci scalar; the lattice samples every 0.0042 H^-1, i.e. 3 samples per 2 m_chi period, so
    the raw R is aliased) and m_s^2 = running mean (width 0.025 = two chi^2 periods) of the box-mean m_s2_mean.
  * t > t_sw: e2e_kinetic's dephased kinetic records <prefix>_dep_rec.txt (its scratch; t >= 3), converted with the
    column definitions of its export_fluid (cycle-mean R = (f_sig + f_chi) rho/M_P^2; m_s: HFa g_s sqrt(v^2 + rho_chi/m_chi^2),
    HFb sqrt(9 - 3R/4)), thinned to <= 1500 rows log-spaced in t.

SPLICE RULE (model choice [A]): t_sw = 3 (the lattice's end, inside both tracks' advised windows).  For t > t_sw the
kinetic universe is used with its scale factor multiplied by lambda = a_lat(t_sw)/a_kin(t_sw) and its own H(t) (so
K_final = lambda^2 K_kin).  The H mismatch at t_sw (|H_lat/H_kin - 1|, about 2%) and the composition mismatch
(f_sig, f_chi, f_rad on both sides) are recorded; the composition-weighted correction to K is estimated [E] in
diagnostics (K ~ (comoving chi energy)^(2/3) for a long-lived matter component).
Spectator: s/s_* integrated on the spliced table (check column only; the engine evolves its own s).
Rescaling rule for H != design: mf_chain hist_build's (V0 free; Gamma_s ~ H_i^2/H_ref).  The lattice block is the one
computed nearest in H (design, or 8.4e10); the late kinetic history is the one computed at that H [A].
Usage: python3 -B splice.py NAME TAG --lat P0GM --kin fx_natw_nq6400 --H 1.3e11 [--tsw 3]
Labels: [C] arithmetic on [C] inputs; the splice and every choice in it are model choices.
"""
import argparse, json, math
import numpy as np
from scipy.integrate import solve_ivp
import e2c_common as E

REC_COLS = ['t', 'a', 'H', 'f_sig', 'f_chi', 'f_N', 'f_hard', 'f_th', 'T_th', 'Tcal2', 'T_SM', 'tau_KL']
for _j in range(4):
    REC_COLS += ['fN%d' % (_j + 1), 'Gband%d' % (_j + 1), 'Gall%d' % (_j + 1), 'nband%d' % (_j + 1), 'A%d' % (_j + 1), 'nN%d' % (_j + 1)]
REC_COLS += ['n_h', 'rho', 'rho_a4']


def lattice_block(name, key):
    d = json.load(open(E.TRK / 'e2e_lattice' / 'results' / ('hist_%s.json' % name)))
    if key == d['primary']:
        L, head, cfg = d['lattice'], d['head'], d['config']
    else:
        A = d['alternates'][key]
        L, head, cfg = A['lattice'], A['head'], A['config']
    g = lambda k: np.array([np.nan if x is None else float(x) for x in L[k]])
    return dict(t=g('t'), a=g('a'), H=g('H'), R=g('R'), R_cyc=g('R_cyc'), ms2=g('m_s2_mean'), ms2_std=g('m_s2_std'),
                f_sig=g('f_sig'), f_chi=g('f_chi'), f_grad=g('f_grad'), f_rad=g('f_rad'), f_SM=g('f_SM'), s_std=g('s_std'),
                eps=g('eps_H'), K=g('K_running'), head=head, cfg=cfg)


def kin_block(prefix, name, H_GeV):
    if str(prefix).endswith('.json'):                    # a kinetic export (bench format), e.g. our moments hand-over
        h = json.load(open(prefix))['history']
        g = lambda k: np.array([np.nan if x is None else float(x) for x in h[k]])
        return dict(t=g('t'), a=g('a'), H=g('H'), R=g('R'), ms=g('m_s'), f_sig=g('f_sig'), f_chi=g('f_chi'), f_rad=g('f_rad'),
                    f_SM=g('f_SM'), f_N=g('f_N'), T_SM=g('T_SM'), rec=None)
    a = np.loadtxt(str(E.KIN_SP / ('%s_dep_rec.txt' % prefix)), ndmin=2)
    r = {c: a[:, i] for i, c in enumerate(REC_COLS)}
    o = np.argsort(r['t'], kind='stable')
    r = {k: v[o] for k, v in r.items()}
    t = r['t']
    tt = np.unique(np.r_[t[0], np.geomspace(t[0], t[-1], 1400), t[-1]])
    idx = np.unique(np.clip(np.searchsorted(t, tt), 0, len(t) - 1))
    r = {k: v[idx] for k, v in r.items()}
    pr = json.load(open(E.TRK / 'mf_hybrid' / 'results' / ('bench_%s.json' % name)))['params']
    mP = E.MP_GEV / H_GeV
    rho = 3 * mP ** 2 * r['H'] ** 2
    R = (r['f_sig'] + r['f_chi']) * rho / mP ** 2
    if pr['smech'] == 'a':
        v = pr['v_MP'] * mP
        ms = (3.0 / v) * np.sqrt(v * v + r['f_chi'] * rho / pr['m_chi_over_Hi'] ** 2)
    else:
        ms = np.sqrt(np.maximum(9.0 - 0.75 * R, 0.0))
    return dict(t=t[idx] if False else r['t'], a=r['a'], H=r['H'], R=R, ms=ms, f_sig=r['f_sig'], f_chi=r['f_chi'],
                f_rad=r['f_N'] + r['f_hard'] + r['f_th'], f_SM=r['f_hard'] + r['f_th'], f_N=r['f_N'], T_SM=r['T_SM'],
                rec=r)


def running_mean(t, y, w):
    out = np.empty_like(y)
    for i in range(len(t)):
        s = (t >= t[i] - w / 2) & (t <= t[i] + w / 2)
        out[i] = np.mean(y[s])
    return out


def build(name, tag, lat='P0GM', kin='fx_natw_nq6400', H_GeV=E.H_REF, t_sw=3.0, write=True, ms2spread=0.0, ms2scale=1.0):
    bench = json.load(open(E.TRK / 'mf_hybrid' / 'results' / ('bench_%s.json' % name)))
    L = lattice_block(name, lat)
    Kb = kin_block(kin, name, H_GeV)
    # ---- lattice part 0 <= t <= t_sw
    m = (L['t'] >= 0) & (L['t'] <= t_sw + 1e-9)
    tl, al, Hl, Rl = L['t'][m], L['a'][m], L['H'][m], L['R_cyc'][m]
    ms2l = running_mean(L['t'], L['ms2'], 0.025)[m]
    ms2l = ms2l * ms2scale                              # uniform rescaling of the lattice-stage m_s^2 (sensitivity check)
    if ms2spread:                                       # separate-universe patch at ms2spread x the lattice's spatial std of m_s^2
        ms2l = ms2l + ms2spread * running_mean(L['t'], L['ms2_std'], 0.025)[m]
    # ---- kinetic part t > t_sw
    mk = Kb['t'] > t_sw + 1e-9
    j0 = int(np.argmin(np.abs(Kb['t'] - t_sw)))
    a_kin_sw, H_kin_sw = float(np.interp(t_sw, Kb['t'], Kb['a'])), float(np.interp(t_sw, Kb['t'], Kb['H']))
    a_lat_sw, H_lat_sw = float(np.interp(t_sw, L['t'], L['a'])), float(np.interp(t_sw, L['t'], L['H']))
    lam = a_lat_sw / a_kin_sw
    tk, ak, Hk, Rk, msk = Kb['t'][mk], Kb['a'][mk] * lam, Kb['H'][mk], Kb['R'][mk], Kb['ms'][mk]
    t = np.r_[tl, tk]; a = np.r_[al, ak]; H = np.r_[Hl, Hk]; R = np.r_[Rl, Rk]; ms2 = np.r_[ms2l, msk ** 2]
    ms = np.sqrt(np.maximum(ms2, 0.0))
    N = np.log(a)
    fsig = np.r_[L['f_sig'][m], Kb['f_sig'][mk]]; fchi = np.r_[L['f_chi'][m] + L['f_grad'][m], Kb['f_chi'][mk]]
    frad = np.r_[L['f_rad'][m], Kb['f_rad'][mk]]; fsm = np.r_[L['f_SM'][m], Kb['f_SM'][mk]]
    eps = -np.gradient(H, t) / H ** 2
    # spectator (check column; exact to t = 60, then envelope ~ a^-3/2)
    Hi_GeV = H_GeV * float(H[0]) if abs(H[0] - 1) < 1e-3 else H_GeV
    fH = lambda q: np.interp(q, t, H); fm = lambda q: np.interp(q, t, ms2)
    tmax = min(60.0, t[-1])
    ss = solve_ivp(lambda q, y: [y[1], -3 * fH(q) * y[1] - fm(q) * y[0]], (0, tmax), [1.0, 0.0], method='DOP853', rtol=1e-9,
                   atol=1e-12, dense_output=True, max_step=0.01)
    c = np.where(t <= tmax, ss.sol(np.minimum(t, tmax))[0], np.nan)
    cd = np.where(t <= tmax, ss.sol(np.minimum(t, tmax))[1], np.nan)
    j60 = np.searchsorted(t, tmax) - 1
    amp60 = math.hypot(c[j60], cd[j60] / ms[j60])
    env = amp60 * (a[j60] / a) ** 1.5
    s_out = np.where(t <= tmax, c, env)
    env_flag = (t > tmax).astype(float)
    mP = E.MP_GEV / H_GeV
    sst = bench['params']['s_star_MP'] * mP
    rho_tot = 3 * mP ** 2 * H ** 2
    E_s = np.where(t <= tmax, 0.5 * cd ** 2 + 0.5 * ms2 * c ** 2, 0.5 * ms2 * env ** 2) * sst ** 2
    f_s = E_s / rho_tot
    # ---- assemble (pre-slice bench rows first)
    h = bench['history']
    tb = np.array(h['t']); pre = tb < 0
    new = dict(t=t, N=N, H=H, eps_H=eps, R=R, sigma_MP=np.full_like(t, np.nan), chi_MP=np.full_like(t, np.nan), s_over_sstar=s_out,
               m_s=ms, f_sig=fsig, f_chi=fchi, f_rad=frad, f_SM=fsm, f_s=f_s, a=a, K_running=a * a * H,
               R_osc_sig=np.where(t > 0, 0.0, np.nan), R_osc_chi=np.where(t > 0, 0.0, np.nan), s_is_envelope=env_flag)
    hist = {}
    for k in new:
        old = [h[k][i] for i in range(len(tb)) if pre[i]] if k in h else [None] * int(pre.sum())
        hist[k] = old + [None if not np.isfinite(x) else float('%.10g' % x) for x in new[k]]
    K = float(a[-1] ** 2 * H[-1])
    # ---- diagnostics of the splice
    ia = int(np.argmin(np.abs(L['t'] - t_sw)))
    lat_sw = dict(t=float(L['t'][ia]), a=a_lat_sw, H=H_lat_sw, f_sig=float(L['f_sig'][ia]), f_chi=float(L['f_chi'][ia]),
                  f_grad=float(L['f_grad'][ia]), f_rad=float(L['f_rad'][ia]), f_SM=float(L['f_SM'][ia]), R_cyc=float(L['R_cyc'][ia]),
                  m_s2=float(ms2l[-1]))
    kin_sw = dict(t=float(Kb['t'][j0]), a=a_kin_sw, H=H_kin_sw, f_sig=float(Kb['f_sig'][j0]), f_chi=float(Kb['f_chi'][j0]),
                  f_rad=float(Kb['f_rad'][j0]), f_SM=float(Kb['f_SM'][j0]), R=float(Kb['R'][j0]), m_s2=float(Kb['ms'][j0] ** 2))
    K_kin = float(Kb['a'][-1] ** 2 * Kb['H'][-1])
    # comoving chi(+gradient, i.e. chi quanta) energy ratio lattice / kinetic at t_sw
    chi_ratio = ((lat_sw['f_chi'] + lat_sw['f_grad']) * H_lat_sw ** 2 * a_lat_sw ** 3) / (kin_sw['f_chi'] * H_kin_sw ** 2 * a_kin_sw ** 3)
    # Friedmann consistency of the table: d ln a / dt vs H, outside the joint
    dlna = np.gradient(np.log(a), t)
    okm = (t > 0.05) & (np.abs(t - t_sw) > 0.05)
    fr_err = float(np.nanmax(np.abs(dlna[okm] / H[okm] - 1))) if okm.any() else None
    diag = dict(t_sw=t_sw, lattice=lat_sw, kinetic=kin_sw, lam=lam, H_jump=H_kin_sw / H_lat_sw - 1, K_kin=K_kin, K_spliced=K,
                K_composition_estimate_E=K * chi_ratio ** (2.0 / 3.0), chi_comoving_ratio=chi_ratio,
                lattice_budget=dict(end=L['head'].get('budget_end'), absmax=L['head'].get('budget_absmax')),
                friedmann_max_rel_err_off_joint=fr_err,
                K_drift_last_decade=float(a[-1] ** 2 * H[-1] / np.interp(t[-1] / 10, t, a * a * H) - 1),
                f_s_at_60=float(np.interp(60.0, t, f_s)), f_s_over_a_late=float(f_s[-50] / a[-50]))
    pr = dict(bench['params'])
    pr.update(H_i_GeV=Hi_GeV, design_H_GeV=H_GeV, s_star_GeV=bench['params']['s_star_MP'] * E.MP_GEV)
    out = dict(name=tag, label='e2e_chain splice: e2e_lattice (%s) to t_sw = %g, e2e_kinetic %s after; model choices throughout' % (lat, t_sw, kin),
               slice='a = 1, H = H_i at sigma = sigma_c', H_i_GeV=Hi_GeV, params=pr,
               spectator=dict(Gamma_s_assumed=1.3e-6 * (H_GeV / E.H_REF)),
               reheating=dict(K=K, K_kin=K_kin, lattice=lat, kinetic=kin, t_sw=t_sw),
               splice=diag, history_columns=bench['history_columns'], history=hist)
    if write:
        path = E.SP / 'hist' / ('bench_%s.json' % tag)
        E.dump(path, out)
        out['path'] = str(path)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('name'); ap.add_argument('tag')
    ap.add_argument('--lat', default='P0GM'); ap.add_argument('--kin', default='fx_natw_nq6400')
    ap.add_argument('--H', type=float, default=E.H_REF); ap.add_argument('--tsw', type=float, default=3.0)
    ap.add_argument('--ms2spread', type=float, default=0.0); ap.add_argument('--ms2scale', type=float, default=1.0)
    a = ap.parse_args()
    o = build(a.name, a.tag, a.lat, a.kin, a.H, a.tsw, ms2spread=a.ms2spread, ms2scale=a.ms2scale)
    print(json.dumps(dict(path=o['path'], H_i_GeV=o['H_i_GeV'], splice=o['splice']), indent=1, default=float))
