#!/usr/bin/env python3
"""mf_ideal (Claude, own code): run the M1r chain (release -> vectors -> S/Q, e_crit, relic, H_i) on a benchmark history
file bench_<name>.json (the mf_hybrid export format; see README.md for the fields used).

Steps
  1. Read history.{t, N, H, R, s_over_sstar, m_s, s_is_envelope, R_osc_sig, R_osc_chi} and params.{H_i_GeV, s_star_GeV}
     (spectator.Gamma_s_assumed if present).  Rows with t >= 0 form the MODE 5 table (a = 1, H = 1 at t = 0 is checked).
  2. R: --R raw (default) uses the bench's R as given (valid when its oscillation is resolved, as in bench_HF*);
     --R trace: in the homogeneous part (R_osc_* null) R is replaced by its cycle mean
     6 Hbar^2 (2 - epsbar), Hbar a running mean of H over --win (H_i^-1), epsbar = -dHbar/dt/Hbar^2; in the fluid part the
     bench's R is already the cycle mean and is used as given (use it when the sampling does not resolve 2 m_chi).
  3. s: --s evolve (default): the engine evolves c = s/s_* from c(0) = 1, c'(0) = 0 with m_s^2(t) from the m_s column
     (MFI_SM2TAB), Gamma_s from the bench (or --gs), back-reaction (H_i/s_*)^2; the bench's own s(t) is compared as a check.
     --s table: c(t) from s_over_sstar (exact part only); the table is then cut where the envelope starts.
  4. Pre-slice exit factor: f(N) = [h_r F_s/F(u_slice, R = 6(2 - eps) h_r^2)]^2 from the bench's t < 0 rows (s frozen),
     reported at 3-4 e-folds before the slice (the vectors themselves use f = 1, exact for constant H).
  5. release_mfi.release (dense + long), vectors_mfi (own background), chain_mfi.analyse.
Usage: python3 -B bench_run.py BENCH.json [--tag NAME] [--R trace|raw] [--win 0.05] [--s evolve|table] [--gs G]
                                          [--tlong 1e4] [--no-chain]
Output: ../results/release/release_<tag>.json, ../results/vectors/vectors_<tag>.json, ../results/chain/chain_<tag>.json,
        ../results/bench/bench_run_<tag>.json (interface diagnostics).  Every choice is a model choice.
"""
import argparse, json, math, os
import numpy as np
import mfi_common as mc


def col(h, k):
    return np.array([np.nan if v is None else float(v) for v in h[k]])


def running_mean(t, y, w):
    out = np.empty_like(y)
    for i in range(len(t)):
        s = (t >= t[i] - w / 2) & (t <= t[i] + w / 2)
        out[i] = np.mean(y[s])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('bench'); ap.add_argument('--tag'); ap.add_argument('--R', default='raw'); ap.add_argument('--win', type=float, default=0.05)
    ap.add_argument('--s', default='evolve'); ap.add_argument('--gs', type=float); ap.add_argument('--tlong', type=float, default=1.0e4)
    ap.add_argument('--no-chain', action='store_true')
    a = ap.parse_args()
    d = json.load(open(a.bench))
    name = d.get('name', os.path.basename(a.bench)[6:-5])
    tag = a.tag or ('bench_%s' % name)
    h = d['history']
    t, N, H, R = col(h, 't'), col(h, 'N'), col(h, 'H'), col(h, 'R')
    c, ms = col(h, 's_over_sstar'), col(h, 'm_s')
    env_flag = col(h, 's_is_envelope') if 's_is_envelope' in h else np.zeros_like(t)
    hom = np.isnan(col(h, 'R_osc_sig')) if 'R_osc_sig' in h else np.ones_like(t, bool)
    Hi_GeV = float(d['params'].get('H_i_GeV', d.get('H_i_GeV', mc.H_M1)))
    sst = float(d['params'].get('s_star_GeV', 1e16))
    gs = a.gs if a.gs is not None else float(d.get('spectator', {}).get('Gamma_s_assumed', 0.0) or 0.0)
    diag = dict(bench=os.path.abspath(a.bench), tag=tag, H_i_GeV=Hi_GeV, s_star_GeV=sst, Gamma_s=gs, R_mode=a.R, s_mode=a.s)
    # slice checks
    i0 = int(np.argmin(np.abs(t)))
    diag['slice_row'] = dict(t=float(t[i0]), N=float(N[i0]), H=float(H[i0]), R=float(R[i0]), c=float(c[i0]))
    post = t >= 0
    tp, Np, Hp, Rp, cp, msp, ep, hp = t[post], N[post], H[post], R[post], c[post], ms[post], env_flag[post], hom[post]
    if tp[0] > 1e-12:                                         # insert an exact t = 0 row (interpolated; c = 1 frozen)
        ins = lambda arr, v: np.r_[v, arr]
        tp, Np, Hp, Rp = ins(tp, 0.0), ins(Np, np.interp(0.0, t, N)), ins(Hp, np.interp(0.0, t, H)), ins(Rp, np.interp(0.0, t, R))
        cp, msp, ep, hp = ins(cp, 1.0), ins(msp, np.interp(0.0, t, ms)), ins(ep, 0.0), ins(hp, True)
    Np = Np - np.interp(0.0, t, N)                            # a = 1 at the slice
    dts = np.diff(tp[hp]) if hp.sum() > 2 else np.array([np.nan])
    diag['homogeneous_part'] = dict(t_end=float(tp[hp][-1]) if hp.any() else None, median_dt=float(np.nanmedian(dts)))
    if a.R == 'trace':
        Hb = running_mean(tp, Hp, a.win)
        eb = -np.gradient(Hb, tp) / Hb ** 2
        Rt = 6 * Hb ** 2 * (2 - eb)
        Rp = np.where(hp, Rt, Rp)
        Rp[0] = R[i0] if abs(t[i0]) < 1e-12 else Rp[0]
        diag['R_trace'] = dict(window=a.win, R_at_0=float(Rp[0]), R_mean_hom=float(np.mean(Rp[hp])) if hp.any() else None)
    # spectator
    sel = np.ones_like(tp, bool)
    if a.s == 'table':
        sel = ep < 0.5
        diag['s_table_cut_at'] = float(tp[sel][-1])
    tp, Np, Hp, Rp, cp, msp = tp[sel], Np[sel], Hp[sel], Rp[sel], cp[sel], msp[sel]
    cd = np.gradient(cp, tp)
    path = mc.SP / ('benchtab_%s.txt' % tag)
    with open(path, 'w') as fh:
        fh.write('%d 7\n' % len(tp))
        for row in zip(tp, Np, Hp, Rp, cp, cd, msp ** 2):
            fh.write(' '.join('%.17g' % x for x in row) + '\n')
    diag['table'] = dict(path=str(path), rows=int(len(tp)), t_end=float(tp[-1]), K_end=float(np.exp(2 * Np[-1]) * Hp[-1]))
    # pre-slice exit factor (s frozen -> u = u_slice; F re-solved at R = 6(2 - eps) h_r^2)
    pre = t < 0
    if pre.sum() > 3:
        st0 = mc.prep(mc.U_PREP, 12.0); Fs = mc.F_of(st0['z'])
        Nb = np.interp(0.0, t, N) - N[pre]; Hb = H[pre]; eb = col(h, 'eps_H')[pre] if 'eps_H' in h else np.zeros(pre.sum())
        fx = []
        for nb in (1.0, 2.0, 3.0, 4.0):
            if nb > Nb.max():
                continue
            j = int(np.argmin(np.abs(Nb - nb)))
            hr = Hb[j] / np.interp(0.0, t, H)
            st = mc.prep(mc.U_PREP, 6 * (2 - eb[j]) * hr * hr, guess=st0['z'])
            fx.append(dict(N_before=float(Nb[j]), h_r=float(hr), eps=float(eb[j]), f=float((hr * Fs / mc.F_of(st['z'])) ** 2)))
        diag['exit_factor'] = fx
    from release_mfi import release
    env_extra = dict(MFI_SM2TAB='1') if a.s == 'evolve' else dict(MFI_STAB='1')
    # release with the table background; spec_env builds MFI_TABLE etc.; extra keys through a thin wrapper
    import release_mfi as rm
    orig = rm.spec_env

    def spec_env(*args, **kw):
        env, gam = orig(*args, **kw)
        env.update(env_extra)
        return env, gam
    rm.spec_env = spec_env
    bcap = max(4.0, float(np.nanmax(msp)))
    info = release(tag, diag['table']['K_end'], float(msp[-1]), 0.0, gs=gs, sstar=sst, bg='table:%s' % path, Tl=a.tlong, bcap=bcap)
    rm.spec_env = orig
    # compare the evolved s with the bench's exact s
    S = np.load(mc.RUNS / (tag + '_dense.npz'))['S']
    ts, cs = S[:, 0], S[:, mc.COLS.index('phihat')]
    ex = (ep < 0.5)[sel] if a.s == 'evolve' else np.ones_like(tp, bool)
    m = ex & (tp <= ts[-1])
    diag['s_check'] = dict(max_abs_diff_vs_bench=float(np.max(np.abs(np.interp(tp[m], ts, cs) - cp[m]))) if m.any() else None,
                           t_range=[float(tp[m][0]), float(tp[m][-1])] if m.any() else None)
    os.makedirs(mc.RES / 'bench', exist_ok=True)
    if not a.no_chain:
        import vectors_mfi, chain_mfi
        vectors_mfi.main(tag)
        o = chain_mfi.analyse(tag)
        json.dump(o, open(mc.RES / 'chain' / ('chain_%s.json' % tag), 'w'), indent=1, default=float)
        diag['chain'] = {k: o[k] for k in ('K', 'eps', 'E_over_Ec', 'N_YL', 'N_YH', 'S_per_Q_e_inf', 'e_crit_bulk', 'e_crit_Q2e5')}
        diag['chain']['e1e-6'] = o['e1e-6']; diag['chain']['e3e-6'] = o['e3e-6']
    json.dump(diag, open(mc.RES / 'bench' / ('bench_run_%s.json' % tag), 'w'), indent=1, default=float)
    print(json.dumps({k: v for k, v in diag.items() if k != 'chain'}, indent=1, default=float)[:3000])
    if 'chain' in diag:
        c6 = diag['chain']['e1e-6']
        print('%s: K=%.4f eps=%.4f E/Ec=%.3f N_YL=%.4g | e=1e-6: R_vec=%.3g R_Y=%.3g S/Q=%.4f %s | e_crit %s / %s' % (
            tag, diag['chain']['K'], diag['chain']['eps'], diag['chain']['E_over_Ec'], diag['chain']['N_YL'], c6['R_vec'], c6['R_Y'],
            c6['S_per_Q'], c6['branch'], diag['chain']['e_crit_bulk'], diag['chain']['e_crit_Q2e5']))


if __name__ == '__main__':
    main()
