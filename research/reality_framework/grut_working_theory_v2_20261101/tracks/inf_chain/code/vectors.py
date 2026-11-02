#!/usr/bin/env python3
"""inf_chain item 2 (Claude, own code): longitudinal hidden vectors on each benchmark's own recorded background, with the
horizon-exit dependence of the benchmark's pre-slice history.

Method [C]:
  1. Background table for t >= 0 (a = 1, H = H_i = 1 at the slice): the dense engine run (dt = 0.005) gives a, H and the
     Higgs norm F^2 = 2(|X|^2 + 9|Y|^2 + 36|Sigma|^2); afterwards
       Gamma_I = 0.5 runs: exact radiation a = (a_d^2 + 2K(t - T_d))^{1/2};
       Gamma_I ~ 0 runs (TB*): own fluid continuation (matter x0 e^{-Gamma t}/a^3 + radiation from the decay);
     F^2 -> F_v^2 + 2 chi^2/a^2 (the X part redshifts as a^-2).
  2. Modes (vmodes.c): exact massless de Sitter Goldstone start at k/aH = 50 with the slice's prepared F_s, through the
     recorded release, read out as the adiabatic number once Omega/H > 100 (checked again at 200) after the release.
     Single-Higgs-direction Proca action, as dq_heat and its skeptic.
  3. Horizon-exit factor (the dq_heat SKEPTIC's frozen-angle method with the benchmark's actual exit-time H and F): the
     frozen Goldstone power is (H/2 pi F)^2 at exit, so each mode with k < 1 is multiplied by
        f(N) = [h_r(N) F_s / F(s(N), h_r(N), eps(N))]^2 ,  k = e^{-N} h_r(N),
     with F re-solved at u = mu + s (10i - mu), R = 6 (2 - eps) h_r^2, s = phi_I/phi_I(cal), h_r = H/H_i from the
     benchmark's background (ch_common.Bench.history); f = 1 for modes that exit after the slice.
  4. N_A = int dlnk k^3 n_k f / 2 pi^2;  deposit per vector m_A <a_A>, <a_A> = (pi K/2 Gamma_A)^{1/2} (decay in radiation),
     Gamma_A = e^2 m_A/48 pi, m_A = e F_v;  R_vec = N_A m_A <a_A>/E_X.
Validation: the M1r reference (MODE 0, (3, 0.5)) against dq_heat's per-mode spectrum (vector_modes.json).
Usage: python3 -B vectors.py KEY [gam_fluid] ; KEY = M1r | HT6t | T1 | HT6t_Send | TB1 | TB2 ...
Output: ../results/vectors_<KEY>[_g<gam>].json
"""
import json, math, os, subprocess, sys
import numpy as np
from scipy.integrate import solve_ivp
import ch_common as ch
from benchmarks_reg import get_bench

fc = ch.fc
VBIN = ch.SP / 'build' / 'vmodes'
ES = tuple(float(x) for x in os.environ.get('CH_ES', '1e-7,1e-6,3e-6,1e-5,1e-4').split(','))
KLO = np.geomspace(1e-4, 1.0, 41)
KHI = np.geomspace(1.0, 30.0, 61)[1:]
KG = np.concatenate([KLO, KHI])


def m1r_dense(Td=100.0):
    """M1r reference background: MODE 0, (m_I, Gamma_I) = (3, 0.5), prep at u = 10i, R = 12 (dq_heat's run)."""
    f = ch.RUNS / 'M1r_dense.npz'
    if f.exists():
        return np.load(f)['S']
    ch.build()
    p = ch.P05
    st = fc.stationary(ch.U_PREP, p, R=12.0)
    z = st['z']
    G, I0, ph0 = fc.trilinear_G(z, p, 3.0)
    g = lambda v: '%.17g' % float(v)
    args = [str(ch.BIN), '0', '3', '0.5', g(p['xi']), g(fc.ETA2), '0', g(p['lx']), g(p['ly']), g(p['ls']), g(p['kap']), g(p['M']),
            g(p['h']), g(p['mu']), '0', '0', g(G.real), g(G.imag), '0', g(Td), '1e-10', '1e-13', str(int(Td / 0.005)), '0', '0', 'F']
    for c in z:
        args += [g(c.real), g(c.imag)]
    args += ['0'] * 6
    r = subprocess.run(args, capture_output=True, text=True)
    S, Cb, END, I0d, W = ch.parse(r.stdout)
    np.savez_compressed(f, S=S)
    json.dump(dict(END=END, I0=I0d, F_s=ch.F_of(z)), open(ch.RUNS / 'M1r_dense_meta.json', 'w'), indent=1, default=float)
    return S


def fluid_late(t0, a0, x0, gam, rho_r0=0.0, tmax=1e10):
    """matter (x0 e^{-gam t}/a^3) + radiation (fed by the decay); 3H^2 = rho (H_i units).  Returns t, ln a, H on a log grid."""
    def rhs(lt, y):
        t = math.exp(lt); la, Rr = y            # Rr = rho_r a^4
        a = math.exp(la)
        rm = x0 * math.exp(-gam * t) / a ** 3
        H = math.sqrt((rm + Rr / a ** 4) / 3)
        return [t * H, t * gam * rm * a ** 4]
    lts = np.linspace(math.log(t0), math.log(tmax), 6000)
    s = solve_ivp(rhs, (lts[0], lts[-1]), [math.log(a0), rho_r0 * a0 ** 4], t_eval=lts, method='LSODA', rtol=1e-11, atol=1e-14)
    t = np.exp(s.t); la = s.y[0]; Rr = s.y[1]
    a = np.exp(la)
    rm = x0 * np.exp(-gam * t) / a ** 3
    H = np.sqrt((rm + Rr / a ** 4) / 3)
    return t, la, H, Rr


def build_table(key, gam_fluid=None, Td=100.0, tmax=1e10):
    """background table for vmodes; returns dict with FS, K, F_v2, chi2, path."""
    if key == 'M1r':
        S = m1r_dense(Td)
        FS = json.load(open(ch.RUNS / 'M1r_dense_meta.json'))['F_s']
        rel = None
    else:
        S = np.load(ch.RUNS / (key + '_dense.npz'))['S']
        rel = json.load(open(ch.RES / ('release_%s.json' % key)))
        FS = rel['prep']['F']
    C = ch.COLS
    t, a, H = S[:, 0], S[:, 1], S[:, 2]
    X, Y, Sv = S[:, C.index('absX')], S[:, C.index('absY')], S[:, C.index('S_over_vev')]
    Sg = Sv * ch.P05['h'] / math.sqrt(2)
    F2 = 2 * (X ** 2 + 9 * Y ** 2 + 36 * Sg ** 2)
    late = t >= 0.75 * Td
    chi2 = float(np.mean(X[late] ** 2 * a[late] ** 2)); Fv2 = float(np.mean(72 * Sg[late] ** 2))
    ad, Hd = float(a[-1]), float(H[-1])
    rows = [(float(tt), math.log(aa), float(hh), float(ff)) for tt, aa, hh, ff in zip(t, a, H, F2)]
    if gam_fluid is None:
        K = ad * ad * Hd
        tl = np.geomspace(Td * 1.0001, tmax, 4000)
        al = np.sqrt(ad ** 2 + 2 * K * (tl - Td))
        for tt, aa in zip(tl, al):
            rows.append((float(tt), math.log(aa), K / aa ** 2, Fv2 + 2 * chi2 / aa ** 2))
        info = dict(K=K, late='radiation')
    else:
        x0 = rel['result']['x0']
        tl, la, Hl, Rr = fluid_late(Td, ad, x0, gam_fluid, tmax=tmax)
        K = float(math.exp(2 * la[-1]) * Hl[-1])
        for tt, l, hh in zip(tl[1:], la[1:], Hl[1:]):
            rows.append((float(tt), float(l), float(hh), Fv2 + 2 * chi2 / math.exp(2 * l)))
        info = dict(K=K, late='fluid', x0=x0, gam=gam_fluid, Rr_end=float(Rr[-1]),
                    K_formula=0.50128 * x0 ** (2 / 3) * gam_fluid ** (-1 / 3))
    tag = key + ('' if gam_fluid is None else '_g%.3e' % gam_fluid)
    path = ch.SP / ('table_%s.txt' % tag)
    with open(path, 'w') as fh:
        fh.write('%d\n' % len(rows))
        for r in rows:
            fh.write('%.17g %.17g %.17g %.17g\n' % r)
    info.update(FS=FS, Fv2=Fv2, Fv=math.sqrt(Fv2), chi2=chi2, path=str(path), tag=tag, a_d=ad, H_d=Hd, Td=Td)
    return info


def run_modes(info, ks, e, TR=None, AD=100.0, nproc=2):
    TR = info['Td'] if TR is None else TR
    chunks = [ks[i::nproc] for i in range(nproc)]
    procs = []
    for c in chunks:
        args = [str(VBIN), info['path'], '%.17g' % info['FS'], '%.17g' % TR, '%.17g' % AD, 'D', '0']
        for k in c:
            args += ['%.17g' % k, '%.6g' % e]
        procs.append(subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True))
    out = {}
    for pr in procs:
        o, er = pr.communicate()
        if pr.returncode != 0:
            raise RuntimeError(er[:300])
        for line in o.split('\n'):
            if line.startswith('M '):
                v = line.split()
                out[float(v[1])] = dict(n1=float(v[3]), n2=float(v[4]), t1=float(v[5]), t2=float(v[6]), steps=int(v[7]))
    return [out[float('%.10g' % k)] if float('%.10g' % k) in out else out[min(out, key=lambda x: abs(x - k))] for k in ks]


def exit_factor(b, Ngrid=np.linspace(0.0, 10.0, 81)):
    """f(N) along the benchmark's pre-slice path (see docstring); F_s = F of the slice preparation."""
    P = b.history(Ngrid)
    Fs = ch.F_of(b.prep()['z'])
    fN, info = [], []
    guess = None
    for (s, hr, eps), N in zip(P, Ngrid):
        u = ch.MU + (ch.U_PREP - ch.MU) * s
        R = 6.0 * (2.0 - eps) * hr * hr
        st = None
        for gg in (guess, None):
            try:
                st = fc.stationary(u, ch.P05, R=R, guess=gg)
                if st['grad'] < 1e-8:
                    break
            except Exception:
                st = None
        guess = st['z']
        F = ch.F_of(st['z'])
        fN.append((hr * Fs / F) ** 2)
        info.append(dict(N=float(N), s=float(s), h_r=float(hr), eps=float(eps), F=F, f=fN[-1],
                         m2min_over_H2=float(sorted(st['eig'])[1] / hr ** 2), amps=[float(x) for x in np.abs(st['z'])]))
    hN = np.array([x[1] for x in P])
    return np.array(Ngrid), np.array(fN), hN, info


def f_of_k_factory(Ng, fN, hN):
    lnk_of_N = -Ng + np.log(hN)

    def f_of_k(k):
        if k >= hN[0]:
            return 1.0
        if math.log(k) < lnk_of_N[-1]:
            return float(fN[-1])
        N = float(np.interp(math.log(k), lnk_of_N[::-1], Ng[::-1]))
        return float(np.interp(N, Ng, fN))
    return f_of_k


def assemble(ks, modes, e, info, f_of_k, EX):
    n = np.array([m['n2'] for m in modes]); n1 = np.array([m['n1'] for m in modes])
    fk = np.array([f_of_k(k) for k in ks])
    dN = ks ** 3 * n / (2 * math.pi ** 2)
    lk = np.log(ks)
    lo = ks <= 1.0
    NA_lo_static = float(np.trapezoid(dN[lo], lk[lo]))
    NA_lo = float(np.trapezoid((dN * fk)[lo], lk[lo]))
    hi = ks >= 1.0
    NA_hi = float(np.trapezoid(dN[hi], lk[hi])) if hi.sum() > 1 else 0.0
    mA = e * info['Fv']; GA = e * e * mA / (48 * math.pi)
    aA = math.sqrt(math.pi * info['K'] / (2 * GA))
    dep = mA * aA
    NA = NA_lo + NA_hi
    return dict(e=e, N_A=NA, N_A_k_le_1=NA_lo, N_A_k_le_1_static=NA_lo_static, f_eff=NA_lo / NA_lo_static, N_A_k_ge_1=NA_hi,
                m_A=mA, Gamma_A=GA, a_mean=aA, deposit_per_vector=dep, R_vec=NA * dep / EX,
                max_readout_change=float(np.max(np.abs(n - n1) / np.maximum(np.abs(n), 1e-300))),
                spectrum=[dict(k=float(k), n_k=float(x), f=float(ff), dN_dlnk=float(d)) for k, x, ff, d in zip(ks, n, fk, dN)])


def main(key, gam_fluid=None):
    info = build_table(key, gam_fluid)
    if key == 'M1r':
        EX = 188706.29275778163
        f_of_k = lambda k: 1.0
        exitinfo = None
    else:
        rel = json.load(open(ch.RES / ('release_%s.json' % key)))
        EX = rel['result']['E_X']
        b = get_bench(key)
        Ng, fN, hN, exitinfo = exit_factor(b)
        f_of_k = f_of_k_factory(Ng, fN, hN)
    out = dict(method=__doc__, key=key, background={k: v for k, v in info.items() if k != 'path'}, E_X=EX, exit_path=exitinfo, by_e=[])
    for e in ES:
        modes = run_modes(info, list(KG), e)
        r = assemble(KG, modes, e, info, f_of_k, EX)
        r['steps_total'] = int(sum(m['steps'] for m in modes))
        out['by_e'].append(r)
        print('%s e=%.0e: N_A=%.5g (k<=1 %.5g, static %.5g, f_eff %.4f; k>=1 %.4g) dep=%.4g R_vec=%.5g K=%.5g Fv=%.3f readout chg %.1e' % (
            info['tag'], e, r['N_A'], r['N_A_k_le_1'], r['N_A_k_le_1_static'], r['f_eff'], r['N_A_k_ge_1'], r['deposit_per_vector'],
            r['R_vec'], info['K'], info['Fv'], r['max_readout_change']), flush=True)
    if key == 'M1r':
        vm = json.load(open(ch.TRK / 'dq_heat' / 'results' / 'vector_modes.json'))
        cmp = []
        for c in vm['cases']:
            if c['variant'] != 'full':
                continue
            mine = [x for x in out['by_e'] if abs(x['e'] - c['e']) < 1e-12 * c['e']]
            if not mine:
                continue
            mine = mine[0]
            ks = np.array([s['k'] for s in mine['spectrum']]); nk = np.array([s['n_k'] for s in mine['spectrum']])
            rat = []
            for s in c['spectrum']:
                if 1e-4 <= s['k'] <= 30:
                    rat.append(dict(k=s['k'], ratio=float(np.exp(np.interp(math.log(s['k']), np.log(ks), np.log(nk)))) / s['n_k']))
            cmp.append(dict(e=c['e'], N_A_k_le_1_dqheat=c['N_A_k_le_1'], N_A_k_le_1_mine=mine['N_A_k_le_1'],
                            ratio_N_A_k_le_1=mine['N_A_k_le_1'] / c['N_A_k_le_1'], N_A_dqheat=c['N_A'], N_A_mine=mine['N_A'],
                            per_mode_ratio=rat))
            rr = [x['ratio'] for x in rat]
            print('  vs dq_heat e=%.0e: N_A(k<=1) %.5f / %.5f = %.4f ; N_A %.5f / %.5f ; per-mode ratio [%.3f, %.3f] (k<=0.3: [%.3f, %.3f])' % (
                c['e'], mine['N_A_k_le_1'], c['N_A_k_le_1'], mine['N_A_k_le_1'] / c['N_A_k_le_1'], mine['N_A'], c['N_A'], min(rr), max(rr),
                min(x['ratio'] for x in rat if x['k'] <= 0.3), max(x['ratio'] for x in rat if x['k'] <= 0.3)), flush=True)
        out['validation_vs_dq_heat'] = cmp
    tag = info['tag']
    json.dump(out, open(ch.RES / ('vectors_%s.json' % tag), 'w'), indent=1, default=float)


if __name__ == '__main__':
    k = sys.argv[1]
    g = float(sys.argv[2]) if len(sys.argv) > 2 else None
    main(k, g)
