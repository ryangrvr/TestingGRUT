#!/usr/bin/env python3
"""inf_potential (Claude, own code): calibrate the benchmark and reference potentials, and tabulate for each
  - the CMB predictions at the pivot k = 0.05/Mpc (exact Mukhanov-Sasaki): A_s, n_s, alpha_s, r, n_t, H_*;
  - N_* for the engine-like history (constant Gamma_I = 0.5 H_i) and for (w, T_RH) histories;
  - the last 0-15 e-folds: H(N), phi(N), holder field phi_I(N) = phi - phi_min, eps_H(N), eta_V(N), and the holder ratio
    s(N) = phi_I(N)/phi_I(slice) with h_r(N) = H(N)/H(slice) for the slices S1 (1 e-fold before eps_H = 1, PRIMARY)
    and Send (eps_H = 1);
  - the end of inflation and the first oscillations (Gamma = 0 and Gamma = 0.5 H_i): m_I/H_i, m_I/H_end, turning points,
    frequency, <w>, a^2 H;
  - the holder coupling g_I re-fixed so that mu_eff = 10i H_i at the slice (|10i - mu| = 12.2987 H_i);
  - the open inflaton -> Y Y channels through the holder itself (thresholds 2 m_L = 44.50, m_L + m_H = 58.35, 2 m_H =
    72.20 H_i at the late vacuum) and the hidden heat they would deposit [E].
Usage: python3 benchmarks.py <name> [<name> ...]   names: T1 E1 HT6t HT6 HT4 HT8 PL2 (each at H_i = 1.29e11 and 1.50e11)
Output: ../results/bench_<name>.json.   Labels: [C] given the model choices; the hidden-heat estimate is [E].
"""
import json, math, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parent))
import inf_common as ic
import calibrate as cal

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
HI_LIST = (1.29e11, 1.50e11)
U_GAP = abs(10j - 7.159456159)          # |10i - mu| in H_i  [Q fx_holder: g_I phi_I(0) = 12.2987 H_i]
ALPHA_I = math.atan2(10.0, -7.159456159)  # arg(10i - mu) = 2.192144074 [Q]
M_L, M_H = 22.249, math.sqrt(900 + 2 * 7.159456159 * 40 / math.sqrt(2))   # light/heavy Y at the late vacuum [Q dq_heat]
SIG = 40 / math.sqrt(2)                 # |Sigma| = h/sqrt2 in H_i
DEPOSIT_YL = 0.6565                     # conformal deposit per Y_L in units of the M1r E_X (K = 1.986) [Q dq_heat]
EX_M1R = 188706.29                      # M1r X release energy (H_i^4, a = 1 at the M1r slice) [Q dq_heat]
NEB = [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 2.7, 3.0, 3.5, 3.8, 4.0, 4.5,
       4.8, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 15.0]


def calib(name, Hi, hist=('gamma', 0.5)):
    if name == 'T1':
        return cal.calibrate_alpha('T', Hi, hist)
    if name == 'E1':
        return cal.calibrate_alpha('E', Hi, hist)
    if name == 'HT6t':
        g = (0.3469 * (Hi / 1.29e11) ** (2 / 3), 1.553e-5 * (Hi / 1.29e11) ** 0.5)
        return cal.calibrate_hilltop_ns(Hi, hist, 6, x0=g)
    if name in ('HT4', 'HT6', 'HT8'):
        p = int(name[2]); vg = {4: 0.7406, 6: 0.4169, 8: 0.3391}[p] * (Hi / 1.29e11) ** (2 / 3)
        return cal.calibrate_hilltop(Hi, hist, p, 0.0, vev_guess=vg)
    if name == 'PL2':
        return cal.calibrate_param('PL', 'mu', 7.3246e-5 * (Hi / 1.29e11) ** 2, Hi, ('inst',), fixed=dict(q=2))
    raise ValueError(name)


def turning_points(t, f, fd):
    idx = np.where(np.sign(fd[1:]) != np.sign(fd[:-1]))[0]
    out = []
    for i in idx:
        # linear interpolation of the velocity zero
        w = fd[i] / (fd[i] - fd[i + 1])
        out.append((t[i] + w * (t[i + 1] - t[i]), f[i] + w * (f[i + 1] - f[i])))
    return out


def oscillation(m, bg, V0, h_slice, gam_Hi, tmax_Hi, npts=None):
    """post-inflation run from eps_H = 1 (a = 1 there); times in H_i^-1 (H_i = H at S1)."""
    sq = math.sqrt(V0)
    g = gam_Hi * h_slice
    tend = tmax_Hi / h_slice                      # in 1/sqrt(V0) units (t_Hi = t_code * h_slice)
    sol = ic.post_inflation(m, bg, g, tend, rtol=1e-10, atol=1e-14, dense=True)
    npts = npts or int(max(200001, 80 * math.sqrt(m.m2_min()) * tend))
    t = np.linspace(0, tend, npts)
    Y = sol.sol(t)
    lna, f, fd, rr = Y
    rho = 0.5 * fd ** 2 + np.array([m.v(x) for x in f]) + rr
    H = np.sqrt(rho / 3)
    pmin = m.phi_min()
    tps = turning_points(t, f, fd)
    tp_t = np.array([x[0] for x in tps]); tp_f = np.array([x[1] for x in tps])
    # half-periods -> frequency omega = pi / (t_{j+1} - t_j)
    om = math.pi / np.diff(tp_t) if len(tp_t) > 1 else np.array([])
    mI = math.sqrt(m.m2_min())
    rows = []
    for j in range(min(len(tp_t), 4000)):
        if j % max(1, len(tp_t) // 60) and j > 12:
            continue
        aj = math.exp(float(sol.sol(tp_t[j])[0]))
        rows.append(dict(t_Hi=float(tp_t[j] * h_slice), a=aj, phi_turn=float(tp_f[j]), phiI_turn=float(tp_f[j] - pmin),
                         V_over_V0=float(m.v(tp_f[j])),
                         omega_over_mI=float(om[j] / mI) if j < len(om) else None,
                         omega_over_Hi=float(om[j] / h_slice) if j < len(om) else None))
    # equation of state averaged over windows of ~10 oscillations
    w_inst = (0.5 * fd ** 2 - np.array([m.v(x) for x in f])) / (0.5 * fd ** 2 + np.array([m.v(x) for x in f]))
    win = []
    edges = np.linspace(0, tend, 41)
    for a0, a1 in zip(edges[:-1], edges[1:]):
        sel = (t >= a0) & (t < a1)
        win.append(dict(t_Hi=float(0.5 * (a0 + a1) * h_slice), w_avg=float(np.mean(w_inst[sel])),
                        H_over_Hi=float(np.mean(H[sel]) / h_slice), a=float(np.exp(np.mean(lna[sel]))),
                        rad_frac=float(np.mean(rr[sel] / rho[sel]))))
    K_end = float(np.exp(2 * lna[-1]) * H[-1] / h_slice)
    return dict(gamma_over_Hi=gam_Hi, t_max_Hi=tmax_Hi, n_turning=int(len(tp_t)), turning=rows, windows=win,
                K_at_tmax_endnorm=K_end, rad_frac_end=float(rr[-1] / rho[-1]),
                omega_first_over_mI=float(om[0] / mI) if len(om) else None,
                omega_min_over_Hi=float(om.min() / h_slice) if len(om) else None,
                omega_max_over_Hi=float(om.max() / h_slice) if len(om) else None)


def yy_channels(mI_Hi, gI):
    """perturbative inflaton -> Y Y widths (H_i units) through -(g_I e^{i alpha_I} phi_I Sigma* Y^2 + h.c.) at the late
    vacuum (u = mu real, Sigma = h/sqrt2): couplings phi a^2: |Sigma| g cos(alpha), phi b^2: -|Sigma| g cos(alpha),
    phi a b: -2 |Sigma| g sin(alpha), with b the light (m_L) and a the heavy (m_H) real component [D]."""
    ca, sa = math.cos(ALPHA_I), math.sin(ALPHA_I)
    def beta2(m1, m2):
        x = 1 - (m1 + m2) ** 2 / mI_Hi ** 2; y = 1 - (m1 - m2) ** 2 / mI_Hi ** 2
        return math.sqrt(x * y) if x > 0 else 0.0
    G_LL = (SIG * gI * ca) ** 2 / (8 * math.pi * mI_Hi) * beta2(M_L, M_L)
    G_HH = (SIG * gI * ca) ** 2 / (8 * math.pi * mI_Hi) * beta2(M_H, M_H)
    G_LH = (SIG * gI * sa) ** 2 / (4 * math.pi * mI_Hi) * beta2(M_L, M_H)
    return dict(G_LL=G_LL, G_LH=G_LH, G_HH=G_HH, G_tot=G_LL + G_LH + G_HH, thresholds=dict(LL=2 * M_L, LH=M_L + M_H, HH=2 * M_H))


def run(name):
    t0 = time.time()
    out = dict(name=name, method=__doc__, by_Hi=[])
    for Hi in HI_LIST:
        m, r = calib(name, Hi)
        bg = r['bg']; V0 = r['V0']; hs = r['h_slice']
        sq = math.sqrt(V0)
        pmin = m.phi_min()
        fS1, fNS1 = ic.bg_at(bg, 1.0); fe, fNe = ic.bg_at(bg, 0.0)
        he, _ = ic.h_of(m, fe, fNe)
        mI = math.sqrt(m.m2_min()) if m.m2_min() > 0 else float('nan')     # sqrt(V0) units
        mI_Hi = mI / hs
        phiI_S1 = fS1 - pmin; phiI_end = fe - pmin
        Hi_MP = Hi / ic.MPL_GEV
        gI_S1 = U_GAP * Hi_MP / abs(phiI_S1); gI_end = U_GAP * Hi_MP / abs(phiI_end)
        # sign convention: mu_eff = mu + G phihat_I with G = (10i - mu)/phihat_I(slice): alpha_I = arg(10i - mu) if
        # phi_I(slice) > 0, arg(10i - mu) - pi if phi_I(slice) < 0 (hilltop: phi_I = phi - v < 0 during inflation)
        aI = ALPHA_I if phiI_S1 > 0 else ALPHA_I - math.pi
        rec = dict(H_i_GeV=Hi, params=m.p, V0_MP4=V0, V0_quarter_GeV=V0 ** 0.25 * ic.MPL_GEV, phi_min_MP=pmin,
                   Nstar=r['Nstar'], As=r['As'], ns=r['ns'], ns_pull_sigma=(r['ns'] - ic.NS_TARGET) / ic.NS_SIG,
                   r=r['r'], alpha_s=r['alpha_s'], n_t=r['n_t'], Hstar_GeV=r['Hstar_GeV'],
                   H_S1_GeV=Hi, H_end_GeV=sq * he * ic.MPL_GEV, Hend_over_Hi=he / hs, Hstar_over_Hi=r['Hstar_GeV'] / Hi,
                   T_inst_GeV=r['ninfo']['T_inst'], rho_end_GeV4=r['ninfo']['rho_end'],
                   m_I_GeV=mI * sq * ic.MPL_GEV, m_I_over_Hi=mI_Hi, m_I_over_Hend=mI / he,
                   phi_S1_MP=fS1, phi_end_MP=fe, phiI_S1_MP=phiI_S1, phiI_end_MP=phiI_end,
                   phidot_S1_over_HiMP=fNS1 * 1.0, phidot_end_over_HiMP=fNe * he / hs,
                   g_I_S1=gI_S1, g_I_end=gI_end, alpha_I_S1=aI,
                   K_end_gamma05=r['cache']['info']['K_end_norm'] if 'info' in r['cache'] else None)
        # slice tables
        tab = []
        for Ne in NEB + [r['Nstar']]:
            f, fN = ic.bg_at(bg, Ne)
            h, eps = ic.h_of(m, f, fN)
            tab.append(dict(N_before_end=Ne, phi_MP=f, phiI_MP=f - pmin, H_GeV=sq * h * ic.MPL_GEV, H_over_Hi=h / hs,
                            eps_H=eps, eta_V=m.d2v(f) / m.v(f), s_S1=(f - pmin) / phiI_S1, s_end=(f - pmin) / phiI_end,
                            hr_S1=h / hs, hr_end=h / he))
        rec['history'] = tab
        # alternative histories for N_*
        alts = []
        for hist in [('inst',), ('gamma', 0.1), ('gamma', 2.0)] + [('wT', w, T) for w in (0.0, 1 / 3) for T in (1e14, 1e12, 1e10, 1e8)]:
            if name == 'PL2' and hist[0] == 'gamma':
                continue
            if hist[0] == 'wT' and hist[1] == 1 / 3 and hist[2] != 1e10:
                continue
            try:
                ev = cal.evaluate(m, Hi, hist, need_spec=True)
                alts.append(dict(history=list(hist), Nstar=ev['Nstar'], As_at_fixed_shape=ev['As'], ns=ev['ns'], r=ev['r'],
                                 alpha_s=ev['alpha_s']))
            except Exception as e:
                alts.append(dict(history=list(hist), error=str(e)))
        rec['Nstar_alternatives'] = alts
        # oscillation (not for PL2: m_I/H ~ 3e4, homogeneous run not useful)
        if name != 'PL2':
            tmax = 12.0 if mI_Hi < 100 else 6.0
            rec['osc_gamma0'] = oscillation(m, bg, V0, hs, 0.0, tmax)
            rec['osc_gamma05'] = oscillation(m, bg, V0, hs, 0.5, tmax)
        # inflaton -> YY through the holder (late vacuum), g_I at S1 and at Send
        yy = {}
        for lab, gI in (('S1', gI_S1), ('Send', gI_end)):
            ch = yy_channels(mI_Hi, gI)
            # hidden heat if the total width is Gamma_tot: N_Y = 2 BR rho_end/m_I (a = 1 at the end), R = N_Y * DEPOSIT [E]
            rho_end_Hi4 = 3 * (he / hs) ** 2 / Hi_MP ** 2
            NI = rho_end_Hi4 / mI_Hi
            ch['N_inflaton_quanta_Hi3'] = NI
            ch['R_YY_vs_Gamma_tot'] = {str(G): 2 * ch['G_tot'] / G * NI * DEPOSIT_YL for G in (0.1, 0.5, 2.0, 10.0)}
            yy[lab] = ch
        rec['inflaton_to_YY'] = yy
        out['by_Hi'].append(rec)
        print('%s H_i=%.3g: params %s N*=%.3f As=%.4e ns=%.5f (%.2f sig) r=%.3e als=%.2e H*/Hi=%.4f Hend/Hi=%.4f mI/Hi=%.2f mI/Hend=%.2f '
              'phiI(S1)=%.4e phiI(end)=%.4e gI(S1)=%.3e  [%.0fs]' % (name, Hi, m.p, r['Nstar'], r['As'], r['ns'], rec['ns_pull_sigma'], r['r'],
              r['alpha_s'], rec['Hstar_over_Hi'], rec['Hend_over_Hi'], mI_Hi, rec['m_I_over_Hend'], phiI_S1, phiI_end, gI_S1, time.time() - t0), flush=True)
    json.dump(out, open(RES / ('bench_%s.json' % name), 'w'), indent=1, default=float)


if __name__ == '__main__':
    for n in sys.argv[1:]:
        run(n)
