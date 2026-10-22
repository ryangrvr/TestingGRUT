"""Part 4: realistic expansion history of QV1 from T = 1e3 GeV to today and into the future.

Two q-masses:
  natural  m_q = sqrt(2/3) M = 1.988e18 GeV (KV shape with chi = 1/M^4, q0 K' = K0): validated reduction R1 only
           (m_q/H >= 1e30 throughout);
  bracket  m_q = 4 meV (lowest unscreened value in the ledger, Eot-Wash alpha = 1/3 [E]): R1 where m_q/H >= 100,
           and the FULL equations from T = 1e3..1e6 GeV through the electroweak crossover (m_q/H ~ 3 at 1 TeV).
Matter: SM plasma with the Saikawa-Shirai EoS (p(T=0) = 0, i.e. normalised as P_q dictates), plus dust
(baryons + CDM, omega_b = 0.0224, omega_c = 0.120: Planck 2018 abstract [Q], borrowed from a LambdaCDM fit).
Under P_q the late universe has no vacuum term. Reduction R1 (validated in s2):
    3 K0 H^2 (1 + (S + S_N)/(3 K0 m^2)) = rho + S^2/(12 K0 m^2),   S = rho - 3p (plasma trace + dust)
q-sector outputs: rho_V = S^2/(12 K0 m^2); the q-correction delta rho_q = [S^2/(12K0) - H^2 (S + S_N)]/m^2;
delta G/G = -Y = -S/(3 K0 m^2); Gdot/G = -H S_N/(3 K0 m^2).
Run: python3 -I s4_realistic_history.py [outdir] [eos_dir]   (about 5-10 min)
"""
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qv1_core import KV_SHAPE, QV1, SMMatter, integrate_full  # noqa: E402
from sm_eos import DEFAULT_DIR, SMEoS  # noqa: E402

M = 2.435e18                      # GeV, reduced Planck mass (brief)
K0 = M ** 2
HBAR = 6.582119569e-25            # GeV s
YR = 3.15576e7                    # s (Julian year)
GYR = 1e9 * YR
KB = 8.617333262e-14              # GeV/K
T0 = 2.7255 * KB                  # GeV; Fixsen 2009 (0911.1955) abstract [Q]: 2.72548 +- 0.00057 K
MPC_M = 3.0856775814913673e22
H100 = 100e3 / MPC_M * HBAR       # GeV
RHO_CRIT100 = 3 * K0 * H100 ** 2  # GeV^4 per h^2
OMEGA_B, OMEGA_C = 0.0224, 0.120  # Planck 2018 VI abstract [Q]
RHO_LAMBDA = 2.52e-47             # GeV^4 (brief; Planck 2018 inputs) - used only for QV1-cal comparisons
M_NAT = math.sqrt(2.0 / 3.0) * M
M_4MEV = 4e-12


class Background:
    """GR + SM plasma + dust history as a function of T (entropy conservation), on a fine grid."""

    def __init__(self, eos, T_hi=1e3, a_future=1e3, n=40000):
        self.e = eos
        s0 = eos.s(T0)
        self.rho_d0 = (OMEGA_B + OMEGA_C) * RHO_CRIT100
        lnT = np.linspace(math.log(T_hi), math.log(T0 / a_future), n)
        T = np.exp(lnT)
        Tc = np.maximum(T, eos.T[0])
        grho = eos.g_rho(Tc)
        gs = eos.g_s(Tc)
        self.T = T
        self.rho_pl = math.pi ** 2 / 30 * grho * T ** 4
        self.s = 2 * math.pi ** 2 / 45 * gs * T ** 3
        self.rho_d = self.rho_d0 * self.s / s0
        self.N = -np.log(T / T0) - np.log(gs / eos.g_s(T0)) / 3.0      # N = ln(a/a0)
        Sfrac = eos.S_over_rho(Tc)
        Sfrac = np.where(T < eos.T[0], 0.0, Sfrac)
        self.S_pl = Sfrac * self.rho_pl
        self.S_err = eos.S_over_rho_err(Tc) * self.rho_pl
        self.S = self.S_pl + self.rho_d
        self.rho = self.rho_pl + self.rho_d
        self.S_N = np.gradient(self.S, self.N)
        self.p = (1 - Sfrac) / 3.0 * self.rho_pl

    def history(self, m, const=0.0, order=1):
        """H(N) for the reduction (order 0: GR + const; 1: + adiabatic q-correction), t(N) by quadrature."""
        rho = self.rho + const
        if order == 0 or m is None:
            H2 = rho / (3 * K0)
        else:
            H2 = (rho + self.S ** 2 / (12 * K0 * m ** 2)) / (3 * K0 * (1 + (self.S + self.S_N) / (3 * K0 * m ** 2)))
        H = np.sqrt(np.maximum(H2, 0.0))
        # t: start at T_hi with t = 1/(2H) (radiation era; error O(dln g/dlnT) < 1e-3)
        dt = np.concatenate([[1.0 / (2 * H[0])], np.diff(self.N) * 0.5 * (1 / H[1:] + 1 / H[:-1])])
        t = np.cumsum(dt)
        return H, t


def age_and_epochs(bg, H, t):
    i0 = np.argmin(np.abs(bg.T - T0))
    return dict(t0_Gyr=float(t[i0] * HBAR / GYR), H0_kms_Mpc=float(H[i0] / H100 * 100))


def _dHH(H1, H0, drho, rho):
    """H_QV1/H_GR - 1; analytic first order delta rho_q/(2 rho) when below double precision."""
    if abs(drho / rho) < 1e-9:
        return drho / (2 * rho)
    return H1 / H0 - 1


def epoch_table(bg, m, label):
    H1, t1 = bg.history(m, order=1)
    H0g, t0g = bg.history(None, order=0)
    HL, tL = bg.history(None, const=RHO_LAMBDA, order=0)
    rows = []
    marks = [('T=1e3 GeV', 1e3), ('EW crossover T=160 GeV', 160.0), ('T=100 GeV', 100.0), ('T=10 GeV', 10.0),
             ('QCD crossover T=155 MeV', 0.155), ('T=10 MeV (nu decoupling scale)', 0.010),
             ('e+e- annihilation T=0.2 MeV', 2e-4), ('BBN (D bottleneck) T=0.07 MeV', 7e-5)]
    # equality, recombination, today, future
    zeq_i = int(np.argmin(np.abs(np.log(bg.rho_d / bg.rho_pl))))
    marks.append(('matter-radiation equality (QV1)', float(bg.T[zeq_i])))
    marks.append(('recombination z=1090', T0 * 1091))
    marks.append(('today T0=2.7255 K', T0))
    marks.append(('future a=10 a0', T0 / 10))
    marks.append(('future a=1000 a0', T0 / 1000 * 1.0001))
    for name, Tm in marks:
        i = int(np.argmin(np.abs(np.log(bg.T / Tm))))
        S, SN, rho, H = bg.S[i], bg.S_N[i], bg.rho[i], H1[i]
        Y = S / (3 * K0 * m ** 2)
        Gdot = -H * SN / (3 * K0 * m ** 2)
        rows.append(dict(epoch=name, T_GeV=float(bg.T[i]), z=float(math.exp(-bg.N[i]) - 1), t_s=float(t1[i] * HBAR),
                         H_GeV=float(H), m_over_H=float(m / H) if H > 0 else None,
                         H_QV1_over_GR_minus_1=float(_dHH(H1[i], H0g[i], (S ** 2 / (12 * K0) - H ** 2 * (S + SN))
                                                         / m ** 2, rho)),
                         H_QV1_over_LCDM=float(H1[i] / HL[i]),
                         S_over_rho=float(S / rho), S_err_over_rho=float(bg.S_err[i] / rho),
                         rhoV_q_GeV4=float(S ** 2 / (12 * K0 * m ** 2)),
                         delta_rho_q_GeV4=float((S ** 2 / (12 * K0) - H ** 2 * (S + SN)) / m ** 2),
                         delta_rho_q_over_rho=float((S ** 2 / (12 * K0) - H ** 2 * (S + SN)) / m ** 2 / rho),
                         rhoV_q_over_rhoLambda=float(S ** 2 / (12 * K0 * m ** 2) / RHO_LAMBDA),
                         dG_over_G=float(-Y), Gdot_over_G_per_yr=float(Gdot / HBAR * YR)))
    # series for figures (200 points)
    idx = np.unique(np.linspace(0, len(bg.T) - 1, 220).astype(int))
    series = dict(T_GeV=bg.T[idx].tolist(), z=(np.exp(-bg.N[idx]) - 1).tolist(), t_s=(t1[idx] * HBAR).tolist(),
                  H_QV1_GeV=H1[idx].tolist(),
                  H_QV1_over_GR_minus_1=[float(_dHH(H1[i], H0g[i], (bg.S[i] ** 2 / (12 * K0) - H1[i] ** 2 *
                                                                     (bg.S[i] + bg.S_N[i])) / m ** 2, bg.rho[i]))
                                         for i in idx],
                  H_QV1_over_LCDM=(H1[idx] / HL[idx]).tolist(), S_over_rho=(bg.S[idx] / bg.rho[idx]).tolist(),
                  rhoV_q_GeV4=(bg.S[idx] ** 2 / (12 * K0 * m ** 2)).tolist(),
                  dG_over_G=(-bg.S[idx] / (3 * K0 * m ** 2)).tolist(),
                  Gdot_over_G_per_yr=(-H1[idx] * bg.S_N[idx] / (3 * K0 * m ** 2) / HBAR * YR).tolist())
    ages = dict(QV1=age_and_epochs(bg, H1, t1), GR_no_vacuum=age_and_epochs(bg, H0g, t0g),
                QV1cal_LCDM=age_and_epochs(bg, HL, tL))
    # asymptotic behaviour: t H at the last point (2/3 for matter domination; QV1-cal -> de Sitter)
    asym = dict(QV1_tH_end=float(t1[-1] * H1[-1]), QV1cal_H_end_over_HLambda=float(HL[-1] / math.sqrt(RHO_LAMBDA /
                                                                                                        (3 * K0))),
                QV1_rhoV_q_end_GeV4=float(bg.S[-1] ** 2 / (12 * K0 * m ** 2)))
    return dict(label=label, m_q_GeV=m, rows=rows, ages=ages, asymptotic=asym, series=series)


# ------------------------------------------------------------------------------------------- full run, 4 meV
def reduced_sm(mat, q, x0, rd0, N0, N1, order=1):
    """R0/R1 reduction on the SM matter path (same w(ln rho) table as the full run); units K0 = 1, m = q.m."""
    xg, wg = mat.xg, mat.wg
    dw = np.gradient(wg, xg)
    m2 = q.m2

    def SN_of(x, rd):
        w = np.interp(x, xg, wg)
        wx = np.interp(x, xg, dw)
        e = math.exp(x)
        xN = -3 * (1 + w)
        return e * (1 - 3 * w), (e * (1 - 3 * w) - 3 * e * wx) * xN - 3 * rd

    def Hof(N, x, rd):
        rho = math.exp(x) + rd
        if order == 0:
            return math.sqrt(rho / 3.0)
        Spl, SN = SN_of(x, rd)
        S = Spl + rd
        corr = (S + SN) / (3 * m2)
        if abs(corr) > 0.3 or S * S / (12 * m2) > 0.3 * rho:     # m/H not large: R1 meaningless, fall back to R0
            return math.sqrt(rho / 3.0)
        return math.sqrt((rho + S * S / (12 * m2)) / (3.0 * (1 + corr)))

    def f(N, u):
        x, rd, t, J = u
        H = Hof(N, x, rd)
        return [-3 * (1 + np.interp(x, xg, wg)), -3 * rd, 1.0 / H, H]
    sol = solve_ivp(f, (N0, N1), [x0, rd0, 0.0, 0.0], method='DOP853', rtol=1e-12, atol=1e-30, dense_output=True,
                    max_step=0.01)
    return sol, Hof, SN_of


def T_at_m_over_H(eos, m, ratio):
    """Temperature at which m/H = ratio on the GR + SM history (radiation era)."""
    lo, hi = math.log(1e-4), math.log(1e17)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        H = math.sqrt(eos.rho(math.exp(mid)) / (3 * K0))
        if m / H > ratio:
            lo = mid
        else:
            hi = mid
    return math.exp(0.5 * (lo + hi))


def full_run_4meV(eos, T_i, mode, T_f=None, rtol=1e-11, n_windows=150, Sfrac_override=None, m=M_4MEV,
                  mH_end=1e5):
    """mode 'adiabatic': Y, V on the quasi-static track (1st order) at T_i; mode 'q0': Y = V = 0 (q at the P_q
    equilibrium q0). Units: K0 = 1, time 1/m, rho in K0 m^2. Default end: m/H = mH_end."""
    t0w = time.time()
    unit = K0 * m ** 2
    if T_f is None:
        T_f = T_at_m_over_H(eos, m, mH_end)
    e = eos
    if Sfrac_override is not None:
        e = Sfrac_override
    xg, wg = e.w_of_lnrho(unit, Tmin=0.3 * T_f, Tmax=3 * T_i)
    mat = SMMatter(xg, wg)
    q = QV1(1.0, **KV_SHAPE)
    rho_i = e.rho(T_i) / unit
    s0 = e.s(T0)
    rd_i = (OMEGA_B + OMEGA_C) * RHO_CRIT100 * e.s(T_i) / s0 / unit
    x_i = math.log(rho_i)
    red1, H1f, SN_of = reduced_sm(mat, q, x_i, rd_i, 0.0, 30.0, order=1)
    red0, H0f, _ = reduced_sm(mat, q, x_i, rd_i, 0.0, 30.0, order=0)
    # N at T_f: solve rho(N) = rho(T_f)
    x_f = math.log(e.rho(T_f) / unit)
    Ns = np.linspace(0, 30, 30001)
    xs = red1.sol(Ns)[0]
    N_f = float(np.interp(-x_f, -xs, Ns))
    Spl, SN = SN_of(x_i, rd_i)
    H_i = H1f(0.0, x_i, rd_i)
    if mode == 'adiabatic':
        Y_i = (Spl + rd_i) / 3.0
        V_i = H_i * SN / 3.0
    else:
        Y_i, V_i = 0.0, 0.0
    H_i = q.H_from_constraint(Y_i, V_i, math.exp(x_i) + rd_i)
    t_i = 1.0 / (2 * H_i)
    t_f = t_i + float(red1.sol(N_f)[2])
    period = 2 * math.pi
    lo = max(t_i + 2 * period, t_i * 1.5)
    centres = np.exp(np.linspace(math.log(lo), math.log(t_f - period), n_windows))
    k_per = 1
    pts = [c + k_per * period * (np.arange(17) / 16 - 0.5) for c in centres]
    keep = [0]
    for i in range(1, len(pts)):
        if pts[i][0] > pts[keep[-1]][-1]:
            keep.append(i)
    pts = np.array([pts[i] for i in keep])
    y0 = [Y_i, V_i, H_i, 0.0, 0.0, x_i, rd_i]
    # error-control scale for Y: the quasi-static value when the q-mode is heavy; when m/H_i < 1 that value is
    # meaningless (Y never reaches it) and the misalignment is ~ (S/rho) x (e-folds) ~ 1e-4..1e-3
    Ysc = max(abs(Y_i), (Spl + rd_i) / 3.0, 1e-30) if 1.0 / H_i > 1.0 else 1e-4
    scale = [Ysc, Ysc * max(H_i, 1.0), H_i, 1.0, 1.0, abs(x_i) + 1.0, rd_i]
    out = integrate_full(q, mat, y0, t_i, pts.ravel(), rtol=rtol, atol=1e-16, scale=scale)
    out = out.reshape(pts.shape[0], pts.shape[1], 7)
    rows = []
    for k in range(pts.shape[0]):
        Y, V, H, N, J, x, rd = out[k].T
        Hbar = (J[-1] - J[0]) / (N[-1] - N[0])
        Hbar0 = (red0.sol(N[-1])[3] - red0.sol(N[0])[3]) / (N[-1] - N[0])
        Hbar1 = (red1.sol(N[-1])[3] - red1.sol(N[0])[3]) / (N[-1] - N[0])
        E = []
        for j in range(16):
            xr = red1.sol(N[j])
            Sp, SNj = SN_of(xr[0], xr[1])
            Hr = H1f(N[j], xr[0], xr[1])
            Yq = (Sp + xr[1]) / 3.0
            Vq = Hr * SNj / 3.0
            E.append(0.75 * ((V[j] - Vq) ** 2 + (Y[j] - Yq) ** 2) * math.exp(3 * N[j]))
        rho_m = math.exp(x[8]) + rd[8]
        C, Crel = q.constraint(Y[8], V[8], H[8], rho_m)
        T_here = float(np.interp(-x[8], -np.log(e.rho(e.T[::-1]) / unit), e.T[::-1]))
        Edust = float(np.mean(E)) * math.exp(-1.5 * (N[0] + N[-1]))
        rows.append(dict(N=float(N[8]), T_GeV=T_here, m_over_H=float(1.0 / Hbar), Y=float(np.mean(Y[:16])),
                         eps0=float(Hbar / Hbar0 - 1), eps1=float(Hbar / Hbar1 - 1),
                         eps2=float(Hbar / math.sqrt(Hbar1 ** 2 + Edust / 3.0) - 1),
                         Eosc_a3=float(np.mean(E)), Eosc_a3_spread=float(np.std(E)), C_rel=float(Crel)))
    # relic: E a^3 at the end (comoving to a(N) normalised with a = 1 at T_i), scaled to today by entropy
    late = [r for r in rows if r['m_over_H'] > 300] or rows[-5:]
    E_end = np.mean([r['Eosc_a3'] for r in late[-20:]])
    E_end_spread = np.std([r['Eosc_a3'] for r in late[-20:]]) / E_end if E_end > 0 else None
    # convert: E(T) = E_a3 e^{-3N(T)} (units K0 m^2); rho_osc,0 = E(T) s0/s(T)
    rowf = late[-1]
    E_Tf = E_end * math.exp(-3 * rowf['N']) * unit
    rho_osc0 = E_Tf * s0 / e.s(rowf['T_GeV'])
    Omega_h2 = rho_osc0 / RHO_CRIT100
    sel = [r for r in rows if 100 < r['m_over_H'] < 3e4]
    return dict(T_i=T_i, mode=mode, m_q_GeV=m, m_over_H_i=1.0 / H_i, Y_i=Y_i, V_i=V_i,
                N_f=N_f, T_f=T_f, E_a3_late=float(E_end), E_a3_late_rel_spread=E_end_spread,
                rho_osc_today_GeV4=float(rho_osc0), Omega_osc_h2=float(Omega_h2),
                Omega_osc_h2_over_Omega_c_h2=float(Omega_h2 / OMEGA_C),
                max_abs_eps1_100_to_3e4=float(max(abs(r['eps1']) for r in sel)) if sel else None,
                max_abs_eps0_100_to_3e4=float(max(abs(r['eps0']) for r in sel)) if sel else None,
                max_abs_eps2_100_to_3e4=float(max(abs(r['eps2']) for r in sel)) if sel else None,
                T_osc_GeV_mH3=float(T_at_m_over_H(e, m, 3.0)),
                max_abs_C_rel=float(max(abs(r['C_rel']) for r in rows)), wall_s=time.time() - t0w,
                rows=rows[::3])


class _EoSVariant:
    """Wrap SMEoS with a modified trace fraction above T_cut (to bracket the high-T uncertainty)."""

    def __init__(self, base, factor, T_cut=300.0, add_sigma=0.0):
        self.__dict__.update(base.__dict__)
        self.base = base
        Sf = base.Sfrac.copy()
        hi = base.T > T_cut
        Sf[hi] = Sf[hi] * factor + add_sigma * base.Serr[hi]
        self.Sfrac = Sf

    def __getattr__(self, k):
        return getattr(self.base, k)

    def w_of_lnrho(self, rho_unit, Tmin=1e-3, Tmax=1e7):
        sel = (self.T >= Tmin) & (self.T <= Tmax)
        T = self.T[sel]
        rho = math.pi ** 2 / 30 * self.grho[sel] * T ** 4 / rho_unit
        w = (1 - self.Sfrac[sel]) / 3.0
        x = np.log(rho)
        o = np.argsort(x)
        return x[o], w[o]


def _full_job(arg):
    eos_dir, j = arg
    j = dict(j)
    eos = SMEoS(eos_dir)
    tag = j.pop('tag')
    var = j.pop('variant', None)
    if var is not None:
        eos_use = _EoSVariant(eos, var[0], add_sigma=var[1])
    else:
        eos_use = eos
    m = j.get('m', M_4MEV)
    if 'mH_i' in j:
        j['T_i'] = T_at_m_over_H(eos, m, j.pop('mH_i'))
    r = full_run_4meV(eos_use, **j)
    r['tag'] = tag
    if var is not None:
        r['eos_variant'] = var
    # decay estimate [E]: Gamma <= 4 m^3/(192 pi M^2) (Higgs doublet, if kinematically open), vs H at T = 1 MeV
    Gam = 4 * m ** 3 / (192 * math.pi * M ** 2)
    r['Gamma_max_GeV'] = Gam
    r['Gamma_max_over_H_at_1MeV'] = Gam / math.sqrt(eos.rho(1e-3) / (3 * K0))
    return r


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                                                                  'results')
    eos_dir = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_DIR
    T00 = time.time()
    eos = SMEoS(eos_dir)
    bg = Background(eos)
    res = {}
    res['natural'] = epoch_table(bg, M_NAT, 'natural m_q = sqrt(2/3) M')
    res['bracket_4meV_reduction'] = epoch_table(bg, M_4MEV, 'm_q = 4 meV (reduction; valid where m/H > ~100)')
    # full-equation runs at 4 meV (and a mass scan), in parallel
    jobs = [dict(T_i=1e3, mode='adiabatic', tag='4meV adiabatic at T=1 TeV'),
            dict(T_i=1e3, mode='q0', tag='4meV q=q0 at T=1 TeV'),
            dict(T_i=1e4, mode='q0', tag='4meV q=q0 at T=10 TeV'),
            dict(T_i=1e5, mode='q0', tag='4meV q=q0 at T=100 TeV'),
            dict(T_i=1e6, mode='q0', tag='4meV q=q0 at T=1000 TeV'),
            dict(T_i=1e5, mode='q0', variant=(0.0, 0.0), tag='4meV q=q0 at 100 TeV, S/rho = 0 above 300 GeV'),
            dict(T_i=1e5, mode='q0', variant=(1.0, 1.0), tag='4meV q=q0 at 100 TeV, S/rho + 1 sigma above 300 GeV'),
            dict(T_i=1e3, mode='adiabatic', rtol=1e-12, tag='4meV adiabatic at 1 TeV, rtol 1e-12 (tolerance check)'),
            dict(T_i=1e3, mode='adiabatic', mH_end=3e5, tag='4meV adiabatic at 1 TeV, run to m/H = 3e5')]
    for mq, lab in [(1e-9, '1 eV'), (1e-6, '1 keV'), (1e-3, '1 MeV'), (1.0, '1 GeV')]:
        jobs.append(dict(m=mq, mH_i=3.0, mode='adiabatic', tag='m=%s adiabatic at m/H=3' % lab))
        jobs.append(dict(m=mq, mH_i=0.01, mode='q0', tag='m=%s q=q0 at m/H=0.01' % lab))
    with ProcessPoolExecutor(max_workers=4) as ex:
        runs = list(ex.map(_full_job, [(eos_dir, j) for j in jobs]))
    for r in runs:
        print(r['tag'], 'Omega h2 %.3e' % r['Omega_osc_h2'], 'eps2 %s' % r['max_abs_eps2_100_to_3e4'],
              'wall %.1f' % r['wall_s'], flush=True)
    res['bracket_4meV_full'] = runs
    res['inputs'] = dict(M_GeV=M, T0_GeV=T0, T0_K=2.7255, omega_b=OMEGA_B, omega_c=OMEGA_C,
                         rho_crit_100_GeV4=RHO_CRIT100, rho_d0_GeV4=bg.rho_d0, rho_Lambda_GeV4=RHO_LAMBDA,
                         m_nat_GeV=M_NAT, m_bracket_GeV=M_4MEV, hbar_GeV_s=HBAR)
    res['meta'] = dict(wall_s=time.time() - T00)
    with open(os.path.join(outdir, 'realistic_history.json'), 'w') as fh:
        json.dump(res, fh, indent=1)
    for lab in ('natural', 'bracket_4meV_reduction'):
        print(lab, res[lab]['ages'])
        for rw in res[lab]['rows']:
            print('  %-40s T=%.3e z=%.3e t=%.3e s m/H=%.2e dH/H=%.2e rhoVq=%.2e dG/G=%.2e Gdot/G=%.2e/yr' % (
                rw['epoch'], rw['T_GeV'], rw['z'], rw['t_s'], rw['m_over_H'], rw['H_QV1_over_GR_minus_1'],
                rw['rhoV_q_GeV4'], rw['dG_over_G'], rw['Gdot_over_G_per_yr']))
    print('total', time.time() - T00)


if __name__ == '__main__':
    main()
