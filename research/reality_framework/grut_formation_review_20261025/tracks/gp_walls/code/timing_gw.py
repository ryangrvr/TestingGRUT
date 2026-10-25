#!/usr/bin/env python3
"""Wall-network annihilation, gravitational waves, energy injection and Q-ball loading
(Claude review, track gp_walls). Own code; Codex's point.py was read but not run.

Standard scaling estimate (Saikawa 1703.02576 eqs. 2.27-2.28, read):
  t_ann = C_ann A sigma / Delta,  radiation era H = 1/(2t),  H = pi sqrt(g*/90) T^2 / M_red.
g*(T), g*s(T) are computed self-consistently (photons, 3 nu, e, mu, pi) for T < 150 MeV.
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

MRED = 2.435e18
HBAR = 6.582119569e-25
T0 = 2.7255 * 8.617333262e-14          # GeV
LAM, V, GX, LAMX = 0.03, 0.924089637375362, 0.1, 0.001
DELTA = 1e-22
SIGMA = 2 * np.sqrt(2 * LAM) * V**3 / 3
MX = GX * V
MS = np.sqrt(2 * LAM) * V
HD = 0.01
MD = np.sqrt(GX**2 + HD**2) * V
BDM = 4.4e-10                           # calibrated rho_DM/s, GeV (Codex input, unchanged)
YQ = 4.8294009076988175e-09             # fixed net charge yield (Codex input)
QREF, MQ = 200000.0, 17064.537454       # reference object (own profile reproduces E)
U_WALL_CORE = 1398.7                    # own wall_core.py, product trial, GeV (attraction)
OMEGA_Q = 0.8611951308                  # own core_profile.py


def _rho_s(m, T, g, fermion):
    """g* (energy) and g*s (entropy) contributions of one species of mass m, internal dof g."""
    y = m / T
    if y > 40:
        return 0.0, 0.0
    sgn = 1.0 if fermion else -1.0
    E = lambda u: np.sqrt(u * u + y * y)
    e = quad(lambda u: u * u * E(u) / (np.exp(E(u)) + sgn), 0, 80)[0]
    sp_ = quad(lambda u: (u * u * E(u) + u**4 / (3 * E(u))) / (np.exp(E(u)) + sgn), 0, 80)[0]
    return g * 15 / np.pi**4 * e, g * 45 / (4 * np.pi**4) * sp_


def gstar(T):
    gs, gss = 2.0 + 5.25, 2.0 + 5.25       # photons + 3 neutrinos (coupled, T >~ MeV)
    for m, g, f in [(0.000511, 4, True), (0.10566, 4, True), (0.13957, 2, False), (0.13498, 1, False)]:
        a, b = _rho_s(m, T, g, f)
        gs += a; gss += b
    return gs, gss


def H_of_T(T):
    g, _ = gstar(T)
    return np.pi * np.sqrt(g / 90) * T * T / MRED


def T_of_H(H):
    return brentq(lambda T: H_of_T(T) - H, 1e-7, 0.149)


def entropy(T):
    return 2 * np.pi**2 / 45 * gstar(T)[1] * T**3


def standard(A, C):
    t = C * A * SIGMA / DELTA
    H = 1 / (2 * t)
    T = T_of_H(H)
    return dict(A=A, C_ann=C, t_ann_s=t * HBAR, H_ann_GeV=H, T_ann_GeV=T, gstar=gstar(T)[0])


def T_fixed_g(A, C, g):
    H = DELTA / (2 * C * A * SIGMA)
    return np.sqrt(H * MRED / (np.pi * np.sqrt(g / 90)))


def gw_today(A, C, eps_gw=0.7):
    r = standard(A, C)
    H, T = r['H_ann_GeV'], r['T_ann_GeV']
    g, gs = gstar(T)
    om_emit = eps_gw * A**2 * SIGMA**2 / (24 * np.pi * H**2 * MRED**4)       # rho_gw = eps G A^2 sigma^2
    Omega_gamma_h2 = 2.473e-5
    gs0 = 3.931
    om0h2 = om_emit * Omega_gamma_h2 * (g / 2) * (gs0 / gs)**(4 / 3)
    f0 = H / (2 * np.pi * HBAR) * 2 * np.pi * (T0 / T) * (gs0 / gs)**(1 / 3)  # f ~ a_ann H_ann / a0 (Saikawa 3.10)
    # Saikawa eq. (3.9)/(3.10) closed forms for cross-checking
    s39 = 7.2e-18 * eps_gw * A**2 * (gs / 10)**(-4 / 3) * (SIGMA / 1e9)**2 * (T / 1e-2)**(-4)
    s310 = 1.1e-9 * (g / 10)**0.5 * (gs / 10)**(-1 / 3) * (T / 1e-2)
    return dict(A=A, C_ann=C, eps_gw=eps_gw, T_ann_GeV=T, Omega_emit=om_emit, Omega_gw_h2_today=om0h2,
                f_peak_Hz=f0, Saikawa_3_9=s39, Saikawa_3_10_Hz=s310)


def babichev_variants(Ti_list=(40.0, 86.0)):
    """Babichev et al. 2504.07902 (read via HTML: eqs. 25-29). Bias eps chi^3, V_bias ~ 2 eps v^3."""
    epsr = DELTA / (2 * LAM * V**4)                  # eps/(lam v)
    rows = []
    for Ck in (0.9, 1.0, 1.4):
        H = 0.1 * np.sqrt(LAM) * V * epsr**(2 / 3) / Ck**2        # their eq. 27 (H_i = sqrt(lam) v)
        rows.append(dict(variant='B1: eq.27 as normalised (formation at H_i = sqrt(lam) v)', C_kcut=Ck,
                         H_ann_GeV=H, T_ann_GeV=float(T_of_H(H)) if H < H_of_T(0.149) else float(np.sqrt(H * MRED / (np.pi * np.sqrt(61.75 / 90))))))
    for Ti in Ti_list:
        gi, gsi = 86.25, 86.25
        Hi = np.pi * np.sqrt(gi / 90) * Ti**2 / MRED
        for A in (0.8,):
            rho_i = A * SIGMA * 2 * Hi                 # rho_wall = A sigma / t
            growth = (rho_i / DELTA)**(1 / 3)          # dust-like rho_wall ~ a^-3 until rho_wall = V_bias
            T = Ti / growth * (gsi / 10.75)**(1 / 3)
            rows.append(dict(variant='B2: dust-like reading from the actual formation time', T_form_GeV=Ti,
                             a_ratio=growth, T_ann_GeV=T))
        for Ck in (1.0,):
            ratio = 3 * Ck * (1 / epsr)**(1 / 3)       # tau_ann/tau_i, their eq. 26, literal transfer
            T = Ti / ratio * (gsi / 3.931)**(1 / 3)
            rows.append(dict(variant='B3: literal eq.26 with tau_i at the actual formation time', T_form_GeV=Ti,
                             tau_ratio=ratio, T_ann_GeV=T))
    return dict(eps_over_lam_v=epsr, simulated_range_eps_over_lam_v=[0.0035, 0.05],
                extrapolation_decades=float(np.log10(0.0035 / epsr)), rows=rows)


def widths(kappa=9.6e-10, Cgamma=10.0):
    bX = np.sqrt(1 - 4 * MX**2 / MS**2)
    bD = np.sqrt(max(1 - 4 * MD**2 / MS**2, 0))
    GX_ = GX**4 * V**2 * bX / (4 * np.pi * MS)
    GD = HD**4 * V**2 * bD / (8 * np.pi * MS)
    vH, mh, mmu, me, alpha = 246.22, 125.25, 0.10566, 0.000511, 1 / 137.036
    th = kappa * V * vH / (mh**2 - MS**2)
    lep = sum(m * m * MS / (8 * np.pi * vH**2) * (1 - 4 * m * m / MS**2)**1.5 for m in (me, mmu) if MS > 2 * m)
    gam = Cgamma**2 * alpha**2 * MS**3 / (256 * np.pi**3 * vH**2)
    Gv = th**2 * (lep + gam)
    tot = GX_ + GD + Gv
    return dict(Gamma_XX_GeV=GX_, Gamma_DD_GeV=GD, Gamma_visible_GeV=Gv, theta=th, lifetime_s=HBAR / tot,
                BR_XX=GX_ / tot, BR_DD=GD / tot, BR_visible=Gv / tot, beta_X=bX,
                X_momentum_from_s_at_rest_GeV=np.sqrt(MS**2 / 4 - MX**2))


def injection(A=0.8, C=5.0, Tform=86.0):
    r = standard(A, C)
    T = r['T_ann_GeV']
    t = C * A * SIGMA / DELTA
    rho_term = A * SIGMA / t + DELTA / 2          # walls + released false-vacuum energy (half the volume)
    s = entropy(T)
    rad = np.pi**2 / 30 * gstar(T)[0] * T**4
    # network radiation integrated from formation to annihilation: d(E a^3)/dt ~ H rho_w a^3, rho_w = A sigma/t
    # pairs per entropy dY = rho_w H dt / (E_s s); E_s >= m_s gives an upper bound. Use ln T as variable.
    def integrand(lnT):
        TT = np.exp(lnT)
        H = H_of_T(TT)
        rho_w = 2 * A * SIGMA * H
        return rho_w / (MS * entropy(TT))           # per d ln a ~ -d ln T
    lo, hi = np.log(T), np.log(min(Tform, 0.149))
    Y_net = quad(integrand, lo, hi, limit=200)[0]
    # above 149 MeV use fixed g*=86.25 scaling ~ 1/T (negligible contribution, bounded analytically)
    Y_term = rho_term / (MS * s)
    Y_pairs = Y_net + Y_term
    budget = 0.01 * BDM / (2 * MX)
    return dict(A=A, C_ann=C, T_ann_GeV=T, rho_terminal_GeV4=rho_term, terminal_over_radiation=rho_term / rad,
                terminal_energy_per_entropy_GeV=rho_term / s, terminal_over_calibrated_DM=rho_term / s / BDM,
                pair_yield_upper_terminal=Y_term, pair_yield_upper_network=Y_net, pair_yield_upper_total=Y_pairs,
                codex_1pct_pair_budget=budget, fraction_of_1pct_budget=Y_pairs / budget,
                pair_rest_mass_over_DM=2 * MX * Y_pairs / BDM)


def loading(A=0.8, vw=0.4):
    """Q-ball population vs wall network. rho_DM = BDM * s (objects dominate late matter)."""
    rows = []
    Tann = standard(0.8, 5.0)['T_ann_GeV']
    for T in (86.0, 40.0, 1.0, 0.1, Tann):
        if T < 0.149:
            H = H_of_T(T); s = entropy(T)
        else:
            g = 86.25 if T > 5 else 61.75
            H = np.pi * np.sqrt(g / 90) * T**2 / MRED; s = 2 * np.pi**2 / 45 * g * T**3
        t = 1 / (2 * H)
        rho_dm = BDM * s
        rho_w = A * SIGMA / t
        nQ = YQ * s / QREF
        rows.append(dict(T_GeV=T, rho_DM_GeV4=rho_dm, rho_wall_scaling_GeV4=rho_w, rho_DM_over_rho_wall=rho_dm / rho_w,
                         swept_DM_per_area_per_Hubble_over_sigma=rho_dm * vw * t / SIGMA,
                         capture_probability_threshold=SIGMA / (rho_dm * vw * t),
                         Qball_number_density_GeV3=nQ, Qball_spacing_GeVinv=nQ**(-1 / 3),
                         Qball_induced_bias_2dE_nQ_over_Delta=2 * 1.738e-17 * nQ / DELTA))
    t_i = 1 / (2 * np.pi * np.sqrt(86.25 / 90) * 86.0**2 / MRED)
    t_a = 0.8 * 5.0 * SIGMA / DELTA
    crossings = A * vw * np.log(t_a / t_i)
    KE = MQ * (1 / np.sqrt(1 - vw**2) - 1)
    vcap = {f'{f}': float(np.sqrt(2 * f * U_WALL_CORE / MQ)) for f in (0.01, 0.1, 1.0)}
    # snowplow limit: fully loaded walls move x ~ t sqrt(Delta/rho_DM); domains of size ~t clear when rho_DM ~ Delta
    s_star = DELTA / BDM
    T_star = (s_star / (2 * np.pi**2 / 45 * 3.931))**(1 / 3)
    # wall-bound X mode (Poschl-Teller): -psi'' + g^2 v^2 tanh^2(z/w) psi = mu^2 psi
    Wr = np.sqrt(2 * GX**2 / LAM)
    nu = (-1 + np.sqrt(1 + 4 * Wr**2)) / 2
    mu0 = np.sqrt(1 - nu**2 / Wr**2)
    # numerical check by finite differences
    zz = np.linspace(-40, 40, 8001); hz = zz[1] - zz[0]
    pot = np.tanh(zz / Wr)**2
    main = 2 / hz**2 + pot; off = -np.ones(len(zz) - 1) / hz**2
    from scipy.linalg import eigh_tridiagonal
    ev = eigh_tridiagonal(main, off, select='i', select_range=(0, 1))[0]
    return dict(rows=rows, wall_crossings_per_Qball_formation_to_annihilation=crossings,
                wall_frame_KE_at_vw_GeV=KE, wall_core_attraction_GeV=U_WALL_CORE,
                capture_speed_for_inelastic_fraction=vcap,
                snowplow_limit_T_GeV=T_star,
                wall_bound_X_gap_over_mX_analytic=mu0, wall_bound_X_gap_over_mX_numeric=float(np.sqrt(ev[0])),
                Qball_Omega=OMEGA_Q, gap_below_Omega=bool(mu0 < OMEGA_Q),
                energy_gain_per_charge_moved_to_wall_mode_GeV=(OMEGA_Q - mu0) * MX,
                bound_mode_count=int(np.floor(nu) + 1))


def main(outp):
    A_list = [0.5, 0.7, 0.8, 0.9, 2.0]
    C_list = [1.0, 2.0, 3.0, 5.0, 10.0]
    grid = [standard(A, C) for A in A_list for C in C_list]
    central = [r for r in grid if r['A'] in (0.7, 0.8, 0.9) and r['C_ann'] in (2.0, 3.0, 5.0)]
    codex_rows = [dict(A=A, C=C, g=g, T=T_fixed_g(A, C, g)) for A in (0.5, 0.8, 2.0) for C in (1.0, 2.0, 10.0) for g in (10.75, 61.75)]
    Hmin = H_of_T(1e-3)
    dmin = {f'A{A}_C{C}': 2 * C * A * SIGMA * Hmin for A, C in [(0.8, 5.0), (2.0, 10.0), (1.0, 1.0)]}
    dom = 4 * 5.0 * 0.8**2 * SIGMA**2 / (3 * MRED**2)       # Saikawa (2.29) with C=5, A=0.8
    gws = [gw_today(A, C, e) for A, C in [(0.8, 5.0), (0.8, 2.0), (2.0, 10.0), (0.5, 1.0)] for e in (0.3, 0.7, 1.1)]
    # PTA: NANOGrav 15-yr, A=2.4e-15 at f=1/yr, gamma=13/3 => Omega_gw h^2(f) = 2 pi^2/(3 H100^2) f^2 h_c^2
    H100 = 100 / 3.0857e19
    fyr = 1 / 3.15576e7
    def om_pta(f):
        hc = 2.4e-15 * (f / fyr)**(-2 / 3)
        return 2 * np.pi**2 / (3 * H100**2) * f**2 * hc**2
    best = gw_today(0.8, 5.0, 0.7)
    out = dict(inputs=dict(sigma_GeV3=SIGMA, Delta_GeV4=DELTA, m_s=MS, m_X=MX, m_D=MD),
               standard_grid=grid,
               standard_central_T_range_GeV=[min(r['T_ann_GeV'] for r in central), max(r['T_ann_GeV'] for r in central)],
               standard_full_T_range_GeV=[min(r['T_ann_GeV'] for r in grid), max(r['T_ann_GeV'] for r in grid)],
               codex_formula_fixed_g_reproduction=dict(min_T_GeV=min(r['T'] for r in codex_rows), rows=codex_rows),
               Delta_min_for_T_ann_1MeV_GeV4=dmin, Delta_over_Delta_min=({k: DELTA / val for k, val in dmin.items()}),
               domination_bound_Saikawa_2_29_GeV4=dom,
               babichev=babichev_variants(),
               gw=gws, pta_Omega_h2_at_f_peak=om_pta(best['f_peak_Hz']), pta_Omega_h2_at_10nHz=om_pta(1e-8),
               lvk_O3_isotropic_bound_flat=5.8e-9,
               s_decays=widths(), injection=[injection(0.8, 5.0), injection(2.0, 10.0), injection(0.5, 1.0)],
               loading=loading())
    Path(outp).write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(json.dumps({k: out[k] for k in ['standard_central_T_range_GeV', 'standard_full_T_range_GeV',
                                          'Delta_over_Delta_min', 'domination_bound_Saikawa_2_29_GeV4',
                                          'pta_Omega_h2_at_f_peak']}, indent=1, default=float))
    print('codex min T (fixed g)', out['codex_formula_fixed_g_reproduction']['min_T_GeV'])
    print(json.dumps(out['babichev'], indent=1, default=float))
    print(json.dumps(gws[:6], indent=1, default=float))
    print(json.dumps(out['s_decays'], indent=1, default=float))
    print(json.dumps(out['injection'], indent=1, default=float))
    print(json.dumps(out['loading'], indent=1, default=float))


if __name__ == '__main__':
    main(sys.argv[1])
