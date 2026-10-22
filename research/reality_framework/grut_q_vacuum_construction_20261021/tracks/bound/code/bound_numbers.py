"""Track `bound`: the class-wide numbers of NOTES section A.4-A.6 and B.5, with every input and its provenance.
Writes results/bound_numbers.json. Pure numpy; runs in < 1 s.  Run: python3 bound_numbers.py

Eot-Wash 2020 digitization. The alpha-lambda exclusion curve of Lee et al., PRL 124, 101101 (2020), arXiv:2002.11761,
Fig. 5 bottom (95% CL), was read from the vector paths of the arXiv source file fig5b1.pdf (kept in the scratchpad, not in
the repository): frame 2e-6..1e-3 m (x, log), 1e-3..1e6 (|alpha|, log). The (lambda, alpha) vertices nearest to the
crossings are embedded below. Calibration check: the same procedure gives the alpha = 1 crossing at 38.35 um (paper: 38.6 um)
and, for the Eot-Wash 2007 curve, 55.2 um (Kapner et al. 2007 abstract: 56 um).
"""
import json
import os
import numpy as np

HBARC_EV_M = 1.973269804e-7          # [R] CODATA
HBAR_EV_S = 6.582119569e-16          # [R] CODATA
MPC_KM = 3.0856775814913673e19       # [R]
M_RED_GEV = 2.435e18                 # BRIEF.md
RHO_L_GEV4 = 2.52e-47                # BRIEF.md (Planck 2018 inputs)
H0_KMSMPC = 67.4                     # Planck 2018 [Q via equations/ledger tracks]
OMEGA_M = 0.315                      # Planck 2018, ledger row 9 [Q]
OMEGA_L = 0.685
R_SUN_M = 6.957e8                    # IAU nominal solar radius [R]

# --- digitized vertices (lambda [m], |alpha|), Lee et al. 2020 Fig. 5 (bottom), 95% CL
EW2020 = {'alpha1': [(3.7210e-05, 1.1208), (3.8484e-05, 0.98945)],
          'alpha13': [(5.1642e-05, 0.37337), (5.3604e-05, 0.33469), (5.5699e-05, 0.30004)]}
EW2007 = {'alpha1': [(5.4697e-05, 1.0401), (5.6307e-05, 0.92339)],
          'alpha13': [(7.2046e-05, 0.34913), (7.4551e-05, 0.30702)]}
HUST2020 = {'alpha1': [(4.7120e-05, 1.0745), (4.8331e-05, 0.91462)],          # black curve labelled HUST 2012/16/20
            'alpha13': [(5.6074e-05, 0.37267), (5.7546e-05, 0.32087)]}


def crossing(pts, alpha):
    """log-log interpolation of lambda at given alpha along the digitized vertices."""
    pts = sorted(pts, key=lambda p: -p[1])
    for (l1, a1), (l2, a2) in zip(pts[:-1], pts[1:]):
        if a1 >= alpha >= a2:
            t = (np.log(alpha) - np.log(a1)) / (np.log(a2) - np.log(a1))
            return float(np.exp(np.log(l1) + t * (np.log(l2) - np.log(l1))))
    raise ValueError


def main():
    out = {}
    H0_s = H0_KMSMPC / MPC_KM
    H0_eV = H0_s * HBAR_EV_S
    out['inputs'] = {
        'H0 [eV]': H0_eV, 'H0 source': 'Planck 2018 H0 = 67.4 km/s/Mpc [Q via ledger]; hbar [R]',
        'Omega_m': OMEGA_M, 'rho_Lambda [GeV^4]': RHO_L_GEV4, 'M reduced [GeV]': M_RED_GEV,
        'alpha (unscreened q-mode)': 1.0 / 3.0, 'alpha source': 'equations track 4.3 [D]: G_eff = G*(1 + e^{-mr}/3)',
    }
    # ----------------------------------------------- Eot-Wash
    l1_2020 = crossing(EW2020['alpha1'], 1.0)
    l13_2020 = crossing(EW2020['alpha13'], 1.0 / 3.0)
    l1_2007 = crossing(EW2007['alpha1'], 1.0)
    l13_2007 = crossing(EW2007['alpha13'], 1.0 / 3.0)
    l1_hust = crossing(HUST2020['alpha1'], 1.0)
    l13_hust = crossing(HUST2020['alpha13'], 1.0 / 3.0)
    cal = 38.6e-6 / l1_2020                       # rescale to the paper's stated alpha = 1 value
    lam13 = l13_2020 * cal
    dig_unc = max(abs(cal - 1), abs(56e-6 / l1_2007 - 1))
    m_min = HBARC_EV_M / lam13
    out['EotWash'] = {
        'source': 'Lee, Adelberger, Cook, Fleischer, Heckel, PRL 124, 101101 (2020), arXiv:2002.11761, Fig. 5 bottom (95% CL)',
        'read level': 'figure: vector data digitized from the arXiv source figure; text: abstract + the constraint paragraph',
        'alpha=1 crossing digitized [um]': l1_2020 * 1e6, 'alpha=1 stated [um] [Q]': 38.6,
        'EW2007 alpha=1 digitized [um]': l1_2007 * 1e6, 'EW2007 alpha=1 stated (Kapner+ 2007 abstract) [um] [Q]': 56.0,
        'alpha=1/3 crossing digitized [um]': l13_2020 * 1e6,
        'alpha=1/3 crossing, rescaled to the stated alpha=1 point [um]': lam13 * 1e6,
        'digitization uncertainty (relative)': dig_unc,
        'EW2007 alpha=1/3 [um]': l13_2007 * 1e6, 'HUST alpha=1 [um]': l1_hust * 1e6, 'HUST alpha=1/3 [um]': l13_hust * 1e6,
        'strongest at alpha = 1/3': 'Eot-Wash 2020',
        'm_q lower bound [meV] = hbar c / lambda(1/3)': m_min * 1e3,
        'm_q lower bound range [meV]': [HBARC_EV_M / (lam13 * (1 + dig_unc)) * 1e3, HBARC_EV_M / (lam13 * (1 - dig_unc)) * 1e3],
        'longer ranges': ('alpha = 1/3 is excluded for every lambda from ~54 um to ~1 cm (Eot-Wash 2004/2007, Irvine: '
                          '|alpha| < 0.1; Adelberger-Heckel-Nelson 2003 Fig. midlam), from ~4 mm to ~3e14 m (laboratory, '
                          'geophysical, LAGEOS, LLR, planetary: |alpha| < 0.1; same review Fig. longlam, figure read), and for '
                          'lambda >> 1 AU by Cassini gamma (below). Read level: figures rendered from the arXiv source.'),
    }
    # ----------------------------------------------- Cassini
    gm1, sg = 2.1e-5, 2.3e-5                         # Bertotti, Iess, Tortora 2003 [Q via ledger row 7]
    lower = gm1 - 2 * sg                              # 2 sigma lower edge
    # gamma(r) - 1 = -(2/3) e^{-x}/(1 + e^{-x}/3) >= lower
    e = (-lower) / (2.0 / 3.0 + (-lower) / 3.0)
    xmin = -np.log(e)
    out['Cassini'] = {'gamma-1 [Q]': [gm1, sg], '2 sigma lower edge': lower, 'required m b': xmin,
                      'm_q > [eV] (b = R_sun)': xmin * HBARC_EV_M / R_SUN_M,
                      'light q (m b << 1): gamma - 1 = -1/2, excess in sigma': 0.5 / sg}
    # ----------------------------------------------- regime (i) relation and its maximum
    def omega_dyn(m_eV):
        return 2.25 * OMEGA_M**2 * (H0_eV / m_eV)**2
    rho_c = RHO_L_GEV4 / OMEGA_L
    out['regime_i'] = {
        'relation': 'Omega_dyn(t0) = (9/4) Omega_m^2 (H0/m_q)^2 [1 + O(c (H0/m_q)^2)], w_dyn = +1 (stiff, ~ a^-6); '
                    'rho_DE(t0) = Lam_cal(mu) + Omega_dyn rho_c, Lam_cal = -q0 (mu - mu0) + O((mu - mu0)^2)',
        'H0/m at Eot-Wash bound': H0_eV / m_min,
        'Omega_dyn max (unscreened, Eot-Wash alpha=1/3 bound)': omega_dyn(m_min),
        'rho_dyn max [GeV^4]': omega_dyn(m_min) * rho_c,
        'Omega_dyn at the Cassini-only bound (2.9e-15 eV)': omega_dyn(xmin * HBARC_EV_M / R_SUN_M),
        'Omega_dyn at natural m_q = sqrt(2/3) M': omega_dyn(np.sqrt(2.0 / 3.0) * M_RED_GEV * 1e9),
        '1 + w_DE(t0) max = 2 Omega_dyn/Omega_L': 2 * omega_dyn(m_min) / OMEGA_L,
        'stage-1 rho_V = rho_m^2/(12 K m^2) vs full effective (3/4) rho_m^2/(K m^2)': 'factor 9 (the -3H Kdot and G-shift terms)',
    }
    # ----------------------------------------------- regime (ii) mass requirement and gap
    H_eq_over_H0 = np.sqrt(2 * OMEGA_M * 3400.0**3)          # H at matter-radiation equality (z_eq ~ 3400) [R]
    out['regime_ii'] = {
        'm_q needed to be Hubble-frozen today': 'm_q <~ H0 = %.3e eV' % H0_eV,
        'm_q below which the dust trace drives K (phi-MDE, a ~ t^(1/2); APT07)': 'm_q <~ H(z_eq) ~ %.2e H0 = %.2e eV' % (H_eq_over_H0, H_eq_over_H0 * H0_eV),
        'gap to the Eot-Wash bound m_min/H0': m_min / H0_eV,
        'gap to the Cassini-only bound': xmin * HBARC_EV_M / R_SUN_M / H0_eV,
        'G_cos/G_lab for light q (= 3/4) vs Lamine+ 2025 G/G_N = 1.011 +- 0.018, sigma': (1.011 - 0.75) / 0.018,
    }
    # ----------------------------------------------- regime (iii) relic estimate [E]
    gstar = 106.75
    T_osc = (90 / (np.pi**2 * gstar))**0.25 * np.sqrt(m_min * 1e-9 * M_RED_GEV)       # GeV, H(T) = m
    alpha_s = 0.085                                                                     # alpha_s(~1.6 TeV) [R]
    kappa_qcd = (2940.0 / 360.0) * alpha_s**2 / (np.pi**2 * gstar / 30)                # leading trace anomaly [E]
    T0 = 2.7255 * 8.617333e-14                                                          # GeV [R]
    growth = (T_osc / T0) * (3.91 / gstar)**(1 / 3)
    Om_r = 9.1e-5
    rows = {}
    for label, L2 in [('drive over 1 e-fold of H (ln^2 = 1)', 1.0), ('drive from T = 1e15 GeV: ln(H_i/m) = 2 ln(T_i/T_osc)', (2 * np.log(1e15 / T_osc))**2)]:
        f_on = 0.1 * kappa_qcd**2 * L2
        rows[label] = {'relic fraction at onset': f_on, 'Omega_relic today (if stable)': f_on * growth * Om_r}
    tau_lim_m = (192 * np.pi * (M_RED_GEV**2) * (HBAR_EV_S * 1e-9))**(1 / 3)          # m with tau = 1 s: m^3/(192 pi M^2) = hbar/1s
    out['regime_iii_estimate'] = {
        'label': '[E] order-of-magnitude; coefficient 0.1 from frw_regimes.json part III (constant-kappa toy)',
        'T_osc at m = m_min [GeV]': T_osc, 'kappa_SM(T_osc) from the QCD trace anomaly [E]': kappa_qcd,
        'growth factor relic/radiation from T_osc to today': growth, 'cases': rows,
        'Omega_c (Planck 2018, h = 0.674) [Q]': 0.120 / 0.674**2,
        'm_q above which the relic decays before 1 s (Gamma ~ m^3/(192 pi M^2), one scalar channel) [E] [GeV]': tau_lim_m,
    }
    # ----------------------------------------------- B: threshold stability of a disappearing-Lambda structure
    D_EW = 1.19e8                                  # GeV^4, equations track 7 (tree level)
    K = M_RED_GEV**2
    R_dS = 4 * RHO_L_GEV4 / K
    out['B_threshold'] = {
        'dR = 2 c/(q0 K\') with natural q0 K\' = K, c = Delta_EW [GeV^2]': 2 * D_EW / K,
        'R_dS = 4 rho_L/K [GeV^2]': R_dS,
        'dR/R_dS (natural linear-K realization)': 2 * D_EW / K / R_dS,
        'required q0 K\'/K across the feature (dR < R_dS)': 2 * D_EW / (K * R_dS),
    }
    # ----------------------------------------------- B: Appleby-Battye Minkowski vacuum, G_Mink/G = 1 + e^{2b}
    R0_planck = 12 * (H0_eV / (M_RED_GEV * 1e9))**2            # R0 = 12 H0^2 in reduced-Planck units
    # (i) as written in 0705.3199: 8 b e^{-2b} << 2 R0  (taken with M = reduced Planck mass [A])
    b1 = 1.0
    for _ in range(50):
        b1 = 0.5 * np.log(8 * b1 / (2 * R0_planck))
    # (ii) physical reading F''(R0) ~ (4b/R0) e^{-2b} << x^2 with x = 1 mm
    x_eVinv = 1e-3 / HBARC_EV_M
    R0_eV2 = 12 * H0_eV**2
    b2 = 1.0
    for _ in range(50):
        b2 = 0.5 * np.log(4 * b2 / (R0_eV2 * x_eVinv**2))
    out['B_AppleByBattye_Minkowski'] = {
        'G_Mink/G = 1/F\'(0) = 1 + e^{2b} (F\'(R_1) ~ 1)': 'derived from 0705.3199 eq (ab3)',
        'b at equality of their condition 8 b e^{-2b} = 2 R0, R0 = 12 H0^2/M^2': b1, 'G_Mink/G (>=)': float(np.exp(2 * b1)),
        'b at equality of F\'\'(R0) = (1 mm)^2': b2, 'G_Mink/G (>=), lab reading': float(np.exp(2 * b2)),
    }
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, '..', 'results', 'bound_numbers.json'), 'w') as fh:
        json.dump(out, fh, indent=1, default=float)
    for k, v in out.items():
        print(k, json.dumps(v, indent=1, default=float))


if __name__ == '__main__':
    main()
