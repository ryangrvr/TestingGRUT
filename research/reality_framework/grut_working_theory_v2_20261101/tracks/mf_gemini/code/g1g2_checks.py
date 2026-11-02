"""mf_gemini: checks of G1 and G2 (Claude's estimates in BRIEF_MULTIFIELD.md), own numbers.

G1: energy-amplitude identity, kinematic blocking, one-loop potential of an O(1) Yukawa.
G2: the proposal as written: v_chi, lam_chi, Coleman-Weinberg vs tree, blocking ratio, the
    tachyonic rate, the Higgs portal route (width, Higgs mass in inflation, tree threshold,
    Landau pole, tachyonic windows), fermionic and bosonic preheating with radiatively
    allowed couplings.

Physical units GeV; reduced M_P.  Inputs read from results/waterfall_bg.json when present
(the post-waterfall oscillation amplitude); otherwise a stated placeholder.

Run: python3 -B g1g2_checks.py -> results/g1g2.json
"""
import json
import os
import numpy as np
from scipy.optimize import minimize_scalar, brentq

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
MP = 2.435e18


def cw_residual(y, lam, v, Nf=1, shifted=False, nchi=4001):
    """max |V_CW - best counterterms| / V0 over chi in [0, 1.5 v].

    V_CW = -(Nf/32 pi^2) m_N^4 [ln(m_N^2/mu^2) - 3/2] for Nf Majorana fermions (2 dof each),
    m_N = y chi (shifted=False) or y (chi - v) (shifted: bare mass tuned so N is light in vacuum).
    Counterterms: {1, chi^2, chi^4} (Z2-even) or {1, chi, chi^2, chi^3, chi^4} (shifted).  The fit is a
    least-squares fit relative to V0 on the grid; mu-dependence is absorbed by the chi^4 term.
    """
    x = np.linspace(0, 1.5, nchi)          # chi/v
    mN = y * (x - 1.0) if shifted else y * x
    with np.errstate(divide='ignore', invalid='ignore'):
        f = np.where(np.abs(mN) > 0, mN**4 * np.log(mN**2), 0.0)   # in units v^4 (mu = v)
    Vcw = -(Nf / (32 * np.pi**2)) * f * v**4
    V0 = lam * v**4 / 4
    pw = [0, 1, 2, 3, 4] if shifted else [0, 2, 4]
    A = np.vstack([x**p for p in pw]).T * v**4
    coef, *_ = np.linalg.lstsq(A / V0, Vcw / V0, rcond=None)
    res = Vcw / V0 - (A / V0) @ coef
    return float(np.max(np.abs(res)))


def main():
    out = {}
    Hs = [1.25e11, 1.30e11, 1.55e11]
    # ---------------- G1 ----------------
    g1 = {}
    for b in [2.0, 3.0, 4.0]:
        g1[f'm={b}H'] = dict(Phi_over_MP=float(np.sqrt(6.0) / b))
    H = 1.3e11
    b = 3.0
    m = b * H
    Phi = np.sqrt(6.0) * H * MP / m
    # Majorana -(y/2) phi NN: Gamma = y^2 m beta^3/16 pi (inf_reheat R4)
    y_half = np.sqrt(16 * np.pi * 0.5 * H / m)
    g1['y_for_Gamma_0.5H_Majorana_Nf1'] = float(y_half)
    g1['y_for_Gamma_0.5H_Dirac_Nf1'] = float(np.sqrt(8 * np.pi * 0.5 * H / m))
    g1['blocking_yPhi_over_m2_at_that_y'] = float(y_half * Phi / (m / 2))
    lam_I = m**2 / Phi**2 * 0  # quadratic oscillator: tree V = m^2 phi^2/2
    Vtree = 0.5 * m**2 * Phi**2
    Vcw = y_half**4 * Phi**4 / (32 * np.pi**2)
    g1['CW_over_tree_at_Phi'] = float(Vcw / Vtree)
    g1['label'] = 'D (identity) / C (ratios)'
    out['G1'] = g1

    # ---------------- G2 ----------------
    g2 = {}
    for Hx in Hs:
        hub = Hx / MP
        v = np.sqrt(8.0 / 3.0) * MP
        lam = 27.0 / 16.0 * hub**2
        mchi = 3 * Hx
        d = dict(v_over_MP=float(v / MP), lam=float(lam), m_chi_GeV=float(mchi),
                 V0_quarter_GeV=float((3 * MP**2 * Hx**2)**0.25))
        for y in [1.5, 2.0]:
            d[f'y={y}'] = dict(
                Gamma_over_H_Majorana_unblocked=float(y**2 * mchi / (16 * np.pi) / Hx),
                Gamma_over_H_Dirac_unblocked=float(y**2 * mchi / (8 * np.pi) / Hx),
                naive_y4_16pi2_over_lam=float(y**4 / (16 * np.pi**2) / lam),
                CW_residual_over_V0_Z2=cw_residual(y, lam, v),
                CW_residual_over_V0_shifted=cw_residual(y, lam, v, shifted=True),
                MN_vac_over_mchi_half=float(y * v / (mchi / 2)))
        # radiative cap: residual <= 1% of V0
        def f(lny):
            return np.log(cw_residual(np.exp(lny), lam, v)) - np.log(0.01)
        ycap = float(np.exp(brentq(f, np.log(1e-6), np.log(1.0))))
        def fs(lny):
            return np.log(cw_residual(np.exp(lny), lam, v, shifted=True)) - np.log(0.01)
        ycap_s = float(np.exp(brentq(fs, np.log(1e-6), np.log(1.0))))
        d['y_cap_CW_1pct'] = ycap
        d['y_cap_CW_1pct_shifted'] = ycap_s
        d['Gamma_over_H_at_ycap_Majorana_unblocked'] = float(ycap**2 * mchi / (16 * np.pi) / Hx)
        d['tachyonic_mass2_at_chi0_over_H2'] = -4.5
        s = (-3 + np.sqrt(9 + 4 * 4.5)) / 2
        d['tachyonic_growth_rate_per_efold_max'] = float(s)
        d['efolds_H_over_2pi_to_v_at_max_rate'] = float(np.log(v / (Hx / (2 * np.pi))) / s)
        g2[f'H={Hx:.3g}'] = d
    out['G2'] = g2

    # post-waterfall oscillation amplitude about v (from waterfall_bg.json, else placeholder)
    amp = None
    try:
        with open(os.path.join(RES, 'waterfall_bg.json')) as fh:
            bg = json.load(fh)
        amp = bg['representative']['osc_amp_over_v_first']
        amp_src = 'results/waterfall_bg.json (first swing after the end of inflation)'
    except Exception:
        amp = 0.5
        amp_src = 'placeholder 0.5 (waterfall_bg.json missing)'
    out['osc_amp_over_v'] = dict(value=amp, source=amp_src)

    # ---------------- preheating with radiatively allowed couplings ----------------
    pre = {}
    H = 1.3e11
    hub = H / MP
    v = np.sqrt(8.0 / 3.0) * MP
    lam = 27.0 / 16.0 * hub**2
    mchi = 3 * H
    Phi = amp * v
    V0 = 3 * MP**2 * H**2
    for lab, y in [('y=2 (proposal)', 2.0), ('y=1.5 (proposal)', 1.5),
                   ('y=cap (shifted, CW<=1%)', g2['H=1.3e+11']['y_cap_CW_1pct_shifted'])]:
        # Greene-Kofman first-crossing energy fraction (inf_reheat), N mass y(chi - v) tuned
        f1 = y**2.5 * (Phi / mchi)**0.5 / (2 * np.pi**3)
        pre[lab] = dict(y=float(y), GK_energy_fraction_per_crossing=float(f1),
                        mN_at_turning_point_over_MP=float(y * Phi / MP),
                        blocking_ratio_yPhi_over_m2=float(y * Phi / (mchi / 2)))
    # bosonic instant preheating, 1/2 g^2 (chi - v)^2 phi^2 (added field phi), CW cap:
    # V_CW = (g^2 (chi-v)^2)^2/(64 pi^2) ln -> residual ~ g^4 v^4/(64 pi^2) * O(0.1) <= 0.01 V0
    def bos_res(g2c):
        x = np.linspace(0, 1.5, 4001)
        mm = g2c * (x - 1.0)**2
        with np.errstate(divide='ignore', invalid='ignore'):
            ff = np.where(mm > 0, mm**2 * np.log(mm), 0.0)
        Vb = ff / (64 * np.pi**2) * v**4
        A = np.vstack([x**p for p in [0, 1, 2, 3, 4]]).T * v**4
        c, *_ = np.linalg.lstsq(A / (lam * v**4 / 4), Vb / (lam * v**4 / 4), rcond=None)
        return float(np.max(np.abs(Vb / (lam * v**4 / 4) - (A / (lam * v**4 / 4)) @ c)))
    g2cap = float(np.exp(brentq(lambda l: np.log(bos_res(np.exp(l))) - np.log(0.01), np.log(1e-12), np.log(1.0))))
    gb = np.sqrt(g2cap)
    q = g2cap * Phi**2 / (4 * mchi**2)
    # first-crossing (Kofman-Linde-Starobinsky) occupation n_k = exp(-pi k^2/(g Phi m)) ->
    # n = (g Phi m)^{3/2}/(8 pi^3); energy at turning point ~ n g Phi
    nb = (gb * Phi * mchi)**1.5 / (8 * np.pi**3)
    fb = nb * gb * Phi / (0.5 * mchi**2 * Phi**2)
    pre['bosonic (chi-v)^2 phi^2, CW cap'] = dict(g2=g2cap, q=float(q), first_crossing_energy_fraction=float(fb),
                                                 note='crossings of chi = v: ~1 per H^-1 at m = 3H')
    out['preheating_caps'] = pre

    # ---------------- Higgs portal route ----------------
    hp = {}
    for Gam in [0.5, 1.0]:
        A = np.sqrt(8 * np.pi * mchi * Gam * H)       # Gamma(chi -> H H^dag) = A^2/(8 pi m), 4 real dof
        lchiH = A / (2 * v)
        thr = lchiH**2 / lam                            # tree threshold: lam_H(full) = lam_SM + lchiH^2/lam
        lamH_full = 0.0 + thr                           # lam_SM(m_chi ~ 4e11 GeV) ~ 0 (near instability scale) [R]
        # one-loop Landau pole from 24 lam^2/16pi^2 (dominant at large lam)
        tLP = 16 * np.pi**2 / (24 * lamH_full)
        mH2_infl = lchiH * v**2                         # |m_H^2| at chi = 0 if the EW vacuum is tuned
        # fraction of V0 released by the Higgs at chi = 0 (lchiH > 0): (lchiH v^2)^2/(4 lam_H) / V0, |H|^4 normalisation
        dV = (lchiH * v**2)**2 / (4 * lamH_full) / V0
        # tachyonic window (lchiH < 0) during chi > v: m_H^2 = -|l|(chi^2 - v^2), chi = v + Phi cos(m t)
        tt = np.linspace(-np.pi / 2, np.pi / 2, 20001)
        chi = v + Phi * np.cos(tt)
        mneg = np.sqrt(np.clip(lchiH * (chi**2 - v**2), 0, None))
        expo = np.trapezoid(mneg, tt) / mchi
        # CW from the Higgs loop on the chi potential: 4 real dof, m^2 = lchiH chi^2
        cwH = 4 * (lchiH)**2 / (64 * np.pi**2) / (lam / 4)
        hp[f'Gamma={Gam}H'] = dict(A_over_H=float(A / H), lam_chiH=float(lchiH), threshold_lchiH2_over_lam=float(thr),
                                   lamH_full=float(lamH_full), LandauPole_efolds_of_scale=float(tLP),
                                   LandauPole_GeV=float(mchi * np.exp(tLP)),
                                   abs_mH2_inflation_over_H2=float(mH2_infl / H**2),
                                   Higgs_energy_release_at_chi0_over_V0=float(dV),
                                   tachyonic_window_exponent_per_half_cycle=float(expo),
                                   Higgs_CW_quartic_over_lam_quarter=float(cwH),
                                   mH_osc_amplitude_over_mchi_half=float(np.sqrt(2 * lchiH * v * Phi) / (mchi / 2)))
    # safe portal: |m_H^2| from the portal below H^2 at chi = 0 (inflation) and in the windows
    l_safe_infl = H**2 / v**2
    l_safe_win = H**2 / (2 * v * Phi)
    for lab, l in [('safe_inflation', l_safe_infl), ('safe_windows', l_safe_win)]:
        A = 2 * l * v
        hp[lab] = dict(lam_chiH=float(l), Gamma_over_H=float(A**2 / (8 * np.pi * mchi) / H))
    hp['note'] = ('lam_SM at 4e11 GeV taken as 0 (it is near its zero there for central M_t); the threshold '
                  'term dominates.  Boundedness of the full potential with lam_chiH < 0 requires lam_SM > 0 '
                  '[D: lam_chiH^2 < lam_chi lam_H(full) <=> lam_SM > 0].')
    out['Higgs_portal'] = hp

    # ---------------- G3 check (heavy waterfall, perturbative) ----------------
    g3 = {}
    for Hx in [1.25e11, 1.3e11, 1.55e11]:
        hub = Hx / MP
        # Gamma = y^2 m/16pi (Majorana) with y^2 <= lam/2 (M_N = y v < m/2): Gamma_max = lam m/32 pi
        # V0 = lam v^4/4 = 3 MP^2 H^2, m = sqrt(2 lam) v -> m/H = sqrt(2) 12^{1/4} lam^{1/4} (MP/H)^{1/2}
        def gam(lam):
            mm = np.sqrt(2) * 12**0.25 * lam**0.25 / np.sqrt(hub)
            return lam * mm / (32 * np.pi)
        lmin = brentq(lambda l: gam(l) - 1.0, 1e-8, 1.0)
        lmin05 = brentq(lambda l: gam(l) - 0.5, 1e-8, 1.0)
        for lab, l in [('Gamma=H', lmin), ('Gamma=0.5H', lmin05)]:
            mm = np.sqrt(2) * 12**0.25 * l**0.25 / np.sqrt(hub)
            vv = (12 / l)**0.25 * np.sqrt(hub)
            y = np.sqrt(l / 2)
            g3[f'H={Hx:.3g},{lab}'] = dict(lam=float(l), m_over_H=float(mm), v_over_MP=float(vv), y_max=float(y),
                                          y4_16pi2_over_lam=float(y**4 / (16 * np.pi**2) / l))
    out['G3_check'] = g3
    with open(os.path.join(RES, 'g1g2.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
