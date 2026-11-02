#!/usr/bin/env python3
"""inf_reheat items 1-2 (Claude, own code): reheating channels for M1r, family-independent numbers.

Units: H_i = 1 (the slice Hubble rate), M_P = M_P/H_i.  GeV conversions use H_i = 1.291e11 GeV unless stated.
Sections (keys of ../results/channels.json):
  widths      decay-width conventions (Dirac trace computed with explicit gamma matrices; Majorana by Wick counting) [D/C]
  optionB     reflection broken by A phi_I |H|^2: the A needed for Gamma/m = 1/6, Higgs mass^2 +-A phi_I, tachyonic growth,
              tree quartic shift, one-loop flatness, the stability bounds on A and what a safe A gives [C/E]
  blocking    singlet-fermion channel at the engine amplitude: y needed, y Phi vs m/2, open fraction, unblocking epoch,
              Greene-Kofman per-crossing production and energy fraction (checked numerically in fermion_modes.py) [C/E]
  seesaw      type-I seesaw N: y_nu^2 = mt M/v^2, Gamma_N = mt M^2/(8 pi v^2) in H_i units; M needed for tau_N < 8 [C]
  hidden      inflaton -> hidden leakage through the holder itself (2-body Y_L Y_L, 4-body Y_L + 3X off-shell, 6-body) and the
              heat it deposits, R_inf = 3 M_P^2 K^2 BR / E_X [C/E]
  higgs_portal  the reflection-even portal (lambda_IH/2) phi_I^2 |H|^2: Higgs held through reheating; window [E]
  thermalisation  Harigaya-Mukaida (1312.3097) condition and Davidson-Sarkar (hep-ph/0009078) timescale [Q/E]
Labels in the JSON: each block carries its label.
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad, dblquad

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rh_common as rc

H = rc.H_CAL_FREEX
MPh = rc.MP / H                       # M_P in H_i units (1.886e7)
EX = rc.EX_REF_3_05


# ------------------------------------------------------------------------------------------- widths
def dirac_trace_check(mphi=3.0, m=0.4):
    """sum_spins |u-bar(p1) v(p2)|^2 = Tr[(p1slash + m)(p2slash - m)] with explicit Dirac matrices; compare 2 mphi^2 beta^2."""
    g0 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]], complex)
    s = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]
    gs = [np.block([[np.zeros((2, 2)), si], [-si, np.zeros((2, 2))]]) for si in s]
    E = mphi / 2; p = math.sqrt(E * E - m * m)
    p1 = [E, 0, 0, p]; p2 = [E, 0, 0, -p]
    sl = lambda q: q[0] * g0 - q[1] * gs[0] - q[2] * gs[1] - q[3] * gs[2]
    tr = np.trace((sl(p1) + m * np.eye(4)) @ (sl(p2) - m * np.eye(4))).real
    beta2 = 1 - 4 * m * m / mphi ** 2
    return dict(trace=tr, expected=2 * mphi ** 2 * beta2, Gamma_Dirac_over_y2m=(p / (8 * math.pi * mphi ** 2)) * tr / mphi,
                analytic=beta2 ** 1.5 / (8 * math.pi))


def widths():
    d = dirac_trace_check()
    return dict(label='D (trace computed)', dirac_check=d,
                conventions={
                    'Dirac, L = -y phi psibar psi': 'Gamma = y^2 m beta^3 / (8 pi)',
                    'Majorana, L = -(y/2) phi Nbar N  [m_N(phi) = M + y phi; adopted]': 'Gamma = y^2 m beta^3 / (16 pi): vertex 2 x (y/2) = y, identical-particle 1/2',
                    'Majorana, L = -y phi N^T C N = -y phi Nbar N': 'Gamma = y^2 m beta^3 / (4 pi): vertex 2y, identical 1/2; here m_N(phi) = M + 2 y phi'},
                Nf_y2_for_Gamma_over_m_1_6=dict(dirac=8 * math.pi / 6, majorana_half=16 * math.pi / 6, majorana_full=4 * math.pi / 6),
                note='WORKING_THEORY "N_f y^2 = 4" is the Dirac normalisation (8 pi/6 = 4.19); with the adopted Majorana normalisation it is 8.38.')


# ------------------------------------------------------------------------------------------- option B
def optionB(b=3.0, Gam=0.5, Phi=math.sqrt(6) / 3, phi_slice=math.sqrt(6) / 3, xiH=0.0878):
    """A phi |H|^2 (4 real components, -(A/2) phi h_i^2 each): Gamma(phi -> H H*) = A^2/(8 pi m)."""
    A = math.sqrt(8 * math.pi * b * Gam)                       # H_i
    PhiH = Phi * MPh; phsH = phi_slice * MPh
    mH2_inf = A * phsH                                         # H_i^2, during inflation at the slice value
    grow = 2.3963 * math.sqrt(A * PhiH) / b                    # int_{-pi/2}^{pi/2} sqrt(cos) = 2.39628
    dlam = A ** 2 / (2 * b ** 2)                               # tree-level shift lambda_low = lambda_high - A^2/(2 m^2)
    V = 3 * MPh ** 2                                           # H_i^4
    cw = 4 * (A * phsH) ** 2 / (64 * math.pi ** 2) / V         # |V_CW| / V (log ~ 1)
    # slope: d/dphi [(A phi)^2 ln] ~ 2 A^2 phi ln / (16 pi^2) against V' = sqrt(2 eps) V / M_P, eps = r/16
    eps = rc.r_of_H(H) / 16
    slope_ratio = (2 * A ** 2 * phsH / (16 * math.pi ** 2)) / (math.sqrt(2 * eps) * V / MPh)
    # bounds on A: (i) negative-branch inflation stability |A| phi_* < 0.5 H^2 (half of 12 xi H^2 - ...);
    # (ii) no tachyonic resonance after inflation q_A = |A| Phi / m^2 < 1 (Enqvist et al. 1608.08848 find |sigma| <~ 1e8 GeV,
    #      i.e. q_A ~ 1, for their chaotic model)
    A_inf = 0.5 / phsH
    A_res = b ** 2 / PhiH
    out = dict(label='C (numbers) / E (criteria)', A_for_Gamma_over_m_1_6=A, A_GeV=A * H, mH2_during_inflation_H2=mH2_inf,
               mH_during_inflation_GeV=math.sqrt(mH2_inf) * H, tachyonic_growth_exponent_per_half_cycle=grow, tree_quartic_shift=dlam,
               CW_over_V_inflation=cw, CW_slope_over_tree_slope=slope_ratio,
               A_max_inflation_branch=A_inf, A_max_inflation_branch_GeV=A_inf * H, A_max_resonance=A_res, A_max_resonance_GeV=A_res * H)
    for nm, Ab in (('at_A_max_resonance', A_res), ('at_A_max_inflation', A_inf)):
        G = Ab ** 2 / (8 * math.pi * b)
        out[nm] = dict(A=Ab, Gamma_hh=G, Gamma_over_m=G / b, T_RH_GeV=rc.T_RH(G * H))
    # SM temperature at tau = 8, instant thermalisation, direct Higgs decays (Gamma_hh << Gamma_I(total) case: the SM share
    # of the radiation is Gamma_hh/Gamma_tot, so rho_SM ~ (2/5)(Gamma_hh/H) rho_I at early times)
    for nm in ('at_A_max_resonance', 'at_A_max_inflation'):
        G = out[nm]['Gamma_hh']
        tau = 8.0; Ht = 1.0 / (1.0 + 1.5 * tau)            # matter era from H = 1 at tau = 0
        rhoI = 3 * MPh ** 2 * Ht ** 2
        rho_sm = 0.4 * G / Ht * rhoI * H ** 4
        out[nm]['T_SM_tau8_GeV_instant'] = rc.T_of_rho(rho_sm)
        out[nm]['T_SM_tau8_over_Hi'] = rc.T_of_rho(rho_sm) / H
    return out


# ------------------------------------------------------------------------------------------- blocking and preheating
def gk_per_crossing(y, Phi, m):
    """Greene-Kofman: n_k = exp(-pi k^2 / (y |phidot_0|)) at a zero crossing, |phidot_0| = Phi m.  Majorana: 2 helicities.
    n = (y Phi m)^{3/2} / (4 pi^3); energy at the next turning point ~ n y Phi; fraction of rho_phi = m^2 Phi^2 / 2."""
    n = (y * Phi * m) ** 1.5 / (4 * math.pi ** 3)
    f1 = n * y * Phi / (0.5 * m * m * Phi * Phi)
    return n, f1


def blocking(b=3.0, Phi=math.sqrt(6) / 3, ys=(1.0, 0.1, 1e-2, 1e-3, 1e-4)):
    PhiH = Phi * MPh
    rows = []
    for y in list(ys) + [math.sqrt(16 * math.pi / 6)]:
        ratio = y * PhiH / (b / 2)
        open_frac = (2 / math.pi) * math.asin(min(1.0, 1 / ratio))
        a_unb = max(1.0, ratio ** (2 / 3))                      # Phi ~ a^{-3/2} in the matter era
        n, f1 = gk_per_crossing(y, PhiH, b)
        pert_half = y * y / 16                                 # Gamma_0 * (pi/m) = y^2/16 per half period
        rows.append(dict(y=y, Gamma0_over_m=y * y / (16 * math.pi), yPhi_over_half_m=ratio, open_fraction=open_frac,
                         a_unblock=a_unb, tau_unblock_matter=(a_unb ** 1.5 - 1) / 1.5, GK_n_per_crossing_H3=n,
                         GK_energy_fraction_per_crossing=f1, perturbative_fraction_per_half_period=pert_half,
                         kF_over_m=math.sqrt(y * PhiH * b / math.pi) / b))
    return dict(label='C (formulae evaluated) / E (Pauli saturation and back-reaction not solved)', b=b, Phi_over_MP=Phi,
                note='y = sqrt(16 pi/6) = 2.89 is the Majorana coupling for Gamma/m = 1/6 (N_f = 1).', rows=rows)


# ------------------------------------------------------------------------------------------- seesaw N
def seesaw(Ms=(0.25, 0.5, 1.0, 1.5, 3.0, 10.0, 30.0, 100.0), mt=rc.MNU, Hg=H):
    rows = []
    for M in Ms:
        MG = M * Hg
        ynu2 = mt * MG / rc.VEW ** 2
        GN = ynu2 * MG / (8 * math.pi)
        rows.append(dict(M_over_Hi=M, M_GeV=MG, ynu2=ynu2, Gamma_N_over_Hi=GN / Hg, tau_N_over_Hi=Hg / GN))
    M8 = math.sqrt(8 * math.pi * rc.VEW ** 2 * Hg / (8.0 * mt)) / Hg      # tau_N = 8 H_i^-1
    return dict(label='C', mt_eV=mt * 1e9, rows=rows, M_for_tau_N_8_over_Hi=M8,
                note='Gamma_N = (y_nu^dag y_nu)_11 M / (8 pi) = mt M^2 / (8 pi v^2), v = 174.1 GeV; at rest (boost E/M lengthens it).')


# ------------------------------------------------------------------------------------------- hidden leakage via the holder
def hidden(b_list=(3.0, 5.0, 10.0, 20.0, 25.0, 30.0, 40.0, 50.0, 100.0), lam=0.5):
    import fx_common as fc
    p = fc.params(lam)
    M, h, mu, kap = p['M'], p['h'], p['mu'], p['kap']
    mL = math.sqrt(M ** 2 - math.sqrt(2) * mu * h); mH = math.sqrt(M ** 2 + math.sqrt(2) * mu * h)
    msig = math.sqrt(2 * lam) * h
    dbeta = 2.192144074179184                                  # arg(G) - arg(mu): alpha_I
    G3 = lambda q: 3 * kap ** 2 * q / (256 * math.pi ** 3)     # Y_L -> 3X at mass q [Q ledger]
    lam_k = lambda a, b_, c: a * a + b_ * b_ + c * c - 2 * a * b_ - 2 * a * c - 2 * b_ * c
    rows = []
    for b in b_list:
        # g_I re-fixed for this b at the static slice: phihat0 = (sqrt(6)/b) M_P/H_i
        phs = math.sqrt(6) / b * MPh
        gI = abs(10j - mu) / phs
        AL = math.sqrt(2) * gI * h * abs(math.cos(dbeta)); ALH = math.sqrt(2) * gI * h * abs(math.sin(dbeta))
        two = AL ** 2 / (32 * math.pi * b) * math.sqrt(max(0.0, 1 - 4 * mL ** 2 / b ** 2)) if b > 2 * mL else 0.0
        lh = (ALH ** 2 / (16 * math.pi * b ** 3)) * math.sqrt(lam_k(b * b, mL * mL, mH * mH)) if b > mL + mH else 0.0
        four = 0.0
        if mL < b < 2 * mL:
            f = lambda q2: (AL ** 2 * math.sqrt(max(0.0, lam_k(b * b, mL * mL, q2))) / (16 * math.pi * b ** 3)) * \
                math.sqrt(q2) * G3(math.sqrt(q2)) / ((q2 - mL ** 2) ** 2) / math.pi
            four = quad(f, 0.0, (b - mL) ** 2, limit=200)[0]
        # six-body: both Y_L off-shell (q1 + q2 < b), Gamma_2 with the identical-particle 1/2
        def f6(q2b, q2a):
            l = lam_k(b * b, q2a, q2b)
            if l <= 0:
                return 0.0
            g2 = AL ** 2 * math.sqrt(l) / (32 * math.pi * b ** 3)
            return g2 * (math.sqrt(q2a) * G3(math.sqrt(q2a)) / (q2a - mL ** 2) ** 2 / math.pi) * \
                (math.sqrt(q2b) * G3(math.sqrt(q2b)) / (q2b - mL ** 2) ** 2 / math.pi)
        six = dblquad(f6, 0.0, b * b, lambda qa: 0.0, lambda qa: (b - math.sqrt(qa)) ** 2)[0] if b < mL else 0.0
        rows.append(dict(b=b, g_I=gI, A_L=AL, Gamma_YLYL=two, Gamma_YLYH=lh, Gamma_YL_3X_offshell=four, Gamma_6X=six,
                         Gamma_hidden=two + lh + four + six))
    return dict(label='C (rates in the narrow off-shell factorisation) / E (O(1) factors)', m_YL=mL, m_YH=mH, m_sigma=msig,
                thresholds=dict(two_body=2 * mL, LH=mL + mH, four_body=mL, sigma_pair=2 * msig), rows=rows,
                R_inf_formula='R_inf = 3 M_P^2 K^2 (Gamma_hid/Gamma_I) / E_X = %.4g K^2 BR (M_P/H_i = %.4g, E_X = %.1f)' % (3 * MPh ** 2 / EX, MPh, EX),
                R_inf_coefficient=3 * MPh ** 2 / EX,
                note='Hidden decay products thermalise into the hidden plasma before onset; R_inf adds to R_total (S/Q closed form). '
                     'The induced portals A_ISigma phi |Sigma|^2 etc. (renormalisation condition) are not included; they obey the same bound.')


# ------------------------------------------------------------------------------------------- even Higgs portal
def higgs_portal(b=3.0, Phi=math.sqrt(6) / 3, phi_star=None, q_max=(1.0, 350.0)):
    """(lambda/2) phi^2 |H|^2 -> Higgs mass^2 lambda phi^2 >= 0 on both sides (no tachyonic branch).  In the matter era
    lambda <phi^2>/H^2 = 3 lambda M_P^2 / m^2 is constant.  Window: lambda <phi^2> >= 10 H^2 and q = lambda Phi^2/(4 m^2) <= q_max."""
    lam_min = 10 * b ** 2 / (3 * MPh ** 2)
    PhiH = Phi * MPh
    out = dict(label='E (criteria) / C (numbers)', lam_min=lam_min, ratio_in_MD_per_lambda=3 * MPh ** 2 / b ** 2)
    for qm in q_max:
        out['lam_max_q%g' % qm] = 4 * qm * b ** 2 / PhiH ** 2
    out['q_at_lam_min'] = lam_min * PhiH ** 2 / (4 * b ** 2)
    out['note'] = ('q_max = 1: no broad resonance (conservative).  q_max = 350: the q of the bound lambda < 3e-8 of Enqvist et al. '
                   '1608.08848 at their m = 1.6e13 GeV, phi_end ~ sqrt(2) M_P (transferred, not recomputed).')
    return out


# ------------------------------------------------------------------------------------------- thermalisation
def thermalisation(b=3.0, Gammas=(1e-2, 1e-4, 1e-6, 1e-8, 1e-10), alpha=1 / 40.0, Hg=H):
    """Harigaya-Mukaida (1312.3097) condition for thermalising before reheating ends: alpha^{8/5} >> (m/M_P)(M_P^2 Gamma/m^3)^{1/5}."""
    rows = []
    m = b * Hg
    for G in Gammas:
        Gg = G * Hg
        rhs = (m / rc.MP) * (rc.MP ** 2 * Gg / m ** 3) ** 0.2
        rows.append(dict(Gamma_over_Hi=G, lhs=alpha ** 1.6, rhs=rhs, margin=alpha ** 1.6 / rhs))
    return dict(label='Q (abstract-level condition) / E', alpha=alpha, rows=rows)


def main():
    res = dict(method=__doc__, H_i_GeV=H, MP_over_Hi=MPh)
    res['widths'] = widths()
    res['optionB'] = optionB()
    res['blocking'] = blocking()
    res['seesaw'] = seesaw()
    res['hidden'] = hidden()
    res['higgs_portal'] = higgs_portal()
    res['thermalisation'] = thermalisation()
    json.dump(res, open(rc.RES / 'channels.json', 'w'), indent=1, default=float)
    print(json.dumps({k: v for k, v in res.items() if k not in ('method',)}, indent=1, default=lambda x: float(x))[:9000])


if __name__ == '__main__':
    main()
