"""Track m1_transport, item 1: tree-level elastic cross sections of the light X sector.

Independent code (own derivation, no Codex code imported). Units m_X = 1 unless stated.

Model (Codex convention, kin_rates NOTES section 1 and bf8ab3b DERIVATION.md):
    U = lambda_phi/4 (phi^2 - v^2)^2 - eps v (phi - v)^3 + g^2 phi^2 |X|^2 + lambda_X |X|^4
    X complex, m_X = g v, m_phi^2 = U''(v) = 2 lambda_phi v^2 (the tilt term has vanishing second
    derivative at v), trilinear h X X* vertex = 2 g^2 v = 2 g m_X.

Feynman rules (iM = sum of diagrams, all momenta incoming, metric +---):
    contact (X* X)^2 with coefficient lambda_X:  iM = -i lambda_X * (2 ways for the two X) * (2 ways
    for the two X*) = -4 i lambda_X, both for XX->XX and for X Xbar -> X Xbar.
    scalar exchange: iM = (-i g3)^2 * i/(q^2 - m_phi^2)  ->  M = g3^2/(m_phi^2 - q^2), g3 = 2 g m_X.

Cross section for 2->2 equal masses: dsigma/dOmega = |M|^2/(64 pi^2 s); identical final particles
carry an extra 1/2 in the total cross section.  At threshold s = 4 m_X^2 and t = u = 0.
"""
import json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
g = 0.1
lam = 1e-3
lam_phi = 0.03
eps = 0.002
v_GeV = 0.924089637375362
mX_GeV = g * v_GeV
mphi2 = 2 * lam_phi * v_GeV**2 / mX_GeV**2       # in m_X^2 units
g3 = 2 * g                                        # trilinear / m_X
Gamma_phi = g3**2 * math.sqrt(1 - 4 / mphi2) / (16 * math.pi * math.sqrt(mphi2))  # m_X units


def sigma_identical(M, s=4.0):
    """Total cross section, identical final particles, |M| constant: (1/2) * |M|^2/(16 pi s)."""
    return 0.5 * M**2 / (16 * math.pi * s)


def sigma_distinct(M, s=4.0):
    return M**2 / (16 * math.pi * s)


def main():
    out = {}
    # ---- XX -> XX at threshold --------------------------------------------------------------
    M_contact = -4 * lam
    M_t = g3**2 / (mphi2 - 0.0)       # t = 0 at threshold
    M_u = g3**2 / (mphi2 - 0.0)       # u = 0 at threshold
    M_XX = M_contact + M_t + M_u
    out['XX_to_XX'] = {
        'convention': 'U contains lambda_X |X|^4 with lambda_X = 1e-3 (Codex); M = -4 lambda_X + 2 g3^2/m_phi^2 at threshold',
        'M_contact': M_contact, 'M_phi_t_plus_u': M_t + M_u, 'M_total': M_XX,
        'ratio_exchange_over_contact_amplitude': (M_t + M_u) / abs(M_contact),
        'sigma_contact_only_mX2': sigma_identical(M_contact),
        'sigma_contact_only_formula': 'lambda_X^2/(8 pi m_X^2)',
        'sigma_contact_only_check': lam**2 / (8 * math.pi),
        'sigma_phi_exchange_only_mX2': sigma_identical(M_t + M_u),
        'sigma_total_mX2': sigma_identical(M_XX),
        'sigma_total_formula': '(4 lambda_X - g_h^2/3)^2/(128 pi m_X^2), g_h = 2 g = 0.2, m_phi^2 = 6 m_X^2',
        'sigma_total_check': (4 * lam - (2 * g)**2 / 3)**2 / (128 * math.pi),
        'codex_sigma_mX2': 2.16627561e-7,
        'ratio_total_over_codex': sigma_identical(M_XX) / 2.16627561e-7,
        'ratio_codex_over_contact_only': 2.16627561e-7 / sigma_identical(M_contact),
    }
    # task-text convention V = lambda |X|^4 / 4 : contact amplitude is -lambda
    out['XX_to_XX_task_convention_V_lambda_X4_over_4'] = {
        'M_contact': -lam, 'sigma_contact_only_mX2': sigma_identical(-lam),
        'sigma_contact_only_formula': 'lambda^2/(128 pi m_X^2)',
        'note': 'Codex and all M1 packets use U = lambda_X|X|^4 (no 1/4); the 1/4 convention would make the contact piece 16x smaller in sigma.'}
    # velocity dependence of the exchange term at v ~ 0.07: |t| <= 4 p^2 = 4 (m v)^2
    for vrel in (0.05, 0.0684, 0.1, 0.2):
        p = vrel / 2  # CM momentum for relative velocity vrel (NR)
        tmax = 4 * p**2
        M_tu_back = g3**2 / (mphi2 + tmax) + g3**2 / mphi2  # backward: t = -4p^2, u = 0
        out.setdefault('XX_velocity_dependence', {})[str(vrel)] = {
            'p_cm_over_mX': p, 'abs_t_max_over_mphi2': tmax / mphi2,
            'M_total_backward': M_contact + M_tu_back, 'M_total_threshold': M_XX,
            'fractional_change_of_sigma_at_backward_angle': (M_contact + M_tu_back)**2 / M_XX**2 - 1}
    # ---- X Xbar -> X Xbar at threshold -----------------------------------------------------
    s = 4.0
    M_s = g3**2 / (mphi2 - s)         # s-channel below the pole (pole at s = m_phi^2 = 6)
    M_t_bar = g3**2 / mphi2
    M_XXbar = M_contact + M_t_bar + M_s
    out['XXbar_to_XXbar'] = {
        'M_contact': M_contact, 'M_t': M_t_bar, 'M_s_threshold': M_s, 'M_total': M_XXbar,
        'sigma_total_mX2_threshold': sigma_distinct(M_XXbar),
        'sigma_without_s_channel_mX2': sigma_distinct(M_contact + M_t_bar),
        'ratio_to_XX_sigma': sigma_distinct(M_XXbar) / sigma_identical(M_XX),
        'resonance': {'sqrt_s_pole_mX': math.sqrt(mphi2), 'Gamma_phi_mX': Gamma_phi,
                      'Gamma_phi_GeV': Gamma_phi * mX_GeV,
                      'kinetic_energy_per_particle_at_pole_mX': math.sqrt(mphi2) / 2 - 1,
                      'lab_speed_each_at_pole_c': math.sqrt(1 - 1 / (math.sqrt(mphi2) / 2)**2),
                      'peak_sigma_unitarity_s_wave_mX2_at_pole': 16 * math.pi / (mphi2 - 4),
                      'note': 'at v ~ 0.07 the gas is at s - 4 = v^2 ~ 0.005 << 2 = m_phi^2 - 4: far below the pole; Codex kinetics used the resonance only at the 20 MeV source clock (k ~ m_X).'},
        'note': 'relevant only if antiparticles are present; the 2 MeV checkpoint has N_minus = -1.4e-13 (solver residual).'}
    # ---- M1 tilt: masses and couplings -----------------------------------------------------
    out['M1_tilt'] = {
        'm_phi_over_mX': math.sqrt(mphi2), 'm_phi_GeV': math.sqrt(mphi2) * mX_GeV,
        'm_phi2_formula': "U''(v) = 2 lambda_phi v^2; the -eps v (phi-v)^3 term has zero second derivative at v",
        'phi_cubic_candidate1_over_v': 6 * lam_phi, 'phi_cubic_M1_over_v': 6 * (lam_phi - eps),
        'phi_cubic_fractional_change': -eps / lam_phi,
        'XX_elastic_tree_amplitude_changed_by_tilt': False,
        'reason': 'tree XX->XX uses only the h X X* vertex 2 g^2 v and the X quartic; the phi cubic enters XX elastic scattering first at one loop (and at tree level only in phi phi -> phi phi, X Xbar -> phi phi).'}
    # ---- sigma/m in SIDM units for completeness --------------------------------------------
    GeV2_to_cm2 = 0.3894e-27
    GeV_to_g = 1.78266e-24
    sig_cm2 = sigma_identical(M_XX) / mX_GeV**2 * GeV2_to_cm2
    out['SIDM_units'] = {'sigma_cm2': sig_cm2, 'sigma_over_m_cm2_per_g': sig_cm2 / (mX_GeV * GeV_to_g)}
    (ROOT / 'results' / 'xsec.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
