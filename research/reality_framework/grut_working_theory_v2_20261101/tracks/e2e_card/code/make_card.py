"""e2e_card: writes CANDIDATE_M2.json, the full parameter card of M2 (M2-X and M2-Q).

Every entry: value, units, status (supplied | fitted: <datum> | derived | computed | estimate | open), source track,
window (where tuned or bounded; null if not scanned). Numbers are copied from the tracks' result files (read-only)
and from this track's spot checks (results/*.json). Nothing here is fixed by closure; nothing is a GRUT prediction.
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
TR = os.path.abspath(os.path.join(HERE, '..', '..'))
RES = os.path.join(HERE, '..', 'results')
m2 = json.load(open(os.path.join(TR, 'e2e_chain', 'results', 'm2_points.json')))
lep = json.load(open(os.path.join(RES, 'lepto_summary.json')))
hig = json.load(open(os.path.join(RES, 'higgs_spot.json')))
MP = 2.435e18
HREF = 1.3e11


def E(v, u, st, src, win=None, note=None):
    d = dict(value=v, units=u, status=st, source=src, window=win)
    if note:
        d['note'] = note
    return d


common = {
    'gravity_and_vacuum': {
        'M_P_reduced': E(MP, 'GeV', 'measured input', 'standard'),
        'rho_vac': E(2.52e-47, 'GeV^4', 'measured input, underived', 'WORKING_THEORY s2'),
        'vacuum_constant': E('V(true vacuum) = rho_vac (constant added)', None, 'supplied', 'mf_hybrid'),
    },
    'inflaton_sector_HFa': {
        'architecture': E('V = (lam/4)(chi^2 - v^2)^2 + (1/2) g^2 sigma^2 chi^2 - c chi + c1 (sigma - sigma_c) - (1/2) m2 (sigma - sigma_c)^2 + V_1loop + V_s + V_N + const',
                          None, 'supplied', 'mf_hybrid NOTES s3'),
        'v': E(0.019596, 'M_P', 'supplied (fixed under the V0 rescaling)', 'mf_hybrid'),
        'lam_design': E(2.3196e-7, '1', 'supplied at H_ref = 1.3e11 GeV; scales as (H_i/H_ref)^2 under the V0 rescaling [A]', 'mf_hybrid; mf_chain cmb_rescale'),
        'g_design': E(1.6347e-3, '1', 'supplied at H_ref; scales as H_i/H_ref', 'mf_hybrid'),
        'sigma_c': E(0.0057735, 'M_P', 'derived (sqrt(lam) v / g)', 'mf_hybrid'),
        'm_chi': E(250.0, 'H_i', 'derived (sqrt(2 lam) v), held fixed in H_i units', 'mf_hybrid'),
        'm_sigma': E(600.0, 'H_i', 'derived (g v)', 'mf_hybrid'),
        'c_bias': E(1e5, 'H_i^3', 'supplied (soft Z2 breaking)', 'mf_hybrid', 'not scanned; walls absent at this value (e2e_lattice)'),
        'c1': E(2.2546e4, 'H_i^3', 'fitted: A_s = 2.1e-9 with m2 (Planck 2018, 1807.06211)', 'mf_hybrid'),
        'm2': E(0.021963, 'H_i^2', 'fitted: n_s = 0.9649 (Planck 2018) [HFa]', 'mf_hybrid'),
        'mu_R_factor_kappa': E(10.0, '1', 'supplied: mu_R^2 = kappa g^2 (2 sigma_c)^2 (one-loop renormalisation condition; loop gives 58% of eta at the pivot)', 'mf_hybrid'),
        'higgs_portals_sigma_chi_s': E(0.0, '1', 'supplied (no tree portals)', 'mf_hybrid'),
        'hidden_portals_sigma_chi': E(0.0, '1', 'supplied (sequestered)', 'mf_hybrid'),
    },
    'inflaton_sector_HFb_alternative': {
        'c1': E(1.75e4, 'H_i^3', 'fitted: A_s', 'mf_hybrid'),
        'm2': E(0.0098, 'H_i^2', 'fitted: n_s = 0.9743 (ACT DR6 P-ACT, 2503.14452)', 'mf_hybrid'),
        'spectator_mass': E('m_s^2 = m_s0^2 + xi_s R, m_s0 = 3 H_i, xi_s = -3/4', None, 'supplied', 'mf_hybrid'),
        'note': 'HFb reproduces HFa to < 1% in the waterfall (e2e_lattice), K (e2e_kinetic, e2e_chain) and S/Q (e2e_chain); only n_s differs',
    },
    'spectator_holder': {
        'mass_mechanism_HFa': E('m_s^2 = g_s^2 chi^2, g_s v = 3 H_i (g_s = 8.17e-6 at H_ref)', None, 'supplied', 'mf_hybrid'),
        'holder_term': E('-(g_h e^{i alpha} s Sigma* Y^2 + h.c.)', None, 'supplied; the only s-odd term [D/A]', 'mf_hybrid'),
        'alpha_plus_sign_s': E('the DM sign', None, 'supplied (CP sign; not fixable by closure)', 'WORKING_THEORY s2/s7'),
        'lambda_s': E(1e-14, '1', 'supplied: required for boundedness, lambda_s >= 7.48 g_h^4 (V_min(|u|) -> -7.48 |u|^4) [C, e2e_chain skeptic]',
                      'e2e_chain/SKEPTIC.md', '>= 3.1e-15 (M2-X card s_*), >= 1.4e-15 (M2-Q card s_*)'),
        'y_s': E(4.6e-3, '1', 'supplied (s Nbar_s N_s; largest y_s keeping s light)', 'mf_hybrid'),
        'Gamma_s': E(1.3e-6, 'H_i', 'derived [E]', 'mf_hybrid'),
        'initial_condition_s_dot': E(0.0, '-', 'supplied (frozen during the observable e-folds; drift 4 H per e-fold, 0.2% of c0 over 50 e-folds)', 'e2e_chain skeptic'),
        'naturalness': E('s_* is a pre-observable initial condition: no de Sitter equilibrium sits in the window (tree V_min(|u(s)|) tachyonic, m2_eff = -6.3e-4 H^2 at M2-X); P_SY(window) <~ exp(-6e7)',
                         None, 'computed [C/D]', 'e2e_chain + skeptic'),
    },
    'singlet_neutrinos': {
        'content': E('four axial Majorana N_j: -(M_j/2) Nbar_j exp(i gamma5 phi_j/f_j) N_j, phi_1 = sigma, phi_2..4 = chi; plus massless N_s', None, 'supplied', 'mf_hybrid; e2e_sm'),
        'M_j': E([250.0, 400.0, 200.0, 120.0], 'H_i', 'supplied', 'mf_hybrid'),
        'f_j': E([1200.0, 2000.0, 2000.0, 2000.0], 'H_i', 'supplied', 'mf_hybrid', 'f/f_edge = 1.4-1.5 (axial EFT edge)'),
        'mtilde_j_S1': E([0.050, 0.09122, 0.050, 0.050], 'eV', 'supplied (scenario S1; N2 fixed by the CP-even sum rule sum eta_j mt_j = sum m_i [D])', 'e2e_sm'),
        'CP_parities_delta0': E([1, -1, 1, 1], '-', 'supplied', 'e2e_sm'),
        'CI_R0': E('3x4 least-norm CP-even solution; column 2 imaginary (sm_card.json CI_R0_re/_im); max|R| = 1.35; cancellation x4.1 in m_nu', None, 'computed', 'e2e_sm'),
        'CI_CP_odd_direction': E('R(delta) = R0 exp(delta A), A = where(eta_j eta_k > 0, i, 1) * A_ODD (A_ODD real antisymmetric in sm_card.json)', None, 'supplied',
                                 'e2e_sm sm_seesaw.family', note='sm_card.json stores only A_ODD; the i pattern is needed to reproduce eps (e2e_card check)'),
        'y2_jj': E([0.0536, 0.1565, 0.0429, 0.0257], '1', 'derived (mt M / v_EW^2) at H_ref; scale as H_i', 'e2e_sm; reproduced e2e_card'),
        'Gamma_rest': E([0.533, 2.49, 0.341, 0.123], 'H_i', 'derived; reproduced e2e_card', 'e2e_sm'),
        'eps_per_unit_delta': E([2.07e-2, 3.17e-3, 9.61e-3, 6.67e-3], '1', 'computed (Covi-Roulet-Vissani, hep-ph/9605319); reproduced by own code to 3e-4', 'e2e_sm; e2e_card'),
        'K_washout': E([46.3, 84.5, 46.3, 46.3], '1', 'derived (mt/m*, m* = 1.08e-3 eV)', 'e2e_sm'),
        'N_s': E('massless dark radiation; no soft mass, no seesaw Yukawa; N_s -> i gamma5 N_s exact', None, 'supplied', 'e2e_sm',
                 note='supersedes e2e_chain spectator.decay "N_s -> L H": the s energy goes to dark radiation, so the chain\'s +0.07%/+0.26% H_req entropy shift is dropped'),
    },
    'axial_UV_completion': {
        'construction': E('Phi_a = (rho_a/sqrt2) e^{i theta_a}, sigma = f1 theta_1, chi = f2 theta_2; M_j = y_j <rho_a>/sqrt2; hybrid potential = monodromy on the angles (0811.1989, 0808.0706)', None, 'supplied [A]', 'e2e_sm'),
        'radial_mass': E('m_rho = 4 pi f (strongly coupled edge)', None, 'supplied [A]; weak coupling fails', 'e2e_sm',
                         'theta-dot/m_rho = 0.6-0.9 (e2e_sm) / 0.92 for sigma with the bench peak theta-dot/f1 = 11.6 (e2e_card); weak coupling needs lambda >= 270 (sigma, e2e_card) / 248 (chi)'),
        'thermal_restoration': E('none: T_max/f1 = 0.27-0.35', None, 'derived', 'e2e_sm; e2e_chain T_max'),
        'theta_WWtilde': E('absent in the mass-phase basis (N gauge singlet: chiral N rotation anomaly-free) [A]; consistent with the field-redefinition bookkeeping of 2006.03148', None, 'argued', 'e2e_sm; e2e_card'),
    },
    'standard_model': {
        'M_t': E(172.57, 'GeV', 'measured input (PDG 2024 [R]), +-0.29', 'mf_higgs/e2e_higgs'),
        'M_h': E(125.2, 'GeV', 'measured input [R]', 'e2e_higgs'),
        'alpha_s_MZ': E(0.1180, '1', 'measured input [R], +-0.0009', 'e2e_higgs'),
        'xi_H_at_Mt': E(0.1, '1', 'supplied', 'WORKING_THEORY s2; e2e_higgs',
                        'lower edge 0.046-0.055 (inflation, Hawking-Moss); upper ~111 [Q mf_higgs skeptic]; own HM estimate 0.039/0.044 (e2e_card)'),
        'light_nu': E('normal ordering, m1 = 0; dm21^2 = 7.49e-5, dm31^2 = 2.513e-3 eV^2', 'eV^2', 'fitted: NuFIT 6.0 IC24 (2410.05380)', 'e2e_sm', 'm1 <~ 4.2 meV (sum < 0.0642 eV, DESI DR2+CMB 2503.14744)'),
        'PMNS': E('th12 33.68, th13 8.56, th23 43.3, delta_CP 212 deg; Majorana phases 0', 'deg', 'fitted (NuFIT 6.0) / supplied', 'e2e_sm'),
    },
    'hidden_sector_M1r': {
        'g_phiX': E(0.1, '1', 'supplied', 'WORKING_THEORY s2'),
        'lambda_X': E(1e-3, '1', 'supplied', 'WORKING_THEORY s2'),
        'lambda_phi': E(0.03, '1', 'supplied', 'WORKING_THEORY s2'),
        'v_phi': E('10 m_X', None, 'supplied', 'WORKING_THEORY s2'),
        'eps_tilt': E(0.002, '1', 'supplied', 'WORKING_THEORY s2', '(0.0017, 0.058) one-vacuum window; coherent-orbit transition sensitive at +-5%'),
        'm_X': E(0.0924, 'GeV', 'derived (g v_phi); v_phi set by the choice m_X = 92.4 MeV [supplied]', 'WORKING_THEORY s2'),
        'm_phi': E('sqrt(6) m_X', None, 'derived', 'WORKING_THEORY s2'),
        'M_Y': E(30.0, 'H_i', 'supplied', 'WORKING_THEORY s2'),
        'h_Sigma': E(40.0, 'H_i', 'supplied', 'WORKING_THEORY s2'),
        'mu_SigmaYY': E(7.16, 'H_i', 'supplied', 'WORKING_THEORY s2'),
        'kappa_YXXX': E(3.46e-3, '1', 'supplied', 'WORKING_THEORY s2'),
        'lambda_Y_lambda_Sigma': E(0.5, '1', 'supplied (perturbative to M_P)', 'WORKING_THEORY s2'),
        'charges': E('X 1, Y 3, Sigma 6', None, 'supplied', 'WORKING_THEORY s2'),
        'e_hidden': E(1e-6, '1', 'supplied (declared); with c0 it selects the branch', 'WORKING_THEORY s2; e2e_chain', 'M2-X needs e < e_crit = 1.7e-6; M2-Q needs e > 0.43e-6 (P0)'),
        'm_A': E('6 e h = 240 e H_i (2.4-2.7e7 GeV at e = 1e-6)', None, 'derived', 'e2e_relic'),
        'xi_X': E(-5.0 / 6.0, '1', 'supplied', 'WORKING_THEORY s2'),
        'c_L_induced_holder_portal': E(0.0, '1', 'supplied: induced holder portals vanish at M_Y (c_L = 1/16pi^2 variant shifts S/Q by +0.05-0.1)', 'e2e_chain', '<~ 1/16 pi^2'),
        'D_target': E(4.4e-10, 'GeV per unit entropy', 'supplied (rounded; Planck 2018 Omega_c h^2 = 0.120 gives 4.37e-10 [R])', 'WORKING_THEORY s2'),
    },
    'initial_conditions': {
        'sigma_at_CMB_exit': E(0.0212, 'M_P', 'on the slow-roll attractor (derived from N_*)', 'mf_hybrid'),
        'chi_during_inflation': E('chi_min = c / M_chi^2(sigma) > 0 everywhere', None, 'derived', 'mf_hybrid'),
        'higgs': E('h ~ 0 with de Sitter fluctuations (12 xi H^2 mass)', None, 'derived', 'e2e_higgs'),
        'hidden_state_at_slice': E('M1r release state prepared by the spectator holder (mf_ideal MODE 5 engine IC): Y at the held minimum, Sigma at h/sqrt2, vectors Bunch-Davies, X empty', None, 'supplied/derived', 'mf_ideal; e2e_chain'),
        'N_and_SM': E('empty at the slice', None, 'supplied', 'e2e_kinetic'),
        'slice': E('a = 1, H = H_i at sigma = sigma_c; eps_H = 1 at t = 0.19 H_i^-1', None, 'convention', 'mf_hybrid'),
    },
    'closures_and_rules (model choices)': {
        'Pauli_blocking': E('GM (box-mean saturation) for M2-X; P0 (local per-site) for M2-Q', None, 'supplied closure [A]; OPEN which is physical', 'e2e_lattice; e2e_chain'),
        'lattice_to_kinetic_handover': E('moments at t_sw = 3 H_i^-1, energy-conserving (degenerate sphere; dephased start)', None, 'supplied', 'e2e_chain'),
        'H_i_interpolation': E('root of ln(H_req/H_i) interpolated in ln H between the 8.4e10 and 1.3e11 GeV final histories', None, 'supplied; +-2% systematic (linear-in-H rule)', 'e2e_chain; skeptic'),
        'relic_rule': E('e2e_relic relic(): free X if S/Q > 3.4716 (M = 1 + 2p, relaxed-start p-map); Q-matter below (qm_micro + skeptic ripening)', None, 'supplied', 'e2e_relic'),
        'seed_rule': E('Q = 2e5 seed/nucleation assumed for 2.93 < S/Q < 3.47 (bulk-only band)', None, 'supplied [A]', 'e2e_chain skeptic'),
        'leptogenesis_source': E('non-thermal N = excess of the kinetic f_Nj over equilibrium; thermal N evolved separately', None, 'supplied [A] (corrects e2e_sm)', 'e2e_card'),
    },
}


def point(key, name, closure, sstar_card, extra):
    p = m2[key]
    H = p['self_consistent']['H_i_GeV']
    s = H / HREF
    iso = p['isocurvature']
    floor = iso['s_star_min_GeV']
    gh = p['holder']['g_h_s_star_over_H_i'] / (sstar_card / H)
    beta = 0.038 * (floor / sstar_card)**2
    fs_dec = p['spectator']['f_s_at_decay_at_s_star_min'] * (sstar_card / floor)**2
    lk = lep['M2X_GM' if key == 'M2X' else 'M2Q_P0']
    d = {
        'closure': E(closure, None, 'supplied closure [A]', 'e2e_chain'),
        'c0': E(p['holder']['c0'], 'units of M1r calibration (mu_eff = mu + c0 (10i - mu))', 'supplied (tuned)', 'e2e_chain; skeptic', extra['c0_window']),
        'holder_phase_delta': E(p['holder']['delta_deg'], 'deg', 'supplied (tuned)', 'e2e_chain; skeptic', extra['phase_window']),
        'mu_eff_slice': E(p['holder']['mu_eff_slice'], 'H_i (re, im)', 'derived', 'e2e_chain m2_points.json (NOTES/digest components were wrong; skeptic)',
                          note='|mu_eff| = %.3f H_i, arg %.2f deg' % (p['holder']['abs_mu_eff'], p['holder']['arg_mu_eff_deg'])),
        'g_h_s_star': E(p['holder']['g_h_s_star_over_H_i'], 'H_i', 'derived from c0', 'e2e_chain'),
        's_star': E(sstar_card, 'GeV', 'supplied initial condition (card value)', 'e2e_card', '>= %.3g GeV (beta_iso < 0.038, Planck 2018 1807.06211; D = %.2f +-10%%)' % (floor, iso['D'])),
        'g_h': E(gh, '1', 'derived (g_h s_* / s_*)', 'e2e_card'),
        'lambda_s_min': E(7.48 * gh**4, '1', 'derived bound', 'e2e_chain skeptic formula'),
        'H_i': E(H, 'GeV', 'fitted: Omega_DM (D = 4.4e-10 GeV) through the self-consistent chain', 'e2e_chain', '+-2% (interpolation rule)'),
        'V0_quarter': E(p['self_consistent']['V0_quarter_GeV'], 'GeV', 'derived', 'e2e_chain'),
        'lam': E(2.3196e-7 * s**2, '1', 'derived (V0 rescaling [A])', 'e2e_card'),
        'g': E(1.6347e-3 * s, '1', 'derived (V0 rescaling [A])', 'e2e_card'),
        'M_N_GeV': E([m * H for m in (250, 400, 200, 120)], 'GeV', 'derived', 'e2e_card'),
        'K': E(p['self_consistent']['K'], '1', 'computed (lattice to t = 3, kinetic network after)', 'e2e_chain'),
        'T_RH': E(extra['T_RH'], 'GeV', 'computed (interpolated in ln H between final histories)', 'e2e_chain histories; e2e_card'),
        'N_star': E(p['cmb']['Nstar'], 'e-folds', 'computed', 'e2e_chain; skeptic'),
        'n_s_HFa': E(p['cmb']['MS']['n_s'], '1', 'fitted target (retuned through m2), not an output', 'e2e_chain'),
        'r': E(p['cmb']['r'], '1', 'computed (MS modes)', 'e2e_chain; skeptic', '+-4-5% from the H_i rule'),
        'alpha_s_running': E(p['cmb']['MS']['alpha_s'], '1', 'computed (MS modes)', 'e2e_chain'),
        'capacity_N_max': E(p['cmb']['capacity_Nmax'], 'e-folds', 'computed', 'e2e_chain'),
        'S_per_Q_1e-6': E(p['S_per_Q']['e1e6'], '1', 'computed', 'e2e_chain; skeptic'),
        'S_per_Q_3e-6': E(p['S_per_Q']['e3e6'], '1', 'computed', 'e2e_chain'),
        'e_crit_bulk': E(p['S_per_Q']['e_crit_bulk'], '1', 'computed', 'e2e_chain'),
        'e_crit_Q2e5': E(p['S_per_Q']['e_crit_Q2e5'], '1', 'computed', 'e2e_chain'),
        'release': E({k: p['release'][k] for k in ('eps', 'E_over_Ec', 'N_YL', 'N_YH', 'N_sigma')}, '-', 'computed', 'e2e_chain'),
        'ledger_1e-6': E({k: p['ledger_1e6'][k] for k in ('R_vec', 'R_Y', 'R_casc', 'R_leak_cL0', 'dR_grav')}, 'heat ratios', 'computed', 'e2e_chain'),
        'relic': E(extra['relic'], '-', 'computed', 'e2e_relic relic(); e2e_chain'),
        'CP_phase_delta_lepton': E(lk['washout_x1']['delta_needed_corrected'], 'rad', 'fitted: eta_B = 6.12e-10 (Planck 2018 [R])', 'e2e_card (corrected e2e_sm source)',
                                   '%.1e - %.1e (washout x0.5 - x2); e2e_sm uncorrected value %.1e' % (lk['washout_x0.5']['delta_needed_corrected'], lk['washout_x2']['delta_needed_corrected'], lk['washout_x1']['delta_needed_as_e2e_sm'])),
        'etaB_per_unit_delta': E(lk['washout_x1']['etaB_per_delta_corrected'], '1', 'computed [C on Q/E rates]', 'e2e_card'),
        'beta_iso_at_card_s_star': E(beta, '1', 'derived', 'e2e_card', '< 0.038'),
        's_energy_fraction_at_decay': E(fs_dec, '1', 'derived (scaled from e2e_chain by (s_*/floor)^2)', 'e2e_card'),
        'DeltaN_eff_Ns': E(2.86 * fs_dec, '1', 'derived [E] (N_s dark radiation)', 'e2e_card'),
        'Mt_crit_grav_only': E(hig['Mt_crit_grav_at_M2'][extra['hk']][extra['hp']], 'GeV', 'computed+interpolated [E] (quadratic in ln H_i through e2e_higgs nodes)', 'e2e_higgs; e2e_card',
                               'alpha_s 0.1171/0.1189: %.2f/%.2f' % (hig['Mt_crit_grav_at_M2']['KIN_alphas_0.1171'][extra['hp']], hig['Mt_crit_grav_at_M2']['KIN_alphas_0.1189'][extra['hp']])),
        'higgs_with_media': E('D < 0 for M_t <= 174.6 GeV on every history (N gas + SM medium)', None, 'computed [C on E]', 'e2e_higgs'),
    }
    d.update(extra.get('more', {}))
    return d


m2x_rel = m2['M2X']['relic']
m2q_rel = m2['M2Q']['relic']
card = {
    'label': 'CANDIDATE_M2: full parameter card of M2 (M1r hidden sector + HFa multi-field inflaton sector). Every entry is a model choice or follows from model choices; nothing is fixed by closure; nothing is a GRUT prediction. Data comparisons are retrospective.',
    'conventions': 'H_i units unless stated; slice a = 1, H = H_i at sigma = sigma_c; c0 = holder strength at the slice in units of M1r calibration mu_eff = 10i H_i; status words: supplied (model choice), fitted (to named datum), derived (algebra), computed (numerical), estimate, open.',
    'common': common,
    'M2-X': point('M2X', 'M2-X', 'GM (box-mean Pauli saturation)', 1.5e16, dict(
        c0_window='self-consistent free X with e_crit at c0 1.45-1.60 (H_i 1.0-1.5e11); e_crit ~ 2e-6 only near 1.55; e_crit 1.3-3.8e-6 over delta -2..+6',
        phase_window='-2..+6 deg (delta = -4 loses e_crit)', T_RH=2.72e12, hk='LAT_GM', hp='M2X_1.123e11',
        relic=dict(branch='free X (asymmetric)', M_over_mX=m2x_rel['M_over_mX'], p=m2x_rel['p'], m_X_GeV=0.0924,
                   n_local_cm3=m2x_rel['observables']['n_local_per_cm3'], sigma_over_m_XX_cm2_g=6.0e-8, sigma_over_m_XXbar_cm2_g=7.07e-7,
                   annihilation='X Xbar -> phi phi closed; no signal'))),
    'M2-Q': point('M2Q', 'M2-Q', 'P0 (local per-site Pauli blocking)', 1.5e16, dict(
        c0_window='1.40 (bistable: free X 1.04e11 and Q-matter 1.17e11) - 1.50; GM + fragment-spread pocket c0 1.485-1.52, delta +2..+3.5 deg (bulk-only)',
        phase_window='at least -4..+4 deg at c0 1.45 (skeptic; e2e_chain had 0..+4)', T_RH=2.40e12, hk='LAT_P0', hp='M2Q_9.89e10',
        relic=dict(branch='Q-matter', M_over_mX=m2q_rel['M_over_mX'], Q_mean=m2q_rel['Q_mean'], f_v=m2q_rel['f_v'],
                   object_mass_ug=9.61, object_mass_over_MPl=0.44, radius_nm=0.730, n_local_m3=4.91e-14,
                   crossing_rate_m2_yr='0.55 (n<v>A, normal) / 0.27 (isotropic flat plane)', impulse_1mg_1mm_kg_m_s=3.6e-27,
                   sigma_over_m_contact_cm2_g=7.0e-9, Pk='CDM-like (zeta ~ 1e-7) [Q qm_transfer fit]'))),
}
json.dump(card, open(os.path.join(HERE, '..', 'CANDIDATE_M2.json'), 'w'), indent=1)
print('written', os.path.getsize(os.path.join(HERE, '..', 'CANDIDATE_M2.json')), 'bytes')
for k in ('M2-X', 'M2-Q'):
    print(k, {kk: card[k][kk]['value'] for kk in ('H_i', 'g_h', 'lambda_s_min', 'beta_iso_at_card_s_star', 'DeltaN_eff_Ns', 'CP_phase_delta_lepton', 'Mt_crit_grav_only', 'lam', 'g')})
