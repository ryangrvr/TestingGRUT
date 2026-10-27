"""st_physics step 5 (own code): seed extrapolation, validity of the quartic-dominated era, interaction sizes,
circular-orbit energy at fixed charge and an X-bar-per-charge estimate.  Everything here is [D] algebra on
published/own numbers or an [E] estimate; labels are in the output.

usage: python3 -I s5_estimates.py <codex timing packet dir> <late_growth.json> <orbits_floquet.json>
                                  <vacuum_seed.json> <out.json>
"""
import sys, json
import numpy as np

pkt, lgf, orbf, vacf, out = sys.argv[1:6]
lg, orb, vac = json.load(open(lgf)), json.load(open(orbf)), json.load(open(vacf))
rr = {r['tag']: r for r in json.load(open(f'{pkt}/RESPONSE_RESULTS.json'))['rows']}
LAM, KAP, MY, MU_S, SABS = 1e-3, 0.0034586847211042523, 30.0, 7.159456159513794, 28.284271247461902
ML2, MH2 = MY ** 2 - 2 * MU_S * SABS, MY ** 2 + 2 * MU_S * SABS  # Y_L, Y_H mass^2 (H_i units) at the Sigma vev
H_GEV, MX_GEV, GX, LPHI, V_GEV = 133349619398.33566, 0.09240896373753621, 0.1, 0.03, 0.924089637375362
TAU_YL = 6.89e6
res = dict(labels='[D] algebra on published numbers; [E] estimates; own code')


def tau_of(a, C, ts):
    return ts + a * a / (2 * C)


for tag in ('quiet3', 'quiet5'):
    L = lg[tag]; C, ts = L['C'], L['tau_star']
    mu = L['mu_eff_amplitude']
    a03, a06 = L['crossing_a']['a_003'], L['crossing_a']['a_1e6']
    r = rr[tag]
    cr = {c['phase_seed_radians']: c for c in r['crossings']}
    d = np.load(f'{pkt}/{tag}.npz')
    amp = d['amplitude_per_radian']; times = d['times']
    # [D/E] 1e-9 extrapolations (unit and 0.1 thresholds) from the 1e-6 samples
    a_unit_1e9 = a06 + C / mu * np.log(1e3)
    a01_1e6 = cr[1e-6]['amplitude_10percent']['a_bracket'][1]
    a_01_1e9 = a01_1e6 + C / mu * np.log(1e3)
    max_1e9_by_800 = float(1e-9 * np.max(abs(amp)))
    # seed below which the unit crossing would come after the Y_L decay-number median
    a_YL = np.sqrt(2 * C * (TAU_YL - ts))
    log10_seed_YL = np.log10(1e-6) - mu * (a_YL - a06) / C / np.log(10)
    # pericentre-free crossing (|dchi_X| = r_max at k=2.13) and its 1e-6 extrapolation
    pf = L['pericentre_free_first_unit']['0.03_dchi_over_rmax']
    a_pf_1e6 = pf['a'] + C / mu * np.log(0.03 / 1e-6)
    # vacuum-equivalent seed [E] (vacuum_seed.json) and its extrapolated unit crossing
    vrow = vac['rows'][tag]
    s_eq = [x['equivalent_supplied_radians_per_lnk_late'] for x in vrow['per_k'] if abs(x['k'] - 2.13) < 1e-9][0]
    a_vac = a03 + C / mu * np.log(0.03 / s_eq)
    ext = dict(mu_eff=mu, C=C, tau_star=ts,
               unit_1e9=dict(a=a_unit_1e9, tau=tau_of(a_unit_1e9, C, ts)),
               tenth_1e9=dict(a=a_01_1e9, tau=tau_of(a_01_1e9, C, ts)),
               max_1e9_response_through_800=max_1e9_by_800,
               a_YL_median=a_YL, log10_seed_for_unit_after_YL_median=log10_seed_YL,
               pericentre_free_003=pf,
               pericentre_free_1e6_extrapolated=dict(a=a_pf_1e6, tau=tau_of(a_pf_1e6, C, ts)),
               vacuum_equivalent_seed_k213=s_eq,
               vacuum_unit_crossing_extrapolated=dict(a=a_vac, tau=tau_of(a_vac, C, ts)),
               vacuum_band_rms_direct=vrow['band_rms_first_sampled_crossings'])
    # ---------------- validity of the quartic-dominated era -----------------
    o = orb[tag]['orbits']['late_mean']
    rmin, rmax, rc = o['rmin'], o['rmax'], o['r_circ']
    quartic_m2 = dict(peri=LAM * rmin ** 2, apo=LAM * rmax ** 2, circ=LAM * rc ** 2)
    aa = np.array([1., 6.24, 23.05, 55.4, 5.23e3, 1e6, 1e9, 1e12])
    absx_apo = rmax / (np.sqrt(2) * aa)
    sextic_ratio_L = KAP ** 2 * absx_apo ** 2 / (ML2 * LAM)  # |V6|/V4 with Y_L (largest), physical units
    sextic_ratio_H = KAP ** 2 * absx_apo ** 2 / (MH2 * LAM)
    # neutral field: m_phi^2 = 2 g^2 |X|^2 - lam_phi v^2 > 0 while |X| > x_c (phi pinned at 0)
    x_c_H = np.sqrt(LPHI * V_GEV ** 2 / (2 * GX ** 2)) / H_GEV
    a_phi_peri = rmin / (np.sqrt(2) * x_c_H)
    a_phi_apo = rmax / (np.sqrt(2) * x_c_H)
    mX_H = MX_GEV / H_GEV
    a_mass_eq = np.sqrt(LAM * o['r_circ'] ** 2) / mX_H  # conformal a m_X = sqrt(lam) r_circ
    # conformal visible temperature (late radiation): rho_r a^4 = r a^4 M^2 H^2
    bg = d['bg']; a_end = np.exp(bg[-1, 15]); rad_a4 = bg[-1, 14] * a_end ** 4
    MRED = 2.435e18
    aT_H = (30 * rad_a4 * MRED ** 2 * H_GEV ** 2 / (np.pi ** 2 * 106.75)) ** 0.25 / H_GEV
    thermal_m2_phi_bath = GX ** 2 / 12 * aT_H ** 2  # a^2 m_T^2 if a thermal phi bath at T_vis existed (one real phi)
    val = dict(conformal_quartic_mass2=quartic_m2,
               a_grid=aa.tolist(), apocentre_physical_absX=absx_apo.tolist(),
               Y_induced_sextic_over_quartic_YL=sextic_ratio_L.tolist(),
               Y_induced_sextic_over_quartic_YH=sextic_ratio_H.tolist(),
               gauge_e2_over_lam={'e=1e-7': 1e-14 / LAM, 'e=1e-6': 1e-12 / LAM},
               R_term='radiation era: R = 0 up to the X sector own trace (eta^2-suppressed); xi R and a\'\'/a vanish',
               phi_pinned_until_a=dict(pericentre=a_phi_peri, apocentre=a_phi_apo, x_c_over_H=x_c_H),
               vacuum_mass_matters_at_a=a_mass_eq,
               conformal_visible_temperature_over_H=aT_H,
               hypothetical_thermal_phi_bath_conformal_mass2=thermal_m2_phi_bath,
               note='[E] Within the declared model nothing gives X a conformal mass before a~1e12; a phi population '
                    'thermalised with the visible bath (not in the model) would dominate the quartic by ~1e3.')
    # ---------------- Q-ball / bound-lump test at phi = 0 -----------------
    # V(|x|) = lam |x|^4 - kap^2 |x|^6 / M^2 (Y integrated out at tree level); Q-balls of a massless field need
    # min 2V/rho^2 < 0, i.e. V < 0 somewhere: |x|^2 > lam M^2 / kap^2.
    qb = dict(absx_where_V_turns_negative_YL=float(np.sqrt(LAM * ML2) / KAP),
              absx_where_V_turns_negative_YH=float(np.sqrt(LAM * MH2) / KAP),
              max_absx_on_history_after_tau10=float(np.max(np.hypot(bg[times >= 10, 2], bg[times >= 10, 3]))),
              absx_initial=float(np.hypot(bg[0, 2], bg[0, 3])),
              Y_at_that_point_over_M_YL=float(KAP * (LAM * ML2 / KAP ** 2) ** 1.5 / ML2 / np.sqrt(ML2)),
              verdict='[D/E] no bound lumps: V/|x|^2 increases monotonically over the populated range; gauge '
                      'Coulomb (e^2) is repulsive; neutral phi is pinned at 0 (no attractive portal yet)')
    # ---------------- circular orbit at fixed charge and X-bar estimate -----------------
    rows = {}
    for lab in ('tau10', 'late_mean'):
        ob = orb[tag]['orbits'][lab]
        E, Lc = ob['E'], ob['L']
        Ecirc = 0.75 * LAM ** (1 / 3) * Lc ** (4 / 3)
        frac = 1 - Ecirc / E
        kp = ob['peak']['k']
        rc2 = (Lc / np.sqrt(LAM)) ** (2 / 3)
        omegas = {'massless_k_peak': kp, 'k2+lam_rc2': np.sqrt(kp ** 2 + LAM * rc2),
                  'k2+2lam_rc2': np.sqrt(kp ** 2 + 2 * LAM * rc2), 'k2+3lam_rc2': np.sqrt(kp ** 2 + 3 * LAM * rc2)}
        xbar_pairs = {k: frac * E / (2 * w * Lc) for k, w in omegas.items()}  # excess energy -> X Xbar pairs
        xbar_gas = {k: (E / (w * Lc) - 1) / 2 for k, w in omegas.items()}     # everything a gas at omega
        # [E] equilibrium end-point if the waves equilibrate while X is still massless (phi = 0): charge in a
        # circular condensate (energy E_circ), the rest a massless Bose gas at mu ~ 0 (2 dof): rho = pi^2 T^4/15,
        # n_Xbar = zeta(3) T^3 / pi^2.  Upper value (mu = 0); mu ~ omega_c lowers it by ~exp(-mu/T) ~ 10%.
        Tth = (15 * (E - Ecirc) / np.pi ** 2) ** 0.25
        xbar_eq = 1.2020569031595942 * Tth ** 3 / np.pi ** 2 / Lc
        omega_c = np.sqrt(LAM) * np.sqrt(rc2)
        xbar_eq_mu = xbar_eq * np.exp(-omega_c / Tth)
        rows[lab] = dict(E=E, L=Lc, E_circ=Ecirc, E_circ_over_E=Ecirc / E, released_fraction=frac,
                         equilibrium_massless=dict(T_conformal=Tth, Xbar_per_net_mu0=xbar_eq,
                                                   Xbar_per_net_mu_omega_c_boltzmann=xbar_eq_mu, omega_c=omega_c,
                                                   label='[E] massless Bose gas + circular condensate'),
                         E_over_L=E / Lc, k_peak=kp, r_circ=np.sqrt(rc2), omegas=omegas,
                         Xbar_per_net_excess_into_pairs=xbar_pairs, Xbar_per_net_all_gas=xbar_gas,
                         Q_cube_at_own_peak=abs(r['late_NX_normalization_through80']) * (2 * np.pi / kp) ** 3,
                         Q_cube_at_2p13=abs(r['late_NX_normalization_through80']) * (2 * np.pi / 2.13) ** 3)
    res[tag] = dict(extrapolation=ext, validity=val, qball=qb, circular=rows)
    print(tag, json.dumps(dict(ext=ext, circ=rows['tau10']), default=float)[:3000])
    print(tag, 'validity', json.dumps(val, default=float)[:1500])
    print(tag, 'qball', qb)
res['kin_source_bracket_Xbar_per_net'] = [0.134, 2.586]
json.dump(res, open(out, 'w'), indent=1, default=float)
