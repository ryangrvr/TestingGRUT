"""fx_observables items 3 and 4a (obligations): late-universe channels of conservative M1, hidden radiation, and the
cosmological price of a detectable Higgs portal.  [C] arithmetic on our own numbers, [E] where marked.

Inputs: results/objects.json, results/dd.json (own); tracks/fx_formation results/eq_history.json (T_h(T_vis), S/Q,
read-only sibling result) [Q own]; REVIEW_REHEATING / rp_portals Higgs source (Y_pair = 2.237e-8 at b = 5e-11) [Q own].
Writes results/late_and_portal.json.
"""
import json, math, os
import numpy as np
from scipy.special import kn
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
FXF = os.path.join(HERE, '..', '..', 'fx_formation', 'results', 'eq_history.json')

MX = 0.0924089637375362
YQ = 5.07941423e-9
D = 4.4e-10
HBARC = 1.973269804e-14
INVMX_CM = HBARC / MX
C_CMS = 2.99792458e10
GEV_G = 1.78266192e-24
T_UNIV_S = 4.35e17
MPL = 1.22089e19
VH, MH = 246.22, 125.20
GAMMA_H_SM = 4.10e-3     # GeV [R], SM Higgs width

obj = json.load(open(os.path.join(RES, 'objects.json')))
rows = {r['Q']: r for r in obj['rows']}
dd = json.load(open(os.path.join(RES, 'dd.json')))


def sigma_geo_cm2(Q):
    return rows[Q]['geometric_obj_obj']['r90']['sigma_cm2']


def main():
    out = {}
    # ---------------- 3a: Xbar burn-out on objects (early universe) ----------------
    # Codex bf8ab3b: charge exchange Xbar + core -> X + core(Q-2) at k = 1 m_X: sigma = 2.76/m_X^2 (1.30 % of events)
    # [Q Codex, confirmed as an implication in m1_transport SKEPTIC].  Exothermic (releases 2 omega = 1.76 m_X), so
    # sigma v -> const at low k (Wigner) [D]; sigma v(k=1) = 2.76 * 0.7071.
    sv = 2.76 * (1 / math.sqrt(2))                 # m_X^-2 (c = 1)
    n_core_2MeV = 1.2141071422639381e-18          # m_X^3 at T_vis = 2 MeV for Q = 2e5 (m1_transport) [Q own]
    H_2MeV = 1.78593e-24 / MX                      # m_X units [Q own]
    G_over_H_2MeV = n_core_2MeV * sv / H_2MeV
    # Gamma/H ~ T in RD for constant sigma v -> e-folds from 2 MeV to today ~ Gamma/H(2 MeV) (geometric series in ln a)
    out['Xbar_burnout'] = dict(sigma_v_mX2=sv, Gamma_over_H_at_2MeV_Q2e5=G_over_H_2MeV,
                               efolds_after_2MeV=G_over_H_2MeV, T_vis_where_Gamma_eq_H_eV=2e6 / G_over_H_2MeV,
                               scaling_with_Q='Gamma/H ~ P_abs R^2/Q ~ Q^-1/3 at fixed absorption probability per geometric hit',
                               efolds_Q1e14_if_Pabs_0p025=G_over_H_2MeV * (1e14 / 2e5)**(-1 / 3),
                               conclusion='p_today = exp(-1e5) for Q = 2e5; X-bar does not survive to halos unless the per-hit annihilation probability is < 1e-5')
    # hypothetical: p = 0.005 survived to today, Milky Way at the Sun
    rho_loc = 0.4                                   # GeV/cm^3 [R]
    M2 = rows[2e5]['mass_GeV']
    n_obj = rho_loc / M2
    G_ann = n_obj * sv * INVMX_CM**2 * C_CMS        # per Xbar per s
    E_rel = 2 * rows[2e5]['Omega'] * MX              # GeV per event, into one fast X
    out['Xbar_halo_if_survived'] = dict(p=0.005, n_obj_cm3=n_obj, rate_per_Xbar_s=G_ann, fraction_annihilated_per_tH=G_ann * T_UNIV_S,
                                        energy_per_event_GeV=E_rel, outgoing_X_speed_c=math.sqrt(1 - 1 / (1 + 2 * rows[2e5]['Omega'])**2),
                                        DM_mass_fraction_converted_per_t_H=0.005 * (E_rel / (rows[2e5]['E_over_Qm'] * MX)) * G_ann * T_UNIV_S,
                                        visibility='invisible: products are free X (no SM portal); v >> v_esc, they leave the halo')
    # ---------------- 3b: object-object mergers today ----------------
    env = {'MW_solar': (0.4, 220.0 * math.sqrt(2)), 'dwarf_core': (3.8, 15.0 * math.sqrt(2)), 'cluster_core': (1.0, 1000.0 * math.sqrt(2)),
           'MW_inner_1kpc': (8.0, 150.0 * math.sqrt(2))}
    mer = {}
    mlist = {m['Q']: m for m in obj['mergers']}
    for Q in (1e5, 2e5, 1e6, 1e7, 1e8):
        M = rows[Q]['mass_GeV']
        sig = sigma_geo_cm2(Q)
        fr = mlist[Q]['frac_of_pair_mass'] if Q in mlist else None
        d = {}
        for k, (rho, vrel) in env.items():
            G = rho / M * sig * vrel * 1e5          # per object per s
            d[k] = dict(collisions_per_object_per_tH=G * T_UNIV_S, mass_fraction_released_per_tH=G * T_UNIV_S * (fr or 0) / 2)
        mer['Q%.0e' % Q] = dict(sigma_cm2=sig, sigma_over_M_cm2_g=rows[Q]['geometric_obj_obj']['r90']['sigma_over_M_cm2_g'],
                                release_fraction_per_merger=fr, env=d,
                                kR_halo=M * 1e-3 * rows[Q]['r90'] / MX)   # p R with p = M v, v = 1e-3 (classical regime if >> 1)
    out['mergers_today'] = mer
    out['mergers_note'] = ('mass released per merger = E(Q)+E(Q)-E(2Q) (own solved family); it goes to hidden quanta only '
                           '(free X with v ~ 0.1-0.9 c, phi -> X Xbar); no massless hidden state exists, so no dark radiation; '
                           'decaying-DM-to-DR bounds allow ~2-4 % of the DM to convert (Poulin et al. 2016 arXiv:1606.02073 [R])')
    # ---------------- 3c: hidden radiation / Delta N_eff ----------------
    # ratio rho_DM / rho_one_nu at T_vis: rho_DM = D s, rho_nu1 = (7/8) 2 (pi^2/30) T_nu^4
    def ratio_DM_nu(T, gs, Tnu_over_T):
        s = 2 * math.pi**2 / 45 * gs * T**3
        rnu = 7 / 8 * 2 * math.pi**2 / 30 * (Tnu_over_T * T)**4
        return D * s / rnu
    hid = {}
    hist = json.load(open(FXF))['histories']['R0_heavy|Y_onset']['rows']
    hr = [r for r in hist if 'H_over_mX' in r]
    Tv = np.array([r['T_vis_GeV'] for r in hr]); Th = np.array([r['T'] for r in hr]); eta = np.array([r['eta'] for r in hr])
    for lab, T, gs, tnu in (('BBN_1MeV', 1e-3, 10.75, 1.0), ('BBN_0.1MeV', 1e-4, 3.909, (4 / 11)**(1 / 3))):
        r_all = ratio_DM_nu(T, gs, tnu)
        th = float(np.interp(math.log(T), np.log(Tv[::-1]), Th[::-1]))
        et = float(np.interp(math.log(T), np.log(Tv[::-1]), eta[::-1]))
        thermal_frac = (1 - et) * 1.5 * th / 0.94     # vapour kinetic energy / DM mass
        hid[lab] = dict(DeltaNeff_if_all_DM_counted=r_all, T_h_over_mX=th, vapour_charge_fraction=1 - et,
                        DeltaNeff_hidden_thermal=r_all * thermal_frac, T_h_over_T_vis=th * MX / T)
    out['hidden_radiation'] = dict(rows=hid, light_states='none: lightest hidden state X (92.4 MeV), then phi (226 MeV); hidden U(1) Higgsed (m_A = 6 e h = 3.2e6-3.2e7 GeV for e = 1e-7-1e-6, decays to X Xbar at t ~ 3e-7-3e-4 s); Y, Sigma ~ 1e12 GeV; no Goldstone (gauged); the X-number violation is the tiny sextic',
                                   hidden_entropy_per_visible=YQ * 1.235,
                                   comparison='Planck 2018 N_eff = 2.99 +- 0.17; BBN N_eff ~ 2.9 +- 0.2; SO/CMB-S4-class sigma ~ 0.03-0.06 [R]')
    # ---------------- 4a: portal obligations ----------------
    th_b = 0.0290978                                # theta / b  [Q own rp_portals]
    def Gamma_hss(b):
        beta = math.sqrt(1 - 4 * 0.22635**2 / MH**2)
        return b * b * VH * VH * beta / (8 * math.pi * MH)
    def Ypair(b):
        return 2.23656e-8 * (b / 5e-11)**2         # [Q own rp_portals, g* = 106.75; +18 % with g*(T)]
    # entropy injected into the hidden sector by h -> ss at T_vis ~ 20-200 GeV, thermalised at T_h = r T_vis
    hv = [(r['T_vis_GeV'], r['T']) for r in hr if r['T_vis_GeV'] > 5]
    rT = np.median([th * MX / tv for tv, th in hv if 10 < tv < 200])
    gstar = 106.75
    def dS_over_Q_linear(b):
        G = Gamma_hss(b)
        def integrand(T):
            nh = MH**2 * T / (2 * math.pi**2) * kn(2, MH / T)
            H = 1.66 * math.sqrt(gstar) * T * T / MPL
            s = 2 * math.pi**2 / 45 * gstar * T**3
            # d(Delta S)/dT = [Gamma m_h n_h / (r T)] / (n_Q) * dt/dT, dt = dT/(H T)
            return G * MH * nh / (rT * T) / (YQ * s) / (H * T)
        return quad(integrand, 2.0, 2000.0, limit=200, points=[20, 50, 125])[0]
    def S_over_Q_saturated(b):
        # injection >> hidden heat: relativistic hidden plasma g_h = 3 holding the injected energy; J = 5 Y_hss
        J = 5 * Ypair(b) / 2
        # rho_inj / s_vis = J T  ->  rho_h = J T s;  s_h = (4/3) (pi^2 g_h/30)^{1/4} rho_h^{3/4};  per charge / (YQ s)
        T = 50.0
        s = 2 * math.pi**2 / 45 * gstar * T**3
        rho = J * T * s
        sh = 4 / 3 * (math.pi**2 * 3 / 30)**0.25 * rho**0.75
        return sh / (YQ * s)
    lin_coef = dS_over_Q_linear(5e-11)
    budget = {'R0_heavy (S/Q = 1.24)': 3.47 - 1.235, 'medium e1e6 (S/Q = 2.07)': 3.47 - 2.07}
    bmax = {k: 5e-11 * math.sqrt(v / lin_coef) for k, v in budget.items()}
    # detectability couplings from dd.json (Q = 2e5, phi^2 portal)
    r2 = [r for r in dd['rows'] if r['Q'] == 2e5 and r['portal'] == 'phi2'][0]
    rX = [r for r in dd['rows'] if r['Q'] == 2e5 and r['portal'] == 'X2'][0]
    b_LZ = r2['LZ2024_equiv']['Eup55']['coupling_limit']
    b_XLZD = r2['projections_10_events']['XLZD_200tyr']['coupling_for_10_events']
    port = dict(theta_over_b=th_b, b_LZ2024_limit=b_LZ, b_XLZD_10_events=b_XLZD,
                lambdaHX_LZ2024_limit=rX['LZ2024_equiv']['Eup55']['coupling_limit'],
                lambdaHX_XLZD_10_events=rX['projections_10_events']['XLZD_200tyr']['coupling_for_10_events'],
                T_h_over_T_vis_used=rT, dS_over_Q_per_b5e11_squared_linear=lin_coef)
    for lab, b in (('b_XLZD_10ev', b_XLZD), ('b_LZ2024', b_LZ), ('b_5e-11', 5e-11), ('b_ATLAS_max', 5.0501e-3)):
        port[lab] = dict(b=b, theta=th_b * b, Gamma_h_ss_GeV=Gamma_hss(b), BR_inv=Gamma_hss(b) / (GAMMA_H_SM + Gamma_hss(b)),
                         Y_pair=Ypair(b), pairs_per_charge=Ypair(b) / YQ, f_rest=2 * MX * Ypair(b) / D,
                         erasure_needed=1 - 0.01 / (2 * MX * Ypair(b) / D) if 2 * MX * Ypair(b) / D > 0.01 else 0.0,
                         dS_over_Q_linear=lin_coef * (b / 5e-11)**2, S_over_Q_saturated=S_over_Q_saturated(b))
    port['b_max_object_survival'] = bmax
    port['rate_ratio_survival_vs_XLZD10'] = {k: (v / b_XLZD)**2 for k, v in bmax.items()}
    port['events_XLZD_at_b_max'] = {k: 10 * (v / b_XLZD)**2 for k, v in bmax.items()}
    # trilinear portal A phi |H|^2 (allowed: the M1 tilt breaks phi -> -phi): no tree h -> phi phi; phi is produced
    # through h-phi mixing theta = A v_H/(m_h^2 - m_phi^2) like a light Higgs-mixed scalar.  [E] estimate of the
    # dominant IR freeze-in (top/charm Compton and annihilation, q g -> q phi, q qbar -> g phi), with an effective
    # Yukawa y = theta m_q/v_H and a rate Gamma ~ c alpha_s y^2 T (c = 0.1-1 brackets the collision integral):
    #   Y_phi ~ (n_eq/s) Gamma/H at T ~ m_q/3,   Delta S/Q ~ (Y_phi/Y_Q) (3 T / T_h),  T_h = r T_vis.
    rA = [r for r in dd['rows'] if r['Q'] == 2e5 and r['portal'] == 'phi1'][0]
    A_X = rA['projections_10_events']['XLZD_200tyr']['coupling_for_10_events']
    A_LZ = rA['LZ2024_equiv']['Eup55']['coupling_limit']
    tri = {}
    for lab, A in (('A_XLZD_10ev_GeV', A_X), ('A_LZ2024_GeV', A_LZ)):
        th = A * VH / (MH**2 - 0.22635**2)
        rr = {}
        for q, mq, gs_ in (('top', 172.7, 100.0), ('bottom', 4.18, 75.0), ('charm', 1.27, 60.0)):
            for c in (0.1, 1.0):
                y = th * mq / VH
                T = mq / 3
                GH = c * 0.12 * y * y * MPL / (1.66 * math.sqrt(gs_) * T)
                Y = 1.2 / (math.pi**2) * T**3 / (2 * math.pi**2 / 45 * gs_ * T**3) * GH
                rr['%s_c%g' % (q, c)] = dict(Y_phi=Y, dS_over_Q=(Y / YQ) * 3 / rT)
        # h-decay channel via mixing x phi cubic: b_eff v_H ~ theta 6 (lambda_phi - eps) v
        b_eff = th * 6 * (0.03 - 0.002) * 0.924089637 / VH
        tri[lab] = dict(A=A, theta=th, channels=rr, b_eff_h_phiphi=b_eff,
                        dS_over_Q_from_h_decay=lin_coef * (b_eff / 5e-11)**2,
                        dS_over_Q_total_range=[sum(v['dS_over_Q'] for k, v in rr.items() if k.endswith('c0.1')),
                                               sum(v['dS_over_Q'] for k, v in rr.items() if k.endswith('c1'))])
    port['trilinear_portal'] = tri
    out['portal_obligations'] = port
    print(json.dumps(out['portal_obligations']['trilinear_portal'], indent=1, default=float))
    with open(os.path.join(RES, 'late_and_portal.json'), 'w') as f:
        json.dump(out, f, indent=1, default=float)


if __name__ == '__main__':
    main()
