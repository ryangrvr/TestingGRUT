"""Parts 1-2: background with a persistent unbiased Z2 wall network at the benchmark tension; T_dom, H0 given T_CMB,
age, acceleration, wall fractions at BBN and recombination, the observable R_H. Writes results/history_benchmark.json.
"""
import json
import math
import os
import time
import core as C
import observables as O

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')

def clean(d):
    return {k: v for k, v in d.items() if k not in ('sol1', 'sol2', 'rho_w2', 'rw1')}

def main():
    t0 = time.time()
    sig = C.SIGMA_BENCH
    out = {'module': 'cosmo_observable', 'part': '1-2 benchmark history',
           'inputs': {
               'sigma_GeV3': sig, 'sigma_provenance': 'ACTION tree planar tension (2 sqrt(2 lam_phi)/3) v^3, v=1000 GeV, lam_phi=0.01',
               'A_central': C.A_CENTRAL, 'A_range': [C.A_LOW, C.A_HIGH],
               'A_provenance': 'LIT Hiramatsu-Kawasaki-Saikawa arXiv:1309.5001 eqs.(4),(29): A = rho_wall t/sigma_wall = 0.8+-0.1',
               'T_CMB_K': C.T_CMB_K, 'omega_b': C.OMEGA_B_H2, 'omega_c': C.OMEGA_C_H2,
               'calibration_provenance': 'Planck 2018 arXiv:1807.06209 Table 2 TT,TE,EE+lowE+lensing; T_CMB Fixsen 2009 arXiv:0911.1955',
               'gstar_provenance': 'Saikawa-Shirai arXiv:1803.01038 App. C fits; g_rho0=%.4f g_s0=%.4f' % (C.G_RHO0, C.G_S0),
               'M_P_GeV': C.M_P, 'T_start_GeV': 100.0,
               'rho_L_base': 0.0, 'rho_L_planck_like_GeV4': C.RHO_L_PLANCK,
               'H0_meas_range_kmsMpc_retrospective': list(C.H0_MEAS_RANGE)}}

    # ---- T_dom (rho_w = rho_rad) for A = 0.7, 0.8, 0.9, with and without backreaction -------------------------
    tdom = {}
    for A in (C.A_LOW, C.A_CENTRAL, C.A_HIGH):
        for br in (True, False):
            h = C.history(sig, A=A, n_post=1.0, f_switch=1.0, backreaction=br, T_end=C.T0 * 1.0001)
            key = 'A=%.1f_backreaction=%s' % (A, br)
            tdom[key] = {'T_dom_GeV': h['T_switch_GeV'], 'T_dom_keV': h['T_switch_GeV'] * 1e6,
                         't_dom_s': h['t_switch_GeV-1'] * C.HBAR_GEV_S, 'z_dom': O.z_of_T(h['T_switch_GeV'])}
        hH = C.history(sig, A=A, n_post=1.0, f_switch=1.0, backreaction=True, T_end=C.T0 * 1.0001, convention='H')
        tdom['A=%.1f_convention=sigmaH_AH=2A_backreaction=True' % A] = {
            'T_dom_GeV': hH['T_switch_GeV'], 'T_dom_keV': hH['T_switch_GeV'] * 1e6,
            't_dom_s': hH['t_switch_GeV-1'] * C.HBAR_GEV_S, 'z_dom': O.z_of_T(hH['T_switch_GeV'])}
        # analytic radiation-era estimate (no backreaction, t = 1/(2H)): H_dom = 2 A sigma / (3 M_P^2)
        Hd = 2 * A * sig / (3 * C.M_P**2)
        tdom['A=%.1f_analytic' % A] = {'H_dom_GeV': Hd, 't_dom_s_Saikawa_eq2.19': 3 * C.M_P**2 / (4 * A * sig) * C.HBAR_GEV_S,
                                       't_dom_s_Saikawa_formula_2.93e3/A/sigmaTeV3': 2.93e3 / A / (sig / 1e9)}
    out['T_dom'] = tdom

    # ---- full histories ----------------------------------------------------------------------------------------
    runs = {}
    for n in (1.0, 2.0, 3.0, 4.0):
        for A in (C.A_LOW, C.A_CENTRAL, C.A_HIGH):
            for fsw in ((0.1, 0.3, 1.0, 3.0) if A == C.A_CENTRAL else (1.0,)):
                for rl_name, rl in (('rhoL=0', 0.0), ('rhoL=planck_like', C.RHO_L_PLANCK)):
                    if rl_name != 'rhoL=0' and not (A == C.A_CENTRAL and fsw == 1.0):
                        continue
                    h = C.history(sig, A=A, n_post=n, f_switch=fsw, rho_L=rl)
                    s = O.summarize(h, with_accel=(A == C.A_CENTRAL and fsw == 1.0))
                    key = 'n=%d_A=%.1f_fswitch=%g_%s' % (n, A, fsw, rl_name)
                    entry = {'n_post': n, 'A': A, 'f_switch': fsw, 'rho_L': rl_name}
                    entry.update(s)
                    if n in (1.0, 2.0) and rl == 0.0 and fsw == 1.0:
                        entry['tuned_negative_vacuum_H0=70'] = O.tuned_negative_vacuum(h, 70.0)
                    if n == 3.0:
                        entry['stress_matter_Omega_h2_equiv'] = s['rho_w0_GeV4'] / C.RHO_C100
                    if n == 4.0:
                        entry['stress_radiation_DeltaNeff_at_T_rec'] = s['fractions']['T=%.4g_GeV' % O.T_REC]['DeltaNeff_equiv']
                    runs[key] = entry
    # convention variant (rho_w = 2 A sigma H) for the observable
    for n in (1.0, 2.0):
        h = C.history(sig, A=C.A_CENTRAL, n_post=n, f_switch=1.0, convention='H')
        s = O.summarize(h, with_accel=False)
        e = {'n_post': n, 'A': C.A_CENTRAL, 'f_switch': 1.0, 'rho_L': 'rhoL=0', 'convention': 'sigmaH'}
        e.update(s)
        runs['n=%d_A=0.8_fswitch=1_rhoL=0_convention=sigmaH' % n] = e
    out['runs'] = runs

    # ---- headline extraction -------------------------------------------------------------------------------------
    def rng(n, field):
        vals = []
        for k, e in runs.items():
            if e['n_post'] == n and e['rho_L'] == 'rhoL=0':
                vals.append(field(e))
        return [min(vals), max(vals)]
    base = {n: runs['n=%d_A=0.8_fswitch=1_rhoL=0' % n] for n in (1, 2, 3, 4)}
    head = {
        'T_dom_keV_central': tdom['A=0.8_backreaction=True']['T_dom_keV'],
        'T_dom_keV_A_range': [tdom['A=0.7_backreaction=True']['T_dom_keV'], tdom['A=0.9_backreaction=True']['T_dom_keV']],
        'z_dom_central': tdom['A=0.8_backreaction=True']['z_dom'],
        'T_dom_keV_all_variants_A_backreaction_convention': [min(e['T_dom_keV'] for e in tdom.values() if 'T_dom_keV' in e),
                                                              max(e['T_dom_keV'] for e in tdom.values() if 'T_dom_keV' in e)],
        't_dom_s_central': tdom['A=0.8_backreaction=True']['t_dom_s'],
    }
    for n in (1, 2):
        b = base[n]
        head['n=%d' % n] = {
            'H0_pred_kmsMpc': b['H0_pred_kmsMpc'],
            'R_H_vs_67_70_74': [b['R_H']['67'], b['R_H']['70'], b['R_H']['74']],
            'R_H_range_over_A_fswitch_vs70': rng(n, lambda e: e['R_H']['70']),
            'age_pred_yr': b['age_pred_yr'],
            'q0': b['q0'],
            'acceleration_intervals_z': [[iv['z_begin'], iv['z_end']] for iv in b['acceleration_intervals']],
            'Omega_w_at_T_rec': b['fractions']['T=%.4g_GeV' % O.T_REC]['rho_w/rho_tot'],
            '1-Omega_w_at_T_rec': b['fractions']['T=%.4g_GeV' % O.T_REC]['1-Omega_w'],
            'tuned_negative_vacuum_H0=70': b['tuned_negative_vacuum_H0=70'],
        }
    b1 = base[1]
    for T in O.T_BBN:
        k = 'T=%.4g_GeV' % T
        head['BBN_' + k] = {'rho_w/rho_tot': b1['fractions'][k]['rho_w/rho_tot'],
                            'DeltaNeff_equiv': b1['fractions'][k]['DeltaNeff_equiv'],
                            'DeltaNeff_equiv_A_range': rng(1, lambda e, k=k: e['fractions'][k]['DeltaNeff_equiv'])}
    head['stress_n=3_Omega_h2_equiv'] = base[3]['stress_matter_Omega_h2_equiv']
    head['stress_n=3_R_H_vs70'] = base[3]['R_H']['70']
    head['stress_n=4_DeltaNeff_at_T_rec'] = base[4]['stress_radiation_DeltaNeff_at_T_rec']
    head['no_walls_reference'] = {}
    for rl_name, rl in (('rhoL=0', 0.0), ('rhoL=planck_like', C.RHO_L_PLANCK)):
        H0 = math.sqrt((C.RHO_R0 + C.RHO_M0 + rl) / (3 * C.M_P**2))
        head['no_walls_reference'][rl_name] = {'H0_kmsMpc': O.H_kmsmpc(H0), 'R_H_vs70': O.H_kmsmpc(H0) / 70.0}
    head['T_rec_eV'] = O.T_REC * 1e9
    # gravitational scales of the benchmark wall (estimates; Vilenkin 1983 / Ipser-Sikivie 1984 via secondary sources)
    hb = C.history(sig, A=C.A_CENTRAL, n_post=1.0, f_switch=1.0)
    td = hb['t_switch_GeV-1']
    Td = hb['T_switch_GeV']
    rho_d = C.rho_rad(Td) + C.rho_matter(Td) + hb['rho_w_switch']
    Hd = math.sqrt(rho_d / (3 * C.M_P**2))
    GEV_G = 1.78266192e-24     # grams per GeV
    MSUN_G = 1.98847e33
    M_hor = 4 * math.pi / 3 * rho_d / Hd**3
    M_wall = 4 * math.pi * sig * td**2
    head['gravitational_scales'] = {
        't_sigma_s = 1/(2 pi G sigma)': 1.0 / (2 * math.pi * C.G_N * sig) * C.HBAR_GEV_S,
        'R_selfgrav_s = M_P^2/sigma (r_s(R) = R for a closed spherical wall)': C.M_P**2 / sig * C.HBAR_GEV_S,
        'R_bh_min_m = M_P/sqrt(sigma m_h) (r_s exceeds wall width 1/m_h)': C.M_P / math.sqrt(sig * C.M_H) * C.HBARC_GEV_M,
        'horizon_mass_at_t_dom_Msun': M_hor * GEV_G / MSUN_G,
        'closed_wall_mass_4pi_sigma_tdom^2_Msun': M_wall * GEV_G / MSUN_G,
        'note': 'estimate: horizon-size closed walls near t_dom are self-gravitating; collapse may form black holes '
                '(stable, matter-like: covered by the n=3 stress test). No PBH formation fraction is computed.'}
    out['headline'] = head
    out['runtime_s'] = time.time() - t0
    with open(os.path.join(RES, 'history_benchmark.json'), 'w') as f:
        json.dump(out, f, indent=1)
    return out

if __name__ == '__main__':
    o = main()
    print(json.dumps(o['headline'], indent=1))
    print('runtime %.1f s' % o['runtime_s'])
