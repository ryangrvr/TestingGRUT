"""Run every part of module cosmo_observable and collect the headline numbers in results/headline.json.

    cd research/reality_framework/grut_issue13_wall_history_20261020/modules/cosmo_observable/code
    OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python3 run_all.py
"""
import json
import os
import time
import core as C
import normalization
import history_benchmark
import scale_family
import bubbles_gw
import checks

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')

def main():
    os.makedirs(RES, exist_ok=True)
    times = {}
    t = time.time(); nm = normalization.main(); times['normalization'] = time.time() - t
    t = time.time(); hb = history_benchmark.main(); times['history_benchmark'] = time.time() - t
    t = time.time(); sf = scale_family.main(); times['scale_family'] = time.time() - t
    t = time.time(); bg = bubbles_gw.main(); times['bubbles_gw'] = time.time() - t
    t = time.time(); ck = checks.main(); times['checks'] = time.time() - t

    H = hb['headline']
    gw = bg['GW_domination_epoch']
    head = {
        'labels': 'D=derived closed form, C=computed here, Q=quoted literature, E=order-of-magnitude estimate',
        'sigma_GeV3 [Q, BRIEF/ACTION]': C.SIGMA_BENCH,
        'T_dom_keV (rho_w=rho_rad, A=0.8, backreaction) [C]': H['T_dom_keV_central'],
        'T_dom_keV_A_0.7_0.9 [C]': H['T_dom_keV_A_range'],
        'T_dom_keV all variants (A, backreaction on/off, sigma/t vs 2A sigma H) [C]': H['T_dom_keV_all_variants_A_backreaction_convention'],
        'z_dom [C]': H['z_dom_central'], 't_dom_s [C]': H['t_dom_s_central'],
        'R_H = H0_pred/H0_meas, n=1 frozen (w=-2/3), vs 67/70/74 [C]': H['n=1']['R_H_vs_67_70_74'],
        'R_H n=1 range over A and switch point, vs 70 [C]': H['n=1']['R_H_range_over_A_fswitch_vs70'],
        'R_H, n=2 (w=-1/3), vs 67/70/74 [C]': H['n=2']['R_H_vs_67_70_74'],
        'R_H n=2 range over A and switch point, vs 70 [C]': H['n=2']['R_H_range_over_A_fswitch_vs70'],
        'H0_pred_kmsMpc n=1, n=2 [C]': [H['n=1']['H0_pred_kmsMpc'], H['n=2']['H0_pred_kmsMpc']],
        'age_pred_yr n=1, n=2 [C]': [H['n=1']['age_pred_yr'], H['n=2']['age_pred_yr']],
        'acceleration z-interval n=1 [C]': H['n=1']['acceleration_intervals_z'],
        'acceleration z-interval n=2 [C]': H['n=2']['acceleration_intervals_z'],
        'tuned rho_L<0 to force H0=70: q0, age_yr (n=1) [C]': [H['n=1']['tuned_negative_vacuum_H0=70']['q0'],
                                                               H['n=1']['tuned_negative_vacuum_H0=70']['age_yr']],
        'tuned rho_L<0 to force H0=70: q0, age_yr (n=2) [C]': [H['n=2']['tuned_negative_vacuum_H0=70']['q0'],
                                                               H['n=2']['tuned_negative_vacuum_H0=70']['age_yr']],
        'Omega_w at T=1 MeV, DeltaNeff_equiv [C]': [H['BBN_T=0.001_GeV']['rho_w/rho_tot'], H['BBN_T=0.001_GeV']['DeltaNeff_equiv']],
        'Omega_w at T=0.1 MeV, DeltaNeff_equiv [C]': [H['BBN_T=0.0001_GeV']['rho_w/rho_tot'], H['BBN_T=0.0001_GeV']['DeltaNeff_equiv']],
        'DeltaNeff_equiv A-range at 1 MeV, 0.1 MeV [C]': [H['BBN_T=0.001_GeV']['DeltaNeff_equiv_A_range'],
                                                          H['BBN_T=0.0001_GeV']['DeltaNeff_equiv_A_range']],
        '1-Omega_w at recombination T=%.4f eV, n=1, n=2 [C]' % H['T_rec_eV']: [H['n=1']['1-Omega_w_at_T_rec'], H['n=2']['1-Omega_w_at_T_rec']],
        'stress n=3 (conversion to stable matter at T_dom): Omega h^2 equiv, R_H [C]': [H['stress_n=3_Omega_h2_equiv'], H['stress_n=3_R_H_vs70']],
        'stress n=4 (conversion to radiation at T_dom): DeltaNeff at T_rec [C]': H['stress_n=4_DeltaNeff_at_T_rec'],
        'no-wall reference R_H (rho_L=0) [C]': H['no_walls_reference']['rhoL=0']['R_H_vs70'],
        'GeV^3 -> kg/m^2 [D]': C.GEV3_KG_M2,
        'LMS 1505.03673 printed bound kg/m^2 (R+M, R+M+L) [Q]': [4.22e-9, 3.85e-9],
        'LMS literal SI -> sigma GeV^3 (R+M+L), cube root keV [C]': [nm['LMS_1505.03673']['R+M+Lambda']['literal_sigma_GeV3'],
                                                                     nm['LMS_1505.03673']['R+M+Lambda']['literal_cuberoot_MeV'] * 1e3],
        'LMS reconciled (kg/cm^2 reading) sigma GeV^3 (R+M+L), cube root MeV [C]': [nm['LMS_1505.03673']['R+M+Lambda']['kg_per_cm2_reading_sigma_GeV3'],
                                                                                     nm['LMS_1505.03673']['R+M+Lambda']['kg_per_cm2_reading_cuberoot_MeV']],
        'ZKO-type G sigma t0<1e-5: sigma GeV^3 [E]': nm['ZKO_Saikawa_eq2.21']['sigma_max_GeV3'],
        'Codex G sigma/H0=1e-5: sigma GeV^3, benchmark ratio [C]': [nm['Codex_G_sigma_over_H0']['sigma_GeV3'],
                                                                    nm['Codex_G_sigma_over_H0']['benchmark_over_this']],
        'benchmark / reconciled LMS bound [C]': nm['bounds_table']['LMS_R+M+Lambda_reconciled']['benchmark_over_bound'],
        'v_max(no domination by today) MeV, range over A, rho_L, freezing [C]': sf['v_max_dom_range_MeV'],
        'v_max(no domination) MeV, A=0.8, Planck-like rho_L [C]': sf['v_max_dom']['A=0.800_rhoL=planck_like']['v_max_MeV'],
        'v_max(CMB, reconciled LMS R+M+L) MeV; m_X max MeV [C]': [nm['bounds_table']['LMS_R+M+Lambda_reconciled']['v_max_MeV'],
                                                                  nm['bounds_table']['LMS_R+M+Lambda_reconciled']['m_X_max_MeV_gX0.1']],
        'v_max(CMB, literal LMS) MeV [C]': nm['bounds_table']['LMS_R+M+Lambda_literal_SI']['v_max_MeV'],
        '|r| allowed by m_X|r|<63.93 GeV at m_X,max(CMB) [C]': nm['bounds_table']['LMS_R+M+Lambda_reconciled']['r_max_from_mX_r<63.93GeV'],
        'v giving T_dom = 1 MeV, 0.1 MeV, T_rec (GeV) [C]': [sf['characteristic_v']['T_dom=1MeV']['v_GeV'],
                                                             sf['characteristic_v']['T_dom=0.1MeV']['v_GeV'],
                                                             sf['characteristic_v']['T_dom=T_rec']['v_GeV']],
        'collapsing bubbles: rho/s max for Omega h^2<=0.12 (GeV) [C]': bg['rho_over_s_max_GeV_for_Oh2_0.12'],
        'collapsing bubbles: min T_end (GeV) for f c_b p_- = 1, 0.1, 0.01, 1e-3 [C,E]': [bg['T_end_min_GeV_for_Oh2<=0.12'][k] for k in
                                                                                         ('f_cb_p=1', 'f_cb_p=0.1', 'f_cb_p=0.01', 'f_cb_p=0.001')],
        'GW at domination: Omega_gw(T_dom) vs total, eps=0.7 A=0.8 [E]': gw['cases']['eps=0.7_A=0.8']['Omega_gw_peak_at_T_dom_vs_rho_tot'],
        'GW: rho_gw0/rho_c100 (Omega_GW h^2 equiv) eps=0.7 A=0.8 [E]': gw['cases']['eps=0.7_A=0.8']['rho_gw0_over_rho_c100_(Omega_GW_h2_equiv)'],
        'GW: peak frequency today Hz [E]': gw['f_peak_today_Hz'],
        'checks pass/total': [ck['n_pass'], ck['n_checks']], 'checks failing': ck['failing'],
        'runtime_s': times,
    }
    with open(os.path.join(RES, 'headline.json'), 'w') as f:
        json.dump(head, f, indent=1)
    return head

if __name__ == '__main__':
    t0 = time.time()
    h = main()
    print(json.dumps(h, indent=1))
    print('total runtime %.1f s' % (time.time() - t0))
