#!/usr/bin/env python3
"""inf_chain (Claude, own code): per-benchmark table -> ../results/summary.json; also slims the vector JSONs (per-mode
spectra kept only at e = 1e-6).  Labels as in NOTES.md."""
import json, math
import ch_common as ch

R = ch.RES


def load(n):
    return json.load(open(R / n))


def slim():
    for f in R.glob('vectors_*.json'):
        d = json.load(open(f))
        for x in d['by_e']:
            if abs(x['e'] - 1e-6) > 1e-12 and 'spectrum' in x:
                del x['spectrum']
        json.dump(d, open(f, 'w'), indent=0, default=float)


def row(key, vtag, pot):
    rel = load('release_%s.json' % key)['result']
    relf = load('release_%s.json' % key)
    c = load('chain_%s.json' % vtag)
    v = load('vectors_%s.json' % vtag)
    v6 = [x for x in v['by_e'] if abs(x['e'] - 1e-6) < 1e-12][0]
    rl = c['relic_e1e6']
    out = dict(key=key, slice='S1 (1 e-fold before eps_H = 1)', holder_cal=pot.get('holder_cal', 'S1'), Gamma_I_over_Hi=pot['gam'],
               m_I_over_Hi=relf['bench']['m_I_over_Hi'], g_I=relf['bench']['g_I'], alpha_I=relf['bench']['alpha_I'],
               H_i_potential_GeV=pot['Hi'], H_star_GeV=pot['Hstar'], n_s=pot['ns'], r=pot['r'], N_star=pot['Nstar'],
               release=dict(N_X=rel['N_X'], N_X_err=rel['N_X_err'], E_X=rel['E_X'], eps=rel['eps'], E_over_Ec=rel['E_over_Ec'], K=c['K'],
                            N_YL=rel['N_YL'], N_YH=rel['N_YH'], N_sigma=rel['N_sigma'], gauss=rel['gauss'], hidden_work=rel['hidden_work'],
                            H_consistency=rel['H_consistency'], std_slope_NX=relf['long']['std_slope_NX'],
                            normalisation='a = 1 at S1 (divide comoving numbers by e^3, E_X by e^4, K by e^2 for a = 1 at eps_H = 1)'),
               prep=dict(amps=relf['prep']['amps'], F=relf['prep']['F'], lightest_m2=relf['prep']['eig'][1]),
               vectors_e1e6=dict(N_A=v6['N_A'], f_eff=v6['f_eff'], N_A_k_ge_1=v6['N_A_k_ge_1'], deposit_over_EX=v6['deposit_per_vector'] / rel['E_X'],
                                 R_vec=v6['R_vec']),
               heat_e1e6=dict(R_vec=c['e1e6']['R_vec'], R_Y=c['e1e6']['R_Y'], R_casc=c['e1e6']['R_casc'], R_YY=c['e1e6']['R_YY'],
                              R_total=c['e1e6']['R_total'], W=c['e1e6']['W']),
               S_per_Q_e1e6=c['e1e6']['S_per_Q'], S_per_Q_e1e6_noYY=c['e1e6']['S_per_Q_noYY'], S_per_Q_e_inf=c['S_per_Q_e_inf'],
               e_crit_bulk=c['e_crit_bulk'], e_crit_Q2e5=c['e_crit_Q2e5'],
               branch_e1e6=rl['branch'], p=rl.get('p'), M=rl['M'], Y_over_YQ_at_HM1=c['Y_over_YQ_at_HM1'],
               H_i_req_GeV=rl['H_i_req_GeV'], H_i_req_pmap_x2_GeV=rl.get('H_i_req_pmap_x2_GeV'),
               H_i_req_noYY_GeV=c['relic_e1e6_noYY']['H_i_req_GeV'], consistency_ratio=rl['H_i_req_GeV'] / pot['Hi'])
    return out


def main():
    slim()
    bp = lambda n: json.load(open(ch.TRK / 'inf_potential' / 'results' / ('bench_%s.json' % n)))['by_Hi'][0]
    H6, T1 = bp('HT6t'), bp('T1')
    rc_h, rc_t = load('recal_HT6t_cal1.json'), load('recal_T1_cal1.json')
    tb = {k: load('calib_%s.json' % k)['steps'][-1] for k in ('TB1', 'TB2')}
    pots = {
        'HT6t': dict(gam=0.5, Hi=H6['H_i_GeV'], Hstar=H6['Hstar_GeV'], ns=H6['ns'], r=H6['r'], Nstar=H6['Nstar']),
        'HT6t_cal1': dict(gam=0.5, Hi=rc_h['H_i_GeV'], Hstar=rc_h['Hstar_GeV'], ns=rc_h['ns'], r=rc_h['r'], Nstar=rc_h['Nstar']),
        'HT6t_Send': dict(gam=0.5, Hi=H6['H_i_GeV'], Hstar=H6['Hstar_GeV'], ns=H6['ns'], r=H6['r'], Nstar=H6['Nstar'], holder_cal='eps_H = 1'),
        'T1': dict(gam=0.5, Hi=T1['H_i_GeV'], Hstar=T1['Hstar_GeV'], ns=T1['ns'], r=T1['r'], Nstar=T1['Nstar']),
        'T1_cal1': dict(gam=0.5, Hi=rc_t['H_i_GeV'], Hstar=rc_t['Hstar_GeV'], ns=rc_t['ns'], r=rc_t['r'], Nstar=rc_t['Nstar']),
        'TB1': dict(gam=tb['TB1']['Gamma_cal_over_Hi'], Hi=tb['TB1']['H_i_new_GeV'], Hstar=tb['TB1']['H_star_GeV'], ns=tb['TB1']['ns'],
                    r=tb['TB1']['r'], Nstar=tb['TB1']['Nstar_new']),
        'TB2': dict(gam=tb['TB2']['Gamma_cal_over_Hi'], Hi=tb['TB2']['H_i_new_GeV'], Hstar=tb['TB2']['H_star_GeV'], ns=tb['TB2']['ns'],
                    r=tb['TB2']['r'], Nstar=tb['TB2']['Nstar_new']),
    }
    vt = {'TB1': 'TB1_g7.512e-09', 'TB2': 'TB2_g2.113e-09'}
    rows = {k: row(k, vt.get(k, k), p) for k, p in pots.items()}
    for k in ('TB1', 'TB2'):
        rows[k]['reheating'] = {kk: tb[k][kk] for kk in ('Gamma_cal_over_Hi', 'Gamma_over_mI', 'T_RH_GeV', 'N_matter_efolds', 'K', 'y_majorana_N1')}
    out = dict(rows=rows,
               iteration=dict(HT6t=[pots['HT6t']['Hi'], rows['HT6t']['H_i_req_GeV'], rows['HT6t_cal1']['H_i_req_GeV']],
                              T1=[pots['T1']['Hi'], rows['T1']['H_i_req_GeV'], rows['T1_cal1']['H_i_req_GeV']]),
               sq_check=load('sq_check.json')['max_abs_rel_converged'],
               vector_validation_M1r=[dict(e=x['e'], ratio=x['ratio_N_A_k_le_1']) for x in load('vectors_M1r.json')['validation_vs_dq_heat']],
               pmap=dict(rows=load('pmap_kin.json')['rows'], self_consistent={k: dict(p_fixed=v['p_fixed_yield'], p_sc=v['p_self_consistent'])
                                                                             for k, v in load('pmap_kin.json')['self_consistent'].items()}))
    json.dump(out, open(R / 'summary.json', 'w'), indent=1, default=float)
    for k, r in rows.items():
        print('%-10s m_I=%7.1f g_I=%.2e | NX=%.4g EX=%.4g eps=%.3f E/Ec=%.2f K=%.4g NYL=%.4g NYH=%.4g | f_eff=%.3f R(1e-6): vec %.3g Y %.3g casc %.3g YY %.3g | S/Q(1e-6)=%.4g (noYY %.4g) inf=%.3g e_crit=%s | p=%.3g M=%.3g Y/YQ=%.4g H_req=%.3e H_pot=%.3e ratio=%.3f H*=%.3e r=%.3e ns=%.4f' % (
            k, r['m_I_over_Hi'], r['g_I'], r['release']['N_X'], r['release']['E_X'], r['release']['eps'], r['release']['E_over_Ec'], r['release']['K'],
            r['release']['N_YL'], r['release']['N_YH'], r['vectors_e1e6']['f_eff'], r['heat_e1e6']['R_vec'], r['heat_e1e6']['R_Y'], r['heat_e1e6']['R_casc'],
            r['heat_e1e6']['R_YY'], r['S_per_Q_e1e6'], r['S_per_Q_e1e6_noYY'], r['S_per_Q_e_inf'], r['e_crit_bulk'], r['p'] or float('nan'), r['M'],
            r['Y_over_YQ_at_HM1'], r['H_i_req_GeV'], r['H_i_potential_GeV'], r['consistency_ratio'], r['H_star_GeV'], r['r'], r['n_s']))


if __name__ == '__main__':
    main()
