"""qm_micro SKEPTIC (Claude, 2 Nov 2026): collects the skeptic results into skeptic/results/sk_summary.json and the
compact coarsening tables into skeptic/results/sk_coarsen.json.  Own code."""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
SKR = os.path.join(HERE, '..', 'results')
QMR = os.path.join(HERE, '..', '..', 'results')
DQH = os.path.join(HERE, '..', '..', '..', 'dq_heat', 'NOTES.md')
SCR = os.environ.get('SK_SCR', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/qm_micro_skeptic')
H_M1 = 1.342842735e11
YRATIO = 0.991023
EQ_OBJ = 0.775866          # charge-weighted E/Q of the objects at <Q> >~ 1e19 (qm_micro; E(Q) fit)
DELTA = 1 - 0.7758647074595397


def Hi(M):
    return H_M1 * (0.9374 / (M * YRATIO)) ** (2 / 3)


def main():
    out = {}
    # ---- item 1: independent scattering vs qm_micro
    sk = json.load(open(os.path.join(SKR, 'sk_scatter.json')))
    qs = json.load(open(os.path.join(QMR, 'scatter1d.json')))
    T = np.array([r['T'] for r in qs['X_normal']]); s = np.array([r['s'] for r in qs['X_normal']])
    db = np.array([r['destroyed'] for r in qs['Xbar_normal']])
    o = qs['oblique']['X']; th = np.array(o['theta_deg']); Tg = np.array(o['T']); tab = np.array(o['one_minus_R'])
    cmp_rows = []
    for r in sk['X']:
        row = dict(T=r['T'], theta=r['theta'], s_skeptic=r['s'], unitarity=r['unitarity'])
        if r['theta'] == 0:
            i = int(np.argmin(abs(T - r['T'])))
            if abs(T[i] - r['T']) < 1e-9:
                row['s_qm_micro'] = float(s[i]); row['Xbar_destroyed_skeptic'] = r['Xbar_destroyed']; row['Xbar_destroyed_qm_micro'] = float(db[i])
        else:
            i = int(np.argmin(abs(Tg - r['T']))); j = int(np.argmin(abs(th - r['theta'])))
            if abs(Tg[i] - r['T']) < 1e-9:
                row['s_qm_micro'] = float(tab[i, j])
        cmp_rows.append(row)
    out['item1_scatter'] = dict(method='collocation BVP profile + collocation linear BVP with Robin radiation conditions (own code)',
                                sigma_planar=sk['profile']['sigma'], Omega_minus_Omega_b=sk['profile']['Omega_minus_Omega_b'],
                                bulk=sk['bulk'], rows=cmp_rows, threshold=sk['threshold'],
                                dispersion=dict(c_s=sk['dispersion']['c_s_smallk'], gaps=sk['dispersion']['gaps_k0']))
    out['item1_thermal_quadrature'] = json.load(open(os.path.join(SKR, 'sk_thermal_check.json')))['rows']
    bel = json.load(open(os.path.join(SKR, 'sk_beliaev.json')))
    out['item1_beliaev'] = dict(u=bel['u'], c_s=bel['c_s'], e_plus_P=bel['e_plus_P'], rows=bel['rows'])
    # ---- non-isothermal factor and Bose factor
    hf = json.load(open(os.path.join(SKR, 'sk_heat.json')))['rows']
    out['F_h'] = [{k: r[k] for k in ('T', 'F_h_sphi1', 'F_h_sphi0', 'F_h_X_only_Boltzmann_s1', 'F_h_closed_form_NR')} for r in hf]
    bo = json.load(open(os.path.join(SKR, 'sk_bose.json')))['rows']
    out['Bose_factor'] = bo[::3]
    wg = json.load(open(os.path.join(SKR, 'sk_wagner.json')))
    out['wagner_selection'] = [dict(shape=r['shape'], alpha_late=r['alpha_late'], spread=r['alpha_spread'], Rmax_over_Rc=r['Rmax_over_Rc_late']) for r in wg['runs']]
    # ---- coarsening variants
    co = {}
    compact = {}
    for SQ in ('1.5', '2.5', '3.0', '3.2'):
        runs = {}
        for f in ('cheat_%s.json' % SQ, 'cheat2_%s.json' % SQ):
            p = os.path.join(SCR, f)
            if os.path.exists(p):
                for r in json.load(open(p))['runs']:
                    fz = r['freeze']
                    runs[r['case']] = dict(T_f=fz['T'], T_vis_f_GeV=fz['T_vis_GeV'], Q_charge_weighted=fz['Q_charge_weighted'],
                                           Q_capture_eff=fz['Q_capture_eff'], n_obj_per_charge=fz['n_obj_per_charge'],
                                           f_v_at_freeze_imposedT=fz['f_v'], f_v_final_imposedT=r['final']['f_v'])
                    compact.setdefault(SQ, {})[r['case']] = r['rows'][::15] + r['rows'][-1:]
        c = runs['control_computed']['Q_charge_weighted']
        for k, v in runs.items():
            v['Q_ratio_to_control'] = v['Q_charge_weighted'] / c
        co[SQ] = runs
    out['coarsen'] = co
    json.dump(compact, open(os.path.join(SKR, 'sk_coarsen.json'), 'w'), indent=None, default=float, separators=(',', ':'))
    # ---- freeze-out with energy balance; M and H_i
    fr = json.load(open(os.path.join(SKR, 'sk_freeze.json')))['runs']
    qm = json.load(open(os.path.join(QMR, 'summary.json')))['item3_4']
    tabM = {}
    for r in fr:
        SQ = '%.1f' % r['S_over_Q']
        e = r['energy']
        fv = e['final_f_v']
        M = EQ_OBJ + DELTA * fv
        q = qm[SQ]['runs']['Q0_4e11_trap_beliaev' if r['beliaev_trap'] else 'Q0_4e11']
        G = e['Gamma_a_eq_H']
        tabM.setdefault(SQ, []).append(dict(objects=r['objects_from'], beliaev=r['beliaev_trap'], f_v_energy=fv,
                                            f_v_qm_like_b=r['qm_b']['final_f_v'], f_v_qm_like_a=r['qm_a']['final_f_v'],
                                            f_v_qm_micro=[q['f_v_final_decoupled'], q['f_v_final_eqT']],
                                            ratio_to_qm_upper=fv / q['f_v_final_eqT'], ratio_to_qm_lower=fv / q['f_v_final_decoupled'],
                                            at_Gamma_a_eq_H=dict(T_eq=G['T_eq'], T_vis_MeV=G['T_vis_GeV'] * 1e3, T_v=G['T_v'], T_o=G['T_o']),
                                            M_per_charge=M, H_i_GeV=Hi(M), M_qm_micro=q['M_per_charge'], H_i_qm_micro=q['H_i_GeV']))
    out['freeze_energy_balance'] = tabM
    # ---- e(S/Q) from dq_heat's computed R(e) table (log-linear interpolation)
    eS = [(1e-6, 4.95), (1.5e-6, 4.01), (2e-6, 3.51), (3e-6, 2.99), (5e-6, 2.55), (1e-5, 2.20)]
    le = np.log([x[0] for x in eS]); sq = np.array([x[1] for x in eS])
    out['e_from_dq_heat'] = {SQ: (float(math.exp(np.interp(float(SQ), sq[::-1], le[::-1]))) if float(SQ) >= sq.min() else None)
                             for SQ in ('1.5', '2.5', '3.0', '3.2')}
    out['e_qm_micro'] = {SQ: [qm[SQ]['e_band'].get('e'), qm[SQ]['e_broad'].get('e')] for SQ in qm}
    json.dump(out, open(os.path.join(SKR, 'sk_summary.json'), 'w'), indent=1, default=float)
    # print the key tables
    print('coarsen <Q> (charge-weighted at freeze) and ratio to control:')
    for SQ, runs in co.items():
        print('  S/Q=%s ' % SQ + '  '.join('%s %.3e (x%.2f)' % (k, v['Q_charge_weighted'], v['Q_ratio_to_control']) for k, v in runs.items()))
    print('freeze-out with energy balance, M, H_i:')
    for SQ, rows in tabM.items():
        for r in rows:
            print('  S/Q=%s %-24s f_v=%.4f (qm_micro %.4f..%.4f; x%.2f..x%.2f) Gamma_a=H at T_eq=%.4f Tvis=%.2f MeV  M=%.5f H_i=%.4e (qm_micro M %.5f..%.5f)' % (
                SQ, r['objects'] + (' [Beliaev]' if r['beliaev'] else ''), r['f_v_energy'], r['f_v_qm_micro'][0], r['f_v_qm_micro'][1],
                r['ratio_to_qm_upper'], r['ratio_to_qm_lower'], r['at_Gamma_a_eq_H']['T_eq'], r['at_Gamma_a_eq_H']['T_vis_MeV'],
                r['M_per_charge'], r['H_i_GeV'], r['M_qm_micro'][0], r['M_qm_micro'][1]))
    print('e from dq_heat:', out['e_from_dq_heat'], ' qm_micro (band, broad):', out['e_qm_micro'])
    print('wagner:', out['wagner_selection'])


if __name__ == '__main__':
    main()
