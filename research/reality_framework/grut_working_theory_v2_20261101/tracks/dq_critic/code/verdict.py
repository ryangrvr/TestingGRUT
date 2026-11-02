#!/usr/bin/env python3
"""dq_critic item 5 (Claude, own code, 2 Nov 2026): the final verdict table and the open-input sensitivities.

Reads ../results/{sq_branch, qm_benchmark, hi_calibration, higgs_check, exit_check}.json and adds one sensitivity not
yet propagated by any track: the radiation constant K = a^2 H of the supplied inflaton decay history.  R_vec carries K^2
at fixed E_X (dq_heat [C]; K = 1.75-2.97 over (m_I, Gamma_I) = (2-4, 0.25-1) H_i, dq_heat SKEPTIC N2 [C]); E_X of
those histories is not computed, so this is [E]: S/Q(e) and e_crit with R_vec x (K/1.9858)^2.
Output: ../results/verdict.json (and summary.json, the headline numbers).
"""
import json, math, os, sys
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
sys.path.insert(0, HERE)
import sq_branch as SB

L = lambda n: json.load(open(os.path.join(RES, n + '.json')))


def K_sensitivity():
    ht, ex, late, summ = SB.load()
    cf = summ['for_critic']['SQ_of_W']
    fits = late['long_run']['fits']
    Nf = abs([f for f in fits['NX'] if f['model'] == 'a-2' and f['a_min'] > 800][0]['N_inf'])
    b = SB.Budget(ht['portals_vanish_at_MY_benchmark'], Nf, 188737.6639)
    inc = 0.5 * sum(SB.INC)
    out = {}
    for K in (1.75, 1.9858, 2.97):
        k2 = (K / 1.9858052754013398) ** 2
        def sq(e, k2=k2):
            p = b.parts(e)
            return SB.SQ_closed(b.W(k2 * (p['A'] + p['B']) + p['C'] + p['Y']), cf) + inc
        ec = {th: math.exp(brentq(lambda x: sq(math.exp(x)) - SB.SIG[th], math.log(1e-8), math.log(1e-3))) for th in ('bulk', 'Q2e5')}
        out['K=%.4g' % K] = dict(Rvec_factor=k2, SQ_e1e6=sq(1e-6), SQ_e3e6=sq(3e-6), e_crit=ec)
    return out


def main():
    sb, qb, hc, hg, ec = L('sq_branch'), L('qm_benchmark'), L('hi_calibration'), L('higgs_check'), L('exit_check')
    Ks = K_sensitivity()
    rows = {r['e']: r for r in sb['rows']}
    hrows = {r['e']: r for r in hc['rows']}
    ecr = sb['e_crit']
    table = []

    def add(e_label, e, relic, observable, status):
        r, h = rows[e], hrows[e]
        table.append(dict(e=e_label, SQ=r['SQ_final'], SQ_band_outer=r['SQ_band_outer'], branch=r['branch'],
                          SQ_MP=r['SQ_MP'], branch_MP=h['MP']['branch'], SQ_rolling=r['SQ_paths'], relic=relic,
                          H_i=h['H_i'], H_i_dephasing=h['H_i_dephasing'], H_i_MP=h['MP']['H_i'], r=h['r'], r_MP=h['MP']['r'],
                          observable=observable, status=status))

    add('1e-7', 1e-7, 'free X, cold asymmetric scalar gas with a large symmetric part (p ~ 0.8 [E])',
        'none (CDM-like; sigma/m ~ 6e-8 cm^2/g)', 'C on the heat; p and H_i [E]')
    add('1e-6 (declared)', 1e-6, 'static slice: free X, p ~ 1e-3; rolling-inflaton embeddings: bulk Q-matter',
        'none in either case', 'branch conditional on the unsupplied pre-slice inflaton history (f_eff >= 0.469 needed)')
    add('1.5e-6', 1.5e-6, 'free X (Q-matter if portals vanish at M_P)', 'none', 'C')
    add('2e-6', 2e-6, 'free X, margin 0.3-2.2% (Q-matter with M_P portals)', 'none', 'C, marginal')
    add('3e-6 (benchmark)', 3e-6,
        'bulk Q-matter: <Q> = 1.2e20 (charge-weighted), vapour 12% (perfect trapping) - 17% (Beliaev) of the charge, '
        'M/m_X = 0.802-0.815; objects 0.7 M_Pl, R ~ 1 nm', 'P/P_LCDM - 1 >= -4e-5 at k <= 100/Mpc; Planck-mass objects '
        '(gravitational direct detection only, proposal stage)', 'C (heat, S/Q), C/E (relic), E (P(k) below the code floor)')
    add('5e-6', 5e-6, 'Q-matter (binding objects survive above 4.99e-6); <Q> ~ 8e19', 'none in P(k); Planck-mass objects', 'C/E')
    add('1e-5', 1e-5, 'Q-matter; <Q> ~ 6e19 (interpolated toward S/Q = 1.5)', 'none in P(k); Planck-mass objects', 'C/E')
    add('1e-4 (-> floor)', 1e-4, 'Q-matter at the mediator-heat floor S/Q = 1.84', 'none in P(k); objects ~0.2 M_Pl', 'C/E')

    hpts = {(round(p['Mt'], 2), round(p['H'] / 1e9)): p for p in hg['points']}
    higgs = dict(
        D60_points=[dict(Mt=p['Mt'], H=p['H'], D60=p['D60'], D60_cH05=p.get('D60_cH=0.5'), D60_cH2=p.get('D60_cH=2.0')) for p in hg['points']],
        Mt_crit={k: v['Mt_crit'] for k, v in hg.items() if k.startswith('mtcrit')},
        Mt_PDG2024=(172.57, 0.29),
        Mt_crit_sigma={k: (v['Mt_crit'] - 172.57) / 0.29 for k, v in hg.items() if k.startswith('mtcrit')},
        criterion_spread_GeV=0.09,
        inflation_xi_min=(0.044, 0.055), grav_lower_edge_central=(0.072, 0.085), grav_lower_edge_plus1sigma=(0.096, 0.112),
        plasma_deadline_tau=(8, 10))
    summary = dict(
        e_crit=ecr, f_crit_e1e6=sb['f_crit_e1e6'], SQ_e1e6=rows[1e-6]['SQ_final'], SQ_e1e6_band=rows[1e-6]['SQ_band_outer'],
        SQ_e2e6=rows[2e-6]['SQ_final'], SQ_e2e6_band=rows[2e-6]['SQ_band_outer'], SQ_e1e6_rolling=rows[1e-6]['SQ_paths'],
        SQ_e1e6_MP=rows[1e-6]['SQ_MP'], K_sensitivity=Ks, e_of_SQ=qb['e_of_SQ'], benchmark=qb['benchmark'],
        relic=qb['relic'], Pk_perfect=qb['Pk_perfect']['at_k'], Pk_beliaev=qb['Pk_beliaev']['at_k'],
        m_eff=qb['m_eff_extrapolated_keV'], probe_reach=qb['probe_reach'], objects=qb['objects'],
        lambda_fs=qb['lambda_fs_check'], H_i={str(r['e']): dict(branch=r['branch'], H_i=r['H_i'], MP=r['MP']['H_i'], r=r['r']) for r in hc['rows']},
        e1e6_rolling_H_i=hc['e1e6_rolling_paths'], higgs=higgs,
        exit_check=[dict(s=r['s'], h_r=r['h_r'], f=r['f'], dq_heat_estimate=r['dq_heat_estimate']) for r in ec['rows']])
    json.dump(dict(method=__doc__, table=table, summary=summary), open(os.path.join(RES, 'verdict.json'), 'w'), indent=1, default=float)
    json.dump(summary, open(os.path.join(RES, 'summary.json'), 'w'), indent=1, default=float)
    for t in table:
        print('%-18s S/Q %.3f [%.3f,%.3f] %-16s | M_P %.3f %-8s | rolling %s | H_i %.4e-%.4e (M_P %.4e-%.4e) r %.2e' % (
            t['e'], t['SQ'], *t['SQ_band_outer'], t['branch'], t['SQ_MP'], t['branch_MP'],
            '/'.join('%.2f' % v for v in t['SQ_rolling'].values()), *t['H_i'], *t['H_i_MP'], t['r'][0]))
    print('K sensitivity:', json.dumps(Ks, indent=None, default=float))
    print('Higgs:', json.dumps(higgs, default=float))


if __name__ == '__main__':
    main()
