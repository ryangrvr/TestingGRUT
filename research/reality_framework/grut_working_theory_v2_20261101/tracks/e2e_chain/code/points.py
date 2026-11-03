#!/usr/bin/env python3
"""e2e_chain items 2-4 (Claude): assemble M2-X and M2-Q and export results/m2_points.json.

Inputs (this track's results, all computed here): scan_<tag>.json (chain runs), fixpoints.json, cmb.json,
handover_<tag>.json (kinetic K from the lattice moments).  Relic and observables: e2e_relic relic() / observables()
(read-only).  Isocurvature: D = d ln Omega / d ln c0 from the c0 (1 +- 0.01) runs at the node nearest the point
(Omega ~ M |N_X| K^-3/2 at fixed H_i; c0 ~ |g_h| s_*), s_* > |D| H_i / (2 pi sqrt(beta A_s/(1 - beta))), beta = 0.038
(Planck 2018 uncorrelated CDM isocurvature, arXiv:1807.06211 [R, as quoted by mf_ideal/mf_chain]).
Naturalness [D/E]: Starobinsky-Yokoyama equilibrium (arXiv:astro-ph/9407016, verified this session) of a light s with
inflationary mass m: <s^2> = 3 H^4/(8 pi^2 m^2), relaxation ~ 3 H^2/m^2 e-folds; the residual masses are mf_hybrid's
[Q]: holder loop 2.5e-7 H^2, HFa waterfall-induced <= 3e-5 H^2; HFb xi_s R term 4.5 eps_H H^2 (eps_H <= 1e-7 at the slice).
"""
import json, math
import numpy as np
from scipy.stats import norm
import e2c_common as E
import relic as RL                                        # e2e_relic, read-only

CL = 'cL0_vanish_at_MY'
BETA, AS = 0.038, 2.1e-9


def run(tag, key):
    return E.load(E.RES / ('scan_%s.json' % tag))['runs'][key]


def iso_D(tag, lo, hi, cl=CL):
    a, b = run(tag, lo), run(tag, hi)
    lnO = lambda o: math.log(o['by_e'][cl]['e1e-6']['M'] * abs(o['N_X']) * o['K'] ** -1.5)
    dl = math.log(b['c0'] / a['c0'])
    return dict(D=(lnO(b) - lnO(a)) / dl, D_charge_only=(math.log(abs(b['N_X'])) - math.log(abs(a['N_X']))) / dl,
                dSQ_dlnc0=(b['by_e'][cl]['e1e-6']['S_per_Q'] - a['by_e'][cl]['e1e-6']['S_per_Q']) / dl, node=tag, c0=[a['c0'], b['c0']])


def sstar_min(D, H):
    return abs(D) * H / (2 * math.pi * math.sqrt(BETA * AS / (1 - BETA)))


def beta_at(D, H, s):
    x = (D * H / (2 * math.pi * s)) ** 2 / AS
    return x / (1 + x)


def point(name, closure, key, fp, cmbkey, iso, nodes, spread):
    H = fp['H_star_GeV']
    lo, hi = run(nodes[0], key), run(nodes[1], key)
    w = fp['w']
    near = hi if w > 0.5 else lo
    b = near['by_e'][CL]
    rel = RL.relic(fp['SQ_1e6'], H_i=H, e=1e-6)
    obs = RL.observables(rel)
    cmb = E.load(E.RES / 'cmb.json')[cmbkey]
    sm = sstar_min(iso['D'], H)
    sp = near['spectator']
    # s energy fraction scales as s_*^2 at fixed history (rho_s ~ m_s^2 s^2 / 2)
    sst_bench = sp['s_star_over_Hi'] * near['H_i_GeV']
    scale = (max(sm, sst_bench) / sst_bench) ** 2
    f_dec = sp['f_s_at_decay'] * scale
    out = dict(
        name=name, closure=closure, holder=dict(c0=near['c0'], delta_deg=near['delta_deg'], mu_eff_slice=near['mu_eff'],
                                                abs_mu_eff=near['abs_mu_eff'], arg_mu_eff_deg=near['arg_mu_eff_deg'],
                                                g_h_s_star_over_H_i=12.30 * near['c0']),
        self_consistent=dict(H_i_GeV=H, interpolation_weight=w, extrapolated=fp['extrapolated'], nodes_H_GeV=[lo['H_i_GeV'], hi['H_i_GeV']],
                             dln_ratio_dlnH=fp['slope_dlnratio_dlnH'], K=fp['K'], V0_quarter_GeV=(3 * E.MP_GEV ** 2 * H ** 2) ** 0.25,
                             rule='mf_chain iterate.py V0 rescaling; root of ln(H_req/H_i) interpolated in ln H between two computed final histories'),
        cmb=dict(Nstar=cmb['Nstar'], r=cmb['r_constH'], capacity_Nmax=cmb['capacity']['Nmax'], capacity_ok=cmb['capacity_ok'],
                 MS=dict(n_s=cmb['MS']['ms']['ns'], A_s=cmb['MS']['ms']['As'], alpha_s=cmb['MS']['ms']['alpha_s'], r=cmb['MS']['ms']['r']),
                 note='n_s is the retuned target (HFa Planck 0.9649), not an output'),
        release=dict(eps=fp['eps'], E_over_Ec=fp['E_over_Ec'], N_YL=fp['N_YL'], N_YH=fp['N_YH'], N_sigma=near['N_sigma'],
                     checks_nearest_node=near['checks']),
        ledger_1e6=dict(R_vec=fp['R_vec'], R_Y=fp['R_Y'], R_casc=fp['R_casc'], R_leak_cL0=fp['R_leak'], dR_grav=fp['dR_grav'],
                        R_leak_cL_1_16pi2_nearest_node=near['by_e']['cL_1/16pi^2']['e1e-6']['R_leak'],
                        SQ_cL_1_16pi2_nearest_node=near['by_e']['cL_1/16pi^2']['e1e-6']['S_per_Q'],
                        c_L_statement='c_L = 0: induced holder portals vanish at M_Y with no threshold piece (declared); 1/16pi^2 shown as a variant'),
        S_per_Q=dict(e1e6=fp['SQ_1e6'], e3e6=fp['SQ_3e6'], e_inf=fp['SQ_inf'], branch_1e6=fp['branch_1e6'], branch_3e6=fp['branch_3e6'],
                     e_crit_bulk=fp['e_crit_bulk'], e_crit_Q2e5=fp['e_crit_Q2e5'], SQ_of_e_nearest_node=b['SQ_of_e']),
        relic=dict(branch=rel['branch'], M_over_mX=rel.get('M_central'), p=rel.get('p'), Q_mean=rel.get('Q_mean'), Q_mean_range=rel.get('Q_mean_range'),
                   f_v=rel.get('f_v'), Y_Q=rel['Y_Q'], m_A_GeV=rel.get('m_A_GeV'), observables=obs),
        heavy_remnants=dict(N_YL=fp['N_YL'], N_YH=fp['N_YH'], N_sigma=near['N_sigma'], R_Y=fp['R_Y'],
                            note='Y_L/Y_H decay deep in radiation (Gamma_L = 1.0e-7 H_i); their heat R_Y is in the ledger; none survives'),
        gravitational=dict(dR_grav=fp['dR_grav'], **near['grav']),
        isocurvature=dict(iso, s_star_min_GeV=sm, beta_at_bench_s_star=beta_at(iso['D'], H, sst_bench), bench_s_star_GeV=sst_bench,
                          s_star_over_MP_min=sm / E.MP_GEV),
        spectator=dict(f_s_over_a_bench=sp['f_s_over_a'], Gamma_s_over_Hi=sp['Gamma_s'], a_decay=sp['a_decay'],
                       f_s_at_decay_at_s_star_min=f_dec, entropy_ratio_at_s_star_min=(1 + 0.8862 * f_dec) ** 0.75,
                       H_req_shift=((1 + 0.8862 * f_dec) ** 0.75) ** (2 / 3) - 1,
                       decay='N_s -> L H (model choice; otherwise N_s is dark radiation)'),
        spectator_spread=spread)
    return out


def natural(H, sstar_min, c0, window):
    rows = {}
    for lab, m2 in (('holder_loop_2.5e-7', 2.5e-7), ('HFa_waterfall_max_3e-5', 3e-5), ('HFb_xiR_4.5e-7', 4.5e-7), ('HFb_xiR_4.5e-8', 4.5e-8)):
        srms = math.sqrt(3 / (8 * math.pi ** 2 * m2))          # in H units
        # (i) g_h fixed by the isocurvature floor at the window centre: c0 = c0_centre s/s_min
        z = sstar_min / (H * srms)
        p_fixed = 2 * (norm.sf(z * window[0] / c0) - norm.sf(z * window[1] / c0))
        # (ii) g_h free: best case over the rms (window probability maximised), but then s ~ s_rms violates isocurvature
        best = max((2 * (norm.sf(window[0] / sc) - norm.sf(window[1] / sc)), sc) for sc in np.geomspace(0.1, 10, 400))
        rows[lab] = dict(m2_over_H2=m2, s_rms_over_H=srms, s_rms_GeV=srms * H, N_relax_efolds=3 / m2,
                         ratio_s_star_min_over_s_rms=z, P_window_gh_fixed_by_isocurvature='exp(-%.3g)' % (0.5 * (z * window[0] / c0) ** 2) if p_fixed == 0 else p_fixed,
                         P_window_best_case_gh_free=best[0], sigma_c0_best=best[1],
                         beta_iso_at_s_rms=beta_at(1.0, H, srms * H))
    return rows


if __name__ == '__main__':
    fps = E.load(E.RES / 'fixpoints.json')[CL]
    # spectator-spread checks [C]
    def sq(tag, key):
        return run(tag, key)['by_e'][CL]['e1e-6']['S_per_Q'], run(tag, key)['by_e'][CL]['e_crit_bulk']
    isoX = iso_D('hHFa_GM', '1.534:+0', '1.565:+0')
    isoQ = iso_D('hH84_P0', '1.435:+0', '1.464:+0')
    lat = E.load(E.TRK / 'e2e_lattice' / 'results' / 'hist_HFa.json')
    ds = 6e-4                                                   # lattice delta s / s_* (max, vacuum noise) [Q e2e_lattice]
    spreadX = dict(physical=dict(delta_s_over_s_star=ds, dSQ=isoX['dSQ_dlnc0'] * ds,
                                 note='m_s^2 inhomogeneity sits at k ~ m_chi >> m_s: s sees the box mean; only delta s enters the holder'),
                   separate_universe_pm1sigma_ms2=dict(zip(('+1', '-1'), (sq('hHFa_GM_sp1', '1.550:+0'), sq('hHFa_GM_sp-1', '1.550:+0')))),
                   mean_ms2_pm3pct=dict(zip(('x1.03', 'x0.97'), (sq('hHFa_GM_ms1.03', '1.550:+0'), sq('hHFa_GM_ms0.97', '1.550:+0')))),
                   base=sq('hHFa_GM', '1.550:+0'))
    spreadQ = dict(physical=dict(delta_s_over_s_star=ds, dSQ=isoQ['dSQ_dlnc0'] * ds),
                   separate_universe_pm1sigma_ms2=dict(zip(('+1', '-1'), (sq('hH84_P0_sp1', '1.450:+0'), sq('hH84_P0_sp-1', '1.450:+0')))),
                   mean_ms2_pm3pct={'x1.03': sq('hH84_P0_ms1.03', '1.450:+0')}, base=sq('hH84_P0', '1.450:+0'))
    MX = point('M2-X', 'GM (box-mean Pauli saturation; agrees with the kinetic homogeneous blocking at t_sw)', '1.550:+0', fps['GM']['1.550:+0'],
               'HFa:1.1229e11:15.66', isoX, ('hH84_GM', 'hHFa_GM'), spreadX)
    MQ = point('M2-Q', 'P0 (local per-site Pauli blocking)', '1.450:+0', fps['P0']['1.450:+0'], 'HFa:9.8864e10:9.05', isoQ, ('hH84_P0', 'hHFa_P0'), spreadQ)
    cmbb = E.load(E.RES / 'cmb.json').get('HFb:1.1229e11:15.66')
    hb = run('hHFb_GM', '1.550:+0')['by_e'][CL]
    MX['HFb'] = dict(note='HFb final history (lattice B_P0GM + kinetic hand-over), node 1.3e11 only; H dependence taken as HFa',
                     K=run('hHFb_GM', '1.550:+0')['K'], SQ_1e6_node=hb['e1e-6']['S_per_Q'], e_crit_node=hb['e_crit_bulk'],
                     H_req_over_H_i_node=hb['e1e-6']['H_req_over_H_i'],
                     HFa_same_node=dict(SQ_1e6=run('hHFa_GM', '1.550:+0')['by_e'][CL]['e1e-6']['S_per_Q'], e_crit=run('hHFa_GM', '1.550:+0')['by_e'][CL]['e_crit_bulk']),
                     cmb=dict(Nstar=cmbb['Nstar'], r=cmbb['r_constH'], Nmax=cmbb['capacity']['Nmax'], MS=cmbb['MS']['ms']) if cmbb else None)
    # windows [C at the scanned points; E between them]
    MX['window'] = dict(c0_selfconsistent_free_X_with_e_crit=[1.45, 1.60], H_star_range_GeV=[fps['GM']['1.600:+0']['H_star_GeV'], fps['GM']['1.450:+0']['H_star_GeV']],
                        e_crit_range=[fps['GM']['1.500:+0']['e_crit_bulk'], fps['GM']['1.550:+0']['e_crit_bulk']],
                        delta_deg=[-2, 6], note='c0 1.45 is an extrapolated root (w = 1.34); phase edges from the K = 17.3 splice scan (delta -4, -6, +8 lose the e_crit) and delta = +4 on the final GM history; e_crit ~ 2e-6 needs c0 ~ 1.55 and moves to 9.6e-6 / 2.3e-6 if the lattice-stage m_s^2 is x1.03 / x0.97',
                        Q_matter_at_1e6_on_GM='none: S/Q(1e-6) floor 3.78 (c0 1.5, delta +4, 1.3e11 node), 3.49 at the 8.4e10 node (not self-consistent there, ratio 1.56); 3.49 with the fragment spread sigma = 0.3 (K 12.6)')
    MQ['window'] = dict(c0=[1.40, 1.50], delta_deg=[0, 4], note='c0 = 1.40 is marginal (3.46 interpolated; nodes disagree: 3.81 free X at 8.4e10, 3.29 Q at 1.3e11); 1.45 and 1.50 Q-matter at the 1.3e11 node; delta +4 Q-matter at c0 1.45 (1.3e11 node), free X at c0 1.40; only on the local-blocking closure',
                        exists_only_on='P0 (local Pauli blocking); not on GM')
    nat = dict(M2X=natural(MX['self_consistent']['H_i_GeV'], MX['isocurvature']['s_star_min_GeV'], 1.55, (1.45, 1.60)),
               M2Q=natural(MQ['self_consistent']['H_i_GeV'], MQ['isocurvature']['s_star_min_GeV'], 1.45, (1.40, 1.50)),
               verdict='de Sitter equilibrium gives s_rms = 0.195 H/m = 36-900 H (4e12-1e14 GeV) and needs 1e4-7e7 e-folds; the isocurvature floor needs s_* >= 1.1-1.4e16 GeV ~ 1e5 H. So s_* cannot come from equilibrium: it is an initial condition set before the observable e-folds (frozen: tadpole drift 10 H over 50 e-folds [Q mf_hybrid]). With g_h fixed by the floor, the window has probability ~exp(-1e4) under equilibrium; with g_h free the best case is ~5% but then beta_iso ~ 1.',
               refs=dict(SY='arXiv:astro-ph/9407016 (verified this session)'))
    K = {t: E.load(E.RES / ('handover_%s.json' % t))['summary'] for t in ('hoHFa_GM', 'hoHFa_P0', 'hoHFa_GM_g03', 'hoHFb_GM', 'hoH84_GM', 'hoH84_P0')}
    out = dict(label='e2e_chain M2 points; [C] given the model choices; nothing here is a prediction or validation',
               M2X=MX, M2Q=MQ, naturalness=nat,
               histories={t: dict(K=v['K'], T_RH_GeV=v['T_RH_GeV'], t_SM_half=v['t_SM_half']) for t, v in K.items()},
               fixpoints_cL0=fps)
    E.dump(E.RES / 'm2_points.json', out)
    for P in (MX, MQ):
        print(P['name'], json.dumps(dict(H=P['self_consistent']['H_i_GeV'], K=P['self_consistent']['K'], SQ=P['S_per_Q']['e1e6'], ecrit=P['S_per_Q']['e_crit_bulk'],
              relic=P['relic']['branch'], M=P['relic']['M_over_mX'], Q=P['relic']['Q_mean'], D=P['isocurvature']['D'], smin=P['isocurvature']['s_star_min_GeV'],
              fsdec=P['spectator']['f_s_at_decay_at_s_star_min'], cmb=P['cmb']), default=float))
    print(json.dumps(nat, default=float)[:2500])
