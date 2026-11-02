#!/usr/bin/env python3
"""dq_critic item 2 (Claude, own code, 2 Nov 2026): the Q-matter branch at a benchmark e, combined end to end.

Steps
 1. e(S/Q) with dq_heat's computed R(e) (sq_branch.Budget, frame-corrected, + item-3 entropy) for the S/Q points of
    qm_micro, the e -> infinity floor, and the benchmark e_b = 3e-6.
 2. Relic at e_b: qm_micro's microphysics with its skeptic's corrections (Bose + latent-heat-limited exchange for <Q>,
    energy-conserving f_v, M/m_X) interpolated to S/Q(e_b) [Q, C].  Two trapping cases: perfect interior trapping and
    Beliaev-only trapping (qm_micro variant; qm_transfer SKEPTIC: zeta_eff = 0.216 s/l for every object << l).
 3. P(k)/P_LCDM at k = 10, 20, 50, 100 /Mpc.  qm_transfer's direct runs at qm_micro's points sit on the code's
    numerical floor (|1 - P| <= 1.2e-3), so the physical deficit is obtained from qm_transfer's own grid
    (transfer_grid.json, family X_fc0.1_Tf0.01: f_v = 0.150 by mass, T_f = 0.010, thermo kinetics, s = 1) by
      - removing the code's floor pattern: d_rel(zeta,k) = 1 - P(zeta)/P(zeta = 1e-7) (the floor is the same to 5e-5 in
        every zeta = 1e-7 model of the grid), and correcting for the zeta = 1e-7 member's own deficit;
      - extrapolating in sigma_p (the probe control parameter of qm_transfer R5) to the benchmark sigma_p with the local
        exponent n(k) = dln d/dln sigma_p measured between zeta = 1e-6 and 1e-5 (upper estimate: the exponent rises toward
        the k^2 sigma^2 law as sigma falls) and with n = 2 from zeta = 1e-5 (lower estimate);
      - multiplying the upper estimate by the K_rsa = 800 / 200 deficit ratio of qm_transfer SKEPTIC's converged runs;
      - scaling linearly with f_v (first order in the warm fraction).  Labels: [C] on the grid, [E] for the extrapolation.
 4. Probe reach: qm_transfer's scaling law m_eff(sigma_p, f_v) (fit range 1.5-80 keV, extrapolated [E]); Viel et al. 2005
    (astro-ph/0501562) thermal-relic transfer at 5.7 and 20 keV for side-by-side deficits [R formula].
 5. Object properties: mass vs the Planck mass, radius, number flux at Earth (rho_DM = 0.3-0.4 GeV/cm^3, v = 230 km/s)
    [D on C].
 6. Quick check of the free-streaming length: qm_micro (0.34-0.61 kpc) vs qm_transfer ("~10 pc").
Output: ../results/qm_benchmark.json
"""
import json, math, os, sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
TR = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
import sq_branch as SB

MX = 0.0924089637                    # GeV
EQ_OBJ = 0.77587                     # charge-weighted <E/Q> of objects at <Q> >~ 1e19 (qm_micro) [C]
NB = 449.31336788454803              # bulk charge density (m_X^3)
HBARC_CM = 1.973269804e-14           # GeV cm
GEV_G = 1.78266192e-24               # g per GeV
MPL = 1.220890e19                    # GeV (Planck mass)
KQ = [10.0, 20.0, 50.0, 100.0]


def e_of_SQ():
    ht, ex, late, summ = SB.load()
    cf = summ['for_critic']['SQ_of_W']
    fits = late['long_run']['fits']
    Nf = abs([f for f in fits['NX'] if f['model'] == 'a-2' and f['a_min'] > 800][0]['N_inf'])
    b = SB.Budget(ht['portals_vanish_at_MY_benchmark'], Nf, 188737.6639)
    inc = 0.5 * sum(SB.INC)
    sq = lambda e: SB.SQ_closed(b.W(b.R(e)), cf) + inc
    out = {}
    for t in (2.5, 3.0, 3.2, 3.4716191256166837):
        out['%g' % t] = math.exp(brentq(lambda x: sq(math.exp(x)) - t, math.log(1e-7), math.log(1e-2)))
    out['SQ_floor_e_inf'] = SB.SQ_closed(b.W(b.Y), cf) + inc
    out['SQ_e3e-6'] = sq(3e-6); out['SQ_e5e-6'] = sq(5e-6)
    return out


def relic(SQb):
    """qm_micro + skeptic corrections, linear in S/Q between the 2.5 / 3.0 / 3.2 points [Q]."""
    S = [2.5, 3.0, 3.2]
    Qc = [(7.67e19, 7.96e19), (1.17e20, 1.27e20), (1.45e20, 1.50e20)]            # Bose + heat (s_phi 0..1)
    fv_p = [(0.0971, 0.0983), (0.1197, 0.1213), (0.1292, 0.1315)]                 # energy-conserving, perfect trapping
    fv_b = [0.141, 0.175, 0.189]                                                 # energy-conserving, Beliaev
    M_p = [0.7974, 0.8023, 0.8048]; M_b = [0.8076, 0.8152, 0.8182]
    Tf_p = [0.0087, 0.0089, 0.0091]; Tf_b = [0.0123, 0.0125, 0.0126]             # T_f / m_X (qm_micro criterion)
    Tvis_kin = (0.54, 0.99); Tvis_kin_b = (10.9, 12.9)                           # MeV, energy-conserving freeze-out
    sbar = [0.465, 0.47, 0.47]
    li = lambda arr: float(np.interp(SQb, S, arr))
    Q = (li([q[0] for q in Qc]), li([q[1] for q in Qc]))
    fvp = (li([f[0] for f in fv_p]), li([f[1] for f in fv_p]))
    r = dict(SQ=SQb, Q_charge_weighted=Q, f_v_charge_perfect=fvp, f_v_charge_beliaev=li(fv_b),
             M_over_mX=(li(M_p), li(M_b)), T_f=(li(Tf_p), li(Tf_b)), T_vis_kin_MeV=dict(perfect=Tvis_kin, beliaev=Tvis_kin_b),
             s_bar_Tf=li(sbar))
    r['f_v_mass_perfect'] = [f / r['M_over_mX'][0] for f in fvp]
    r['f_v_mass_beliaev'] = r['f_v_charge_beliaev'] / r['M_over_mX'][1]
    # zeta: s_eff <Q>^{-1/3}; capture-effective Q 1.2-1.4x smaller than charge-weighted (qm_micro SKEPTIC issue 10)
    r['zeta_perfect'] = [r['s_bar_Tf'] * (Q[1] / 1.2) ** (-1 / 3), r['s_bar_Tf'] * (Q[0] / 1.4) ** (-1 / 3)]
    r['zeta_beliaev'] = [3.6e-9, 5.6e-9]                                          # qm_transfer SKEPTIC, Q-independent
    return r


def sigma_bench():
    """sigma_p (qm_transfer convention, thermo) at the benchmark: interpolate qm_transfer's evaluated qm_micro points
    to S/Q = 3.0 and rescale zeta ~ Q^{-1/3} to the corrected <Q>; Beliaev from qm_transfer's Beliaev runs."""
    ev = {m['name']: m for m in json.load(open(os.path.join(TR, 'qm_transfer/results/qm_micro_eval.json')))['models']}
    a, b = ev['M_base_SQ2.5'], ev['M_base_SQ3.2']
    t = (3.0 - 2.5) / 0.7
    lQ = (1 - t) * math.log(1.2e20) + t * math.log(2.7e20)
    sig = math.exp((1 - t) * math.log(a['sigma_p']) + t * math.log(b['sigma_p']))
    lo = sig * (math.exp(lQ) / 1.27e20) ** (1 / 3)
    hi = sig * (math.exp(lQ) / 1.17e20) ** (1 / 3) * 1.12
    bel = [ev['M_bel_SQ3.2']['sigma_p'], ev['M_bel_SQ1.5']['sigma_p']]
    return dict(perfect=(lo, hi), beliaev=(min(bel), max(bel)), interp_Q=math.exp(lQ), interp_sigma=sig)


def grid_deficits(sig_b, fv_scale, fam='X_fc0.1_Tf0.01'):
    G = json.load(open(os.path.join(TR, 'qm_transfer/results/transfer_grid.json')))
    ks = np.array(G['k_Mpc'])
    M = {m['name']: m for m in G['models']}
    z = lambda s: M['%s_z%s' % (fam, s)]
    f7, f6, f5, f35 = z('1e-07'), z('1e-06'), z('1e-05'), z('3e-05')
    P7 = np.array(f7['P_z0'])
    d = {k: 1 - np.array(m['P_z0']) / P7 for k, m in (('6', f6), ('5', f5), ('35', f35))}
    s6, s5, s35 = f6['sigma_p'], f5['sigma_p'], f35['sigma_p']
    noise = 1e-5                     # floor reproducibility between zeta = 1e-7 models: <= 2e-5
    # K_rsa correction (qm_transfer SKEPTIC converged runs): deficit ratio K800/K200 for P_SQ3.4_z3e-04
    kr = json.load(open(os.path.join(TR, 'qm_transfer/skeptic/results/krsa_raw_runs.json')))['gridK_a.json']['P_SQ3.4_z3e-04_K800']
    p800 = np.array([kr['%s' % k] for k in ks]); p200 = np.array(M['P_SQ3.4_z3e-04']['P_z0'])
    kfac = np.maximum((1 - p800) / (1 - p200), 1.0)
    rows = []
    up, lo = [], []
    for i, k in enumerate(ks):
        good6 = d['6'][i] > 3 * noise   # d6 above 3e-5
        n65 = math.log(d['5'][i] / d['6'][i]) / math.log(s5 / s6) if good6 else None
        n53 = math.log(d['35'][i] / d['5'][i]) / math.log(s35 / s5)
        n_up = n65 if n65 is not None else n53
        # true deficit of the zeta = 1e-6 / 1e-5 members (add the 1e-7 member's own deficit)
        d6t = d['6'][i] / (1 - (f7['sigma_p'] / s6) ** n_up) if good6 else None
        d5t = d['5'][i] + (d6t * (f7['sigma_p'] / s6) ** n_up if good6 else d['5'][i] * (f7['sigma_p'] / s5) ** n_up)
        if good6:
            u = d6t * (sig_b[1] / s6) ** n_up
        else:
            u = d5t * (sig_b[1] / s5) ** n_up
        l = d5t * (sig_b[0] / s5) ** 2
        u *= kfac[i] * fv_scale; l *= fv_scale
        up.append(u); lo.append(l)
        rows.append(dict(k=float(k), d_rel_1e6=float(d['6'][i]), d_rel_1e5=float(d['5'][i]), n_65=n65, n_53=n53,
                         Krsa_factor=float(kfac[i]), deficit_upper=float(u), deficit_lower=float(l)))
    lk = np.log(ks)
    interp = lambda arr, kq: float(math.exp(np.interp(math.log(kq), lk, np.log(np.maximum(arr, 1e-300)))))
    at = {('%g' % kq): dict(P_over_PLCDM_min=1 - interp(up, kq), P_over_PLCDM_max=1 - interp(lo, kq),
                           deficit_range=[interp(lo, kq), interp(up, kq)]) for kq in KQ}
    return dict(family=fam, f_v_mass_family=f5['f_v_mass'], T_f_family=f5.get('T_f_eq'), sigma_bench=sig_b,
                fv_scale=fv_scale, per_k=rows, at_k=at)


def viel(mkeV, k, h=0.6736, om=0.1200 / 0.6736 ** 2):
    nu = 1.12
    al = 0.049 * mkeV ** (-1.11) * (om / 0.25) ** 0.11 * (h / 0.7) ** 1.22 / h      # Mpc
    return (1 + (al * k) ** (2 * nu)) ** (-10 / nu)                                  # P ratio = T^2


def m_eff(sig, fv):
    return dict(lya=19.3 * (sig / 1e-8) ** -0.49 * (fv / 0.1) ** -0.45,
                sat=19.5 * (sig / 1e-8) ** -0.44 * (fv / 0.1) ** -0.43,
                lens=23.2 * (sig / 1e-8) ** -0.36 * (fv / 0.1) ** -0.38)


def lam_fs(sig_1d, a_f, h=0.6736, Om=0.3153, Or_h2=4.18e-5):
    """comoving free-streaming length (kpc) of a NR gas with 1D comoving dispersion sig_1d today (units of c)."""
    Or = Or_h2 / h ** 2; OL = 1 - Om - Or
    cH0 = 2997.92458 / h                                                        # Mpc
    Hn = lambda a: math.sqrt(Or / a ** 4 + Om / a ** 3 + OL)
    v = math.sqrt(3) * sig_1d
    I = quad(lambda la: v / math.exp(2 * la) / Hn(math.exp(la)), math.log(a_f), 0.0, limit=200)[0]   # dl = v/(a^2 H) dln a
    return I * cH0 * 1e3


def main():
    out = dict(method=__doc__)
    es = e_of_SQ(); out['e_of_SQ'] = es
    eb = 3e-6; SQb = es['SQ_e3e-6']
    out['benchmark'] = dict(e=eb, SQ=SQb, why='S/Q in 2.5-3: direct heat gives 2.993 at 3e-6 (bulk Q-matter; Q = 2e5 objects not)')
    r = relic(SQb); out['relic'] = r
    sb = sigma_bench(); out['sigma_p'] = sb
    # P(k): perfect trapping (f_v,mass ~0.149) and Beliaev (0.215), scaled to the family's 0.150
    G = json.load(open(os.path.join(TR, 'qm_transfer/results/transfer_grid.json')))
    fvfam = [m for m in G['models'] if m['name'] == 'X_fc0.1_Tf0.01_z1e-05'][0]['f_v_mass']
    out['Pk_perfect'] = grid_deficits(sb['perfect'], max(r['f_v_mass_perfect']) / fvfam)
    out['Pk_beliaev'] = grid_deficits(sb['beliaev'], r['f_v_mass_beliaev'] / fvfam)
    # direct runs of qm_transfer at qm_micro's points (numerical floor)
    ev = {m['name']: m for m in json.load(open(os.path.join(TR, 'qm_transfer/results/qm_micro_eval.json')))['models']}
    out['direct_runs_floor'] = {n: max(abs(1 - p) for p in ev[n]['P_z0']) for n in ev}
    # references
    out['viel_reference'] = {('%g keV' % m): {('%g' % k): viel(m, k) for k in KQ} for m in (5.7, 6.5, 10.0, 20.0)}
    out['m_eff_extrapolated_keV'] = dict(perfect=m_eff(sb['perfect'][1], max(r['f_v_mass_perfect'])),
                                         beliaev=m_eff(sb['beliaev'][1], r['f_v_mass_beliaev']))
    # what a 20 keV-class probe needs (qm_transfer allowed_region, f_v = 10-20 %, T_f = 0.009-0.010)
    need = dict(zeta_20keV_fv10=2.3e-5, zeta_20keV_fv20=1.2e-5)
    need['factor_below'] = [need['zeta_20keV_fv20'] / max(r['zeta_perfect']), need['zeta_20keV_fv10'] / max(r['zeta_perfect'])]
    need['Q_needed_perfect_trapping'] = [(r['s_bar_Tf'] / need['zeta_20keV_fv10']) ** 3, (r['s_bar_Tf'] / need['zeta_20keV_fv20']) ** 3]
    need['note'] = 'Beliaev-only trapping: zeta_eff = 0.216 s/l ~ 4-6e-9 for every Q << 1e24 -> no object size reaches the band'
    out['probe_reach'] = need
    # objects
    obj = {}
    for lab, Q in (('Q_lo', min(r['Q_charge_weighted'])), ('Q_hi', max(r['Q_charge_weighted'])), ('Q_beliaev_hi', 1.6e21 * 0.7),
                   ('Q_SQ1.84_lo', 3.6e19)):
        Mg = Q * EQ_OBJ * MX
        Rm = (3 * Q / (4 * math.pi * NB)) ** (1 / 3) / MX * HBARC_CM                # cm
        flux = [rho / Mg * 2.3e7 * 1e4 * 3.156e7 for rho in (0.3, 0.4)]           # per m^2 per yr
        obj[lab] = dict(Q=Q, M_GeV=Mg, M_g=Mg * GEV_G, M_over_MPl=Mg / MPL, R_cm=Rm,
                        sigma_over_m_cm2_g=math.pi * Rm ** 2 / (Mg * GEV_G), flux_per_m2_yr=flux)
    out['objects'] = obj
    # free-streaming length check
    lf = {}
    for n in ('M_base_SQ2.5', 'M_base_SQ3.2', 'M_bel_SQ3.2'):
        m = ev[n]; lf[n] = dict(sigma_p=m['sigma_p'], a_f=m['a_f'], lam_fs_kpc=lam_fs(m['sigma_p'], m['a_f']))
    # qm_micro: v_rms(T_f) = 0.164 c (3D) at T_vis = 0.59 MeV -> sigma_1d today = v/sqrt3 * a_f
    af_m = 2.35e-13 / 0.59e-3 * (3.91 / 10.75) ** (1 / 3)                          # T0/T_vis with g*s(0.59 MeV) = 10.75
    lf['qm_micro_SQ3.0'] = dict(sigma_p=0.164 / math.sqrt(3) * af_m, a_f=af_m, lam_fs_kpc=lam_fs(0.164 / math.sqrt(3) * af_m, af_m))
    out['lambda_fs_check'] = lf
    json.dump(out, open(os.path.join(HERE, '..', 'results', 'qm_benchmark.json'), 'w'), indent=1, default=float)
    # print digest
    print('e(S/Q):', {k: '%.3e' % v if k[0].isdigit() else '%.4f' % v for k, v in es.items()})
    print('benchmark e=%.0e S/Q=%.4f' % (eb, SQb))
    print('relic:', json.dumps(r, default=float))
    print('sigma_p:', sb)
    for lab in ('Pk_perfect', 'Pk_beliaev'):
        print(lab, {k: ['%.2e' % x for x in v['deficit_range']] for k, v in out[lab]['at_k'].items()})
        for row in out[lab]['per_k']:
            print('   k=%6.1f d6=%.2e d5=%.2e n65=%s n53=%.2f Kf=%.3f  d=[%.2e, %.2e]' % (
                row['k'], row['d_rel_1e6'], row['d_rel_1e5'], '%.2f' % row['n_65'] if row['n_65'] else '  - ', row['n_53'],
                row['Krsa_factor'], row['deficit_lower'], row['deficit_upper']))
    print('direct-run floor:', {k: '%.1e' % v for k, v in out['direct_runs_floor'].items()})
    print('Viel refs:', {m: {k: '%.4f' % v for k, v in d.items()} for m, d in out['viel_reference'].items()})
    print('m_eff (extrapolated):', out['m_eff_extrapolated_keV'])
    print('probe reach:', need)
    for k, v in obj.items():
        print('object %s: M = %.2e GeV = %.2e g = %.2f M_Pl, R = %.1e cm, sigma/m = %.1e cm2/g, flux %.2f-%.2f /m2/yr' % (
            k, v['M_GeV'], v['M_g'], v['M_over_MPl'], v['R_cm'], v['sigma_over_m_cm2_g'], *v['flux_per_m2_yr']))
    print('lambda_fs:', lf)


if __name__ == '__main__':
    main()
