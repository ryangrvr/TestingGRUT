#!/usr/bin/env python3
"""Track sp_vector (Claude, own code): longitudinal-vector ledger at the 373d37e leading point.

Inputs (own results): background.json (vmodel.py), proca_modes.json (proca_modes.py),
phase_transfer.json (phase_transfer.py + vrel.c), floquet_vector.json (floquet_vector.py, floquet_band.py),
and the nominal vrel history (runs/hist_nominal.txt, theta_S(a)).

Pieces:
  A  inflationary (Bunch-Davies) Goldstone: frozen gauge angle alpha = H/(2 pi F_inf) -> late B_L = k alpha/e.
     Computed with the longitudinal Proca mode and evolving mass (profiles P1, P4); control P0 = GMR, vacuum mass.
  B1 release-induced, large scale: the separate-universe map from the source phase psi2 to the late Sigma phase,
     smoothed over small-scale noise, times the inflationary psi2 power at k ~ k_*; converted to vectors with the
     P0 transfer, which corresponds to a Goldstone-angle power (H_i/(2 pi F_vac))^2 = 2.974^2 rad^2 per ln k.
  B2 release-induced, Hubble scale: patch-to-patch Sigma-charge white noise set at a_s ~ 100, accumulated as Sigma
     phase until horizon entry (order-of-magnitude only).
Pairs per net charge: each vector decays to one X Xbar pair (only open channel). Normalized to the actual net
comoving charge |N_X| (both comoving numbers are fixed after release), and also to Y_Q (the convention of our
earlier gp_leftovers/gp_critic estimate, which corresponds to removing the 66-fold yield without diluting vectors).
Usage: python3 vector_ledger.py ../results ../results/vector_ledger.json <hist_nominal.txt>
"""
import json
import sys
from pathlib import Path

import numpy as np

ALLOW = 0.0049296413205678555
NX_NOM = 238.64868
YQ_FACTOR = 65.98227759048056   # raw |Y_X|/Y_Q at the point


def main(rdir, out, hist):
    r = Path(rdir)
    bg = json.loads((r / 'background.json').read_text())
    pm = json.loads((r / 'proca_modes.json').read_text())
    pt = json.loads((r / 'phase_transfer.json').read_text())
    fq = json.loads((r / 'floquet_vector.json').read_text())
    P = pm['profiles']
    res = dict(inputs=dict(F_inf_GeV=bg['F_g_GeV'], m_A_inf_GeV=bg['m_A_inf_GeV'], m_A_vac_GeV=bg['m_A_vac_GeV'],
                           mass_ratio_sq=bg['mass_ratio_sq_vac_over_inf'], alpha_BD=bg['alpha_BD_rad'],
                           a_star=bg['a_star'], k_star=bg['k_star_over_Hi'], allowance=ALLOW))
    # ---- decay of the vector
    e = 1e-7; mA = bg['m_A_vac_GeV']; mX = 0.092409
    beta = np.sqrt(1 - 4 * mX**2 / mA**2)
    width = e**2 * mA * beta**3 / (48 * np.pi)
    g = 106.75; Mred = 2.435e18
    T_HG = (90 / (np.pi**2 * g))**0.25 * np.sqrt(width * Mred)
    T_mean = T_HG / (np.sqrt(np.pi) / 2 * np.sqrt(2))    # <a> weighting of an exponential decay in radiation
    res['decay'] = dict(width_GeV=width, lifetime_s=6.582e-25 / width, T_vis_at_H_eq_Gamma_GeV=T_HG,
                        T_vis_decay_weighted_GeV=T_mean,
                        channels='A -> X Xbar only (Y, Sigma radial, Goldstone closed/eaten; A -> phi phi forbidden by '
                                 'Bose symmetry; kinetic mixing forbidden by dark C in the 373d37e EFT audit [Q])',
                        pairs_per_vector=1, energy_per_vector_in_E_X_at_decay=15.2,
                        note='energy multiple from sp_yledger (same mass and width) [Q]')
    # ---- piece A
    A = {}
    for p in ('P0', 'P1', 'P4', 'P2', 'P5'):
        q = P[p]
        A[p] = dict(N_A=q['N_A'], N_A_superhorizon_k_lt_1_3000=q['N_A_k_below_1_over_3000'],
                    ratio_to_P0=q['ratio_to_P0'],
                    pairs_actual=q['N_A'] / NX_NOM, pairs_actual_superhorizon=q['N_A_k_below_1_over_3000'] / NX_NOM,
                    pairs_YQnorm=q['N_A'] / NX_NOM * YQ_FACTOR,
                    Omega_over_DM_instant_reheat=q['physical']['Omega_A_over_DM'])
    A['P0_vs_GMR_formula'] = pm['P0_over_GMR_formula']
    A['ratio_test_low_k'] = 'n_k(P1)/n_k(P0) = 2.357e-6 for k <= 1e-4 vs (m_vac/m_inf)^2 = 2.3607e-6'
    res['piece_A'] = A
    # ---- transfer from P0 (dN/dlnk per unit Goldstone-angle power)
    ks = np.array([m['k'] for m in P['P0']['modes']])
    dN0 = ks**3 * np.array([m['n_k'] for m in P['P0']['modes']]) / (2 * np.pi**2)
    PG = (1 / (2 * np.pi * 6 * 1e11 / bg['H_i_GeV']))**2   # = 2.974^2
    lk = np.log(ks)
    sel = ks <= 1 / 3000.0
    nu = bg['light_nu']; amp_exit = bg['dpsi2_per_lnk_at_exit']
    Ppsi = amp_exit**2 * ks**(2 * (1.5 - nu))           # psi2 power per ln k at the end of inflation
    # ---- piece B1
    lin = pt['linear_phase']
    dthS_lin = float(np.median([x['dthS_inf'] for x in lin if x['dd'] <= 1e-7]))
    sm = pt['smoothed_response']
    scan = json.loads((r / 'phase_scan_smoothed.json').read_text()) if (r / 'phase_scan_smoothed.json').exists() else None
    B1 = dict(linear_dthetaS_dpsi=dthS_lin, linear_dthetaG_dpsi=dthS_lin / 6,
              linear_range_rad=1e-6, smoothed_crn=sm, smoothed_scan=scan)
    nat = pt['ensembles']['0.489']
    NX_nat = abs(nat['NX']['mean'])
    slopes = []
    for s_ in sm:
        slopes.append((f"crn dpsi={s_['dpsi']}", s_['dthS_dpsi'], s_['dthS_dpsi_sem']))
    if scan:
        slopes.append(('scan +-0.6', scan['slope_thS'], scan['slope_thS_se']))
    rows = []
    for lab, sl, se in slopes:
        for which, val in (('central', abs(sl)), ('2sigma_upper', abs(sl) + 2 * se)):
            slopeG = val / 6
            NA = float(np.trapezoid(dN0[sel] * slopeG**2 * Ppsi[sel] / PG, lk[sel]))
            rows.append(dict(source=lab, which=which, dthetaS_dpsi=val, dthetaG_at_kstar=slopeG * bg['dpsi2_per_lnk_kstar_end'],
                             N_A=NA, pairs_per_net_natural=NA / NX_nat, ratio_to_allowance_natural=NA / NX_nat / ALLOW))
    B1['estimates_natural_spread'] = rows
    B1['NX_natural_net'] = NX_nat
    B1['criterion_dthetaG_kstar_max_rad'] = float(np.sqrt(ALLOW / (P['P0']['N_A'] / NX_NOM)) * np.sqrt(PG))
    B1['criterion_dpsi_kstar_max_if_linear'] = B1['criterion_dthetaG_kstar_max_rad'] / (abs(dthS_lin) / 6)
    linr = pt['linear_radial_Sigma']
    dthS_rad = float(np.median([x['dthS_inf'] for x in linr if x['ds'] <= 1e-7]))
    B1['linear_dthetaS_dlnSigma'] = dthS_rad
    B1['criterion_dlnSigma_kstar_max_if_linear'] = B1['criterion_dthetaG_kstar_max_rad'] / (abs(dthS_rad) / 6)
    rt = json.loads((r / 'radial_transfer.json').read_text()) if (r / 'radial_transfer.json').exists() else None
    if rt:
        B1['radial_modes_linear'] = [dict(m2_over_H2=md['m2_over_H2'], u_XYS=md['u_canonical_XYS'],
                                          dthetaS_dc=md['response']['1e-07']['dthS_dc'],
                                          dc_per_lnk_kstar=md['dc_per_lnk_kstar'], dthetaG_kstar_linear=md['dthetaG_kstar'],
                                          ratio_to_criterion_linear=md['ratio_to_criterion']) for md in rt['modes']]
        B1['radial_noise_ensemble'] = rt['radial_noise_ensemble']
    res['piece_B1'] = B1
    # ---- piece B2 (order of magnitude)
    a_h, th_h = [], []
    for line in open(hist):
        if line.startswith('H '):
            f = line.split(); a_h.append(float(f[1])); th_h.append(float(f[2]))
    a_h = np.array(a_h); th_h = np.array(th_h)
    hhat = 1e11 / bg['H_i_GeV']
    NS_nom = 38.75
    B2 = {}
    cases = [('natural_0.489', pt['ensembles']['0.489']['NS']['std']**2, abs(pt['ensembles']['0.489']['NX']['mean'])),
             ('codex_0.03', pt['ensembles']['0.03']['NS']['std']**2, abs(pt['ensembles']['0.03']['NX']['mean'])),
             ('sigma_1e-3', pt['ensembles']['0.001']['NS']['std']**2, abs(pt['ensembles']['0.001']['NX']['mean']))]
    if rt:
        cases.append(('radial_noise_only_phase_fixed', rt['radial_noise_ensemble']['NS_std']**2,
                      abs(rt['radial_noise_ensemble']['NX_mean'])))
    for case, varN, NXc in cases:
        out_c = {}
        for a_s in (50.0, 100.0):
            th_s = np.interp(np.log(a_s), np.log(a_h), th_h)
            ak = np.geomspace(a_s, min(a_h[-1], 1000.0), 200)
            dth = (np.interp(np.log(ak), np.log(a_h), th_h) - th_s) / NS_nom
            for cw in (0.1, 0.8):
                integ = hhat**2 / 2 * ak * cw * varN * dth**2 * (a_s / ak)**3
                NA = float(np.trapezoid(integ, np.log(ak)))
                out_c[f'a_s={int(a_s)},c_w={cw}'] = dict(N_A=NA, pairs_per_net=NA / NXc)
        vals = [v['pairs_per_net'] for v in out_c.values()]
        B2[case] = dict(var_NS=varN, net_NX=NXc, cases=out_c, range_pairs_per_net=[min(vals), max(vals)],
                        range_ratio_to_allowance=[min(vals) / ALLOW, max(vals) / ALLOW])
    res['piece_B2'] = B2
    # ---- ensembles summary
    res['ensembles_summary'] = {k: dict(eps_mean=v['eps']['mean'], eps_std=v['eps']['std'], NX_mean=v['NX']['mean'],
                                        NS_std=v['NS']['std'], thS_inf_std=v['thS_inf']['std'],
                                        minS_min=v['minS_min']) for k, v in pt['ensembles'].items()}
    # ---- remedy: larger e
    res['remedy_e'] = dict(e_heavy_inflation=bg['e_heavy_during_inflation'],
                           m_A_vac_at_e_heavy_GeV=6 * bg['e_heavy_during_inflation'] * 1e11,
                           floquet=fq['scans'], band_search=fq.get('band_search'),
                           sigma_to_AA_closes_e=7.07e10 / 6e11, YH_to_YL_A_closes_e=7.08e10 / 6e11)
    Path(out).write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(dict(A={k: (v if not isinstance(v, dict) else dict(pairs_actual=v['pairs_actual'],
                                                                          pairs_superh=v['pairs_actual_superhorizon'],
                                                                          pairs_YQ=v['pairs_YQnorm']))
                             for k, v in A.items()},
                          B1=dict(lin=dthS_lin, crit=B1['criterion_dthetaG_kstar_max_rad'],
                                  crit_dpsi=B1['criterion_dpsi_kstar_max_if_linear'],
                                  crit_dlnS=B1['criterion_dlnSigma_kstar_max_if_linear'], rows=rows),
                          B2={k: v['range_pairs_per_net'] for k, v in B2.items()}, decay=res['decay']), indent=1,
                     default=float))


if __name__ == '__main__':
    main(*sys.argv[1:])
