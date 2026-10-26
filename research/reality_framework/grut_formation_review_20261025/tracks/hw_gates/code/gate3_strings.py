#!/usr/bin/env python3
"""hw_gates gate 3 (Claude, own code): Z6 strings (U(1) -> Z6 by Sigma, charge 6) at h = 5.33e12 GeV.

Tension: mu = pi h^2 eps(beta), beta = m_s^2/m_v^2 = 2 lambda_Sigma/(36 e^2) (independent of h). eps(beta) is re-derived
with our gp_stability Nielsen-Olesen solver in s = ln(m_v r) (own code, reused) at beta = 1e4, 1e6, 1e8 and the
large-beta law eps = (1/2) ln beta + c is refitted. G mu = mu/M_Pl^2.
Formation: (i) hidden temperature if all X energy thermalised at release, versus h; (ii) visible reheating temperature
and the absence of a Sigma-visible coupling; (iii) |Sigma| during inflation (aligned minimum) and its radial mass;
(iv) min |Sigma|/vev along the release (own hwgrel run: phase slips?); (v) Sigma winding during the release.
Limits are quoted inputs (read levels in the JSON), retrospective only.
Usage: python3 -B gate3_strings.py <out.json>
"""
import json, sys, time
from pathlib import Path
import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / 'gp_stability' / 'code'))
import spotmodel as S
import z6_strings as Z     # our gp_stability solver (own code, reused)

P = S.params(); H = S.H_SPOT
h_GeV = P['h'] * H


def main(out):
    t0 = time.time()
    rows = []
    for beta in (1e2, 1e4, 1e6, 1e8):
        eps, parts, nodes, _ = Z.solve_string(beta)
        rows.append(dict(beta=beta, eps=eps, eps_minus_half_ln_beta=eps - 0.5 * np.log(beta), nodes=nodes))
        print(rows[-1], time.time() - t0, flush=True)
    c = rows[-1]['eps_minus_half_ln_beta']; c_unc = abs(rows[-1]['eps_minus_half_ln_beta'] - rows[-2]['eps_minus_half_ln_beta'])
    res = dict(method=__doc__, NO_scan=rows, asymptotic_constant=c, constant_uncertainty=c_unc,
               gp_stability_constant_Q=0.15023983484400993)
    def Gmu(e, hG, lamS=1.0):
        beta = 2 * lamS / (36 * e * e); eps = 0.5 * np.log(beta) + c
        return dict(e=e, h_GeV=hG, beta=beta, eps=eps, mu_GeV2=np.pi * hG**2 * eps, Gmu=np.pi * hG**2 * eps / S.MPL**2,
                    m_v_GeV=6 * e * hG, m_s_GeV=np.sqrt(2 * lamS) * hG)
    res['spot'] = Gmu(P['e'], h_GeV)
    res['stated_point_control'] = Gmu(1e-7, 1e11)
    res['e_scan_at_spot_h'] = [Gmu(e, h_GeV) for e in (1e-7, 1e-5, 1e-3, 5e-3, 1e-2)]
    # formation conditions
    EX = 38398.22      # own release window mean (comoving a^4 rho_X at late times); at a = 1 rho_X = V_inf ... use the initial energy
    xa, ya, sa = S.aligned(P)
    rho_rel = S.V(np.array([xa, -ya, -sa]), P, cH=0.0)   # physical energy density right after the Hubble mass is removed (H^4 units)
    gh = 10.0
    T_h_max = (30 * rho_rel / (np.pi**2 * gh))**0.25 * H
    s_RH, T_RH = S.entropy_RH(H)
    za = np.array([xa, -ya, -sa])
    sp = S.spectrum(za, P)
    res['formation'] = dict(
        rho_release_H4=rho_rel, T_hidden_max_if_instantly_thermalised_GeV=T_h_max, T_hidden_max_over_h=T_h_max / h_GeV,
        T_RH_visible_GeV=T_RH, T_RH_over_h=T_RH / h_GeV,
        sigma_visible_coupling='none at tree level (Sigma couples to Y via mu and to the hidden vector; the Higgs portal is on phi, kappa_portal <~ 1e-9) [Q 7737487/gp_walls]',
        abs_Sigma_inflation_over_vev=sa / (P['h'] / np.sqrt(2)), Sigma_radial_mass_over_H=float(np.sqrt(2 * P['ls']) * P['h']),
        Sigma_radial_mode_m2_hessian=sorted(m['m2'] for m in sp['modes'] if not m['gauge'])[-1],
        note='U(1) broken by Sigma (and X) throughout inflation; no thermal restoration: T_h,max < h and the visible bath does not couple to Sigma')
    # release: min |Sigma|/vev and Sigma winding (relaxed preparation, delta = pi/2)
    z0, info = S.prep('relaxed', np.pi / 2, P)
    r = S.run(z0, P, a_end=60.0, rtol=1e-10, atol=1e-13, nsamp=30, win0=30.0)
    sm = r['samples']
    res['release'] = dict(min_absSigma_over_vev=r['end']['minS_over_vev'], thetaS_release_rad=sm[0]['thS'], thetaS_a60_rad=sm[-1]['thS'],
                          thetaS_tail_rad=-sm[-1]['NS'] / (P['h']**2 * 60.0), total_winding_rad=sm[-1]['thS'] - sm[0]['thS'] - sm[-1]['NS'] / (P['h']**2 * 60.0),
                          NS_late=sm[-1]['NS'], windings=abs(sm[-1]['thS'] - sm[0]['thS']) / (2 * np.pi))
    res['limits_quoted'] = dict(
        LVK_O3_2101_12248=dict(model_A_BOS_excluded_above=[9.6e-9, 1e-6], model_B_LRS_excluded_above=[4.0e-15, 6.3e-15],
                                model_C1_excluded_above=[2.1e-15, 4.5e-15], model_C2_excluded_above=[4.2e-15, 7.0e-15],
                                read_level='[Q] as recorded by our gp_stability track from the arXiv abstract/body; re-read of the abstract attempted this session (see NOTES)'),
        PTA=dict(NANOGrav_12p5yr_NG_strings_Gmu_upper='~1e-10 (Blanco-Pillado, Olum, Wachter 2021) [R]',
                 NANOGrav_15yr_new_physics='stable strings fit the signal for G mu ~ 1e-11..1e-10 with specific loop models; not an upper limit [R, see NOTES]',
                 PPTA_2016='G mu < 1.5e-11 for one loop model [R]'),
        Planck_2013_XXV=dict(NG=1.5e-7, AH=3.2e-7, read_level='[Q via gp_stability]'))
    g = res['spot']['Gmu']
    res['comparison'] = dict(Gmu=g, over_LVK_modelB_low=g / 4.0e-15, over_LVK_modelC1_low=g / 2.1e-15, over_LVK_modelA_low=g / 9.6e-9,
                             over_PTA_1em10=g / 1e-10, over_PPTA_1p5em11=g / 1.5e-11, over_Planck_NG=g / 1.5e-7)
    Path(out).write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: v for k, v in res.items() if k not in ('method', 'NO_scan', 'e_scan_at_spot_h')}, indent=1, default=float))
    print('done', time.time() - t0)


if __name__ == '__main__':
    main(sys.argv[1])
