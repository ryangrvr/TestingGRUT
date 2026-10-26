#!/usr/bin/env python3
"""hw_envelope (Claude, own code): Codex 2908e91's envelope bound evaluated at the ls_scaling heavy ladder and spot point.

Closed form (own evaluation): alpha6 = r kappa^2/(2 M_Y^2 (1-r^2)) = mu kappa^2 Sigma/(M^4 - 4 mu^2 Sigma^2);
R = r_i (1+2/cH)^(1/4) with (a) the leading radial relation r_i = sqrt(cH/(2 lam)) H_i, and (b) the actual prepared
|X_i| of the Y-relaxed minimum; bound |a^3 n_X| <= 12 |alpha6| R^6 / H_i; yield over the prompt-radiation entropy.
The simulated a^3 n_X, E_X, |X_i| and yields are read from our ls_scaling track's results/scan_B.json [Q, own track,
not Codex].  Also computed: Codex's matching diagnostics at each ladder point, and the factor decomposition of the
bound between Codex's fixed-coupling row (r_i = M_Y, cH = 10) and the spot point.
Usage: python3 -B spot_bound.py <scan_B.json> <out.json>
"""
import json, math, sys
from pathlib import Path
import numpy as np

MRED, GSTAR, YQ, LAM, CH = 2.435e18, 106.75, 4.8294009076988175e-9, 1e-3, 10.0


def entropy(H):
    T = (90 * MRED**2 * H**2 / (math.pi**2 * GSTAR))**0.25
    return 2 * math.pi**2 / 45 * GSTAR * T**3


def alpha6_H(m, r, kap):
    """alpha6 in H_i^-2 units for M_Y = m H_i."""
    return r * kap**2 / (2 * m**2 * (1 - r**2))


def evaluate(row):
    cpl = row['couplings_over_H']; m, hH, muH = cpl['M'], cpl['h'], cpl['mu']
    kap = row['kappa']; H = row['H_GeV']
    r = math.sqrt(2) * muH * hH / m**2
    a6 = alpha6_H(m, r, kap)
    Sig = hH / math.sqrt(2); D = m**4 - 4 * muH**2 * Sig**2
    a6_direct = muH * kap**2 * Sig / D
    x_lead = math.sqrt(CH / (2 * LAM))
    x_act = row['prep_info']['x']
    out = dict(name=row['name'], delta=row['delta'], prep=row['prep'], H_GeV=H, m=m, h_over_H=hH, mu_over_H=muH, r=r,
               kappa=kap, alpha6_Hunits=a6, alpha6_direct_rel_dev=abs(a6_direct / a6 - 1),
               alpha6_GeVm2=a6 / H**2, c6=a6 / H**2 * MRED**2, eta=a6 * x_lead**2 / LAM, eta_edge=r / 12,
               x_lead=x_lead, x_actual=x_act, x_actual_over_lead=x_act / x_lead)
    NX, EX = abs(row['NX']), row['EX']
    for tag, x0 in (('lead', x_lead), ('actual', x_act)):
        # envelope from the conformal energy at a = 1 with zero physical velocity: R^4 = x0^4 + x0^2/lam
        R = (x0**4 + x0**2 / LAM)**0.25
        N = 12 * a6 * R**6
        out[f'R_{tag}'] = R
        out[f'bound_a3nX_{tag}'] = N
        out[f'margin_{tag}'] = N / NX
        out[f'bound_yield_over_YQ_{tag}'] = N * H**3 / entropy(H) / YQ
        Ncirc = 2 * math.sqrt(2 * LAM) * (LAM * R**4 / (3 * LAM))**0.75
        out[f'eps_bound_{tag}'] = N / Ncirc
    out['R_lead_check_codex_formula'] = x_lead * (1 + 2 / CH)**0.25
    out['late_amplitude_envelope_from_EX'] = (EX / LAM)**0.25      # a|X| envelope implied by the late E_X
    out['envelope_violation_vs_lead'] = (EX / LAM)**0.25 / out['R_lead']
    out['envelope_violation_vs_actual'] = (EX / LAM)**0.25 / out['R_actual']
    out['sim_a3nX'] = NX; out['sim_EX'] = EX; out['sim_eps'] = row['eps']
    out['sim_yield_over_YQ'] = row['yield_over_YQ']
    out['own_yield_over_YQ'] = NX * H**3 / entropy(H) / YQ
    out['saturation_fraction_lead'] = NX / out['bound_a3nX_lead']
    out['eps_over_eta'] = abs(row['eps']) / out['eta']
    # Codex matching diagnostics at R_lead
    R = out['R_lead']; gap2 = m**2 - 2 * muH * Sig; msig = math.sqrt(2) * hH
    Ymax = kap * R**3 / gap2
    out['diag'] = dict(kappaR2_over_M2=kap * R**2 / m**2, omegaX_over_gap=math.sqrt(2 * LAM) * R / math.sqrt(gap2),
                       omegaX_over_msigma=math.sqrt(2 * LAM) * R / msig, Y_quartic_force_ratio=2 * row.get('ly', 1.0) * Ymax**2 / gap2,
                       Sigma_forced_shift=muH * Ymax**2 / (msig**2 * Sig),
                       lamY_Yprep2_over_M2=((row['prep_info']['Y'][0]**2 + row['prep_info']['Y'][1]**2) if 'Y' in row['prep_info']
                                            else row['prep_info']['y']**2) / m**2,
                       mL_over_omegaX=math.sqrt(gap2) / (1.198 * math.sqrt(LAM) * x_lead * math.sqrt(1 + 1 / (LAM * x_lead**2))))
    return out


def fixed_coupling_row():
    """Codex's r_i = M_Y, cH = 10 row in the same H-unit language."""
    M, hs, mu, kap = 1e12, 1e11, 5e11, 1e-3
    r = math.sqrt(2) * mu * hs / M**2
    x_lead = math.sqrt(CH / (2 * LAM))
    H = M / x_lead                            # r_i = M_Y  <=>  M_Y/H_i = sqrt(cH/2lam)
    m = M / H
    a6 = alpha6_H(m, r, kap)
    R = x_lead * (1 + 2 / CH)**0.25
    N = 12 * a6 * R**6
    Ncirc = 2 * math.sqrt(2 * LAM) * (R**4 / 3)**0.75
    return dict(H_GeV=H, m=m, r=r, kappa=kap, alpha6_Hunits=a6, alpha6_GeVm2=a6 / H**2, R=R, bound_a3nX=N,
                eps_bound=N / Ncirc, eta=a6 * x_lead**2 / LAM, eta_edge=r / 12,
                bound_yield_over_YQ=N * H**3 / entropy(H) / YQ)


def fixed_coupling_vs_m():
    """Codex's bound at the FIXED published couplings as a function of m = M_Y/H_i (cH = 10), with the matching
    diagnostics: the quiet regime at fixed couplings versus the heavy window m = 13-43 of ls_scaling."""
    M, hs, mu, kap, ly = 1e12, 1e11, 5e11, 1e-3, 1.0
    r = math.sqrt(2) * mu * hs / M**2; Sig = hs / math.sqrt(2)
    x_lead = math.sqrt(CH / (2 * LAM)); R = x_lead * (1 + 2 / CH)**0.25
    Ncirc = 2 * math.sqrt(2 * LAM) * (R**4 / 3)**0.75
    rows = []
    for m in (70.71, 43.0, 30.0, 13.0, 10.0, 5.957, 2.2995):
        H = M / m; muH, hH = mu / H, hs / H; SigH = Sig / H
        a6 = alpha6_H(m, r, kap); N = 12 * a6 * R**6
        gap2 = m**2 - 2 * muH * SigH; msig = math.sqrt(2) * hH; Ymax = kap * R**3 / gap2
        rows.append(dict(m=m, H_GeV=H, ri_over_M=x_lead / m, alpha6_Hunits=a6, eta=a6 * x_lead**2 / LAM,
                         eps_bound=N / Ncirc, eps_estimate_7p8eta=7.8 * a6 * x_lead**2 / LAM,
                         bound_a3nX=N, bound_yield_over_YQ=N * H**3 / entropy(H) / YQ,
                         estimate_yield_over_YQ_0p38=0.38 * N * H**3 / entropy(H) / YQ,
                         diag=dict(kappaR2_over_M2=kap * R**2 / m**2, omegaX_over_gap=math.sqrt(2 * LAM) * R / math.sqrt(gap2),
                                   omegaX_over_msigma=math.sqrt(2 * LAM) * R / msig, Y_quartic_force_ratio=2 * ly * Ymax**2 / gap2,
                                   Sigma_forced_shift=muH * Ymax**2 / (msig**2 * SigH))))
    return rows


if __name__ == '__main__':
    rows = json.load(open(sys.argv[1]))['rows']
    want = ['B_m10', 'B_m15', 'B_m20', 'B_m30', 'B_m40', 'B_base_Hscaled', 'B_base_codexprep', 'B_r0.0707', 'B_r0.2',
            'B_r0.7', 'B_frac0.25', 'B_frac0.75', 'B_frac0.9']
    res = [evaluate(r) for r in rows if 'eps' in r and r['name'] in want and abs(r['delta'] - math.pi / 2) < 1e-6]
    fx = fixed_coupling_row()
    spot = next(r for r in res if r['name'] == 'B_base_Hscaled')
    # decomposition of bound ratio spot/fixed: alpha6_H ratio (r-factor, kappa^2, m^-2) x R^6 ratio (=1) x H^1.5 ratio
    rf = (spot['r'] / (1 - spot['r']**2)) / (fx['r'] / (1 - fx['r']**2))
    kf = (spot['kappa'] / fx['kappa'])**2
    mf = (fx['m'] / spot['m'])**2
    Hf = (spot['H_GeV'] / fx['H_GeV'])**1.5
    decomp = dict(alpha6_H_ratio=spot['alpha6_Hunits'] / fx['alpha6_Hunits'], r_factor=rf, kappa2_factor=kf, m_factor=mf,
                  product_rkm=rf * kf * mf, R6_ratio=(spot['R_lead'] / fx['R'])**6, H_pow_1p5_ratio=Hf,
                  predicted_bound_ratio=rf * kf * mf * Hf,
                  actual_bound_ratio=spot['bound_yield_over_YQ_lead'] / fx['bound_yield_over_YQ'],
                  GeV_view=dict(alpha6_GeV_ratio=spot['alpha6_GeVm2'] / fx['alpha6_GeVm2'],
                                ri_over_M_spot=spot['x_lead'] / spot['m'], ri_over_M_fixed=1.0,
                                ri_ratio_pow_3p5=(spot['x_lead'] / spot['m'])**3.5,
                                scale_ratio_MY=4.0, scale_pow_3p5=4.0**3.5))
    fvm = fixed_coupling_vs_m()
    out = dict(method=__doc__, ladder=res, codex_fixed_row_ri_eq_MY_cH10=fx, decomposition_spot_over_fixed=decomp,
               fixed_couplings_vs_m=fvm)
    Path(sys.argv[2]).write_text(json.dumps(out, indent=1, default=float))
    for f in fvm:
        print('fixed m=%6.2f H=%.3e ri/M=%.3f a6_H=%.3e eta=%.3e eps_bound=%.4f eps~7.8eta=%.4f yield_bound=%.3e (x0.38: %.3e) diag %s'
              % (f['m'], f['H_GeV'], f['ri_over_M'], f['alpha6_Hunits'], f['eta'], f['eps_bound'], f['eps_estimate_7p8eta'],
                 f['bound_yield_over_YQ'], f['estimate_yield_over_YQ_0p38'], {k: '%.3g' % v for k, v in f['diag'].items()}))
    for r in res:
        print(r['name'], 'own yield from scan NX = %.4f vs scan %.4f' % (r['own_yield_over_YQ'], r['sim_yield_over_YQ']))
    for r in res:
        print('%-18s m=%4.0f r=%.3f kap=%.4g a6_H=%.3e eta=%.4f x0=%.2f (lead %.2f) R_lead=%.2f R_act=%.2f '
              'bound=%.0f/%.0f sim=%.0f margin=%.2f/%.2f yield_bound=%.3f/%.3f sim=%.3f eps=%.4f eps/eta=%.1f env_late/R=%.3f/%.3f'
              % (r['name'], r['m'], r['r'], r['kappa'], r['alpha6_Hunits'], r['eta'], r['x_actual'], r['x_lead'],
                 r['R_lead'], r['R_actual'], r['bound_a3nX_lead'], r['bound_a3nX_actual'], r['sim_a3nX'],
                 r['margin_lead'], r['margin_actual'], r['bound_yield_over_YQ_lead'], r['bound_yield_over_YQ_actual'],
                 r['sim_yield_over_YQ'], r['sim_eps'], r['eps_over_eta'], r['envelope_violation_vs_lead'],
                 r['envelope_violation_vs_actual']))
        print('    diag', {k: '%.3g' % v for k, v in r['diag'].items()})
    print('fixed', json.dumps(fx, indent=1))
    print('decomp', json.dumps(decomp, indent=1))
