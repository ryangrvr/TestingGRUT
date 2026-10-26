#!/usr/bin/env python3
"""ls_scaling (Claude, own code): closed-form scalings [D] checked against our own release scans [C].

Reads results/scan_B.json (and scan_A.json if present); writes results/analytic.json.
Usage: python3 -B analytic.py <results dir>
"""
import json, sys
from pathlib import Path
import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import scal as S  # noqa: E402

R = Path(sys.argv[1])
B = {(r['name'], round(r['delta'], 3)): r for r in json.load(open(R / 'scan_B.json'))['rows'] if 'eps' in r}
out = {}
H, CH, LAM, MRED = S.H_CAND, 10.0, 1e-3, S.MRED
cp0 = dict(S.STATED)

# 1. closed-form ledger at candidate 10 (Codex DECAY_RESULTS inputs [Q]) -> R_Y, R_Sigma, R_cascade
q = dict(NYH=16.780301052357483, NYL=22.184313542703066, NSIG=24.976728038939, EX=33077.302400370805)
L = S.ledger(cp0, H, q['NYH'], q['NYL'], q['NSIG'], q['EX'])
out['closed_form_ledger_cand10'] = dict(ours=dict(R_Y=L['R_Y'], R_Sigma=L['R_Sigma'], R_cascade=L['R_cascade'],
                                                  closed_R_Y=L['closed_form_R_Y'], closed_R_Sigma=L['closed_form_R_Sigma']),
                                        codex=dict(R_Y=113.0689908402781, R_Sigma=86.335210249924,
                                                   R_cascade=29.00161336492347),
                                        deposit_per_YL_Hunits=L['deposit_per_YL'], deposit_per_A_Hunits=L['deposit_per_A'])

# 2. dimensionless ratios at the stated point (H_i of the candidates) [C]
m = S.minima(cp0, H, CH)[0]
x0 = m['x']; omegaX = 1.198 * np.sqrt(LAM) * x0 * np.sqrt(1 + 1 / (LAM * x0**2))  # quartic-oscillator frequency incl.
p = S.pars(cp0, H); ms = S.masses(cp0)                                            # the initial conformal kinetic term
out['stated_ratios'] = dict(x0=x0, omegaX_over_H=omegaX, M_over_H=p['M'], mL_over_H=ms['mYL'] / H,
                            msig_over_H=ms['msig'] / H, mL_over_omegaX=ms['mYL'] / H / omegaX,
                            msig_over_omegaX=ms['msig'] / H / omegaX, lamY_y2_over_M2=m['lamY_y2_over_M2'],
                            kapx2_over_M2=m['kapx2_over_M2'], r=ms['r'], kappa_max=S.kappa_max(cp0),
                            c6=S.c6(cp0)['c6'])

# 3. heavy-limit kink estimate N_Y ~ 9 cH^2 kap^2 x0^6 / m_L^7 versus the m-scan (delta = pi/2) [D vs C]
rows = []
for mm in (10, 15, 20, 30, 40):
    r = B[(f'B_m{mm}', round(np.pi / 2, 3))]
    x0b = r['prep_info']['x']; mL = r['masses']['mYL'] / H
    rows.append(dict(m=mm, kappa=r['kappa'], eps=r['eps'], NYL=r['NYL'], NYH=r['NYH'], NY=r['NYL'] + r['NYH'],
                     kink_estimate=9 * CH**2 * r['kappa']**2 * x0b**6 / mL**7, R_Y=r['R_Y'],
                     R_central=r['R_central'], R_L=r['R_L'], mL_over_omegaX=mL / omegaX))
for a, b in zip(rows, rows[1:]):
    b['slope_NY'] = -np.log(b['NY'] / a['NY']) / np.log(b['m'] / a['m'])
    b['slope_R_Y'] = -np.log(b['R_Y'] / a['R_Y']) / np.log(b['m'] / a['m'])
out['heavy_m_scan'] = dict(rows=rows, predicted_slopes=dict(NY=5.0, R_Y=5.5),
                           note='kappa tied to m at fixed eta/eta_edge = 0.5, r = 0.45, h = 40 H_i')

# 4. vacuum-stability bound and the resulting upper end of the heavy window [D]
eps_t, r_t = 0.15, 0.45
kmax = S.kappa_max(cp0)
out['window_upper_end'] = dict(kappa_max=kmax, m_max_formula='kappa_max/lam_X * sqrt(2 r cH/(eps (1-r^2)))',
                               m_max=kmax / LAM * np.sqrt(2 * r_t * CH / (eps_t * (1 - r_t**2))),
                               m_max_frac_half_edge=kmax / (LAM * np.sqrt(0.5 * 4 * (1 - r_t**2) / (12 * CH))))

# 5. survival coefficient fixed by the charge requirement [D]: c6 = 2 eta lam^2 M_red^2/(cH H^2)
eta = 0.5 * r_t / 12
c6f = lambda Hs: 2 * eta * LAM**2 * MRED**2 / (CH * Hs**2)
out['survival'] = dict(eta=eta, c6_formula_at_Hi=c6f(H), c6_scan=B[('B_m30', round(np.pi / 2, 3))]['c6'],
                       c6_Hscaled=B[('B_base_Hscaled', round(np.pi / 2, 3))]['c6'],
                       H_min_for_ceiling_GeV=H * np.sqrt(c6f(H) / 1.486e13), c6_stated=S.c6(cp0)['c6'])

# 6. scaling symmetry (M, h, mu, H) -> s(M, h, mu, H): check identical eps and R, yield ~ s^1.5, c6 ~ s^-2 [D, C]
a, b = B[('B_m30', round(np.pi / 2, 3))], B[('B_base_Hscaled', round(np.pi / 2, 3))]
s = b['H_GeV'] / a['H_GeV']
out['scaling_symmetry'] = dict(s=s, eps=[a['eps'], b['eps']], R_central=[a['R_central'], b['R_central']],
                               yield_ratio=b['yield_over_YQ'] / a['yield_over_YQ'], s_pow_1p5=s**1.5,
                               c6_ratio=b['c6'] / a['c6'], s_pow_m2=s**-2)

# 7. relative-phase mass during inflation in the heavy direction [D]: m^2 = 18 eta cH (H_i units)
out['phase_mass_over_H'] = float(np.sqrt(18 * eta * CH))

# 8. string tension scaling (conditional on string formation) [Q 3.1e-15 at h = 1e11 GeV, scaled as h^2]
hs = {'stated': 1e11, 'h40_at_Hi': 40 * H, 'h40_at_Hscaled': 40 * b['H_GeV'], 'h20_at_Hscaled': 20 * b['H_GeV']}
out['Gmu_scaled'] = {k: 3.1e-15 * (v / 1e11)**2 for k, v in hs.items()}

# 9. value of one vector quantum per comoving H_i^3 in R units at the heavy spot point [D]
out['R_per_vector_quantum_heavy'] = a['R_cascade'] / a['NYH']
out['vector_allowance_per_net_charge'] = (a['R_L'] - a['R_central']) / (out['R_per_vector_quantum_heavy'] * abs(a['NX']))

# 10. stated-regime one-at-a-time scan summary if present
fa = R / 'scan_A.json'
if fa.exists():
    A = json.load(open(fa))['rows']
    summ = {}
    for r in A:
        if 'eps' not in r:
            continue
        k = r['name'] + '|' + r['prep']
        best = summ.get(k)
        ratio = r['R_lower'] / r['R_L'] if r['R_L'] else None
        cand = dict(delta=r['delta'], eps=r['eps'], R_lower=r['R_lower'], R_central=r['R_central'], R_L=r['R_L'],
                    lower_over_RL=ratio, yield_over_YQ=r['yield_over_YQ'], NY=r['NYH'] + r['NYL'], Nsig=r['Nsig'],
                    c6=r['c6'])
        if ratio is not None and (best is None or best['lower_over_RL'] is None or ratio < best['lower_over_RL']):
            summ[k] = cand
        elif best is None:
            summ[k] = cand
    out['stated_regime_best_phase'] = summ

(R / 'analytic.json').write_text(json.dumps(out, indent=1, default=float))
print(json.dumps(out, indent=1, default=float))
