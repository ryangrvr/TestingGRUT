"""P_q applied to the GRUT hidden scalar sector; radiative costs of the Sigma window; inflation scale vs barrier.

Track cf_vacuum (closure-forcing pass, 27 Oct 2026).  Run: python3 -B c4_hidden_and_costs.py  (~10 s).
Writes ../results/hidden_and_costs.json.  Uses ../results/threshold.json and ../results/sm_relation.json.

Heavy reference (Codex THEORY_DATA.json, 04c53f5): H_i = 1.33350e11, M_Y = 4.00049e12, h = 5.33398e12,
mu = 9.54711e11 GeV, kappa = 3.4587e-3, lam_X = 1e-3, lam_Y = lam_Sigma = 1, e = 1e-7; lam_SigmaY = 0 at matching
(ACTION.md: 'other coefficients were set to zero at a matching reference').
Potential used (aligned phase, X = 0):  V = lam_S (s^2 - h^2/2)^2 + (M_Y^2 + lam_SY s^2 - 2 mu s)|Y|^2 + lam_Y |Y|^4,
s = |Sigma|.  Where m_eff^2(s) = M_Y^2 + lam_SY s^2 - 2 mu s < 0, Y condenses: V_min(s) = lam_S (s^2 - h^2/2)^2
- m_eff^4 / (4 lam_Y).
"""
import json
import os

import numpy as np
from scipy.optimize import brentq, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
H_I = 133349619398.33566
M_Y = 4000488581950.07
HV = 5333984775933.427
MU = 954710753970.2344
out = {'inputs': {'H_i': H_I, 'M_Y': M_Y, 'h': HV, 'mu': MU, 'r = sqrt2 mu h / M_Y^2': np.sqrt(2) * MU * HV / M_Y ** 2}}


def vmin_over_s(lamS, lamY, lamSY, smax_fac=200.0):
    """min over s of V(s, Y*) - V(ours); ours: s0 = h/sqrt2, Y = 0 (V = 0). Returns (min, s_at_min)."""
    s0 = HV / np.sqrt(2)
    s = s0 * np.logspace(-3, np.log10(smax_fac), 20000)
    meff2 = M_Y ** 2 + lamSY * s ** 2 - 2 * MU * s
    V = lamS * (s ** 2 - s0 ** 2) ** 2 - np.where(meff2 < 0, meff2 ** 2 / (4 * lamY), 0.0)
    k = int(np.argmin(V))
    return float(V[k]), float(s[k] / s0)


def lamS_crit(lamY, lamSY):
    f = lambda lg: vmin_over_s(10 ** lg, lamY, lamSY)[0]
    if f(-30) >= 0:
        return 0.0
    if f(2) < 0:
        return np.inf
    return float(10 ** brentq(lambda lg: 1.0 if f(lg) >= 0 else -1.0, -30, 2, xtol=1e-6))


# check: our vacuum is a local minimum (Y mass^2 > 0 there)
s0 = HV / np.sqrt(2)
out['Y_mass2_at_our_vacuum_over_MY2'] = (M_Y ** 2 - 2 * MU * s0) / M_Y ** 2
out['s_where_Y_turns_tachyonic_over_s0'] = M_Y ** 2 / (2 * MU) / s0
crit = {}
for lamY in (1.0, 2.0, 4 * np.pi):
    for lamSY in (0.0, 0.01, 0.03, 0.05, (MU / M_Y) ** 2):
        crit[f'lamY={lamY:.3g}, lamSY={lamSY:.4g}'] = lamS_crit(lamY, lamSY)
out['lamSigma_crit_for_absolute_stability_tree'] = crit
out['analytic_lamSY_removing_the_Y_branch (mu^2/M_Y^2)'] = (MU / M_Y) ** 2
out['reference_point_lamSigma=1_lamY=1: min V, s/s0'] = vmin_over_s(1.0, 1.0, 0.0)

# minimal lam_SY for the P_q Sigma window values of lam_Sigma
win = json.load(open(os.path.join(RES, 'threshold.json')))['C2_GRUT_Sigma_fixed_vev']
need = {}
for key in ('1e+06', '1e+08', '1e+09', '1e+10', '3e+10', '4e+10'):
    lamS = win[key]['lam_Sigma']
    f = lambda lsy: vmin_over_s(lamS, 1.0, lsy)[0]
    lsy = brentq(lambda x: 1.0 if f(x) >= 0 else -1.0, 0.0, 0.2, xtol=1e-7) if f(0.0) < 0 else 0.0
    need[f'M_Sigma={key}'] = {'lam_Sigma': lamS, 'min_lam_SigmaY (lam_Y = 1)': lsy,
                              'depth_without_lamSY_GeV4': vmin_over_s(lamS, 1.0, 0.0)[0],
                              'depth/rho_Lambda': vmin_over_s(lamS, 1.0, 0.0)[0] / 2.52e-47}
out['Sigma_window_needs_lam_SigmaY'] = need

# one-loop pieces feeding lam_Sigma and m_Sigma^2 from the Y sector [D: Coleman-Weinberg expansion, Y masses
# M_Y^2 + lam_SY s^2 +- 2 mu s]
K16 = 16 * np.pi ** 2
dl_mu4 = -MU ** 4 / (24 * np.pi ** 2 * M_Y ** 4)
costs = {'delta_lamSigma_from_mu^4 (threshold at M_Y)': dl_mu4}
for lsy in (0.0, (MU / M_Y) ** 2):
    for key in ('1e+06', '1e+09', '3e+10'):
        MS = float(key)
        lamS = win[key]['lam_Sigma']
        loop = abs(dl_mu4 + lsy * MU ** 2 / (8 * np.pi ** 2 * M_Y ** 2) + lsy ** 2 * np.log(M_Y ** 2 / MS ** 2) / (32 * np.pi ** 2))
        dm2 = MU ** 2 * np.log(M_Y ** 2 / MS ** 2) / (8 * np.pi ** 2) + lsy * M_Y ** 2 / K16 * np.log(M_Y ** 2 / MS ** 2)
        costs[f'lamSY={lsy:.3g}, M_Sigma={key}'] = {'|loop delta lam_Sigma|': loop, 'lam_Sigma needed': lamS,
                                                    'quartic tuning lam/loop': lamS / loop,
                                                    'delta m_Sigma^2 (log) GeV^2': dm2, 'm_Sigma^2': MS ** 2,
                                                    'mass tuning': MS ** 2 / dm2}
out['radiative_costs'] = costs
# Higgs mass parameter shift from the portal needed in the window: lam_HSigma h^2 / 2 against m_h^2 / 2
hig = {}
for key in ('1e+03', '1e+06', '1e+09', '3e+10'):
    r = win[key]['sign=+1']
    for lab in ('at_delta_min', 'at_delta_max'):
        if lab in r:
            lHS = abs(r[lab]['lam_HSigma_GRUTnorm'])
            hig[f'M_Sigma={key} {lab}'] = {'lam_HSigma': lHS, 'theta': r[lab]['theta_mixing'],
                                          'delta m_H^2 = lam_HSigma h^2/2 (GeV^2)': lHS * HV ** 2 / 2,
                                          'tuning m_h^2/(2 delta m_H^2)': 125.13 ** 2 / (lHS * HV ** 2)}
out['Higgs_mass_tuning_from_Sigma_portal'] = hig

# inflation scale against the SM barrier (standard-physics obligation, independent of P_q)
sc = json.load(open(os.path.join(RES, 'sm_relation.json')))['scales']
inf = {}
for k, v in sc.items():
    inf[k] = {'h_max_2loop': v.get('h_max_2loop_GeV'), 'H_i/h_max': H_I / v['h_max_2loop_GeV'] if v.get('h_max_2loop_GeV') else None,
              'Lambda_lambda': v.get('Lambda_lambda_GeV'), 'H_i/Lambda_lambda': H_I / v['Lambda_lambda_GeV']}
out['inflation_vs_barrier'] = inf
out['inflation_note'] = ('random-walk Higgs during N = 60 e-folds: h_rms ~ sqrt(N) H/(2 pi) = %.3g GeV; with H_i ~ 0.5 h_max '
                         'the SM direction is not safe without a Hubble-induced Higgs mass (xi_H R |H|^2 with xi_H >~ 0.05-0.07, '
                         'Mantziris-Markkanen-Rajantie 2011.03763 [Q abstract]) or a stabiliser') % (np.sqrt(60) * H_I / (2 * np.pi))
json.dump(out, open(os.path.join(RES, 'hidden_and_costs.json'), 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
