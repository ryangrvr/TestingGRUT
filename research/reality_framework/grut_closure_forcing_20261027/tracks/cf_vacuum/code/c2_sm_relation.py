"""Freeze the SM absolute-stability (criticality) relation with PDG 2026 inputs; tension; instability scales.

Track cf_vacuum (closure-forcing pass, 27 Oct 2026).  Run: python3 -B c2_sm_relation.py  (~6 min).
Inputs: ../results/c1_fit_grid.json (SMDR fits, c1_fit_grid.py).  Writes ../results/sm_relation.json and
../results/frozen_relation.json.

Method (same pipeline as verifier v_vacua's matching (S), which reproduced Hiller et al. 2401.08811 to 0.01 GeV):
  * MSbar couplings at Q = M_t from SMDR's full fit (Martin-Robertson 1907.02500), interpolated: quadratic in M_t
    (7 fits), linear in M_h and alpha_s (finite differences at M_t = 171.0).
  * RG running with SMDR at loop order 3 (and 5 = 3-loop + leading 4/5-loop pieces) up to M_Pl = 1.2209e19 GeV.
  * Criterion 'msbar': min_{Q <= M_Pl} lambda(Q) = 0 <=> lambda(mu*) = beta_lambda(mu*) = 0 (Bednyakov et al.
    1507.08833): degenerate second minimum = boundary of absolute stability.
  * Instability scales and barrier: own 3-loop RGEs (sm_rge.py, track `vacua`, transcribed from Buttazzo et al.
    1307.3536 App. A/B; verified against SMDR by v_vacua) started from the SMDR-matched couplings; Landau-gauge
    lambda_eff at 1 and 2 loops; V_eff = e^{4 Gamma} lambda_eff h^4 / 4.
"""
import json
import os

import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm

import smdr_common as S
import sm_rge as R

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
MPL = 1.2209e19
MH0, SMH = 125.13, 0.11
AS0, SAS = 0.1180, 0.0009
MT_DIRECT, SMT_DIRECT = 172.60, 0.27
MT_XSEC, SMT_XSEC = 172.1, 0.6
MT_ATLCMS, SMT_ATLCMS = 172.52, 0.33          # ATLAS+CMS 7+8 TeV combination, arXiv:2402.08713
S_MATCH = 0.113                               # track vacua: delta lambda = 3e-4, delta y_t = 5e-4 in quadrature
S_MATCH_EXTRA = 0.088                         # v_vacua C3: 2013-fit vs SMDR matching shift, as an added error
S_CRIT = 0.032                                # half of ('eff' - 'msbar')

grid = json.load(open(os.path.join(RES, 'c1_fit_grid.json')))['fits']
KEYS = ['g3', 'g', 'gp', 'yt', 'yb', 'ytau', 'k']


def build_matching(fits, mh0, as0, mt_ref):
    base = [f for f in fits if abs(f['inputs'][1] - mh0) < 1e-9 and abs(f['inputs'][2] - as0) < 1e-9]
    mts = np.array([f['inputs'][0] for f in base])
    coef = {kk: np.polyfit(mts, [f[kk] for f in base], 2) for kk in KEYS}
    get = lambda p: [f for f in fits if np.allclose(f['inputs'], p)][0]
    dmh = {kk: (get((mt_ref, mh0 + 0.5, as0))[kk] - get((mt_ref, mh0 - 0.5, as0))[kk]) / 1.0 for kk in KEYS}
    das = {kk: (get((mt_ref, mh0, as0 + 0.0009))[kk] - get((mt_ref, mh0, as0 - 0.0009))[kk]) / 0.0018 for kk in KEYS}

    def m(Mt, Mh, als):
        return {kk: float(np.polyval(coef[kk], Mt) + dmh[kk] * (Mh - mh0) + das[kk] * (als - as0)) for kk in KEYS}
    resid = {str(f['inputs']): {kk: m(*f['inputs'])[kk] - f[kk] for kk in ('g3', 'yt', 'k')} for f in fits}
    return m, resid


match, interp_resid = build_matching(grid, MH0, AS0, 171.0)


def mtcrit(Mh, als, L=3, m=match):
    f = lambda Mt: S.min_lambda(L, Mt, m(Mt, Mh, als), Qcut=MPL)[0]
    Mt = brentq(f, 168.0, 175.5, xtol=1e-6)
    return Mt, S.min_lambda(L, Mt, m(Mt, Mh, als), Qcut=MPL)[1]


out = {'inputs': {'PDG2026': {'Mh': [MH0, SMH], 'Mt_direct': [MT_DIRECT, SMT_DIRECT], 'Mt_xsec_pole': [MT_XSEC, SMT_XSEC],
                              'alpha_s': [AS0, SAS]},
                  'ATLAS+CMS 2402.08713': [MT_ATLCMS, SMT_ATLCMS], 'FLAG 2024 alpha_s': [0.1183, 0.0007]},
       'interpolation_residuals': interp_resid}

# ------------------------------------------------------------------ 1. the line
Mc3, mu3 = mtcrit(MH0, AS0, 3)
Mc5, mu5 = mtcrit(MH0, AS0, 5)
Mc2, _ = mtcrit(MH0, AS0, 2)
dMh = (mtcrit(MH0 + 0.5, AS0)[0] - mtcrit(MH0 - 0.5, AS0)[0]) / 1.0
das = (mtcrit(MH0, AS0 + 0.0009)[0] - mtcrit(MH0, AS0 - 0.0009)[0]) / 2.0
out['line'] = {'Mt_crit_loop3': Mc3, 'mu_star_loop3': mu3, 'Mt_crit_loop5': Mc5, 'mu_star_loop5': mu5,
               'Mt_crit_loop2': Mc2, 'loop2_minus_loop3': Mc2 - Mc3, 'loop5_minus_loop3': Mc5 - Mc3,
               'dMtcrit_dMh': dMh, 'dMtcrit_per_0.0009_as': das,
               'formula': f'Mt <= {Mc5:.3f} + {dMh:.3f} (Mh - {MH0}) + {das:.3f} (alpha_s - {AS0})/0.0009 GeV '
                          f'(+-{S_MATCH} th matching, +-{S_CRIT} criterion)'}
# cross-check: v_vacua's grid at the PDG 2025 point (read-only)
VG = os.path.join(HERE, '..', '..', '..', '..', 'grut_q_vacuum_construction_20261021', 'tracks', 'v_vacua', 'results',
                  'v3a_smdr_fit_grid.json')
if os.path.exists(VG):
    vfits = json.load(open(VG))['fits']
    base25 = [f for f in vfits if abs(f['inputs'][1] - 125.20) < 1e-9 and abs(f['inputs'][2] - 0.1180) < 1e-9]
    m25, _ = build_matching(vfits, 125.20, 0.1180, 171.2)
    out['crosscheck_PDG2025_point'] = {'this code on v_vacua grid (125.20, 0.1180), loop 3': mtcrit(125.20, 0.1180, 3, m25)[0],
                                       'v_vacua VERIFY C1 (loop 3 / 5)': [171.026, 171.031],
                                       'this grid extrapolated to (125.20, 0.1180)': mtcrit(125.20, 0.1180, 3)[0],
                                       'Hiller et al. 2401.08811 (PDG 2024 = 125.20, 0.1180)': 171.04}

# Buttazzo et al. 1307.3536 linear formula (TeX eq. stability) and the formula quoted in the brief
butt = lambda Mh, als: 173.34 + (Mh - 129.6) / 2.0 + 0.5 / 2.0 * (als - 0.1184) / 0.0007
brief = lambda Mh, als: 171.53 + 0.15 * (Mh - 125.15) - 0.34 * (als - 0.1184) / 0.0007
out['literature_formulas'] = {
    'Buttazzo 1307.3536 eq. (stability) inverted: Mt < 173.34 + (Mh-129.6)/2.0 + 0.25 (as-0.1184)/0.0007 (+-0.15 th)':
        butt(MH0, AS0),
    'brief formula 171.53 + 0.15 (Mh-125.15) - 0.34 (as-0.1184)/0.0007 (not in the paper; sign of alpha_s term wrong)':
        brief(MH0, AS0),
    'Buttazzo quoted 171.53 +- 0.42 at (125.15, 0.1184); own formula there gives': butt(125.15, 0.1184),
    'Bednyakov et al. 1507.08833 at (125.7, 0.1185)': 171.44, 'this code at (125.7, 0.1185)': mtcrit(125.7, 0.1185)[0]}


# ------------------------------------------------------------------ 2. tension
def ztab(Mtc):
    rows = {}
    s_par = np.hypot(dMh * SMH, das * SAS / 0.0009)
    sets = {'direct PDG2026 172.60+-0.27': (MT_DIRECT, SMT_DIRECT), 'ATLAS+CMS 172.52+-0.33': (MT_ATLCMS, SMT_ATLCMS),
            'xsec pole PDG2026 172.1+-0.6': (MT_XSEC, SMT_XSEC)}
    budgets = {'exp+param+matching': 0.0, '+criterion': S_CRIT, '+criterion+0.3 MC': np.hypot(S_CRIT, 0.3),
               '+criterion+0.6 MC': np.hypot(S_CRIT, 0.6), '+criterion+0.088 matching shift': np.hypot(S_CRIT, S_MATCH_EXTRA)}
    for lab, (Mt, sMt) in sets.items():
        for blab, extra in budgets.items():
            if lab.startswith('xsec') and 'MC' in blab:
                continue
            s = np.sqrt(sMt ** 2 + s_par ** 2 + S_MATCH ** 2 + extra ** 2)
            z = (Mt - Mtc) / s
            rows[f'{lab} | {blab}'] = {'Delta_GeV': Mt - Mtc, 'sigma_tot': s, 'z': z, 'p_one_sided': float(norm.sf(z))}
    return rows


out['tension_PDG2026'] = ztab(Mc5)
# MC-to-pole offset scan: M_t,pole = M_t,MC + d
s_par = np.hypot(dMh * SMH, das * SAS / 0.0009)
out['MC_offset_scan'] = {f'{d:+.1f}': (MT_DIRECT + d - Mc5) / np.sqrt(SMT_DIRECT ** 2 + s_par ** 2 + S_MATCH ** 2 + S_CRIT ** 2)
                         for d in np.arange(-1.6, 0.95, 0.2)}
out['MC_offset_for_z0'] = Mc5 - MT_DIRECT
# FLAG alpha_s variant
McF = Mc5 + das * (0.1183 - AS0) / 0.0009
sF = np.sqrt(SMT_DIRECT ** 2 + (dMh * SMH) ** 2 + (das * 0.0007 / 0.0009) ** 2 + S_MATCH ** 2 + S_CRIT ** 2)
out['FLAG_alpha_s_variant'] = {'Mt_crit': McF, 'z_direct': (MT_DIRECT - McF) / sF, 'sigma_tot': sF}
# needed alpha_s or Mh for stability at central Mt
out['alpha_s_needed_at_Mt_direct'] = AS0 + (MT_DIRECT - Mc5) / das * 0.0009
out['Mh_needed_at_Mt_direct'] = brentq(lambda Mh: S.min_lambda(5, MT_DIRECT, match(MT_DIRECT, Mh, AS0))[0], 120, 140,
                                       xtol=1e-6)
# held-out (data do not exist): FCC-ee and HL-LHC precisions as in track vacua / v_vacua
def heldout(sMt, sas, sMh_, sm=S_MATCH, Mt=MT_DIRECT):
    s = np.sqrt(sMt ** 2 + (dMh * sMh_) ** 2 + (das * sas / 0.0009) ** 2 + sm ** 2 + S_CRIT ** 2)
    return {'sigma_tot': s, 'z_if_central_values_stay': (Mt - Mc5) / s}


out['heldout_illustrative'] = {'FCC-ee (40 MeV, 0.00012, 4 MeV)': heldout(0.04, 0.00012, 0.004),
                               'FCC-ee, matching halved': heldout(0.04, 0.00012, 0.004, sm=S_MATCH / 2),
                               'HL-LHC + lattice (0.2 [R], 0.00045, 0.02)': heldout(0.2, 0.00045, 0.02),
                               'FCC-ee, xsec central 172.1': heldout(0.04, 0.00012, 0.004, Mt=MT_XSEC)}


# ------------------------------------------------------------------ 3. instability scales and barrier
def scales(Mt, Mh=MH0, als=AS0):
    c = match(Mt, Mh, als)
    c0 = dict(g1=np.sqrt(5 / 3) * c['gp'], g2=c['g'], g3=c['g3'], yt=c['yt'], yb=c['yb'], ytau=c['ytau'], lam=c['k'])
    run = R.run(c0, Mt, mu_max=1.3e19, loops=3, n=1600)
    lnmu = np.log(run['mu'])
    lam = run['Y'][6]
    # MSbar lambda zero
    i = np.where(lam < 0)[0]
    res = {'Mt': Mt}
    if len(i):
        j = i[0]
        res['Lambda_lambda_GeV'] = float(np.exp(brentq(lambda x: run['sol'].sol(x - np.log(Mt))[6], lnmu[j - 1], lnmu[j])))
    le1 = np.array([R.lambda_eff(run['Y'][:, k], two_loop=False) for k in range(len(lnmu))])
    le2 = np.array([R.lambda_eff(run['Y'][:, k], two_loop=True) for k in range(len(lnmu))])
    for lab, le in (('1loop', le1), ('2loop', le2)):
        neg = np.where(le < 0)[0]
        if len(neg):
            j = neg[0]
            res[f'Lambda_I_{lab}_GeV'] = float(np.exp(np.interp(0.0, [le[j], le[j - 1]], [lnmu[j], lnmu[j - 1]])))
        V = np.exp(4 * run['Y'][7]) * le * run['mu'] ** 4 / 4
        k = int(np.argmax(np.where(run['mu'] < 1e16, V, -np.inf)))
        if 0 < k < len(V) - 1 and V[k] > 0 and len(neg):
            res[f'h_max_{lab}_GeV'] = float(run['mu'][k])
            res[f'Lambda_V_{lab}_GeV'] = float(V[k] ** 0.25)
            res[f'log10_Lambda_V_{lab}'] = float(np.log10(V[k] ** 0.25))
        res[f'min_lambda_eff_{lab}'] = float(le.min())
    res['min_lambda_msbar'] = float(lam[run['mu'] <= MPL].min())
    res['V_eff_2loop_at_2.435e18_GeV4'] = float(np.interp(np.log(2.435e18), lnmu,
                                                          np.exp(4 * run['Y'][7]) * le2 * run['mu'] ** 4 / 4))
    return res


out['scales'] = {f'Mt={Mt}': scales(Mt) for Mt in (MT_DIRECT - SMT_DIRECT, MT_DIRECT, MT_DIRECT + SMT_DIRECT, MT_XSEC,
                                                    171.5)}
# validation of the scale code against Buttazzo's fit log10 Lambda_V = 9.5 + ... at their reference inputs
bref = scales(173.34, 125.15, 0.1184)
out['scales_validation_1307.3536_point'] = {'this': bref, 'Buttazzo log10 Lambda_V': 9.5,
                                            'Buttazzo Lambda_lambda ~ 2 Lambda_V, Lambda_I ~ 13 Lambda_V': True}
butt_LV = lambda Mh, Mt, als: 9.5 + 0.7 * (Mh - 125.15) - 1.0 * (Mt - 173.34) + 0.3 * (als - 0.1184) / 0.0007
out['Buttazzo_fit_log10_Lambda_V_at_PDG2026_direct'] = butt_LV(MH0, MT_DIRECT, AS0)

json.dump(out, open(os.path.join(RES, 'sm_relation.json'), 'w'), indent=1, default=float)

frozen = {
    'frozen_on': '2026-10-27', 'by': 'Claude track cf_vacuum (closure-forcing pass, issue #13); git read-only, no commit',
    'status': 'retrospective for all data listed in inputs; held out only for measurements published after this date',
    'postulate': 'P_q reading (a)/(c)-worldline: the q-theory final-state equilibrium names the global minimum of the full '
                 'renormalized effective potential; with the SM valid to M_Pl this is absolute stability',
    'relation': {'Mt_crit_GeV': Mc5, 'at': {'Mh': MH0, 'alpha_s': AS0}, 'dMt_dMh': dMh, 'dMt_per_0.0009_alpha_s': das,
                 'theory_error_matching_GeV': S_MATCH, 'criterion_error_GeV': S_CRIT, 'mu_star_GeV': mu5,
                 'criterion': 'min_{mu<=1.2209e19 GeV} lambda_MSbar(mu) = 0 (SMDR v1.3 fit + loop order 5 running)',
                 'inequality': 'M_t(pole) <= Mt_crit + dMt_dMh (Mh - 125.13) + dMt_per (alpha_s - 0.1180)/0.0009',
                 'equality_variant_MPP': 'M_t(pole) = Mt_crit (reading (d), Froggatt-Nielsen; added assumption)'},
    'discriminator': 'z = (M_t,pole - Mt_crit(Mh, alpha_s)) / sigma_tot with a short-distance top mass (threshold scan '
                     'MSR/PS mass converted to pole/MSbar at N3LO) and alpha_s(M_Z) from Z-pole or lattice',
    'failure_condition': 'z > 5 with a theoretically clean short-distance top mass falsifies P_q(a) + SM-to-M_Pl; P_q(a) '
                         'then requires new physics coupled to the Higgs at or below the instability scale',
    'inputs_used': out['inputs']}
json.dump(frozen, open(os.path.join(RES, 'frozen_relation.json'), 'w'), indent=1, default=float)
print(json.dumps({k: out[k] for k in ('line', 'literature_formulas', 'tension_PDG2026', 'MC_offset_for_z0',
                                      'FLAG_alpha_s_variant', 'alpha_s_needed_at_Mt_direct', 'Mh_needed_at_Mt_direct',
                                      'heldout_illustrative', 'scales', 'scales_validation_1307.3536_point',
                                      'Buttazzo_fit_log10_Lambda_V_at_PDG2026_direct')} | {'x': out.get('crosscheck_PDG2025_point')},
                 indent=1, default=float))
