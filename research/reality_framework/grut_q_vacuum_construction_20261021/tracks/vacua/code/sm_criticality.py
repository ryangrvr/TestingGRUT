"""SM stability bound and multiple-point (criticality) relation, validation, retrospective tension, held-out test.

Track `vacua` (NOTES.md section 5). Run:  python3 sm_criticality.py   -> ../results/sm_criticality.json

Criteria (both used; their difference is part of the theory error):
  'msbar': lambda(mu*) = beta_lambda(mu*) = 0 at some mu* <= mu_cut, i.e. min_{mu<=mu_cut} lambda(mu) = 0
           (Bednyakov-Kniehl-Pikelner-Veretin, arXiv:1507.08833, eq. (crit); gauge independent).
  'eff'  : min_{h<=mu_cut} lambda_eff(h) = 0 with the Landau-gauge 2-loop effective potential
           (Buttazzo et al. arXiv:1307.3536, sec. 4.2 and App. B; Degrassi et al. arXiv:1205.6497).
mu_cut = M_Pl = 1.2209e19 GeV (as in 1307.3536 and 1507.08833) unless stated; reduced 2.435e18 GeV also given.
"""
import json, os, time
import numpy as np
from scipy.optimize import brentq, minimize_scalar, fsolve
import sm_rge as S

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'results', 'sm_criticality.json')
MPL = 1.2209e19
MPL_RED = 2.435e18

# ---------------------------------------------------------------- inputs with provenance
PDG = {
    'Mh': (125.20, 0.11, 'PDG 2025 (pdg-2025-v0.2.3.sqlite, S126M, AC, S=1.36), read level: database value'),
    'Mt_direct': (172.56, 0.31, 'PDG 2025 Q007TP "t-Quark Mass (Direct Measurements)", 172.56+-0.31 (S=1.6); '
                                'MC-template mass, interpreted as pole mass with an extra theory error'),
    'Mt_xsec_pole': (172.4, 0.7, 'PDG 2025 Q007TP4 "t-Quark Pole Mass from Cross-Section Measurements"'),
    'alpha_s': (0.1180, 0.0009, 'PDG 2025 QCD review (rev. Aug 2025) eq. (9.25), section read'),
    'MW': (80.3692, 0.0133, 'PDG 2025 S043M displayed average'),
}
FUTURE = {
    'Mt': (0.040, 'FCC-ee ttbar threshold: 6.8 MeV exp (PS mass) + ~35 MeV N3LO theory, Defranchis et al. '
                  'arXiv:2503.18713 abstract; FCC FSR vol.1 arXiv:2505.00272 sec. on theory: O(50) MeV'),
    'alpha_s': (0.00012, 'FCC-ee Z-pole fit, 0.1% (FCC FSR vol.1, arXiv:2505.00272, "Measurements of the strong '
                         'coupling constant" paragraph); lattice: factor ~2 improvement in a decade (same source)'),
    'Mh': (0.004, 'FCC-ee recoil, O(4 MeV) (FCC FSR vol.1, arXiv:2505.00272); HL-LHC ~10-20 MeV [R]'),
}
TH = {'dlam': 0.00030, 'dyt': 0.00050,
      'src': 'Buttazzo et al. 1307.3536 sec. 3: theory errors of the NNLO matching fits'}


def run_to(Mt, Mh, a, MW=80.384, mu_max=2.0e19, loops=3, **kw):
    return S.run(S.matching(Mt, Mh, a, MW, **kw), Mt, mu_max=mu_max, loops=loops, n=60)


def min_lambda(Mt, Mh, a, crit='msbar', mu_cut=MPL, MW=80.384, loops=3, **kw):
    r = run_to(Mt, Mh, a, MW, mu_max=max(2.0e19, 1.5 * mu_cut), loops=loops, **kw)
    sol = r['sol'].sol
    lo, hi = np.log(1.0e7 / Mt), np.log(mu_cut / Mt)
    if crit == 'msbar':
        f = lambda t: sol(t)[6]
    else:
        f = lambda t: S.lambda_eff(sol(t))
    ts = np.linspace(lo, hi, 61)
    vals = np.array([f(t) for t in ts])
    i = int(np.argmin(vals))
    if 0 < i < len(ts) - 1:
        m = minimize_scalar(f, bounds=(ts[i - 1], ts[i + 1]), method='bounded', options={'xatol': 1e-6})
        tmin, vmin = m.x, m.fun
    else:
        tmin, vmin = ts[i], vals[i]
    return float(vmin), float(Mt * np.exp(tmin)), bool(0 < i < len(ts) - 1)


def Mt_crit(Mh, a, crit='msbar', mu_cut=MPL, MW=80.384, loops=3, **kw):
    g = lambda Mt: min_lambda(Mt, Mh, a, crit, mu_cut, MW, loops, **kw)[0]
    Mt = brentq(g, 168.0, 176.0, xtol=1e-5)
    v, mu_star, interior = min_lambda(Mt, Mh, a, crit, mu_cut, MW, loops, **kw)
    return float(Mt), mu_star, interior


def Mh_crit(Mt, a, crit='msbar', mu_cut=MPL, MW=80.384, **kw):
    g = lambda Mh: min_lambda(Mt, Mh, a, crit, mu_cut, MW, **kw)[0]
    return float(brentq(g, 118.0, 142.0, xtol=1e-5))


def fn_double(a, mu_star, MW=80.384):
    """Froggatt-Nielsen 1995 form: lambda(mu*) = beta_lambda(mu*) = 0 at a fixed mu* -> (Mt, Mh)."""
    def eqs(x):
        Mt, Mh = x
        r = run_to(Mt, Mh, a, MW, mu_max=max(2.0e19, 1.5 * mu_star))
        Y = r['sol'].sol(np.log(mu_star / Mt))
        return [Y[6] * 100, S.beta_lambda_lnmu(Y) * 1e3]
    x = fsolve(eqs, [171.0, 128.0], xtol=1e-10)
    return [float(x[0]), float(x[1])], [float(e) for e in eqs(x)]


def main():
    t0 = time.time()
    out = {'inputs': PDG, 'future_precision': FUTURE, 'matching_theory_errors': TH,
           'conventions': 'V = -m^2|H|^2/2 + lambda|H|^4; MSbar; M_Pl = 1.2209e19 GeV unless stated'}

    # ---------------- validation against Buttazzo et al. (2013) Planck-scale values (their inputs)
    ref = (173.34, 125.15, 0.1184)
    r = run_to(*ref)
    cP = S.couplings_at(r, MPL)
    val = {'mine': {k: round(float(v), 5) for k, v in cP.items()},
           'Buttazzo_eq_sysPl': {'g1': 0.6154, 'g2': 0.5055, 'g3': 0.4873, 'yt': 0.3825, 'lam': -0.0143},
           'source': '1307.3536 sec. 4.1, values at M_Pl for Mt=173.34, Mh=125.15, as=0.1184, MW=80.384'}
    dl = {}
    b = cP['lam']
    dl['dlam_dMt_perGeV'] = float(S.couplings_at(run_to(174.34, 125.15, 0.1184), MPL)['lam'] - b)
    dl['dlam_das_per0.0007'] = float(S.couplings_at(run_to(173.34, 125.15, 0.1191), MPL)['lam'] - b)
    dl['dlam_dMh_perGeV'] = float(S.couplings_at(run_to(173.34, 126.15, 0.1184), MPL)['lam'] - b)
    dl['Buttazzo'] = {'dlam_dMt_perGeV': -0.0066, 'dlam_das_per0.0007': 0.0018, 'dlam_dMh_perGeV': 0.0029}
    val['derivatives'] = dl
    lam = r['Y'][6]
    j = np.where(np.diff(np.sign(lam)))[0][0]
    val['Lambda_lambda_GeV(lambda=0 crossing)'] = float(r['mu'][j])
    val['Buttazzo_Lambda_lambda'] = 'Lambda_V = 10^9.5 GeV and Lambda_lambda ~ 2 Lambda_V ~ 6e9 GeV (sec. 4.2)'
    loops = {}
    for L in (1, 2, 3):
        loops[str(L)] = float(S.couplings_at(run_to(*ref, loops=L), MPL)['lam'])
    val['lam_MPl_by_RGE_loop_order'] = loops
    # effective-potential scales (Buttazzo sec. 4.2: log10 Lambda_V = 9.5, Lambda_I ~ 13 Lambda_V, Landau gauge)
    sol = r['sol'].sol
    fle = lambda lm: S.lambda_eff(sol(np.log(10 ** lm / ref[0])))
    zI = brentq(fle, 8, 13)
    lms = np.linspace(6, 12, 600)
    Vv = [np.exp(4 * sol(np.log(10 ** lm / ref[0]))[7]) * fle(lm) * (10 ** lm) ** 4 / 4 for lm in lms]
    lgLV = np.log10(max(Vv)) / 4
    val['log10_Lambda_V'] = float(lgLV)
    val['Lambda_I_over_Lambda_V'] = float(10 ** zI / 10 ** lgLV)
    val['Buttazzo_scales'] = 'log10 Lambda_V = 9.5; Lambda_I ~ 13 Lambda_V (Landau gauge); Lambda_lambda ~ 2 Lambda_V'
    val['note_on_quoted_171.53'] = ('Buttazzo eq. (stability) M_h > 129.6 + 2.0 (M_t - 173.34) - 0.5 (as-0.1184)/0.0007 '
                                    'gives M_t < 171.12 at M_h = 125.15, not the 171.53 quoted in the next equation; '
                                    'this code reproduces 171.07-171.13 (and Bednyakov 171.44 at their inputs to 0.06)')
    # inverse relation coefficients are computed in the relation block
    val['Mt_crit_at_ref'] = {}
    for crit in ('msbar', 'eff'):
        Mt, mus, inter = Mt_crit(125.15, 0.1184, crit)
        val['Mt_crit_at_ref'][crit] = {'Mt_crit': Mt, 'mu_star': mus, 'mu_star_interior': inter}
    val['Mt_crit_at_ref']['Buttazzo_quoted'] = '171.53 +- 0.15 (th) (+-0.23 as, +-0.15 Mh); eff. potential'
    v_b = Mt_crit(125.7, 0.1185, 'msbar')
    val['Mt_crit_at_Bednyakov_inputs(125.7,0.1185)'] = {'mine_msbar': v_b[0], 'mu_star': v_b[1],
                                                        'Bednyakov_quoted': '171.44 +-0.30 (+0.17 -0.36) th; '
                                                        'eff-type criterion 171.64'}
    out['validation'] = val

    # ---------------- stability bound / criticality relation near the PDG 2025 point
    Mh0, sMh, _ = PDG['Mh']
    a0, sa, _ = PDG['alpha_s']
    MW0 = PDG['MW'][0]
    rel = {}
    for crit in ('msbar', 'eff'):
        Mt0, mus0, inter0 = Mt_crit(Mh0, a0, crit, MW=MW0)
        dMh = (Mt_crit(Mh0 + 0.5, a0, crit, MW=MW0)[0] - Mt_crit(Mh0 - 0.5, a0, crit, MW=MW0)[0]) / 1.0
        da = (Mt_crit(Mh0, a0 + 0.0009, crit, MW=MW0)[0] - Mt_crit(Mh0, a0 - 0.0009, crit, MW=MW0)[0]) / 2.0
        rel[crit] = {'Mt_crit_GeV': Mt0, 'mu_star_GeV': mus0, 'mu_star_below_cut': inter0,
                     'dMtcrit_dMh': dMh, 'dMtcrit_per_0.0009_as': da,
                     'fit': f'Mt_crit = {Mt0:.3f} + {dMh:.3f} (Mh-125.20) + {da:.3f} (as-0.1180)/0.0009  GeV'}
        # inverse relation Mh_crit(Mt)
        rel[crit]['Mh_crit_at_Mt_direct'] = Mh_crit(PDG['Mt_direct'][0], a0, crit, MW=MW0)
        rel[crit]['Mh_crit_at_Mt_xsec'] = Mh_crit(PDG['Mt_xsec_pole'][0], a0, crit, MW=MW0)
        # reduced-Planck cut
        rel[crit]['Mt_crit_cut_reduced_MPl'] = Mt_crit(Mh0, a0, crit, mu_cut=MPL_RED, MW=MW0)[0]
    # theory systematics on Mt_crit (msbar)
    base = rel['msbar']['Mt_crit_GeV']
    sy = {}
    sy['dlam_+0.0003'] = Mt_crit(Mh0, a0, 'msbar', MW=MW0, dlam=TH['dlam'])[0] - base
    sy['dyt_+0.0005'] = Mt_crit(Mh0, a0, 'msbar', MW=MW0, dyt=TH['dyt'])[0] - base
    sy['RGE_2loop_minus_3loop'] = Mt_crit(Mh0, a0, 'msbar', MW=MW0, loops=2)[0] - base
    sy['no_yb_ytau'] = Mt_crit(Mh0, a0, 'msbar', MW=MW0, yb=0.0, ytau=0.0)[0] - base
    sy['MW_+1sigma'] = Mt_crit(Mh0, a0, 'msbar', MW=MW0 + 0.0133)[0] - base
    sy['criterion_eff_minus_msbar'] = rel['eff']['Mt_crit_GeV'] - base
    sy = {k: float(v) for k, v in sy.items()}
    s_match = float(np.hypot(sy['dlam_+0.0003'], sy['dyt_+0.0005']))
    s_crit = abs(sy['criterion_eff_minus_msbar']) / 2.0
    rel['theory_shifts_GeV'] = sy
    rel['sigma_th_matching_GeV'] = s_match
    rel['sigma_th_criterion_GeV'] = s_crit
    out['relation'] = rel

    # ---------------- Froggatt-Nielsen double condition at fixed mu*
    fn = {}
    for name, mus in (('MPl_1.22e19', MPL), ('reduced_2.435e18', MPL_RED), ('1e18', 1e18)):
        x, res = fn_double(a0, mus, MW0)
        fn[name] = {'Mt': x[0], 'Mh': x[1], 'residuals(100lam,1e3beta)': res}
    fn['FN1995_quoted'] = 'Mt = 173 +- 5, Mh = 135 +- 9 GeV (2-loop, mu* ~ 1e19 GeV), hep-ph/9511371'
    fn['caveat'] = ('matching fits of 1307.3536 are linear around Mh=125.15; Mh = 128.6 (132.0) is an extrapolation '
                    'by 3.5 (6.9) GeV. Tree-level curvature d2lambda/dMh2 ~ 2 lambda/Mh^2 ~ 1.6e-5 GeV^-2 gives '
                    'lambda errors ~1e-4 (4e-4), i.e. ~0.05 (0.2) GeV in the predicted Mh [E]')
    out['FN_double_criticality'] = fn

    # ---------------- retrospective tension (PDG 2025) and held-out projection
    ten = {}
    for crit in ('msbar', 'eff'):
        R = rel[crit]
        s_par = np.hypot(R['dMtcrit_dMh'] * sMh, R['dMtcrit_per_0.0009_as'] * (sa / 0.0009))
        for lab, (Mt, sMt, _) in (('direct_MC', PDG['Mt_direct']), ('xsec_pole', PDG['Mt_xsec_pole'])):
            for budget, extra in (('exp+param+matching', 0.0), ('+criterion', s_crit),
                                  ('+criterion+MCpole0.3', np.hypot(s_crit, 0.3))):
                if lab == 'xsec_pole' and 'MCpole' in budget:
                    continue
                stot = float(np.sqrt(sMt ** 2 + s_par ** 2 + s_match ** 2 + extra ** 2))
                ten[f'{crit}|{lab}|{budget}'] = {'Delta_GeV': float(Mt - R['Mt_crit_GeV']), 'sigma_GeV': stot,
                                                 'z': float((Mt - R['Mt_crit_GeV']) / stot)}
    out['tension_retrospective'] = ten
    # future
    R = rel['msbar']
    s_par_f = np.hypot(R['dMtcrit_dMh'] * FUTURE['Mh'][0], R['dMtcrit_per_0.0009_as'] * FUTURE['alpha_s'][0] / 0.0009)
    fut = {}
    for lab, sm in (('matching errors as today', s_match), ('matching errors halved', s_match / 2)):
        st = float(np.sqrt(FUTURE['Mt'][0] ** 2 + s_par_f ** 2 + sm ** 2 + s_crit ** 2))
        fut[lab] = {'sigma_GeV': st,
                    'z_if_central_values_stay(direct)': float((PDG['Mt_direct'][0] - R['Mt_crit_GeV']) / st),
                    'z_if_central_values_stay(xsec)': float((PDG['Mt_xsec_pole'][0] - R['Mt_crit_GeV']) / st)}
    # intermediate scenario: HL-LHC top mass ~0.2 GeV [R, order of magnitude] + lattice alpha_s halved (FCC FSR)
    s_par_h = np.hypot(R['dMtcrit_dMh'] * 0.02, R['dMtcrit_per_0.0009_as'] * 0.00045 / 0.0009)
    st = float(np.sqrt(0.2 ** 2 + s_par_h ** 2 + s_match ** 2 + s_crit ** 2))
    fut['HL-LHC Mt 0.2 GeV [R] + lattice as 0.00045 + Mh 0.02'] = {
        'sigma_GeV': st, 'z_if_central_values_stay(direct)': float((PDG['Mt_direct'][0] - R['Mt_crit_GeV']) / st)}
    fut['note'] = ('held out = data that do not yet exist; the z values assume unchanged central values and are '
                   'illustrations of reach, not predictions of outcome')
    out['held_out_projection'] = fut
    from scipy.stats import norm
    for k, v in out['tension_retrospective'].items():
        v['one_sided_p'] = float(norm.sf(v['z']))

    # ---------------- numbers for the vacuum readings (sec. 3 of NOTES)
    r = run_to(PDG['Mt_direct'][0], Mh0, a0, MW0)
    sol = r['sol'].sol
    nums = {}
    for mu in (1e10, 1e12, 1e15, 1e17, MPL_RED, MPL):
        Y = sol(np.log(mu / r['Mt']))
        le = S.lambda_eff(Y)
        nums[f'{mu:.3e}'] = {'lambda': float(Y[6]), 'lambda_eff': float(le),
                             'V_eff_GeV4': float(np.exp(4 * Y[7]) * le * mu ** 4 / 4)}
    lamr = r['Y'][6]
    j = np.where(np.diff(np.sign(lamr)))[0]
    nums['Lambda_lambda_GeV'] = float(r['mu'][j[0]]) if len(j) else None
    out['potential_at_PDG_central'] = nums

    # ---------------- gravitational (CDL) quench estimate near criticality [E]
    out['CDL_quench_estimate'] = cdl_quench(Mh0, a0, MW0, rel['msbar']['Mt_crit_GeV'])
    out['runtime_s'] = time.time() - t0
    with open(OUT, 'w') as f:
        json.dump(out, f, indent=1, default=float)
    return out


def cdl_quench(Mh, a, MW, Mtc, M=2.435e18):
    """Thin-wall estimate [E]: a Minkowski false vacuum decaying to AdS is stable if eps <= (3/4) sigma^2/M^2
    (Coleman-De Luccia 1980 [R]). Potential: RG-improved tree level V = e^{4Gamma} lambda(h) h^4/4 (h >> v).
    For Mt slightly above Mt_crit: depth eps of the deepest point (h <= M_Pl), wall tension
    sigma = int sqrt(2V) dh over the barrier region (V > 0, h < h_min). Reports eps / threshold."""
    rows = []
    for d in (0.001, 0.01, 0.05, 0.1, 0.4, 1.6):
        Mt = Mtc + d
        r = run_to(Mt, Mh, a, MW)
        sol = r['sol'].sol
        hs = np.exp(np.linspace(np.log(1e3), np.log(MPL), 4000))
        Y = np.array([sol(np.log(h / Mt)) for h in hs]).T
        V = np.exp(4 * Y[7]) * Y[6] * hs ** 4 / 4
        i = int(np.argmin(V))
        eps = float(-V[i])
        mask = (V > 0) & (hs < hs[i])
        sigma = float(np.trapezoid(np.sqrt(2 * V[mask]), hs[mask]))
        thr = 0.75 * sigma ** 2 / M ** 2
        rows.append({'dMt_above_crit_GeV': d, 'h_min_GeV': float(hs[i]), 'eps_GeV4': eps, 'sigma_GeV3': sigma,
                     'threshold_GeV4': float(thr), 'eps_over_threshold': float(eps / thr)})
    return {'scan': rows,
            'conclusion': 'eps exceeds the CDL quench threshold already 0.001 GeV above Mt_crit, so gravitational '
                          'quenching moves the stability boundary by < 1e-3 GeV (negligible)',
            'status': '[E] thin wall, RG-improved tree potential, Planckian field values, Einstein gravity'}


if __name__ == '__main__':
    o = main()
    print(json.dumps({k: o[k] for k in ('validation', 'relation', 'FN_double_criticality',
                                         'tension_retrospective', 'held_out_projection',
                                         'potential_at_PDG_central', 'CDL_quench_estimate', 'runtime_s')},
                     indent=1, default=float))
