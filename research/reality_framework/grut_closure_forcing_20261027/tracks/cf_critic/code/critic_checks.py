#!/usr/bin/env python3
"""cf_critic quick checks (Claude, synthesis critic, 27 October 2026). Own code; reads sibling results read-only.

A. Near-criticality line: reconcile cf_vacuum (P3a, SMDR line), cf_scale (P4-SM, Buttazzo eq. 61e fit) and the
   q-vacuum HANDOFF figure on one error budget.
B. Final-vacuum sign required by each capacity/vacuum strengthening (analytic bookkeeping, pairwise compatibility).
C. The P3a Sigma corner after reheating (new): thermal response of a light, portal-coupled Sigma to the visible bath
   (T_RH from hw_gates gate 3). + portal: hidden U(1) restoration temperature, condensate tracking, Kibble Z6 strings
   and their G mu (gate-3 law for type II; gp_stability Nielsen-Olesen solver for type I). - portal: thermal push of
   |Sigma| and the Higgs thermal-mass shift.
D. Portal gap: P3a loop stabiliser vs invisible width vs the dark-matter chain's detector-level portal.
E. P4b band vs the P3a hidden-sector bound and the Higgs-threshold window.
F. TCC at the heavy window (one-line reproduction).
Usage: python3 -B critic_checks.py <out.json>
Labels: [C] computed here, [D] derived here, [Q] quoted from the sibling file named, [E] estimate.
"""
import json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TRK = HERE.parents[1]                                   # .../tracks
FR = TRK.parents[1] / 'grut_formation_review_20261025' / 'tracks'
sys.path.insert(0, str(FR / 'gp_stability' / 'code'))

M_PL = 1.220890e19          # GeV, G = 1/M_PL^2 (gate-3 convention)
M_RED = 2.435e18
H_I = 1.3335e11
h = 40 * H_I                # Sigma vev scale, |Sigma|_vac = h/sqrt2
T_RH = 3.08e14              # gate 3 visible reheating temperature [Q hw_gates]
E_GAUGE = 1e-7
GSTAR = 106.75


def part_A():
    # PDG 2026 inputs as used by both tracks [Q cf_vacuum, cf_scale]
    mt, smt = 172.60, 0.27
    mh, smh = 125.13, 0.11
    als, sals = 0.1180, 0.0009
    lam = lambda mt, mh, a: -0.0143 - 0.0066 * (mt - 173.34) + 0.0018 * (a - 0.1184) / 0.0007 + 0.0029 * (mh - 125.15)
    l0 = lam(mt, mh, als)
    mt_crit_B = mt - l0 / (-0.0066)             # lambda(M_Pl) = 0 (Buttazzo fit, equality form)
    # SMDR line (cf_vacuum frozen_relation) [Q]
    smdr = lambda mh, a: 170.997 + 0.494 * (mh - 125.13) + 0.355 * (a - 0.1180) / 0.0009
    mt_crit_S = smdr(mh, als)
    out = dict(lambda_MPl_Buttazzo_fit=l0, mt_crit_Buttazzo_fit=mt_crit_B, mt_crit_SMDR=mt_crit_S,
               difference_GeV=mt_crit_B - mt_crit_S)
    budgets = {}
    for label, mc in (('no MC-pole term', 0.0), ('0.3 GeV MC-pole', 0.3), ('0.6 GeV MC-pole', 0.6)):
        for th_label, th in (('SMDR theory 0.113+0.032', np.hypot(0.113, 0.032)),
                             ('Buttazzo-fit theory 0.0003/0.0066', 0.0003 / 0.0066),
                             ('Bednyakov upward 0.36', 0.36)):
            s = np.sqrt(smt**2 + mc**2 + (0.494 * smh)**2 + (0.355 * sals / 0.0009)**2 + th**2)
            budgets[f'{label}; {th_label}'] = dict(sigma_GeV=s, z_SMDR=(mt - mt_crit_S) / s, z_Buttazzo=(mt - mt_crit_B) / s)
    out['z_by_budget'] = budgets
    out['xsec_pole_172.1+-0.6'] = dict(z_SMDR=(172.1 - mt_crit_S) / np.sqrt(0.6**2 + 0.054**2 + 0.355**2 + 0.117**2))
    out['EW_subfield_alpha_s_0.1221+-0.0027'] = dict(
        mt_crit_SMDR=smdr(mh, 0.1221),
        z=(mt - smdr(mh, 0.1221)) / np.sqrt(smt**2 + (0.494 * smh)**2 + (0.355 * 0.0027 / 0.0009)**2 + 0.117**2))
    out['reading'] = ('Same line to 0.01 GeV. cf_vacuum 3.45 sigma (no MC-pole term) and cf_scale 2.95 sigma (0.3 GeV '
                      'MC-pole term included) differ only by error budget; HANDOFF 3.1 sits between. P3a is the one-sided '
                      'form (M_t <= line); P4-SM/MPP are the equality form; the equality comparator is A_lambda-dependent '
                      '(cf_scale skeptic: lambda_crit(M_P) = +0.017 at A = 0.0237 with AF g2, g3).')
    return out


def part_B():
    rows = {
        'P2-BF (Banks-Fischler capacity)': 'Lambda_inf > 0, eternal de Sitter, finite S_Lambda',
        'P_q (q-theory final-state equilibrium, cf_vacuum P3a)': 'Lambda_inf = 0 (Minkowski end state); today a non-equilibrium residual',
        'P2-TCC (late-time)': 'no eternal dS: acceleration ends within ~2.4e12 yr; Lambda_inf <= 0 or a rolling field',
        'P1-time (finite four-volume)': 'Lambda_inf < 0 / recollapse',
    }
    keys = list(rows)
    sign = {'P2-BF (Banks-Fischler capacity)': {'+'},
            'P_q (q-theory final-state equilibrium, cf_vacuum P3a)': {'0'},
            'P2-TCC (late-time)': {'0', '-', 'rolling'},
            'P1-time (finite four-volume)': {'-'}}
    compat = {}
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            compat[f'{a} & {b}'] = 'compatible' if sign[a] & sign[b] else 'INCOMPATIBLE (opposite final-vacuum sign)'
    return dict(final_vacuum_required=rows, pairwise=compat,
                candidate1_premise=('Candidate 1 H0 exclusion of the 100 GeV r=1 preparation needs total vacuum >= -(2.28 meV)^4 '
                                    '[Q cf_space 3.7]. Supplied [A] by BF (Lambda_inf>0 plus monotone roll), by P_q only if the '
                                    'present residual is positive (calibrated, i.e. supplied), and reversed by finite four-volume '
                                    '(Luu-Qiu-Tye best fit -6.0e-47 GeV^4 < -(2.28 meV)^4 = -2.7e-47 GeV^4 [C here]).'),
                minus_2p28meV_4_GeV4=-(2.28e-12)**4)


def eps_type2(beta):
    return 0.5 * np.log(beta) + 0.15024     # gate-3 large-beta law [Q hw_gates gate 3]


def part_C(do_solver=True):
    th = json.load(open(TRK / 'cf_vacuum' / 'results' / 'threshold.json'))['C2_GRUT_Sigma_fixed_vev']
    Hub = lambda T: 1.66 * np.sqrt(GSTAR) * T**2 / M_PL
    # Higgs thermal-mass coefficient at high scale [E]: c_H = 3g^2/16 + g'^2/16 + y_t^2/4 + lam/2, couplings at ~1e12-1e14 GeV
    g2, gp, yt, lH = 0.52, 0.45, 0.45, 0.0
    c_H = 3 * g2**2 / 16 + gp**2 / 16 + yt**2 / 4 + lH / 2
    solver_eps = {}
    if do_solver:
        import z6_strings as Z
        for beta in (10.0, 1.0, 1e-2, 1e-4, 1e-6, 1e-8):
            t0 = time.time()
            try:
                xmax = 40.0 / np.sqrt(beta)
                eps, parts, nodes, _ = Z.solve_string(beta, xmin=1e-5, xmax=xmax, n=6000)
                solver_eps[beta] = dict(eps=float(eps), nodes=int(nodes), seconds=time.time() - t0)
            except Exception as ex:      # record, do not hide
                solver_eps[beta] = dict(error=str(ex)[:200])
            print('solver', beta, solver_eps[beta], flush=True)

    def eps_of(beta):
        if beta >= 1e2:
            return eps_type2(beta), 'gate-3 large-beta law'
        # interpolate in ln beta over solver points (beta = 1e-4 ... 1) plus gate-3 anchors; extrapolate linearly in
        # 1/ln(1/beta) below 1e-4 [E]
        pts = sorted((b, d['eps']) for b, d in solver_eps.items() if 'eps' in d)
        pts += [(1e2, eps_type2(1e2) + (0.2378 - 0.15024))]    # gate-3 tabulated eps at 1e2 [Q]
        pts = sorted(pts)
        if beta >= pts[0][0]:
            lb = np.log([p[0] for p in pts]); ev = [p[1] for p in pts]
            return float(np.interp(np.log(beta), lb, ev)), 'solver interpolation'
        b0, e0 = pts[0]
        return float(e0 * np.log(1 / b0) / np.log(1 / beta)) if beta < 1 else e0, '[E] extrapolation eps ~ 1/ln(1/beta)'

    rows = []
    for Mkey, v in th.items():
        if not isinstance(v, dict) or v.get('sign=+1') is None:
            continue
        lamS = v['lam_Sigma']
        Ms = float(Mkey)
        wp, wm = v['sign=+1'], v['sign=-1']
        if not isinstance(wp, dict) or wp.get('delta_min') is None:
            rows.append(dict(M_Sigma_GeV=Ms, lam_Sigma=lamS, window='none'))
            continue
        beta = 2 * lamS / (36 * E_GAUGE**2)
        eps, how = eps_of(beta)
        Gmu = np.pi * h**2 * eps / M_PL**2
        plus = []
        for d in (wp['delta_min'], wp['delta_max']):
            lHS = 2 * np.sqrt(d * lamS)
            Tc = h * (9 * lamS / d) ** 0.25                    # lamHS T^2/6 = lamS h^2
            mT_over_H = np.sqrt(lHS / 6) * Tc / Hub(Tc)
            # [E] HH -> Sigma Sigma* rate n<sigma v> ~ 5e-3 lamHS^2 T; Gamma/H = 1 at T_th (Sigma quanta thermalise)
            T_th = 5e-3 * lHS**2 * M_PL / (1.66 * np.sqrt(GSTAR))
            plus.append(dict(delta=d, lam_HSigma=lHS, T_c_GeV=Tc, T_c_over_h=Tc / h, T_RH_over_T_c=T_RH / Tc,
                             thermal_mass_over_H_at_T_c=mT_over_H, lamHS_over_lamS=lHS / lamS,
                             radiation_fraction_at_T_c=(Tc / T_RH)**4,
                             T_Sigma_quanta_thermalise_E_GeV=T_th, T_th_over_T_c=T_th / Tc,
                             coherent_amplitude_at_T_c_over_vev=Tc / T_RH,
                             vacuum_fluct_over_coherent_at_T_c=(np.sqrt(lHS / 6) * Tc / (2 * np.pi)) / ((h / np.sqrt(2)) * Tc / T_RH)))
        minus = []
        for d in (wm['delta_min'], wm['delta_max']):
            lHS = -2 * np.sqrt(d * lamS)
            # |Sigma|^2(T) = h^2/2 + |lamHS| T^2/(12 lamS) at T >> vev scale (tree + Higgs thermal mass) [D]
            s2 = h**2 / 2 + abs(lHS) * T_RH**2 / (12 * lamS)
            # Higgs effective m^2/T^2 along the thermal valley: c_H + lamHS |Sigma|^2/T^2 -> c_H - delta/3 [D]
            T_star = (h / np.sqrt(2)) * np.sqrt(6) * (lamS / d) ** 0.25      # thermal minimum = vev
            minus.append(dict(delta=d, lam_HSigma=lHS, Sigma_over_vev_at_T_RH=np.sqrt(s2) / (h / np.sqrt(2)),
                              T_star_GeV=T_star, thermal_mass_over_H_at_T_RH=np.sqrt(abs(lHS) / 6) * T_RH / Hub(T_RH),
                              Higgs_m2_over_T2_at_T_RH=c_H + lHS * s2 / T_RH**2, c_H_E=c_H))
        rows.append(dict(M_Sigma_GeV=Ms, lam_Sigma=lamS, beta=beta, string_type='I' if beta < 1 else 'II',
                         eps=eps, eps_method=how, Gmu=Gmu, plus_portal=plus, minus_portal=minus))
    Gs = [r['Gmu'] for r in rows if 'Gmu' in r]
    limits = {'LVK O3 model B (LRS)': 4.0e-15, 'LVK O3 model C-1': 2.1e-15, 'LVK O3 model A (BOS)': 9.6e-9,
              'NANOGrav 15-yr most stringent (log10 < -10.10)': 10**-10.10, 'reference (lam_Sigma=1) gate 3': 8.89e-12}
    return dict(solver_eps=solver_eps, rows=rows, Gmu_range=[min(Gs), max(Gs)], limits_Q_hw_gates_gate3=limits,
                Gmu_over_LVK_B=[min(Gs) / 4.0e-15, max(Gs) / 4.0e-15], Gmu_over_NANOGrav=[min(Gs) / 10**-10.10, max(Gs) / 10**-10.10],
                T_RH_GeV=T_RH, h_GeV=h, c_H_estimate=c_H)


def part_D():
    need_b, need_lHX = 0.117, 0.175          # [Q cf_vacuum D_light_singlets]
    inv_b, inv_lHX = 5.05e-3, 7.14e-3        # [Q rp_portals]
    det_b = (1.0e-10, 1.5e-10)               # detector-level target, f_rest 37-87 at onset [Q REVIEW_REHEATING]
    return dict(need_over_invisible=[need_b / inv_b, need_lHX / inv_lHX],
                need_over_detector_level=[need_b / det_b[1], need_b / det_b[0]],
                note='P3a loop stabilisation by phi needs b >= 0.117; the DM chain already strains at b ~ 1e-10 '
                     '(>= 97-99% erasure of onset pairs uncomputed). Gap ~1e9.')


def part_E():
    lo, hi = 0.241, 0.566                    # [Q cf_scale P4b]
    return dict(min_lamS_lamY_in_band=lo * lo, P3a_hidden_bound=9.05e-4, margin=lo * lo / 9.05e-4,
                M_Sigma_band_GeV=[np.sqrt(2 * lo) * h, np.sqrt(2 * hi) * h], Lambda_lambda_GeV=4.94e10,
                M_Sigma_over_Lambda_lambda=[np.sqrt(2 * lo) * h / 4.94e10, np.sqrt(2 * hi) * h / 4.94e10],
                reading='P4b repair keeps the hidden-sector P3a pass but leaves no Higgs escape; the Sigma corner '
                        '(lam_Sigma <~ 3e-5) is a different branch, outside the common-lambda band.')


def part_F():
    Nmax = np.log(M_RED / H_I)
    return dict(N_max_TCC=Nmax, N_min_heavy=58.4, deficit=58.4 - Nmax)


def main(out):
    t0 = time.time()
    res = dict(track='cf_critic', date='2026-10-27', method=__doc__)
    res['A_near_criticality'] = part_A()
    res['B_final_vacuum_sign'] = part_B()
    res['C_sigma_corner_after_reheating'] = part_C()
    res['D_portal_gap'] = part_D()
    res['E_P4b_vs_P3a'] = part_E()
    res['F_TCC'] = part_F()
    res['runtime_s'] = time.time() - t0
    Path(out).write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: v for k, v in res.items() if k not in ('method',)}, indent=1, default=float)[:12000])


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else str(HERE.parent / 'results' / 'critic_checks.json'))
