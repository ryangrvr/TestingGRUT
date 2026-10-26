#!/usr/bin/env python3
"""kin_thermo (Claude, own code): the verification checks of this track, written to results/checks.json.

Each check records residuals separately from its pass/fail threshold.  'reproduction' items run our earlier
track code (ls_threshold/lsthermo.py) for comparison and are labelled so; everything else is this track's code.
Usage: python3 -B checks.py
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

import kt_common as C
import thermo as TH
import expansion as EX
import kernels as KN
import formation as FM

RES = Path(__file__).resolve().parents[1] / 'results'
MX = C.MX
CODEX_RESERVOIR = {0.3: (0.023049850882506346, 0.021215613259849126), 1.0: (0.45073872791952774, 0.23584943677240486),
                   3.0: (0.8455630737220086, 0.31894780956024205), 12.0: (0.9810698062732441, 0.332365630563511),
                   30.0: (0.9959283980621779, 0.3331773016593269)}
CODEX_STATES = {False: dict(T_over_m=8.11769060417811, condensed=0.004404564224805499, pinv=0.04451797806645164,
                            n_plus=0.05900792776420369, n_minus=0.04451797806645164, rho_plus=0.11714516339538246,
                            P_plus=0.03889436087140684),
                True: dict(T_over_m=7.343904194606052, condensed=0.007072205615206939, pinv=0.06850712020374237,
                           n_plus=0.044288446795821744, n_minus=0.03246613848847116, n_s=0.03604098171527122,
                           rho_s=0.06882200139669796, P_s=0.022454589776708417, scalar_fraction_total=0.327647213576238)}


def check(name, residual, threshold, **extra):
    r = dict(name=name, residual=float(residual), threshold=threshold, passed=bool(residual <= threshold), **extra)
    print('%-58s %s  res=%.3e (thr %.0e)' % (name, 'PASS' if r['passed'] else 'FAIL', residual, threshold), flush=True)
    return r


def main():
    t0 = time.time()
    out = []
    # 1. series vs direct quadrature (independent representation), seven quantities, 28 states
    worst = 0.0
    for t in (0.3, 1.0, 3.0, 12.0):
        for r in (1.0, TH.RS):
            for u in (-3.0, 0.0, 0.5, 0.99):
                a = TH.species(t, r, u); b = TH.quad_species(t, r, u, derivs=True)
                worst = max(worst, max(abs(a[k] / b[k] - 1) for k in ('n', 'P', 'rho', 'dn_du', 'drho_du', 'dn_dt', 'drho_dt')))
    out.append(check('Bessel series vs direct quadrature (7 quantities, 28 states)', worst, 1e-9))
    # 2. Gibbs-Duhem / Maxwell identities by finite differences
    t, u, h = 2.0, 0.4, 1e-5
    a = TH.species(t, 1.0, u)
    dP_du = (TH.species(t, 1.0, u + h)['P'] - TH.species(t, 1.0, u - h)['P']) / (2 * h)
    dP_dt = (TH.species(t + h, 1.0, u)['P'] - TH.species(t - h, 1.0, u)['P']) / (2 * h)
    dn_dt = (TH.species(t + h, 1.0, u)['n'] - TH.species(t - h, 1.0, u)['n']) / (2 * h)
    out.append(check('dP/dmu = n, dP/dT = s, analytic dn/dT (finite differences)',
                     max(abs(dP_du / a['n'] - 1), abs(dP_dt / a['s'] - 1), abs(dn_dt / a['dn_dt'] - 1)), 1e-8))
    # 3. Codex reservoir control rows
    worst = 0.0; rows = []
    for tt, (ns_nm, frac) in CODEX_RESERVOIR.items():
        g = TH.gas(tt, 0.0, 0.0, 0.0)
        rows.append(dict(T_over_m=tt, n_s_over_nminus=g['n_s'] / g['n_minus'], scalar_fraction=g['scalar_energy_fraction'],
                         codex=(ns_nm, frac)))
        worst = max(worst, abs(g['n_s'] / g['n_minus'] / ns_nm - 1), abs(g['scalar_energy_fraction'] / frac - 1))
    out.append(check('Codex THERMODYNAMIC_CONTROL reservoir rows (5 T, two ratios)', worst, 1e-9, rows=rows))
    # 4. Codex matched-energy states at R = 40.5626
    rho = C.rho_onset(40.56262288525937) / MX**4; nq = C.N_Q0 / MX**3
    worst = 0.0; st = {}
    for with_s in (False, True):
        g = TH.state_equilibrium(nq, rho, C.OMEGA_REF, with_s=with_s)
        cx = CODEX_STATES[with_s]
        mine = dict(T_over_m=g['t'], condensed=g['n_cond'] * MX**3, pinv=g['pair_inventory'] * MX**3,
                    n_plus=g['n_plus'] * MX**3, n_minus=g['n_minus'] * MX**3, rho_plus=g['X']['rho'] * MX**4,
                    P_plus=g['X']['P'] * MX**4)
        if with_s:
            mine.update(n_s=g['n_s'] * MX**3, rho_s=g['S']['rho'] * MX**4, P_s=g['S']['P'] * MX**4,
                        scalar_fraction_total=g['S']['rho'] / g['rho_total'], scalar_fraction_gas=g['scalar_energy_fraction'])
        worst = max(worst, max(abs(mine[k] / cx[k] - 1) for k in cx))
        st[str(with_s)] = dict(mine=mine, codex=cx)
    out.append(check('Codex matched-energy states (T/m 8.118 -> 7.344, condensed, species)', worst, 1e-8, states=st))
    # 5. two-potential solver reduces to the equilibrium state when given the equilibrium inventory
    g1 = TH.state_equilibrium(nq, rho, C.OMEGA_REF, with_s=True)
    g2 = TH.state_two_potential(nq, g1['pair_inventory'], rho, C.OMEGA_REF, with_s=True, tguess=g1['t'])
    out.append(check('two-potential solver at the equilibrium inventory = equilibrium solver',
                     max(abs(g2['t'] / g1['t'] - 1), abs(g2['um'] + g1['up']), abs(g2['n_cond'] / g1['n_cond'] - 1)), 1e-9,
                     t=(g1['t'], g2['t']), um=(g2['um'], -g1['up'])))
    # 6. three-potential Newton solver round trip
    g0 = TH.gas(4.0, 0.2, -0.5, 0.3, fast_s=False)
    g3 = TH.state_three_potential(g0['n_Q'], g0['n_minus'], g0['n_s'], g0['rho'], t0=3.0)
    out.append(check('three-potential Newton round trip (n_Q, n_-, n_s, rho)',
                     max(abs(g3['t'] / 4 - 1), abs(g3['up'] - 0.2), abs(g3['um'] + 0.5), abs(g3['us'] - 0.3)), 1e-9,
                     iterations=g3['newton_iterations']))
    # 7. pair-inventory susceptibility: Bose exact vs Codex's MB fixed-T formula 1 + K(n_+ + n_-)
    rows = []
    for tt, uu in ((0.3, 0.0), (1.0, 0.0), (3.0, 0.3), (7.34, 0.86), (12.0, 0.0)):
        S = TH.susceptibility(tt, uu, -uu, fast_s=True)
        g = S['state']
        Kb = TH.K_bose(tt); Km = TH.K_mb(tt)
        rows.append(dict(T_over_m=tt, u=uu, bose_exact=S['dPinv_dnminus_fixed_t_nQ'],
                         MB_formula_with_K_MB=1 + Km * (g['n_plus'] + g['n_minus']),
                         MB_formula_with_K_Bose=1 + Kb * (g['n_plus'] + g['n_minus']), K_bose=Kb, K_mb=Km))
    dev = max(abs(r['bose_exact'] / r['MB_formula_with_K_MB'] - 1) for r in rows)
    out.append(check('dP_inv/dn_- |T,nQ: Bose exact vs MB 1+K(n_+ + n_-) (reported, no threshold)', dev, np.inf, rows=rows))
    # 8. isentropic expansion with the Jacobian conserves comoving entropy and charge
    iso = EX.check_isentropic()
    out.append(check('isentropic expansion via susceptibility Jacobian: s a^3 and n_Q a^3 drift',
                     max(iso['max_rel_entropy_drift'], iso['max_rel_charge_drift']), 1e-9))
    # 9. parent ledger
    for key in ('e1e-7_band', 'e1e-6_band'):
        pl = EX.check_parent_ledger(key)
        out.append(check('parent decay ledger %s: daughters and R vs closed forms' % key, max(abs(pl['rel_err_number']), abs(pl['rel_err_R'])), 1e-10, **pl))
    # 10. objects: reference values, dE/dQ identity, virial
    ob = json.loads((RES / 'objects_EQ.json').read_text())
    ref = ob['reference_Q200000']
    out.append(check('Q=200000: omega/m vs archived 0.8611951308217225', abs(ref['Omega'] / C.OMEGA_REF - 1), 1e-10, mine=ref['Omega']))
    out.append(check('Q=200000: E/(Q m) vs archived 0.9233161353726057', abs(ref['E_over_QmX'] / C.E_OVER_QM_REF - 1), 1e-10, mine=ref['E_over_QmX']))
    out.append(check('Q=200000: E - Q omega vs Codex 1148.1075313763777 GeV', abs(ref['E_minus_Q_omega_GeV'] / 1148.1075313763777 - 1), 1e-8, mine=ref['E_minus_Q_omega_GeV']))
    fd = ob['dEdQ_finite_difference_at_reference'][-1]
    out.append(check('dE/dQ = Omega at the reference (finite difference, dOmega=5e-5)', abs(fd['dEdQ_over_mX'] / ref['Omega'] - 1), 1e-6, **fd))
    out.append(check('max |virial| (Derrick identity) over the table', max(abs(r['virial']) for r in ob['rows']), 1e-9))
    if 'dEdQ_identity_check' in ob:
        out.append(check('dE/dQ = Omega along the table (splines, away from Q_min)', ob['dEdQ_identity_check']['max_abs_dEdQ_over_mX_minus_Omega_away_from_Qmin'], 2e-3))
    wk = EX.check_object_work()
    out.append(check('object-number work identity along the table (path integral vs exact dE)', max(abs(w['rel_err']) for w in wk), 1e-3,
                     linear_comparator_errors=[w['linear_rel_err'] for w in wk]))
    # 11. conversion kernel: table reproduces the two bb211d1 spot values; detailed balance of forward/reverse
    out.append(check('conversion table: cold sigma v vs Codex 1.13039e-16 cm3/s', abs(KN.sv_forward(1e-5) * KN.GEV2_TO_CM3S / 1.13039e-16 - 1), 1e-4))
    out.append(check('conversion table: daughter sigma v (p=0.7071) vs Codex 8.1225e-15 cm3/s', abs(KN.sv_forward(0.7071067811865476) * KN.GEV2_TO_CM3S / 8.1225e-15 - 1), 1e-3))
    db = [KN.detailed_balance_residual(tt, uu) for tt, uu in ((0.3, -1.0), (1.0, -0.8), (3.0, -0.75), (0.1, -0.9))]
    out.append(check('conversion forward = reverse at mu_+ - mu_- = 2 omega (4 states, Bose factors)', max(abs(d['relative_residual']) for d in db), 1e-7, states=db))
    # 12. surface exchange vanishes at mu_+ = omega
    se = KN.surface_exchange(1.0, C.OMEGA_REF)
    out.append(check('surface exchange zero at mu_+ = omega', abs(se['dQ_dt_GeV']) / se['absorbed_number_rate_GeV'], 1e-12))
    # 13. formation: U_eff continuity and dispersion relation sanity (alpha_max <= k_max^2/(2 omega) small-k bound)
    r = FM.rule_fragmentation(absX0_GeV=C.XC)
    out.append(check('U_eff continuity at x_c', abs(FM.U_eff(C.XC * (1 - 1e-9)) / FM.U_eff(C.XC * (1 + 1e-9)) - 1), 1e-6))
    out.append(check('fragmentation: k at alpha_max inside the band (0 < k < k_max)', 0.0 if 0 < r['k_at_alpha_max_over_mX'] < r['k_max_over_mX'] else 1.0, 0.5,
                     k_max_over_m=r['k_max_over_mX'], alpha_max_over_m=r['alpha_max_over_mX']))
    summary = dict(passed=sum(c['passed'] for c in out), total=len(out), seconds=time.time() - t0, checks=out)
    (RES / 'checks.json').write_text(json.dumps(summary, indent=1, default=float) + '\n')
    print('passed %d / %d in %.0f s' % (summary['passed'], summary['total'], summary['seconds']))
    # reproduction (labelled): our ls_threshold/lsthermo.py series against this track's series
    try:
        sys.dont_write_bytecode = True
        sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'ls_threshold' / 'code'))
        import lsthermo as LT  # noqa
        worst = 0.0
        for tt in (0.3, 1.0, 3.0, 12.0):
            for uu in (-2.0, 0.0, 0.5, 0.9):
                B = LT.bessel(tt, 0.95)
                n, P, rho, s = LT.species(tt, uu, B)
                a = TH.species(tt, 1.0, uu)
                worst = max(worst, abs(a['n'] / n - 1), abs(a['P'] / P - 1), abs(a['rho'] / rho - 1), abs(a['s'] / s - 1))
        summary['reproduction_lsthermo_vs_this_track'] = dict(max_rel_diff=worst, label='reproduction with our earlier code, not independent')
        print('reproduction lsthermo vs here: %.2e' % worst)
        (RES / 'checks.json').write_text(json.dumps(summary, indent=1, default=float) + '\n')
    except Exception as exc:
        print('lsthermo comparison skipped:', exc)


if __name__ == '__main__':
    main()
