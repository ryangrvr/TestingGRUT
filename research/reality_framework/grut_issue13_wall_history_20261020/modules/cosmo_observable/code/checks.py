"""Numerical and consistency checks. Each check records value, reference, tolerance and PASS/FAIL; a failing check
is reported as failing (it is not removed). Writes results/checks.json. Requires the other scripts to have run.
"""
import json
import math
import os
import core as C
import observables as O

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')

def load(name):
    with open(os.path.join(RES, name)) as f:
        return json.load(f)

def main():
    checks = []
    def chk(name, value, ref, rel_tol=None, abs_tol=None, note=''):
        if rel_tol is not None:
            ok = abs(value - ref) <= rel_tol * abs(ref)
        else:
            ok = abs(value - ref) <= abs_tol
        checks.append({'name': name, 'value': value, 'reference': ref, 'rel_tol': rel_tol, 'abs_tol': abs_tol,
                       'pass': bool(ok), 'note': note})

    hb, nm, sf, bg = load('history_benchmark.json'), load('normalization.json'), load('scale_family.json'), load('bubbles_gw.json')

    # g_* fits against Saikawa-Shirai eqs. (3.18)-(3.19) and continuity at the 120 MeV junction
    chk('g_rho(T0) vs 3.383 (1803.01038 eq.3.19)', C.G_RHO0, 3.383, abs_tol=0.003)
    chk('g_s(T0) vs 3.931 (1803.01038 eq.3.18)', C.G_S0, 3.931, abs_tol=0.004)
    chk('g_rho continuity at 120 MeV', C.gstar_rho(0.12 * (1 + 1e-9)), C.gstar_rho(0.12 * (1 - 1e-9)), rel_tol=2e-3)
    chk('g_s continuity at 120 MeV', C.gstar_s(0.12 * (1 + 1e-9)), C.gstar_s(0.12 * (1 - 1e-9)), rel_tol=2e-3)
    chk('g_rho(10 MeV) ~ 10.75 (e, nu, gamma)', C.gstar_rho(0.01), 10.75, rel_tol=0.01)
    chk('Omega_gamma h^2 = 2.473e-5 standard', (math.pi**2 / 15 * C.T0**4) / C.RHO_C100, 2.4728e-5, rel_tol=2e-3)
    chk('T_nu/T at 1 keV = (4/11)^(1/3)', C.Tnu_over_T(1e-6), (4 / 11.)**(1 / 3.), rel_tol=1e-6)

    # tension and conversions
    chk('sigma_bench = 9.428090416e7 GeV^3 (BRIEF)', C.SIGMA_BENCH, 9.428090416e7, rel_tol=1e-9)
    chk('sigma^(1/3) = 455.136 GeV (BRIEF)', C.SIGMA_BENCH**(1 / 3.), 455.136, abs_tol=5e-4)
    # independent conversion route: MeV^3 -> MeV/fm^2 -> J/m^2
    MeV_J = 1.602176634e-13
    hbarc_MeV_fm = 197.3269804
    one_GeV3_MeV_fm2 = 1e9 / hbarc_MeV_fm**2               # (1000 MeV)^3/(hbar c)^2 in MeV/fm^2
    one_GeV3_J_m2 = one_GeV3_MeV_fm2 * MeV_J / 1e-30
    chk('1 GeV^3 in J/m^2 via MeV/fm^2 route', C.GEV3_J_M2, one_GeV3_J_m2, rel_tol=1e-12)
    chk('Codex G sigma/H0 = 1e-5 reproduces 2.23e-9 GeV^3', nm['Codex_G_sigma_over_H0']['sigma_GeV3'], 2.23e-9, rel_tol=3e-3)
    chk('Codex benchmark ratio 4.24e16', nm['Codex_G_sigma_over_H0']['benchmark_over_this'], 4.24e16, rel_tol=2e-3)
    chk('LMS kg/cm^2 reading reproduces 0.96 MeV label', nm['LMS_1505.03673']['R+M']['kg_per_cm2_reading_cuberoot_MeV'], 0.96, rel_tol=0.02,
        note='reconciliation check; the literal kg/m^2 reading fails this by a factor 21')
    chk('LMS kg/cm^2 reading reproduces 0.93 MeV label', nm['LMS_1505.03673']['R+M+Lambda']['kg_per_cm2_reading_cuberoot_MeV'], 0.93, rel_tol=0.02)
    chk('LMS literal SI reading reproduces 0.96 MeV label', nm['LMS_1505.03673']['R+M']['literal_cuberoot_MeV'], 0.96, rel_tol=0.02,
        note='EXPECTED TO FAIL: documents the unit inconsistency in arXiv:1505.03673v1')

    # domination time
    td = hb['T_dom']
    chk('t_dom (no backreaction) vs Saikawa eq.(2.19) 3 M_P^2/(4 A sigma)', td['A=0.8_backreaction=False']['t_dom_s'],
        td['A=0.8_analytic']['t_dom_s_Saikawa_eq2.19'], rel_tol=2e-3)
    Hd = td['A=0.8_analytic']['H_dom_GeV']
    rho_target = 3 * C.M_P**2 * Hd**2
    from scipy.optimize import brentq
    T_an = math.exp(brentq(lambda x: math.log(C.rho_rad(math.exp(x)) / rho_target), math.log(1e-7), math.log(1e-4)))
    chk('T_dom (no backreaction) vs closed form rho_rad(T)=3 M_P^2 (2 A sigma/3 M_P^2)^2', td['A=0.8_backreaction=False']['T_dom_GeV'],
        T_an, rel_tol=5e-3)
    h100 = C.history(C.SIGMA_BENCH, T_start=100.0, T_end=C.T0 * 1.001)
    h1000 = C.history(C.SIGMA_BENCH, T_start=1000.0, T_end=C.T0 * 1.001)
    chk('T_dom insensitive to T_start (100 vs 1000 GeV)', h1000['T_switch_GeV'], h100['T_switch_GeV'], rel_tol=1e-6)
    hlo = C.history(C.SIGMA_BENCH, rtol=1e-7)
    hhi = C.history(C.SIGMA_BENCH, rtol=1e-11)
    chk('T_dom ODE tolerance 1e-7 vs 1e-11', hlo['T_switch_GeV'], hhi['T_switch_GeV'], rel_tol=1e-6,
        note='FAILED as first written: rtol=1e-7 differs by 1.4e-6; kept as a record. Production uses rtol=1e-10 (next check)')
    hprod = C.history(C.SIGMA_BENCH, rtol=1e-10)
    htight = C.history(C.SIGMA_BENCH, rtol=1e-12)
    chk('T_dom production rtol 1e-10 vs 1e-12', hprod['T_switch_GeV'], htight['T_switch_GeV'], rel_tol=1e-8)
    chk('rho_w0 production rtol 1e-10 vs 1e-12', hprod['rho_w_end'], htight['rho_w_end'], rel_tol=1e-8)
    chk('H0 ODE-independent: rho_w0 tolerance 1e-7 vs 1e-11', hlo['rho_w_end'], hhi['rho_w_end'], rel_tol=1e-5)

    # late-time analytic limits
    b1 = hb['runs']['n=1_A=0.8_fswitch=1_rhoL=0']
    b2 = hb['runs']['n=2_A=0.8_fswitch=1_rhoL=0']
    H1 = C.kmsmpc_to_gev(b1['H0_pred_kmsMpc']) / C.HBAR_GEV_S
    H2 = C.kmsmpc_to_gev(b2['H0_pred_kmsMpc']) / C.HBAR_GEV_S
    chk('n=1 wall-dominated: q0 = -1/2', b1['q0'], -0.5, abs_tol=1e-6)
    chk('n=2 wall-dominated: q0 = 0', b2['q0'], 0.0, abs_tol=1e-6)
    chk('n=1 age = 2/H0 (a ~ t^2)', b1['age_pred_s'], 2.0 / H1, rel_tol=1e-3)
    chk('n=2 age = 1/H0 (a ~ t)', b2['age_pred_s'], 1.0 / H2, rel_tol=2e-3)
    rw_dom = b1['fractions']  # placeholder to keep structure
    hb1 = C.history(C.SIGMA_BENCH, n_post=1.0)
    rho_w0_an = hb1['rho_w_switch'] * hb1['a_switch_over_a0']
    chk('n=1 rho_w0 = rho_w,dom (a_dom/a0)', b1['rho_w0_GeV4'], rho_w0_an, rel_tol=1e-9)
    tn1 = b1['tuned_negative_vacuum_H0=70']
    tn2 = b2['tuned_negative_vacuum_H0=70']
    H70 = C.kmsmpc_to_gev(70.0) / C.HBAR_GEV_S
    R1 = b1['R_H']['70']
    R2 = b2['R_H']['70']
    chk('tuned rho_L<0, n=1: age = pi/(R_H H0)', tn1['age_s'], math.pi / (R1 * H70), rel_tol=1e-2)
    chk('tuned rho_L<0, n=2: age = (pi/2)/(R_H H0)', tn2['age_s'], (math.pi / 2) / (R2 * H70), rel_tol=1e-2)
    chk('tuned rho_L<0, n=1: q0 = R_H^2/2', tn1['q0'], R1**2 / 2, rel_tol=1e-3)
    b4 = hb['runs']['n=4_A=0.8_fswitch=1_rhoL=0']
    chk('n=4 stress test: R_H equals no-wall value (radiation negligible today)', b4['R_H']['70'],
        hb['headline']['no_walls_reference']['rhoL=0']['R_H_vs70'], rel_tol=1e-3)

    # scale family cross-check: two code paths
    chk('v_max(dom, rhoL=0) root vs scan characteristic v (T_dom=T0)', sf['v_max_dom']['A=0.800_rhoL=0']['v_max_MeV'],
        sf['characteristic_v']['T_dom=T_CMB_today_rhoL0']['v_GeV'] * 1e3, rel_tol=1e-3)
    chk('sigma_max(dom, Planck-like rhoL) vs analytic 3 M_P^2 H0^2 t0/A (order of magnitude)',
        sf['v_max_dom']['A=0.800_rhoL=planck_like']['sigma_max_GeV3'], sf['v_max_dom']['analytic_rho_w0_eq_rho_c0_A0.8']['sigma_GeV3'],
        rel_tol=0.5, note='analytic omits backreaction and the definition rho_w = rho_others; agreement within 50% expected')

    # relic analytic and GW cross-checks
    ac = bg['analytic_check_T_end_1MeV']
    chk('continuous relic rho/s vs 2 sigma/(t s) at T_end = 1 MeV', ac['numeric_rho_over_s'], ac['analytic_2sigma_over_t_s'], rel_tol=0.03)
    g = bg['GW_domination_epoch']['cases']['eps=0.7_A=0.8']
    chk('GW today: direct a^4 redshift vs HKS eq.(35)', g['rho_gw0_over_rho_c100_(Omega_GW_h2_equiv)'],
        g['HKS_eq35_Omega_GW_h2_using_Omega_vs_rho_rad'], rel_tol=0.03)

    n_pass = sum(c['pass'] for c in checks)
    out = {'n_checks': len(checks), 'n_pass': n_pass, 'n_fail': len(checks) - n_pass,
           'failing': [c['name'] for c in checks if not c['pass']], 'checks': checks}
    with open(os.path.join(RES, 'checks.json'), 'w') as f:
        json.dump(out, f, indent=1)
    return out

if __name__ == '__main__':
    o = main()
    for c in o['checks']:
        print('%-4s %s: %.6g vs %.6g' % ('PASS' if c['pass'] else 'FAIL', c['name'], c['value'], c['reference']))
    print('%d/%d pass; failing: %s' % (o['n_pass'], o['n_checks'], o['failing']))
