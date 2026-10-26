#!/usr/bin/env python3
"""kin_thermo (Claude, own code): assemble results/thermo_interface.json for the kinetic solver.

Sections: meta (conventions, import recipe), thermodynamics (species tables, closures, K tables, susceptibility
examples, Codex control, equilibrium comparators on the R grid, differences from the old code, thermal-mass
estimate), expansion (background, sources, energy equation), objects (E(Q), omega(Q) table), kernels (conversion
forward/reverse tables and thermal averages, surface exchange parameters, collision bracket), formation (rules and
parameters), old_proxies (reference f(R) curves and crossings), checks (summary).
Usage: python3 -B build_interface.py
"""
import json
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
R_GRID = [0.0, 1.0, 1.14445, 2.0, 3.0, 5.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 40.5626, 41.707,
          48.499, 49.643, 60.0, 75.0, 90.0]


def f(x):
    return float(x)


def species_tables():
    tg = np.geomspace(0.01, 200.0, 87)
    uX = np.array([-30, -20, -10, -6, -4, -3, -2, -1.5, -1, -0.7, -0.5, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
                   0.6, 0.7, 0.8, 0.85, 0.8611951308217225, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999])
    uS = np.array([-30, -20, -10, -6, -4, -3, -2, -1, -0.5, 0.0, 0.5, 1.0, 1.5, 1.7223902616434450, 2.0, 2.2, 2.3, 2.4, 2.44])
    keys = ('n', 'P', 'rho', 's', 'dn_du', 'drho_du', 'dn_dt', 'drho_dt')
    out = dict(grid_T_over_mX=tg.tolist(), units='m_X = 1: n in m_X^3, P and rho in m_X^4, s in m_X^3; derivatives per m_X',
               note='X and Xbar share the X table (evaluate at mu_+ and mu_-); s at mu_s (mass sqrt6 m_X). Prefer importing thermo.py.')
    for name, r, ug in (('X', 1.0, uX), ('s', TH.RS, uS)):
        tab = {k: np.zeros((len(tg), len(ug))) for k in keys}
        for i, t in enumerate(tg):
            for j, u in enumerate(ug):
                sp = TH.species(float(t), r, float(u))
                for k in keys:
                    tab[k][i, j] = sp[k]
        out[name] = dict(mass_over_mX=r, grid_mu_over_mX=ug.tolist(), **{k: tab[k].tolist() for k in keys})
    return out


def equilibrium_comparators():
    nq = C.N_Q0 / MX**3
    rows = []
    for R in R_GRID:
        for with_s in (False, True):
            for v0 in (True, False):
                rho = C.rho_onset(R, include_V0=v0) / MX**4
                try:
                    g = TH.state_equilibrium(nq, rho, C.OMEGA_REF, with_s=with_s)
                    rows.append(dict(R=R, include_s=with_s, include_V0=v0, T_over_mX=g['t'], mu_over_mX=g['up'],
                                     condensed_charge_GeV3=g['n_cond'] * MX**3, condensed_fraction=g['n_cond'] / nq,
                                     n_plus_GeV3=g['n_plus'] * MX**3, n_minus_GeV3=g['n_minus'] * MX**3, n_s_GeV3=g['n_s'] * MX**3,
                                     Xbar_per_net_charge=g['n_minus'] / nq, s_per_net_charge=g['n_s'] / nq,
                                     pair_inventory_per_net_charge=g['pair_inventory'] / nq,
                                     entropy_per_net_charge=g['s'] / nq, scalar_energy_fraction_gas=g['scalar_energy_fraction'],
                                     scalar_energy_fraction_total=g['S']['rho'] / g['rho_total'], rho_GeV4=g['rho_total'] * MX**4,
                                     P_GeV4=g['P'] * MX**4))
                except Exception as exc:
                    rows.append(dict(R=R, include_s=with_s, include_V0=v0, status='failure %r' % exc))
    return rows


def K_tables():
    tg = np.geomspace(0.03, 100.0, 60)
    return dict(T_over_mX=tg.tolist(), K_bose_mX3=[f(TH.K_bose(t)) for t in tg], K_MB_mX3=[f(TH.K_mb(t)) for t in tg],
                definition='K = n_s^eq/(n_+^eq n_-^eq) at zero chemical potentials; MB version is lsrelic.K_s')


def thermal_mass_estimate():
    """High-T one-loop estimates [E] (exact one-loop with J_B functions is kin_rates' job)."""
    g, lp, lx, v = C.GX, C.LP, C.LX, C.V
    Tc = np.sqrt(lp * v * v / (lp / 4 + g * g / 6))
    rows = []
    for t in (0.3, 1.0, 3.0, 7.34, 12.0, 15.0, 18.0):
        T = t * MX
        vT2 = v * v - T * T * (lp / 4 + g * g / 6) / lp
        vT = np.sqrt(max(vT2, 0.0))
        mX2 = g * g * vT2 + T * T * (lx / 3 + g * g / 12)
        ms2 = 2 * lp * vT2 if vT2 > 0 else (lp / 4 + g * g / 6) * T * T - lp * v * v
        rows.append(dict(T_over_mX=t, vT_over_v=vT / v, mX_T_over_mX=np.sqrt(max(mX2, 0)) / MX, ms_T_over_ms=np.sqrt(max(ms2, 0)) / C.MS,
                         broken=bool(vT2 > 0)))
    return dict(T_c_GeV=f(Tc), T_c_over_mX=f(Tc / MX), rows=rows,
                formulas='v_T^2 = v^2 - T^2 (lP/4 + g^2/6)/lP; m_X^2(T) = g^2 v_T^2 + T^2 (lX/3 + g^2/12); m_s^2 = 2 lP v_T^2 (broken) [E high-T]',
                mean_field_work='at the self-consistent mean field the Gibbs-Duhem identities hold with explicit (T, mu) derivatives; rho includes V_eff(v_T, T); tadpoles from thermo.species dP/dm^2 if needed',
                objects_require='broken phase (T < T_c): no Q-balls above T_c')


def kernel_section():
    p = KN.P_TAB
    svf = KN.SV_TAB
    pp = np.r_[np.linspace(2.533, 6.0, 200), np.linspace(6.1, 40.0, 100)]
    svr = KN.sv_reverse(pp)
    tgrid = [0.03, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 7.34, 12.0]
    avg = []
    for t in tgrid:
        for label, up, um in (('mu=0', 0.0, 0.0), ('cap: mu_+=omega, mu_-=-omega', C.OMEGA_REF, -C.OMEGA_REF)):
            for ex in ('hold', 'geometric'):
                r = KN.conversion_rates(t, up, um, extrap=ex)
                r['label'] = label
                avg.append({k: (f(v) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})
    surf = []
    for t in (0.1, 0.3, 1.0, 3.0, 7.34, 12.0):
        for up in (0.0, 0.5, C.OMEGA_REF, 0.95):
            s = KN.surface_exchange(t, up)
            surf.append({k: f(v) for k, v in s.items()})
    geo = KN.GEO
    return dict(
        conversion_forward=dict(process='Xbar(p_-) + core(Q=200000) -> X(p_+) + core(Q-2); E_+ = E_- + 2 omega m_X',
                                p_minus_over_mX=p.tolist(), sigma_v_GeV2=svf.tolist(), sigma_v_cm3_s=(svf * KN.GEV2_TO_CM3S).tolist(),
                                source='tracks/conversion results/rate_scan.json dense scan + scatter_spot cold value [Q own]; cold profile; '
                                       'frozen-profile weak-wave rate; above p_-=0.708 two labelled extrapolations (hold, geometric)',
                                hard_X_energy='E_out = E_in + 2 omega m_X: cold incoming -> E_out = 2.7224 m_X = 251.57 MeV, K_out = 1.7224 m_X, p_out = 2.532 m_X'),
        conversion_reverse=dict(process='X(p_+) + core(Q-2) -> Xbar + core(Q), reciprocity (sigma v)_r = (sigma v)_f p_- E_-/(p_+ E_+)',
                                p_plus_over_mX=pp.tolist(), sigma_v_GeV2=svr.tolist(), threshold_p_plus_over_mX=f(np.sqrt((1 + 2 * C.OMEGA_REF)**2 - 1))),
        thermal_averages=avg,
        detailed_balance='forward = reverse when mu_+ - mu_- = E(Q) - E(Q-2) = 2 omega m_X; with mu_- = -mu_+ this is mu_+ = omega (checked in checks.json)',
        surface_exchange=dict(model='Langmuir flux with detailed balance at (T, mu = omega); eta in [0,1] bracket (0 for a cold core at linear order: two-body X absorption closed)',
                              radius_GeV=f(geo['r90']), radius_over_mX=f(geo['r90'] * MX), area_GeV2=f(4 * np.pi * geo['r90']**2),
                              hard_sphere_radius_from_sigma_el_over_mX=5.2, examples=surf),
        object_collision_bracket=dict(sigma_GeV2=f(np.pi * (2 * geo['r90'])**2), formula='rate per object = n_Lbar sigma v_rel; not computed'),
        object_charge_conjugate='X on anti-objects: same tables with X <-> Xbar')


def formation_section():
    ob = json.loads((RES / 'objects_EQ.json').read_text())
    frag_onset = {s: FM.rule_fragmentation(absX0_GeV=C.XC, shape=s) for s in (0.1, 0.3, 1.0)}
    frag_gas = []
    for fc in (0.03, 0.1, 0.3, 0.7, 1.0):
        for s in (0.1, 0.3, 1.0):
            r = FM.rule_fragmentation(n_cond_GeV3=fc * C.N_Q0, shape=s)
            r['condensed_fraction_of_net_charge'] = fc
            frag_gas.append(r)
    H_ref = C.hubble(C.T_vis_of_nQ(C.N_Q0))   # consistent g* = g_s = 106.75 background (T_vis = 43.72 GeV, H = 2.685e-15 GeV)
    return dict(
        rules=dict(
            gas_only=FM.rule_gas_only(),
            instant_condensation=dict(callable='formation.rule_instant(t_hint, n_Q, P_inv, rho, variant)', variants=dict(
                omega='e_cap = omega = 0.861195 m_X (consistent linear comparator)', E_over_Q='e_cap = E/Q = 0.923316 m_X (benchmark booking)',
                finite_Q='omega d(n_L Q) + (E - Q omega) dn_L with E(Q_ref) = %.6f GeV, E - Q omega = %.4f GeV' % (ob['reference_Q200000']['E_GeV'], ob['reference_Q200000']['E_minus_Q_omega_GeV'])),
                supplies='instant equilibrium condensation at mu_+ = omega(Q_ref) and dissolution; anti-objects at mu_- = omega; no formation time, no Q distribution'),
            fragmentation=dict(callable='formation.rule_fragmentation(n_cond_GeV3 | absX0_GeV, shape, Ne, f_cap, f_anti, T_GeV)',
                               lambda_eff=FM.LAM_EFF, U_eff='m_X^2|X|^2 - lambda_eff|X|^4 (|X| <= x_c); lambda_X|X|^4 + V0 beyond',
                               dispersion='(alpha^2+k^2)(alpha^2+k^2+U\'\'-omega^2) + 4 omega^2 alpha^2 = 0; band k^2 < 4 lambda_eff |X0|^2',
                               at_onset_release_amplitude_xc=frag_onset, from_gas_condensate=frag_gas,
                               brackets=dict(shape=[0.1, 1.0], Ne=[20, 40], f_cap=[0.5, 1.0], f_anti=[0.0, 1.0]),
                               H_at_onset_GeV=f(H_ref), alpha_max_over_H_at_onset=f(frag_onset[0.3]['alpha_max_GeV'] / H_ref),
                               supplies='k_max, alpha_max, Q_lump(shape), t_frag, object density, energy per lump from the E(Q) table; no lattice')),
        constraints=dict(Q_min_stationary=ob['Q_min'], E_over_Q_equals_m_at=ob.get('E_over_QmX_crosses_1_at'),
                         thin_wall_omega_over_mX=ob.get('thin_wall_omega_over_mX_from_Ueff'),
                         T_c_high_T_estimate_over_mX=f(FM.TC_HIGH_T / MX),
                         statement='fragments with Q < Q_min are not stationary objects (must merge or dissolve); objects need T < T_c'),
        nucleation_from_gas=dict(barrier='F_crit(mu) = max_Q [E(Q) - mu Q] = (E - Q omega) on the thick-wall branch at omega = mu; rate ~ T^4 e^{-F_crit/T} [E]',
                                 table=FM.nucleation_barrier(),
                                 statement='F_crit/T >= 43 even at mu = 0.998 m and T = 12 m: homogeneous nucleation of finite objects from an undersaturated gas is negligible; the gas first saturates (mu -> m), then the condensate fragments'),
        anti_objects='for an elliptic release (epsilon = 0.146) the literature reports Q-balls of both signs [R]; f_anti is a bracket, and anti-object charge must be reported by the solver',
        literature=[
            'Kusenko & Shaposhnikov, PLB 418 (1998) 46, hep-ph/9709492 [abstract read; fragmentation analysis recalled]',
            'Enqvist & McDonald, PLB 425 (1998) 309, hep-ph/9711514 [abstract read]',
            'Kasuya & Kawasaki, PRD 61 (2000) 041301, hep-ph/9909509 [abstract read]',
            'Kasuya & Kawasaki, PRD 62 (2000) 023512, hep-ph/0002285 [abstract read]',
            'Hiramatsu, Kawasaki & Takahashi, JCAP 06 (2010) 008, 1003.1779 [abstract read: 3D peak charge 1.9e-2 (|Phi_in|/m)^2]',
            'Laine & Shaposhnikov, NPB 532 (1998) 376, hep-ph/9804237 [abstract read]',
            'Griest & Kolb 1989; Postma 2002 (solitosynthesis) [recalled]'])


def old_proxies_section():
    p = RES / 'old_proxies.json'
    if not p.exists():
        return dict(status='not yet available')
    d = json.loads(p.read_text())
    rows = d['rows']

    def curve(kind, s_ch, mu_c=None):
        sel = [r for r in rows if r.get('kind') == kind and r.get('s_channel') == s_ch and 'f' in r and (mu_c is None or abs(r.get('mu_c', 0) - mu_c) < 1e-6)]
        sel.sort(key=lambda r: r['R'])
        return [dict(R=r['R'], f=r['f'], p_final=r.get('p_final'), root_fallbacks=r.get('root_fallbacks'),
                     energy_work_cumulative_residual=r.get('energy_work_cumulative_residual'), Th0_over_m=r.get('Th0_over_m'),
                     sigma0=r.get('sigma0'), seconds=r.get('seconds')) for r in sel]

    def crossing(c, target=0.01):
        pts = [(x['R'], x['f']) for x in c if x['f'] is not None]
        for (r0, f0), (r1, f1) in zip(pts, pts[1:]):
            if f0 <= target < f1:
                if f0 <= 0:
                    return dict(bracket=[r0, r1], R_L=None, note='f=0 below (erased)')
                tt = (np.log(target) - np.log(f0)) / (np.log(f1) - np.log(f0))
                return dict(bracket=[r0, r1], R_L=float(np.exp(np.log(r0) + tt * (np.log(r1) - np.log(r0)))))
        return dict(bracket=None, R_L=None)
    out = dict(eps=d['eps'], ymult=d['ymult'], R_grid=d['R_grid'], curves={}, crossings={}, proxy_thresholds=[])
    for name, kind, s_ch, mu_c in (('proxy_isentropic', 'proxy', False, None), ('proxy_isentropic_s_route', 'proxy', True, None),
                                   ('ke_energy_conserving', 'ke', False, 0.8611951308217225), ('ke_energy_conserving_s_route', 'ke', True, 0.8611951308217225),
                                   ('ke_energy_conserving_muc_EoverQ', 'ke', False, 0.9233161353726057)):
        c = curve(kind, s_ch, mu_c)
        out['curves'][name] = c
        out['crossings'][name] = crossing(c)
    out['proxy_thresholds'] = [dict(s_channel=r['s_channel'], bracket=r['bracket']) for r in rows if r.get('kind') == 'proxy_threshold']
    for name, s_ch in (('proxy_isentropic', False), ('proxy_isentropic_s_route', True)):   # the grid is too coarse there
        thr = [r['bracket'] for r in rows if r.get('kind') == 'proxy_threshold' and bool(r.get('s_channel')) == s_ch]
        if thr:
            out['crossings'][name].update(grid_interpolation_R_L=out['crossings'][name].get('R_L'), bisection_bracket=thr[0],
                                          R_L=float(0.5 * (thr[0][0] + thr[0][1])),
                                          note='R_L is the bisection midpoint (N = 400); the coarse-grid log interpolation between R = 75 and 90 is kept only as grid_interpolation_R_L')
    out['failures'] = [r for r in rows if 'status' in r]
    out['label'] = 'reference evaluations with our ls_threshold code (reuse stated); equilibrium proxies, not kinetics'
    return out


def main():
    t0 = time.time()
    ob = json.loads((RES / 'objects_EQ.json').read_text())
    checks = json.loads((RES / 'checks.json').read_text()) if (RES / 'checks.json').exists() else None
    iface = dict(
        meta=dict(track='kin_thermo', date='2026-10-26', author='Claude (own code; reuse stated per section)',
                  conventions=dict(units='GeV; hidden thermodynamics internally in m_X = 1 (m_X = %.10f GeV)' % MX,
                                   species='X (+1, mu_+), Xbar (-1, mu_-), s (neutral, mu_s; m_s = sqrt6 m_X), one bosonic degree each, Bose statistics',
                                   closures=dict(pairs='mu_+ + mu_- = 0', fast_s='mu_s = mu_+ + mu_-', full='both: mu_s = 0',
                                                 caps='mu_+ <= omega (objects), mu_- <= omega (anti-objects), mu_s < m_s'),
                                   pair_inventory='P_inv = n_- + n_s (+ anti-condensate): conserved by s <-> X Xbar, lowered by ss -> X Xbar and 4 -> 2',
                                   f_sym='f = 2 m_X Y_Q (pairs per net charge)/(4.4e-10 GeV) = 2.028545 x pairs per net charge'),
                  import_recipe='sys.path.insert(0, ".../tracks/kin_thermo/code"); import thermo, objects, kernels, formation, expansion',
                  functions=dict(thermo=['species(t,r,u)', 'gas(t,up,um,us,fast_s,with_s)', 'susceptibility(...)', 'state_equilibrium(nq,rho,u_cap,e_cap,with_s)',
                                         'state_two_potential(nq,pinv,rho,u_cap,e_cap,anti_cap,with_s,fast_s)', 'state_three_potential(nq,nm,ns,rho)', 'K_bose(t)', 'K_mb(t)'],
                                 expansion=['parent_source(a,key)', 'hubble_of_a(a)', 'object_work(...)', 'drho_gas_dlna(...)'],
                                 kernels=['sv_forward(p,extrap)', 'sv_reverse(p_plus,extrap)', 'conversion_rates(t,up,um,extrap)', 'surface_exchange(t,up,eta)', 'object_collision_bracket(...)'],
                                 formation=['rule_gas_only()', 'rule_instant(t_hint,nq,pinv,rho,variant)', 'rule_fragmentation(...)', 'nucleation_barrier()']),
                  fixed_point=dict(mX_GeV=MX, ms_GeV=C.MS, gX=C.GX, lambda_X=C.LX, lambda_phi=C.LP, v_GeV=C.V, YQ=C.YQ, DM_per_entropy_GeV=C.DM_PER_S,
                                   allowance_f=C.F_ALLOW, eps=C.EPS, N_X=C.N_X, E_X=C.E_X, H_i_GeV=C.H_I, Q_ref=C.Q_REF,
                                   omega_ref_over_mX=C.OMEGA_REF, E_over_QmX_ref=C.E_OVER_QM_REF, x_c_GeV=C.XC, n_c_GeV3=C.N_C,
                                   rho_c_GeV4=C.RHO_C, V0_GeV4=C.V0, n_Q0_GeV3=C.N_Q0, commits=dict(parent='225028b', kinetic='9877518', conversion='bb211d1')),
                  labels='[D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled',
                  scope='state functions, object sector and formation rules for one supplied model; no Boltzmann solution, no f_sym, no viability'),
        thermodynamics=dict(
            formulas=TH.__doc__,
            species_tables=species_tables(),
            K_tables=K_tables(),
            codex_control=checks['checks'][2:4] if checks else None,
            equilibrium_comparators=dict(definition='full chemical equilibrium (pairs + fast s) at the onset charge n_Q0 and energy rho_c(1+R) (+V0), objects at mu = omega with linear energy omega n_c; comparators, not populations',
                                         rows=equilibrium_comparators()),
            pair_inventory_susceptibility=checks['checks'][6]['rows'] if checks else None,
            thermal_mass_estimate=thermal_mass_estimate(),
            differences_from_old_code=dict(
                lsthermo_thermal_relic='our ls_threshold/lsthermo.py and sp_yledger/thermal_relic.py (gp_leftovers equilibrium.py) carry X and Xbar only; s enters lsrelic only through the rate factor K_MB(T) n_+ n_- (MB statistics); no s number, energy, pressure or entropy in the state; no susceptibilities; single-species caps only',
                here='explicit s species with Bose statistics and its own chemical potential; three closures; analytic susceptibility matrix; pair inventory with s; three state solvers; mass-dependent tadpoles available',
                numbers_at_R_40_5626='T/m 8.11769 -> 7.34390; condensed 0.004405 -> 0.007072 GeV^3; pair inventory 0.04452 -> 0.06851 GeV^3 (s holds 0.03604); scalar energy 32.76% of total (32.85% of gas); see equilibrium_comparators for every R')),
        expansion=dict(background=dict(gstar=C.GSTAR_CONST, hubble='sqrt(pi^2 g*/90) T_v^2/M_Pl, M_Pl = 2.435e18 GeV', a_ref=EX.A_REF,
                                       T_vis_ref_GeV=f(C.T_vis_of_nQ(C.N_Q0)), H_ref_GeV=f(C.hubble(C.T_vis_of_nQ(C.N_Q0))),
                                       H_ref_table_gstar_GeV=f(C.hubble(C.T_vis_of_nQ(C.N_Q0, table=True), table=True)),
                                       comoving_charge='n_Q = Y_Q s_vis(T_v)'),
                       parents=EX.PARENTS, energy_equation=EX.__doc__,
                       checks=dict(isentropic=EX.check_isentropic(), parent_ledger=[EX.check_parent_ledger(k) for k in ('e1e-7_band', 'e1e-6_band')],
                                   object_work=EX.check_object_work())),
        objects=dict(reference_Q200000=ob['reference_Q200000'], Q_min=ob['Q_min'], thin_wall_omega_over_mX=ob.get('thin_wall_omega_over_mX_from_Ueff'),
                     E_over_Q_equals_m_at=ob.get('E_over_QmX_crosses_1_at'), identity_check=ob.get('dEdQ_identity_check'),
                     dEdQ_finite_difference=ob['dEdQ_finite_difference_at_reference'],
                     table_columns='Omega (= omega/m_X = dE/dQ / m_X), Q, E_GeV, E_over_QmX, E_minus_Q_omega_GeV, F0, S0, r50/r90/r99 over m_X, virial, dQ_dOmega, branch',
                     table=[{k: r[k] for k in ('Omega', 'Q', 'E_GeV', 'E_over_QmX', 'E_minus_Q_omega_GeV', 'F0', 'S0', 'r50_over_mX', 'r90_over_mX',
                                               'r99_over_mX', 'virial', 'dQ_dOmega', 'branch', 'dEdQ_over_mX_table')} for r in ob['rows']],
                     work_term='d rho_L = omega d(n_L Q) + (E - Q omega) dn_L; linear comparator: rho_L = omega Q n_L (E_c = mu_c = omega)'),
        kernels=kernel_section(),
        formation=formation_section(),
        old_proxies=old_proxies_section(),
        checks=dict(passed=checks['passed'], total=checks['total'], names=[(c['name'], c['passed'], c['residual']) for c in checks['checks']]) if checks else None,
        build_seconds=None)
    iface['build_seconds'] = time.time() - t0
    (RES / 'thermo_interface.json').write_text(json.dumps(iface, indent=1, default=float) + '\n')
    print('written', (RES / 'thermo_interface.json').stat().st_size / 1e6, 'MB in %.0f s' % iface['build_seconds'])


if __name__ == '__main__':
    main()
