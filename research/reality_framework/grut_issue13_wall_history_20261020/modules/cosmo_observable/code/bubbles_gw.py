"""Part 5: (a) collapsing closed walls of a biased / non-percolating minority population, injected energy and the
stable relic it leaves in the unchanged action; (b) gravitational waves of a scaling network at the domination epoch
(Hiramatsu-Kawasaki-Saikawa arXiv:1309.5001 eqs. (32)-(39); estimate only). Writes results/bubbles_gw.json.

Population model (stated order of magnitude, as instructed): minority domains of size ~ t at time t with number
density n_b ~ p_- / t^3; each releases E ~ sigma t^2 when it collapses within ~ one Hubble time, so
    d rho_inj / d ln t = c_b p_- sigma / t,     c_b in [1, 3] (3 = sphere area/volume with R = t; 1 = n_b p_-/t^3).
A rest-mass fraction f of the injected energy ends in stable hidden quanta (free input, bracket [0.01, 1]).
"""
import json
import math
import os
import time
import numpy as np
from scipy.integrate import solve_ivp, quad
import core as C
import observables as O

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')

def t_of_T_nowalls(T_start=1e3, T_end=None, rho_L=0.0):
    """Standard background (SM radiation with g_*, matter) t(T) as a dense solution in x = ln T."""
    if T_end is None:
        T_end = C.T0
    def rhs(x, y):
        T = math.exp(x)
        H = math.sqrt((C.rho_rad(T) + C.rho_matter(T) + rho_L) / (3 * C.M_P**2))
        return [(-1.0 - C.dlngs_dlnT(T) / 3.0) / (H * math.exp(y[0]))]
    H0 = math.sqrt((C.rho_rad(T_start) + C.rho_matter(T_start)) / (3 * C.M_P**2))
    sol = solve_ivp(rhs, (math.log(T_start), math.log(T_end)), [math.log(1 / (2 * H0))], method='DOP853',
                    rtol=1e-10, atol=1e-12, dense_output=True)
    return sol

def main():
    t0 = time.time()
    out = {'module': 'cosmo_observable', 'part': '5 collapsing bubbles and GW'}
    sig = C.SIGMA_BENCH
    sol = t_of_T_nowalls()
    tT = lambda T: math.exp(sol.sol(math.log(T))[0])

    # ---------- spectrum / stability bookkeeping in the unchanged action -------------------------------------------
    mh, mX, mD = C.M_H, C.M_X, C.M_D
    out['products'] = {
        'm_h_GeV': mh, 'm_X_GeV': mX, 'm_D_GeV': mD,
        'h->XXbar_allowed': mh > 2 * mX, 'h->DD_allowed': mh > 2 * mD,
        'hh->XXbar_allowed_at_rest': 2 * mh > 2 * mX, 'hh->DD_allowed_at_rest': 2 * mh > 2 * mD,
        'XXbar->DD_allowed_at_rest': 2 * mX > 2 * mD, 'DD->XXbar_allowed_at_rest': 2 * mD > 2 * mX,
        'stable': 'X (lightest U(1)_X-charged and lightest hidden state), D (Z2_D odd), h (no allowed two-body '
                  'hidden decay; couples to the SM only through gravity in the unchanged action, so it is stable '
                  'except gravitationally). No hidden state is lighter than X and the action has no portal, so the '
                  'hidden rest mass cannot annihilate into lighter states except via gravity.',
        'rest_mass_retained_if_hh_to_XXbar': 2 * mX / (2 * mh),
    }

    # ---------- injected energy vs time, as a fraction of the total, per unit c_b p_- ----------------------------------
    Ts = np.logspace(3, math.log10(C.T0) + 3, 61)   # 1e3 GeV ... 1e3 T0
    inj = []
    for T in Ts:
        t = tT(T)
        rt = C.rho_rad(T) + C.rho_matter(T)
        inj.append({'T_GeV': T, 't_s': t * C.HBAR_GEV_S,
                    'drho_inj_dlnt_over_rho_tot_per_cb_pminus': sig / t / rt})
    out['injection_fraction_per_Hubble_time'] = inj
    # where the injection per Hubble time equals the total density for c_b p_- = 1, 0.1, 0.01
    def T_equal(cbp):
        g = lambda x: math.log(cbp * sig / tT(math.exp(x)) / (C.rho_rad(math.exp(x)) + C.rho_matter(math.exp(x))))
        from scipy.optimize import brentq
        return math.exp(brentq(g, math.log(1e3), math.log(C.T0 * 1.0001)))
    out['T_where_injection_per_Hubble_time_equals_rho_tot'] = {('cb_p=%g' % c): T_equal(c) for c in (1.0, 0.1, 0.01, 1e-3)}

    # ---------- relic: continuous scale-invariant population until T_end ----------------------------------------------
    # rho_relic(T0) = f c_b p_- int_{t_i}^{t_end} sigma / t^2 (a(t)/a0)^3 dt ; integrate in x = ln T
    def relic_Oh2_per_fcbp(T_end, T_i=1e3):
        def integrand(x):
            T = math.exp(x)
            t = tT(T)
            H = math.sqrt((C.rho_rad(T) + C.rho_matter(T)) / (3 * C.M_P**2))
            dtdx = (-1.0 - C.dlngs_dlnT(T) / 3.0) / H
            return sig / t**2 * C.a_ratio(T)**3 * (-dtdx)
        val, _ = quad(integrand, math.log(T_end), math.log(T_i), limit=400, epsrel=1e-8)
        return val / C.RHO_C100
    def single_epoch_Oh2_per_fp(Tc):
        # minority domains all collapse at t_c (radius ~ t_c): rho_inj = 3 p_- sigma / t_c, then non-relativistic
        return 3 * sig / tT(Tc) * C.a_ratio(Tc)**3 / C.RHO_C100
    rel = []
    for Te in list(np.logspace(1, -9, 41)) + [C.T0 * 10, C.T0 * 1.01]:
        y = relic_Oh2_per_fcbp(Te)
        y1 = single_epoch_Oh2_per_fp(Te)
        rel.append({'T_end_GeV': Te, 'Omega_h2_per_f_cb_pminus_continuous': y,
                    'max_f_cb_pminus_for_Oh2<=0.12_continuous': C.OMEGA_DM_H2_MAX / y,
                    'Omega_h2_per_f_pminus_single_epoch_cb3': y1,
                    'max_f_pminus_single_epoch': C.OMEGA_DM_H2_MAX / y1,
                    'rho_over_s_GeV_per_f_cb_pminus': y * C.RHO_C100 / C.s_ent(C.T0)})
    out['relic_vs_T_end'] = rel
    # T_end at which f c_b p_- = 1, 0.1, 0.01, 1e-3, 1e-6 saturates Omega h^2 = 0.12
    from scipy.optimize import brentq
    sat = {}
    for fp in (1.0, 0.1, 0.01, 1e-3, 1e-6, 1e-10):
        g = lambda x: math.log(fp * relic_Oh2_per_fcbp(math.exp(x)) / C.OMEGA_DM_H2_MAX)
        try:
            sat['f_cb_p=%g' % fp] = math.exp(brentq(g, math.log(C.T0 * 1.01), math.log(9e2)))
        except ValueError:
            sat['f_cb_p=%g' % fp] = 'no saturation for T_end >= 1.01 T0 (relic stays below 0.12 even if collapse continues to today)'
    out['T_end_min_GeV_for_Oh2<=0.12'] = sat
    out['rho_over_s_max_GeV_for_Oh2_0.12'] = C.OMEGA_DM_H2_MAX * C.RHO_C100 / C.s_ent(C.T0)
    # analytic radiation-era check: rho/s = 2 c_b f p_- sigma / (t_end s(T_end))  (continuous)
    Te = 1e-3
    out['analytic_check_T_end_1MeV'] = {
        'numeric_rho_over_s': relic_Oh2_per_fcbp(Te) * C.RHO_C100 / C.s_ent(C.T0),
        'analytic_2sigma_over_t_s': 2 * sig / tT(Te) / C.s_ent(Te),
        'note': 'agreement expected up to g_*s(T) variation of a^3 s'}

    # ---------- GW from the scaling network at the domination epoch (HKS 1309.5001) ----------------------------------
    h = C.history(sig, A=C.A_CENTRAL, n_post=1.0, f_switch=1.0)
    Td = h['T_switch_GeV']
    rw = h['rho_w_switch']
    rtot = C.rho_rad(Td) + C.rho_matter(Td) + rw
    Hd = math.sqrt(rtot / (3 * C.M_P**2))
    a_d = C.a_ratio(Td)
    gw = {}
    for eps_name, eps in (('eps=0.7', 0.7), ('eps=0.3', 0.3), ('eps=1.1', 1.1)):
        for A in (C.A_LOW, C.A_CENTRAL, C.A_HIGH):
            drho = eps * C.G_N * A**2 * sig**2              # HKS eq. (32): (d rho_gw / d ln k)_peak
            Om_d = drho / rtot                              # HKS eq. (33) with rho_c = total at T_dom
            rho0 = drho * a_d**4
            gw['%s_A=%.1f' % (eps_name, A)] = {
                'drho_gw_dlnk_peak_GeV4': drho, 'Omega_gw_peak_at_T_dom_vs_rho_tot': Om_d,
                'Omega_gw_peak_at_T_dom_vs_rho_rad': drho / C.rho_rad(Td),
                'rho_gw0_over_rho_c100_(Omega_GW_h2_equiv)': rho0 / C.RHO_C100,
                'HKS_eq35_Omega_GW_h2_using_Omega_vs_rho_rad': 1.34e-5 * (100.0 / C.gstar_rho(Td))**(1 / 3.) * drho / C.rho_rad(Td),
                'DeltaNeff_equiv_peak_per_lnk_at_T_rec': drho * (a_d / C.a_ratio(O.T_REC))**4 / C.rho_one_nu(O.T_REC),
                'Omega_GW_today_vs_wall_dominated_total_n1': rho0 / O.rho_total(h, C.T0)}
    out['GW_domination_epoch'] = {
        'T_dom_GeV': Td, 'H_dom_GeV': Hd, 'a_dom_over_a0': a_d,
        'f_peak_today_Hz': a_d * Hd / C.HBAR_GEV_S,
        'f_peak_formula': 'HKS eq.(38): f = (a_dom/a0) H_dom (k/a = 2 pi H)',
        'f_wall_width_today_Hz': a_d * C.M_H / (2 * math.pi) / C.HBAR_GEV_S,
        'eps_tilde_provenance': 'HKS eq.(30): 0.7 +- 0.4', 'cases': gw,
        'caveat': 'HKS amplitudes are measured in radiation-era scaling; T_dom is the edge of validity. '
                  'After domination the n=1 network freezes and emission is not modelled. Estimate only.'}
    out['runtime_s'] = time.time() - t0
    with open(os.path.join(RES, 'bubbles_gw.json'), 'w') as f:
        json.dump(out, f, indent=1)
    return out

if __name__ == '__main__':
    o = main()
    print(json.dumps({k: o[k] for k in ('products', 'T_where_injection_per_Hubble_time_equals_rho_tot',
                                        'T_end_min_GeV_for_Oh2<=0.12', 'rho_over_s_max_GeV_for_Oh2_0.12',
                                        'analytic_check_T_end_1MeV', 'runtime_s')}, indent=1))
    g = o['GW_domination_epoch']
    print('f_peak', g['f_peak_today_Hz'], 'f_width', g['f_wall_width_today_Hz'])
    print(json.dumps(g['cases']['eps=0.7_A=0.8'], indent=1))
    for r in o['relic_vs_T_end'][::4]:
        print(r)
