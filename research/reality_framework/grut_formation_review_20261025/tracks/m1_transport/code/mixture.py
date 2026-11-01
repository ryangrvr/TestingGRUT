"""Track m1_transport, item 3/4: two-fluid estimate of the cold + warm-gas power ratio.

Independent code [E]: sub-horizon Newtonian two-fluid linear perturbations in the standard background
of transport_check.py.  Cold cores (fraction 1-f) are pressureless; the gas (fraction f) is a fluid
with sound speed c_s(a) = c_s,i a_i/a.  This is an order-of-magnitude cross-check of Codex's BGK
hierarchy values (P_mix/P_cold = 0.893 at k = 1e5/Mpc etc.), not a reproduction of it: the gas is
treated as a fluid (no Boltzmann hierarchy), radiation only through the background, Phi from
Poisson.  Variables: delta_c, delta_g and their d/dlna.

    d2 delta/dN2 + (1 + dln(aH)/dlna) d delta/dN + (k c_s/aH)^2 delta = (3/2) Omega_m(a) [(1-f) delta_c + f delta_g]

Two gas initial states bracket the matching at 2 MeV: 'clustered' (delta_g = delta_c, same velocity)
and 'erased' (delta_g = 0, already phase-mixed).  A 'never clusters' run (delta_g held at 0) gives the
growth-rate plateau.  An optional gas-core momentum-exchange (drag) term with rate Gamma_d(a) is
included to show the effect of the elastic X-core scattering that Codex's transport omits.
"""
import json, math, os, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import transport_check as tc

ROOT = Path(__file__).resolve().parents[1]
bg = tc.Background()
a_i = bg.a_i
H_i = float(bg.H_GeV(a_i))
fgas = 0.058767411685782034          # Codex's test-model gas fraction of all matter [Q]
v_rms = 0.0684
Mpc_per_GeVinv = tc.c_km_s / (tc.GeV_to_s_inv * tc.Mpc_km)  # Mpc per GeV^-1 ... (c/H in Mpc for H in GeV)


def aH_Mpc(a):
    """comoving aH in Mpc^-1."""
    return a * bg.H_GeV(a) * tc.GeV_to_s_inv * tc.Mpc_km / tc.c_km_s


def dlnH_dlna(a, eps=1e-4):
    return (math.log(bg.H_GeV(a * (1 + eps))) - math.log(bg.H_GeV(a * (1 - eps)))) / (math.log(1 + eps) - math.log(1 - eps))


def run(k, cs2_i, gas_init='erased', never_clusters=False, drag_i=0.0, drag_power=-4.0, f=fgas, N_end=None, with_gas=True):
    """Return (delta_m_mix, delta_c_only) at a = 1 (or a_eq) for wavenumber k [Mpc^-1]."""
    N_end = math.log(1 / a_i) if N_end is None else N_end
    aHi = aH_Mpc(a_i)
    # cold Meszaros growing mode at a_i: delta ~ ln(a/a_k) + const, take dlnδ/dN = 1/(ln(a_i/a_k)+1) with a_k: k = aH
    a_k = a_i * aHi / k if k > aHi else a_i
    slope = 1.0 / (max(math.log(a_i / a_k), 0.0) + 1.0)

    def rhs(N, y, gas_on):
        a = a_i * math.exp(N)
        aH = aH_Mpc(a)
        fric = 1 + 1 + dlnH_dlna(a)           # 1 + dln(aH)/dlna
        Om_a = tc.Om / a**3 / float(bg.E2(a))
        dc, dcp, dg, dgp = y
        if never_clusters and gas_on:
            dm = (1 - f) * dc
        else:
            dm = (1 - f) * dc + f * dg if gas_on else dc
        src = 1.5 * Om_a * dm
        cs2 = cs2_i * (a_i / a)**2
        drag = drag_i * (a_i / a)**(-drag_power) * H_i / float(bg.H_GeV(a)) if drag_i else 0.0  # Gamma(a)/H(a), Gamma ~ n_core v ~ a^-4
        ddc = -fric * dcp + src
        ddg = -fric * dgp + src - (k * math.sqrt(cs2) / aH)**2 * dg
        if drag_i:
            ddg += -drag * (dgp - dcp)
            ddc += -drag * f / (1 - f) * (dcp - dgp)
        if never_clusters and gas_on:
            ddg = 0.0; dgp = 0.0
        return [dcp, ddc, dgp, ddg]

    y0_c = [1.0, slope, 0.0, 0.0]
    if gas_init == 'clustered':
        y0_mix = [1.0, slope, 1.0, slope]
    else:
        y0_mix = [1.0, slope, 0.0, 0.0]
    N_eq = math.log(tc.a_eq / a_i)
    meth = 'Radau' if drag_i else 'DOP853'
    sol_mix = solve_ivp(rhs, (0, N_end), y0_mix, args=(True,), method=meth, rtol=1e-8, atol=1e-12, dense_output=True)
    sol_c = solve_ivp(rhs, (0, N_end), y0_c, args=(False,), method=meth, rtol=1e-8, atol=1e-12, dense_output=True)
    def dm(sol, N, gas):
        y = sol.sol(N)
        return (1 - f) * y[0] + f * y[2] if gas else y[0]
    ratio_0 = (dm(sol_mix, N_end, True) / dm(sol_c, N_end, False))**2
    ratio_eq = (dm(sol_mix, N_eq, True) / dm(sol_c, N_eq, False))**2
    yg = sol_mix.sol(N_end)
    return {'k': k, 'P_mix_over_P_cold_z0': ratio_0, 'P_mix_over_P_cold_eq': ratio_eq,
            'delta_g_over_delta_c_z0': yg[2] / yg[0], 'gas_init': gas_init, 'cs2_i': cs2_i, 'drag_i_over_H': drag_i}


def run_two_stage(k, cs2_i, a_dec, f=fgas):
    """Tight coupling (single fluid, c_eff^2 = f c_s^2) until a_dec, then two fluids with the gas clustered."""
    N_end = math.log(1 / a_i); N_dec = math.log(min(a_dec, 1.0) / a_i)
    aHi = aH_Mpc(a_i); a_k = a_i * aHi / k if k > aHi else a_i
    slope = 1.0 / (max(math.log(a_i / a_k), 0.0) + 1.0)

    def rhs1(N, y, feff):
        a = a_i * math.exp(N); aH = aH_Mpc(a); fric = 2 + dlnH_dlna(a)
        Om_a = tc.Om / a**3 / float(bg.E2(a)); d, dp = y
        cs2 = feff * cs2_i * (a_i / a)**2
        return [dp, -fric * dp + 1.5 * Om_a * d - (k * math.sqrt(cs2) / aH)**2 * d]

    def rhs2(N, y, gas_on):
        a = a_i * math.exp(N); aH = aH_Mpc(a); fric = 2 + dlnH_dlna(a)
        Om_a = tc.Om / a**3 / float(bg.E2(a)); dc, dcp, dg, dgp = y
        dm = (1 - f) * dc + f * dg if gas_on else dc
        src = 1.5 * Om_a * dm; cs2 = cs2_i * (a_i / a)**2
        return [dcp, -fric * dcp + src, dgp, -fric * dgp + src - (k * math.sqrt(cs2) / aH)**2 * dg]
    s1 = solve_ivp(rhs1, (0, N_dec), [1.0, slope], args=(f,), method='DOP853', rtol=1e-8, atol=1e-12)
    d, dp = s1.y[0, -1], s1.y[1, -1]
    s2 = solve_ivp(rhs2, (N_dec, N_end), [d, dp, d, dp], args=(True,), method='DOP853', rtol=1e-8, atol=1e-12)
    sc = solve_ivp(rhs2, (0, N_end), [1.0, slope, 0.0, 0.0], args=(False,), method='DOP853', rtol=1e-8, atol=1e-12)
    y = s2.y[:, -1]; dm = (1 - f) * y[0] + f * y[2]
    return {'k': k, 'a_dec': a_dec, 'P_mix_over_P_cold_z0': (dm / sc.y[0, -1])**2, 'delta_g_over_delta_c_z0': y[2] / y[0]}


def main():
    out = {'scope': __doc__.strip().splitlines()[0], 'a_i': a_i, 'aH_at_2MeV_Mpc^-1': aH_Mpc(a_i), 'f_gas': fgas}
    cs2_adiabatic = 5 / 9 * v_rms**2
    cs2_isothermal = v_rms**2 / 3
    out['sound_speeds'] = {'c_s_adiabatic_c': math.sqrt(cs2_adiabatic), 'c_s_isothermal_c': math.sqrt(cs2_isothermal),
                           'k_sound_2MeV_adiabatic_Mpc^-1': aH_Mpc(a_i) / math.sqrt(cs2_adiabatic),
                           'k_sound_2MeV_isothermal_Mpc^-1': aH_Mpc(a_i) / math.sqrt(cs2_isothermal)}
    # gas Jeans wavenumber history (fluid): k_J = aH sqrt(3 Om(a)/2)/c_s
    kJ = {}
    for a in (a_i, a_i * 10, 1e-6, 1e-5, 1e-4, tc.a_eq, 1e-3, 1e-2, 1e-1, 1.0):
        Om_a = tc.Om / a**3 / float(bg.E2(a))
        kJ[f'{a:.3e}'] = aH_Mpc(a) * math.sqrt(1.5 * Om_a) / (math.sqrt(cs2_isothermal) * a_i / a)
    out['gas_Jeans_k_isothermal_Mpc^-1'] = kJ
    rows = []
    ks = [3e4, 1e5, 3e5, 1e6, 1e7]
    for k in ks:
        for cs2, lab in ((cs2_adiabatic, 'adiabatic'), (cs2_isothermal, 'isothermal')):
            for init in ('erased', 'clustered'):
                r = run(k, cs2, init); r['cs_label'] = lab; rows.append(r)
                print(lab, init, k, round(r['P_mix_over_P_cold_z0'], 4), round(r['P_mix_over_P_cold_eq'], 4), flush=True)
    out['rows'] = rows
    plateau = run(1e6, cs2_isothermal, 'erased', never_clusters=True)
    out['never_clusters_plateau'] = plateau
    out['analytic_plateau_MD_only'] = {'(1-f)^2': (1 - fgas)**2,
                                       'times_(a_eq)^(6f/5)': (1 - fgas)**2 * tc.a_eq**(6 * fgas / 5),
                                       'growth_exponent_p': (-1 + math.sqrt(1 + 24 * (1 - fgas))) / 4}
    # gas-core elastic scattering (omitted by Codex): tight momentum coupling until Gamma_d/H = 1 (E)
    n_core_mX3 = 1.2141071422639381e-18   # at 2 MeV (Codex clock) [Q]
    drag_rows = []
    for sig in (0.01, 1., 100., 1000.):
        GH = n_core_mX3 * sig * v_rms * tc.mX_GeV / H_i
        a_dec = a_i * math.sqrt(GH)
        for k in (1e5, 3e5, 1e6):
            r = run_two_stage(k, cs2_isothermal, a_dec)
            r['sigma_Xcore_mX2'] = sig; r['Gamma_over_H_2MeV'] = GH; drag_rows.append(r)
            print('two-stage', sig, round(GH, 1), k, round(r['P_mix_over_P_cold_z0'], 4), flush=True)
    out['drag_two_stage'] = {'rows': drag_rows, 'note': 'Gamma_d/H at 2 MeV = n_core sigma v_rms/H with n_core(2 MeV) = 1.214e-18 m_X^3 (Codex clock), Gamma_d/H ~ a^-2 in RD; tight coupling modelled as one fluid with c_eff^2 = f c_s^2 until a_dec, then two fluids with the gas clustered; isothermal c_s. Order-of-magnitude only.'}
    def I_unit(a_start, a_stop, N=20000):
        aa = np.geomspace(a_start, a_stop, N)
        return float(np.trapezoid(a_i / (aa**2 * bg.H_s(aa)), np.log(aa)) * tc.c_km_s / tc.Mpc_km)
    mean_v = math.sqrt(8 / (3 * math.pi)) * v_rms
    fsd = {}
    for sig in (0.01, 0.1, 1., 10., 100., 1000., 2686.):
        GH = n_core_mX3 * sig * v_rms * tc.mX_GeV / H_i
        xd = math.sqrt(GH) if GH > 1 else 1.0
        ad = min(a_i * xd, tc.a_eq)
        fsd[f'sigma_{sig:g}'] = {'Gamma_over_H_2MeV': GH, 'x_dec': xd, 'T_dec_keV': float(bg.T(ad)) * 1e6, 'a_dec': ad,
                                 'lambda_fs_mean_pc_after_dec': mean_v * I_unit(ad, tc.a_eq) * 1e6 if ad < tc.a_eq else 0.0,
                                 'k_sound_mixture_2MeV_Mpc^-1': aH_Mpc(a_i) / math.sqrt(fgas * cs2_isothermal)}
    out['drag_decoupling_and_free_streaming'] = fsd
    out['codex_BGK'] = {'1e5': 0.893210, '3e5': 0.849806, '1e6': 0.764169, 'collisionless_1e5': 0.890969}
    (ROOT / 'results' / 'mixture_estimate.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(json.dumps({k: v for k, v in out.items() if k not in ('rows', 'drag_scan')}, indent=1, default=float))


if __name__ == '__main__':
    main()
