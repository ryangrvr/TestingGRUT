#!/usr/bin/env python3
"""v_heating check 3/4/6: freeze-in history in radiation domination, formal thermal bound, nonthermal Pi
criterion, thermalisation estimate, and direct-branching bound.  Independent verifier code.

Inputs: beta_E, per-dof pair rate and <T/E> come from vh2_rates.json (this verifier's own calculation).
"""
import json, os
import numpy as np
from scipy import integrate, optimize

HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, '..', 'results')
MP = 2.435323e18
gstar = 106.75
gh = 4
v = 1000.0
lp, gX2, lX, hD2, lD = 0.01, 0.01, 0.001, 1e-4, 0.001
Tc = 1.547646465 * v
cH = np.sqrt(np.pi**2 * gstar / 90)

r2 = json.load(open(os.path.join(OUT, 'vh2_rates.json')))
betaE = r2['beta_E_4dof']
G1 = r2['pair_rate_per_dof']['total']          # pairs per dof: dN_pairs/dtdV = G1 T^8/MP^4
TinvE = r2['mc_mean_T_over_E']


def TRH_max_analytic(Tc_, betaE_=betaE):
    C7 = (3 / 7) * (4 / 7)**(4 / 3)
    # T_h,max^4 = 10 beta C7 T_RH^7 /(pi^2 g_h c_H MP^3)
    return (Tc_**4 * np.pi**2 * gh * cH * MP**3 / (10 * betaE_ * C7))**(1 / 7)


def Th_of_T_ode(TRH):
    """integrate d(rho_h)/dt + 4 H rho_h = beta T^9/MP^4 in RD in ln a; return max T_h."""
    def rhs(N, y):
        T = TRH * np.exp(-N)
        H = cH * T**2 / MP
        rho = y[0]
        return [betaE * T**9 / MP**4 / H - 4 * rho]
    sol = integrate.solve_ivp(rhs, [0, 6], [0.0], rtol=1e-11, atol=1e-30, dense_output=True,
                              max_step=0.01)
    N = np.linspace(0, 6, 60001)
    rho = sol.sol(N)[0]
    Th = (30 * rho / (np.pi**2 * gh))**0.25
    i = np.argmax(Th)
    return Th[i], TRH * np.exp(-N[i]), Th[-1] / (TRH * np.exp(-6))


def main():
    res = {'inputs': {'betaE': betaE, 'pair_rate_per_dof': G1, 'T_over_E': TinvE, 'MP': MP, 'gstar': gstar,
                      'gh': gh, 'Tc': Tc, 'cH': cH}}
    # ---- (3) formal thermal bound, instantaneous reheating
    TR = TRH_max_analytic(Tc)
    res['TRH_max'] = TR
    res['TRH_max_Tc_pm30'] = [TRH_max_analytic(0.7 * Tc), TRH_max_analytic(1.3 * Tc)]
    res['TRH_max_vs_v'] = {f'{vv}TeV': TRH_max_analytic(1.547646465 * vv * 1e3) / 1 for vv in (0.1, 1, 10, 100)}
    res['scaling_exponent_v'] = float(np.log(TRH_max_analytic(1.547646465e4) / TR) / np.log(10))
    Thmax, Tat, xiT = Th_of_T_ode(TR)
    res['ode_Thmax_over_Tc'] = Thmax / Tc
    res['ode_T_at_max_over_TRH'] = Tat / TR
    res['analytic_T_at_max_over_TRH'] = (4 / 7)**(1 / 3)
    res['xi_T_late'] = xiT
    res['H_inst_reheat_at_TRHmax'] = cH * TR**2 / MP
    res['TRH_max_sensitivity_to_betaE_x1/0.7242'] = TRH_max_analytic(Tc, betaE / 0.72417) / TR
    res['Th_sensitivity_to_betaE_x1/0.7242'] = (1 / 0.72417)**0.25
    # ---- (4) nonthermal Pi from freeze-in gas
    Cavg = (3 * lp + 2 * gX2 + hD2) / 4
    gam_chi = 2 * G1 * TinvE                     # particles per dof times <T/E>
    A = 4 * Cavg * gam_chi / (3 * cH * MP**3)    # Pi(T) = A T^2 (T_RH^3 - T^3)
    res['Cavg'] = Cavg; res['gamma_chi'] = gam_chi
    # cross-check Pi formula by direct integration of sum 1/E over produced particles
    def Pi_direct(T, TRH):
        f = lambda Tp: (8 * G1 * Tp**8 / MP**4) * (TinvE / Tp) * (T / Tp)**2 / (cH * Tp**2 / MP * Tp)
        val = integrate.quad(f, T, TRH, epsrel=1e-12)[0]
        return Cavg * val          # Cavg * (sum over all 4 dof of n<1/E>)  [equal sharing]
    TRt = 1e13
    res['Pi_direct_vs_formula'] = Pi_direct(0.5 * TRt, TRt) / (A * (0.5 * TRt)**2 * (TRt**3 - (0.5 * TRt)**3))
    PiMaxCoef = (3 / 5) * (2 / 5)**(2 / 3)
    lv2 = lp * v**2
    TRH_Pi = (lv2 / (A * PiMaxCoef))**0.2
    res['TRH_Pi_threshold'] = TRH_Pi
    # Hubble-damped: exists T: A T^2 (TRH^3 - T^3) - lv2 - (cH T^2/MP)^2 > 0
    def margin(TRH):
        x = np.logspace(-12, 0, 240001) * TRH
        return np.max(A * x**2 * (TRH**3 - x**3) - lv2 - (cH * x**2 / MP)**2)
    TRH_PiH = optimize.brentq(margin, 1e13, 1e17, xtol=1e6, rtol=1e-10)
    res['TRH_Pi_hubble_damped'] = TRH_PiH
    # closed form for T << TRH:  max_T [A TRH^3 T^2 - cH^2 T^4/MP^2] = A^2 TRH^6 MP^2/(4 cH^2) = lv2
    res['TRH_Pi_hubble_damped_closed'] = (4 * cH**2 * lv2 / (A**2 * MP**2))**(1 / 6)
    res['T_at_hubble_damped_optimum'] = np.sqrt(A * TRH_PiH**3 * MP**2 / (2 * cH**2))
    # first-version bug (linear grid starting at T = 1e-4 TRH) missed this optimum and gave 3.9e16 GeV
    # Pi on the formal thermal boundary
    res['Pi_over_lv2_at_TRHmax'] = A * PiMaxCoef * TR**5 / lv2
    # ---- Gamma22/H for the freeze-in gas near TRH_max
    T = (4 / 7)**(1 / 3) * TR
    n_h = 8 * G1 * (TR**3 - T**3) * T**3 / (3 * cH * MP**3)       # all 4 dof, particles
    E = r2['mc_mean_E_over_T'] * T
    s = 2 * E**2
    M2 = 36 * lp**2
    sig = M2 / (16 * np.pi * s) / 2      # sigma = |M|^2/(16 pi s) x 1/2 (identical final), massless
    Gam = n_h * sig
    H = cH * T**2 / MP
    res['n_h_at_Thmax'] = n_h
    res['Gamma22_over_H_at_Thmax'] = Gam / H
    # ---- (6) direct branching
    def Bmax(Tmax, Tc_=Tc):
        x = (gh / gstar) * (Tc_ / Tmax)**4
        return x / (1 + x)
    res['Bh_max'] = {f'{Tm:.0e}': Bmax(Tm) for Tm in (1e4, 1e6, 1e8, 1e10, 1e12)}
    res['Bh_band_factors'] = [0.7**4, 1.3**4]
    # ODE check: inflaton decay with branching B into hidden radiation
    def bh_ode(B, Gam=1e-3):
        def rhs(t, y):
            rI, rR, rh, a = y
            Hh = np.sqrt((rI + rR + rh) / 3)
            return [-3 * Hh * rI - Gam * rI, -4 * Hh * rR + (1 - B) * Gam * rI, -4 * Hh * rh + B * Gam * rI,
                    Hh * a]
        sol = integrate.solve_ivp(rhs, [0, 5e4], [1.0, 0, 0, 1.0], rtol=1e-10, atol=1e-40, method='LSODA')
        return sol.y[2, -1] / sol.y[1, -1]
    res['Bh_ode_rho_ratio_over_B/(1-B)'] = bh_ode(1e-12) / (1e-12 / (1 - 1e-12))
    # kinetic-equilibrium estimate for a thermal hidden bath at T_h = T_c (order of magnitude)
    nth = gh * 1.2020569 / np.pi**2 * Tc**3
    s_th = 2 * (2.701 * Tc)**2 * 1.0       # <s> ~ 2 <E>^2 (1 - <cos>), <cos> = 0
    Gam_th = nth * M2 / (16 * np.pi * s_th) / 2
    # Gamma/H = 1 with H = cH (Tc/xi)^2/MP  -> xi
    xi_kin = np.sqrt(cH * Tc**2 / MP / Gam_th)
    res['xi_T_kinetic_estimate'] = xi_kin
    res['Tmax_selfconsistent_kinetic'] = Tc / xi_kin
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 'vh3_history.json'), 'w') as f:
        json.dump(res, f, indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()
