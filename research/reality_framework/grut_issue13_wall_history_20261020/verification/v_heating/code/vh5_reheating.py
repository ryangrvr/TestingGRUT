#!/usr/bin/env python3
"""v_heating check 3b: inflaton-dominated (perturbative, matter-like) reheating, formal thermal hidden bound.

Background: d rho_I/dN = -3 rho_I - (Gam/H) rho_I ; d rho_R/dN = -4 rho_R + (Gam/H) rho_I (N = ln a, a_e = 1),
rho_I(0) = 3 MP^2 He^2, rho_R(0) = 0, SM bath instantly thermal with g* = 106.75.
Hidden: d rho_h/dN = -4 rho_h + betaE T^9/(MP^4 H).  Find Gam such that max_N T_h = T_c (T_h^4 = 30 rho_h/(pi^2 g_h)).
Report T_max (max of T over N) and T_RH (T at H = Gam).  Independent verifier code.
"""
import json, os
import numpy as np
from scipy import integrate, optimize

HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, '..', 'results')
MP = 2.435323e18; gstar = 106.75; gh = 4
Tc0 = 1547.646465
betaE = json.load(open(os.path.join(OUT, 'vh2_rates.json')))['beta_E_4dof']


def run(He, Gam, Tc=Tc0, Nmax=None):
    rhoe = 3 * MP**2 * He**2
    sc_h = Tc**4
    def rhs(N, y):
        lI, uR, uh = y
        uI = np.exp(lI)
        rI, rR = uI * rhoe, max(uR, 0.0) * rhoe
        H = np.sqrt((rI + rR) / 3) / MP
        T = (30 * rR / (np.pi**2 * gstar))**0.25
        return [-3 - Gam / H, -4 * uR + Gam / H * uI, -4 * uh + betaE * T**9 / (MP**4 * H) / sc_h]
    if Nmax is None:
        Nmax = (2 / 3) * np.log(He / Gam) + 8
    sol = integrate.solve_ivp(rhs, [0, Nmax], [0.0, 0.0, 0.0], method='LSODA', rtol=1e-10, atol=[1e-10, 1e-40, 1e-30],
                              dense_output=True)
    N = np.linspace(0, Nmax, 40001)
    lI, uR, uh = sol.sol(N)
    uI = np.exp(lI)
    rR = np.maximum(uR, 0) * rhoe
    T = (30 * rR / (np.pi**2 * gstar))**0.25
    Th = (30 * np.maximum(uh, 0) * sc_h / (np.pi**2 * gh))**0.25
    H = np.sqrt((uI + np.maximum(uR, 0)) * rhoe / 3) / MP
    iH = np.argmin(np.abs(np.log(H / Gam)))
    return dict(Tmax=float(T.max()), Thmax=float(Th.max()), TRH_H_eq_Gam=float(T[iH]),
                N_Tmax=float(N[np.argmax(T)]), N_Thmax=float(N[np.argmax(Th)]))


def solve(He, Tc=Tc0):
    f = lambda lg: np.log(run(He, 10**lg, Tc)['Thmax'] / Tc)
    lo, hi = -12.0, np.log10(He)
    lg = optimize.brentq(f, lo, hi, xtol=1e-6)
    r = run(He, 10**lg, Tc); r['Gamma'] = 10**lg; r['He'] = He
    return r


def main():
    rows = [solve(He) for He in (1e3, 1e5, 1e9, 1e13)]
    res = {'rows': rows}
    # fitted power of H_e for Tmax
    lH = np.log10([r['He'] for r in rows]); lT = np.log10([r['Tmax'] for r in rows])
    res['dlogTmax_dlogHe'] = float(np.polyfit(lH, lT, 1)[0])
    # v scaling at He = 1e9 (v = 10 TeV: Tc x 10; but hidden scales are only in Tc here)
    r10 = solve(1e9, Tc=10 * Tc0)
    res['Tmax_v10_over_v1_at_He1e9'] = r10['Tmax'] / rows[2]['Tmax']
    res['expected_10^(4/9)'] = 10**(4 / 9)
    with open(os.path.join(OUT, 'vh5_reheating.json'), 'w') as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
