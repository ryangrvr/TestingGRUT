#!/usr/bin/env python3
"""hw_gates gate 5 (Claude, own code): consequences of h = 5.33e12 GeV (and the spot couplings) for the light sector.

Object: our gp_stability fixed-charge profile solver (own code, reused) at g_X = 0.1, lambda_X = 1e-3, lambda_phi = 0.03,
v = 0.924089637 GeV, Q = 200000 (none of which changes): Omega, E, binding margin, int n^2.
Gauge sector at the spot point: m_A = 6 e h, screening length vs core radius, local screened energy int n^2/(72 h^2),
unscreened Coulomb (depends on e only; quoted from gp_stability [C, own earlier]), m_A vs E, Gamma(A -> X Xbar),
decay epoch, Z6 arithmetic (unchanged).
Survival: c6 = alpha6 M_red^2 at the spot couplings, Gamma t0 with the unit-coefficient X^6 rates [Q], and the
gauge-invariant Planck-suppressed operators that could shift the late alpha6 (Sigma* X^6/M^3, Sigma* Y X^3/M, ...).
Usage: python3 -B gate5_light.py <out.json>
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import simpson

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / 'gp_stability' / 'code'))
import spotmodel as S
import profile_fixed_charge as PF   # our gp_stability solver (own code, reused)

P = S.params(); H = S.H_SPOT
YEAR = 3.15576e7; T0 = 13.8e9 * YEAR
G_COH, G_PAIR = 3.91958726382564e-50, 1.0399701293581432e-46    # unit-coefficient d = 6 core-loss rates /s [Q xcharge, Codex 8986d34]


def main(out):
    t0 = time.time()
    sol = PF.solve(70.0, 3001)
    obs = PF.observables(sol, 70.0)
    Om = obs['Omega']; MX = PF.MX; V = PF.V; G = PF.G
    r = np.linspace(0, 70.0, 140001)
    F, Fp, Sx, Sp, q = PF.fields(sol, r)
    n_shell = 4 * np.pi * r * r * F * F
    Qc = np.cumsum(n_shell) * (r[1] - r[0]); Qc /= Qc[-1]
    r50 = float(np.interp(0.5, Qc, r)) / MX
    I4 = 4 * np.pi * simpson(r * r * F**4, x=r)
    N2 = Om**2 * V**4 * I4 / MX                   # int n^2 d^3x in GeV^3
    res = dict(method=__doc__, object=dict(Omega=Om, E_GeV=obs['E_GeV'], binding_margin_GeV=obs['binding_margin_GeV'],
                                          half_charge_radius_GeVinv=r50, int_n2_GeV3=N2, Q=200000, mod6=200000 % 6,
                                          gp_stability_values_Q=dict(Omega=0.86119513082, E=17064.5374544, margin=1417.2552931, int_n2=47241.833)))
    rows = {}
    for lab, hG, e in (('stated_h1e11', 1e11, 1e-7), ('spot_h5.33e12', P['h'] * H, 1e-7), ('spot_e5e-3', P['h'] * H, 5e-3), ('spot_e1e-2', P['h'] * H, 1e-2)):
        mA = 6 * e * hG
        GA = e**2 * mA / (48 * np.pi) * (1 - 4 * MX**2 / mA**2)**1.5
        Tdec = (90 / (np.pi**2 * S.GSTAR))**0.25 * np.sqrt(GA * S.MRED)
        rows[lab] = dict(h_GeV=hG, e=e, m_A_GeV=mA, screening_length_GeVinv=1 / mA, mA_times_r50=mA * r50,
                         dE_local_screened_GeV=N2 / (72 * hG**2), dE_local_over_margin=N2 / (72 * hG**2) / obs['binding_margin_GeV'],
                         unscreened_coulomb_GeV_scaled=3.4005e-7 * (e / 1e-7)**2, unscreened_over_margin=3.4005e-7 * (e / 1e-7)**2 / obs['binding_margin_GeV'],
                         mA_over_E=mA / obs['E_GeV'], Gamma_A_GeV=GA, tau_A_s=6.582e-25 / GA, T_vis_at_H_eq_Gamma_GeV=Tdec,
                         beta_string=2 / (36 * e * e), sigma_to_AA_open=np.sqrt(2) * hG > 2 * mA, mA_over_msigma=mA / (np.sqrt(2) * hG))
    res['gauge_sector'] = rows
    # survival at the spot couplings
    a6 = S.alpha6(P, H)
    rate = a6['c6']**2 * (G_COH + G_PAIR)
    Sig = P['h'] * H / np.sqrt(2); M, mu, kap = P['M'] * H, P['mu'] * H, P['kap']
    D = M**4 - 4 * (mu * Sig)**2
    X_core = 1.71
    res['survival'] = dict(alpha6_GeVm2=a6['alpha6_GeVm2'], c6=a6['c6'], beta6_GeVm2=a6['beta6_GeVm2'], r=a6['r'],
                           rate_per_s=rate, rate_times_age=rate * T0, c6_ceiling_1pct=np.sqrt(0.01 / (T0 * (G_COH + G_PAIR))),
                           margin_to_ceiling=np.sqrt(0.01 / (T0 * (G_COH + G_PAIR))) / a6['c6'],
                           ls_scaling_c6_Q=1250373.981068544,
                           planck_shifts_relative_alpha6={
                               'Sigma* X^6 / M_red^3 (dim 7, gauge invariant)': (Sig / S.MRED**3) / a6['alpha6_GeVm2'],
                               'Sigma* Y X^3 / M_red (dim 5; Y integrated out, mixes with kappa)': D / (mu * kap * M**2 * S.MRED),
                               'Y* X^3 |X|^2 / M_red^2 at core (relative kappa shift)': X_core**2 / S.MRED**2,
                               'Sigma* Y^2 |X|^2 / M_red (relative mu shift at core)': X_core**2 / (mu * S.MRED)},
                           charge_selection='U(1) gauge: Delta Q_X = 3a + 6b only; X^5-type operators absent; lowest X-only operator after Y, Sigma elimination is X^6 (Z6)',
                           companion_in_core=dict(beta6_X6_over_lamX_X4=abs(a6['beta6_GeVm2']) * X_core**2 / 1e-3,
                                                  alpha6_X6_over_lamX_X4=a6['alpha6_GeVm2'] * X_core**2 / 1e-3))
    Path(out).write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: v for k, v in res.items() if k != 'method'}, indent=1, default=float))
    print('done', time.time() - t0)


if __name__ == '__main__':
    main(sys.argv[1])
