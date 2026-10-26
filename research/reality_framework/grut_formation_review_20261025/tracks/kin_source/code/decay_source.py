#!/usr/bin/env python3
"""kin_source (Claude, own code): the vector-decay source A -> X Xbar as an ODE in ln a, with its own spectrum binning.

Decay law [D]: dN_A/dt = -Gamma_A N_A with H(a) = H_i/a^2 (radiation, constant g*), so in ln a
    dN_A/dln a = -gamma a^2 N_A,   gamma = Gamma_A/H_i,   t H_i = a^2/2.
Each decay gives one X and one Xbar with rest-frame momentum p* = sqrt(m_A^2/4 - m_X^2); the parent is at rest to
better than p/m_A < 1e-6 at decay (checked), so each daughter is born with comoving momentum q = a p* and keeps it.
Hence the time profile and the comoving spectrum are the same function: dN_sign/dln q at q = a p* equals
dN_sign/dln a = gamma a^2 N_A(a).  State vector y = [N_A, N_sign, sum q, sum q^2] (q in H_i units).
Exact analytic reference: N_A(a) = N_A0 exp(-gamma (a^2 - 1)/2); the finite start a = 1 contributes exp(gamma/2) - 1
~ 1e-21 and is kept.  Cross-checks: an independent integration in cosmic time, and the Rayleigh moments.

Usage: python3 -B decay_source.py ../results/decay_source_own.json
"""
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kincommon as K

X_GRID = np.geomspace(1e-4, 10.0, 241)          # x = a/a_unit = q/q0
X_CODEX = np.geomspace(1e-4, 10.0, 161)         # Codex's sample grid (for the comparison only)
QUANTILES = (0.01, 0.1, 0.5, 0.9, 0.99)


def integrate_lna(vp, NA0, H, x_end=10.0):
    gam, pstar = vp['gamma'], vp['p_star_GeV']
    qa = pstar / H                                 # q/H per unit a

    def rhs(l, y):
        a = math.exp(l)
        r = gam * a * a * y[0]
        return [-r, r, r * a * qa, r * (a * qa)**2]

    l_end = math.log(x_end * vp['a_unit'])
    sol = solve_ivp(rhs, (0.0, l_end), [NA0, 0.0, 0.0, 0.0], method='DOP853', rtol=1e-12, atol=1e-30,
                    dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def integrate_time(vp, NA0, H, x_end=10.0):
    """Independent cross-check in cosmic time tau = t H_i (a = sqrt(2 tau)), from a = 1."""
    gam, pstar = vp['gamma'], vp['p_star_GeV']
    qa = pstar / H

    def rhs(tau, y):
        a = math.sqrt(2 * tau)
        r = gam * y[0]
        return [-r, r, r * a * qa]

    tau_end = 0.5 * (x_end * vp['a_unit'])**2
    sol = solve_ivp(rhs, (0.5, tau_end), [NA0, 0.0, 0.0], method='DOP853', rtol=1e-12, atol=1e-30)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[:, -1], sol.nfev


def build_case(label, e, NA0, norm, note):
    H, EX, NX = norm['H'], norm['E_X'], abs(norm['N_X'])
    vp = K.vector_params(e, H=H)
    a_unit, q0, gam, pstar = vp['a_unit'], vp['q0_over_H'], vp['gamma'], vp['p_star_GeV']
    sol = integrate_lna(vp, NA0, H)
    NA_end, Ns_end, Sq, Sq2 = sol.y[:, -1]
    y_t, nfev_t = integrate_time(vp, NA0, H)

    def profile(x):
        a = x * a_unit
        y = sol.sol(np.log(a))
        NA = y[0]
        rate_lna = gam * a * a * NA                 # dN_sign/dln a = dN_sign/dln q
        return a, NA, rate_lna, y[1]

    a_g, NA_g, dNdlna_g, cum_g = profile(X_GRID)
    q_g = a_g * pstar / H
    f_q = 2 * math.pi**2 * dNdlna_g / q_g**3        # comoving phase-space occupation per degree
    exact_NA = NA0 * np.exp(-gam * (a_g**2 - 1) / 2)
    rayleigh = NA0 * X_GRID**2 * np.exp(-X_GRID**2 / 2)
    # moments
    mean_q = Sq / Ns_end
    mean_q2 = Sq2 / Ns_end
    R = 2 * Sq / EX
    mean_a = mean_q * H / pstar
    # quantiles of the decay (number fraction decayed)
    quant = []
    for fq in QUANTILES:
        g = lambda l: sol.sol(l)[1] / NA0 - fq
        l_sol = brentq(g, 1e-9, math.log(10 * a_unit), xtol=1e-14)
        aq = math.exp(l_sol)
        quant.append(dict(number_fraction=fq, a_decay=aq, q_over_H=aq * pstar / H,
                          p_at_reference_over_m=aq * pstar / (K.A_REF * K.MX),
                          analytic_a=math.sqrt(1 - 2 * math.log1p(-fq) / gam)))
    # reference-coordinate numbers (collisionless redshift to a_ref)
    mean_p_ref = mean_q * H / K.A_REF
    # exact conformal energy of one daughter population at a_ref (massive) vs its massless sum
    E_exact_ref = float(np.trapezoid(dNdlna_g * np.sqrt(q_g**2 + (K.A_REF * K.MX / H)**2), np.log(q_g)))
    log_surv_ref = -gam * (K.A_REF**2 - 1) / 2
    pairs_per_net = NA0 / NX
    out = dict(
        label=label, note=note, e=e, normalisation=norm, N_A_over_H3=NA0,
        decay=dict(mA_GeV=vp['mA_GeV'], Gamma_A_GeV=vp['Gamma_A_GeV'], gamma=gam, p_star_GeV=pstar,
                   two_pstar_over_mA_minus_1=2 * pstar / vp['mA_GeV'] - 1, a_unit=a_unit, a_mean_decay=mean_a,
                   a_mean_analytic=vp['a_mean'], q0_over_H=q0, t_unit_over_Hinv=0.5 * a_unit**2,
                   T_vis_at_a_mean_GeV=K.T_vis_of_a(vp['a_mean'], H), T_vis_at_a_ref_GeV=K.T_vis_of_a(K.A_REF, H)),
        moments=dict(N_sign_end=Ns_end, N_A_remaining_end=NA_end, number_closure=(NA_end + Ns_end) / NA0 - 1,
                     mean_q_over_H=mean_q, mean_q_analytic=q0 * math.sqrt(math.pi / 2),
                     mean_q2_over_H2=mean_q2, mean_q2_analytic=2 * q0**2,
                     sum_q_per_sign=Sq, conformal_energy_both_signs=2 * Sq, R_from_energy_moment=R,
                     R_parent_mass_approx=NA0 * vp['mA_over_H'] * vp['a_mean'] / EX,
                     R_time_ode_crosscheck=2 * y_t[2] / EX, time_ode_nfev=nfev_t, lna_ode_nfev=sol.nfev,
                     exact_conformal_energy_one_sign_at_a_ref=E_exact_ref,
                     massive_correction_at_a_ref=E_exact_ref / Sq - 1,
                     pairs_per_net_charge=pairs_per_net, daughters_per_sign_per_net_charge=pairs_per_net,
                     pairs_per_visible_entropy=NA0 * K.YQ / NX,
                     f_rest_if_unprocessed=2 * K.MX * K.YQ * pairs_per_net / K.DM_PER_S,
                     mean_p_at_reference_GeV=mean_p_ref, mean_p_at_reference_over_m=mean_p_ref / K.MX,
                     energy_per_net_charge_at_ref_in_mX=2 * pairs_per_net * mean_p_ref / K.MX,
                     parent_log_survival_at_reference=log_surv_ref,
                     max_occupation_f=float(f_q.max())),
        quantiles=quant,
        table=dict(x=X_GRID, a=a_g, t_over_Hinv=0.5 * a_g**2, q_over_H=q_g, N_A_remaining=NA_g,
                   dN_sign_dlna=dNdlna_g, dN_sign_dt_over_H=gam * NA_g, dN_sign_dlnq=dNdlna_g,
                   dE_conf_sign_dlnq=dNdlna_g * q_g, dR_dlnq_both_signs=2 * dNdlna_g * q_g / EX,
                   f_q=f_q, cumulative_fraction=cum_g / NA0, p_over_mX_at_a_ref=q_g * H / (K.A_REF * K.MX),
                   check_NA_vs_exact=NA_g / exact_NA - 1, check_dNdlnq_vs_rayleigh=dNdlna_g / rayleigh - 1),
    )
    # comparison on Codex's grid
    a_c, NA_c, dN_c, cum_c = profile(X_CODEX)
    out['codex_grid'] = dict(x=X_CODEX, q_over_H=a_c * pstar / H, dN_sign_dlnq=dN_c, dR_dlnq=2 * dN_c * a_c * pstar / H / EX)
    return out


def main(out_path):
    codex_norm = dict(name='Codex 225028b primary window a=30-40', H=K.H_GEV, E_X=K.E_X, N_X=K.N_X, EPS=K.EPS)
    hw_norm = dict(name='our hw_gates release, window a=30-60', H=K.H_GEV_HWGATES, E_X=K.HW['E_X'], N_X=K.HW['N_X'], EPS=K.HW['EPS'])
    cases = [
        build_case('e1e-7_band', 1e-7, K.CODEX_BAND[1e-7], codex_norm,
                   "Codex finite band k/k*=0.2-5 (N_A quoted from 225028b VECTOR_RESULTS 'actual_fine'); conservative partial source."),
        build_case('e1e-7_broad', 1e-7, K.HW['N_A_P1'], codex_norm,
                   'Our hw_gates broader grid k=1e-5..1 H_i, profile P1 (own modes); R against the fixed E_X of the point.'),
        build_case('e1e-7_broad_hwnorm', 1e-7, K.HW['N_A_P1'], hw_norm,
                   'Same broad source against our own hw_gates E_X, N_X, H (the normalisation gate 1 used); sensitivity only.'),
        build_case('e1e-6_band', 1e-6, K.CODEX_BAND[1e-6], codex_norm,
                   "Codex finite band at e=1e-6 (N_A from 'e_sensitivity'); its own m_A, Gamma_A, decay epoch and q0."),
        build_case('e1e-6_broad', 1e-6, K.HW['N_A_P1_e1e6'], codex_norm,
                   'Our hw_gates broader grid at e=1e-6 (P1); same decay law at e=1e-6.'),
        build_case('cascade_vectors_e1e-7', 1e-7, K.N_YH, codex_norm,
                   'Vectors from Y_H -> Y_L + A (one per heavy quantum, N_YH late count [Q]); cold long before decay; same law.'),
        build_case('cascade_vectors_e1e-6', 1e-6, K.N_YH, codex_norm,
                   'Same cascade vectors at e=1e-6.'),
        build_case('sigma_vectors_e1e-7', 1e-7, 2 * K.N_SIG, codex_norm,
                   'Two vectors per Sigma radial quantum (sigma -> A A); negligible, carried for completeness.'),
    ]
    res = dict(method=__doc__, conventions=K.translation_table(), cases=cases)
    K.save(out_path, res)
    for c in cases:
        m = c['moments']
        print('%-26s e=%.0e N_A=%.7f pairs/net=%.7e R=%.5f (time-ODE %.5f) <q>=%.3f <p/m>_ref=%.2f f_rest=%.7f closure=%.1e'
              % (c['label'], c['e'], c['N_A_over_H3'], m['pairs_per_net_charge'], m['R_from_energy_moment'],
                 m['R_time_ode_crosscheck'], m['mean_q_over_H'], m['mean_p_at_reference_over_m'],
                 m['f_rest_if_unprocessed'], m['number_closure']))


if __name__ == '__main__':
    main(sys.argv[1])
