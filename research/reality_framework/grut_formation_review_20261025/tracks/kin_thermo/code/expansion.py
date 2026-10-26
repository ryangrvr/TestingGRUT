#!/usr/bin/env python3
"""kin_thermo (Claude, own code): expansion, sources and the hidden-sector energy equation with the source work.

Background: visible radiation, g* = 106.75 (task statement; the coarse table is available for the old-proxy
comparison), H = sqrt(pi^2 g*/90) T_v^2/M_Pl, a ∝ 1/T_v, the hidden sector decoupled (portal negligible [Q]).
Hidden charge is comoving: n_Q(a) = n_Q0 (a0/a)^3 = Y_Q s_vis(T_v).

Energy equation (per unit physical volume, GeV^5 on the right before dividing by H):
    d rho_gas/d ln a = -3 (rho_gas + P_gas) + [ S_parent + S_med - W_obj - W_mf ] / H
    S_parent = Gamma_A m_A n_A(a)              (vector decays; each decay gives one X and one Xbar of energy m_A/2)
    d n_A/d ln a = -3 n_A - Gamma_A n_A / H    (nonrelativistic parent, zero pressure)
    W_obj   = omega(Q) d(n_L Q)/dt + (E(Q) - Q omega) d n_L/dt     (object sector work; linear comparator drops the
                                                                    second term and uses omega = const)
    W_mf    = mean-field work when masses depend on T or a background field: at a self-consistent mean field
              (dOmega_grand/dphi_bar = 0) the Gibbs-Duhem identities hold with explicit (T, mu) derivatives only and
              the field energy V_eff(phi_bar, T) is part of rho; the tadpoles dP/dm^2 are supplied by thermo.species.
Total conservation: d[(rho_gas + rho_A + rho_L) a^3] = -P_gas d(a^3)   (objects and parents pressureless).
Checks here: (i) isentropic expansion of the equilibrium gas integrated with the susceptibility Jacobian conserves
s a^3 to 1e-9; (ii) the parent/daughter energy ledger closes; (iii) the object work identity against the E(Q) table.
"""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

import kt_common as C
import thermo as TH

RES = Path(__file__).resolve().parents[1] / 'results'
MX = C.MX

# Vector parent (Codex DECAY_SOURCE.json 9877518 [Q]; CONFIRMED by kin_source)
PARENTS = {
    'e1e-7_band': dict(e=1e-7, mA_GeV=3200390.8655600557, Gamma_A_GeV=2.1223251086667907e-10, N_A=2.0650977334710046),
    'e1e-7_broad': dict(e=1e-7, mA_GeV=3200390.8655600557, Gamma_A_GeV=2.1223251086667907e-10, N_A=2.4691556381365065),
    'e1e-6_band': dict(e=1e-6, mA_GeV=32003908.655600563, Gamma_A_GeV=2.1223251086668018e-07, N_A=0.6526449887956494),
    'e1e-6_broad': dict(e=1e-6, mA_GeV=32003908.655600563, Gamma_A_GeV=2.1223251086668018e-07, N_A=0.7712575811986991),
}
A_REF = 7046713943374.485        # case-L density coordinate (DECAY_SOURCE reference_a_from_charge)


def parent_source(a, key='e1e-7_band'):
    """Physical parent density n_A (GeV^3), decay rate per volume per sign (GeV^4) and energy injection (GeV^5)."""
    p = PARENTS[key]
    gam = p['Gamma_A_GeV'] / C.H_I
    nA = p['N_A'] * C.H_I**3 / a**3 * np.exp(-gam * (a * a - 1) / 2)
    return dict(n_A_GeV3=nA, dn_pm_dt_GeV4=p['Gamma_A_GeV'] * nA, dE_dt_GeV5=p['Gamma_A_GeV'] * p['mA_GeV'] * nA,
                p_star_GeV=np.sqrt(p['mA_GeV']**2 / 4 - MX**2), survival=np.exp(-gam * (a * a - 1) / 2))


def hubble_of_a(a, table=False):
    """H(a) with a = 1 at the end of inflation and instant reheating: T_v = T_RH/a for constant g*."""
    T_RH = (90 * C.MRED**2 * C.H_I**2 / (np.pi**2 * C.GSTAR_CONST))**0.25
    return C.hubble(T_RH / a, table)


def object_work(dnLQ_dt, dnL_dt, Q, omega_over_m, E_GeV):
    """W_obj (GeV^5) for given rates of change of the condensed charge n_L Q and object number n_L."""
    return omega_over_m * MX * dnLQ_dt + (E_GeV - Q * omega_over_m * MX) * dnL_dt


def drho_gas_dlna(rho_gas, P_gas, H, S_parent=0.0, S_med=0.0, W_obj=0.0, W_mf=0.0):
    return -3 * (rho_gas + P_gas) + (S_parent + S_med - W_obj - W_mf) / H


# -------------------------------------------------------------------------- checks
def check_isentropic(t0=12.0, u0=0.3, lna_span=3.0, with_s=True):
    """Integrate the equilibrium gas (pairs + fast s) through expansion at fixed comoving charge and energy
    equation d(rho a^3) = -P d(a^3), using the 2x2 Jacobian in (ln t, u); report the drift of s a^3."""
    def rhs(x, y):
        lt, u = y
        t = np.exp(lt)
        g = TH.gas(t, u, -u, 0.0, with_s=with_s)
        S = TH.susceptibility(t, u, -u, 0.0, fast_s=True, with_s=with_s)
        J = S['J3']
        # unknowns (dlnt, du) with du_- = -du_+ : columns t and (u_+ - u_-)
        Jq = np.array([J[0, 0] * t, J[0, 1] - J[0, 2]])
        Jr = np.array([J[3, 0] * t, J[3, 1] - J[3, 2]])
        rhs_vec = np.array([-3 * g['n_Q'], -3 * (g['rho'] + g['P'])])
        return np.linalg.solve(np.vstack([Jq, Jr]), rhs_vec)
    sol = solve_ivp(rhs, (0, lna_span), [np.log(t0), u0], method='DOP853', rtol=1e-11, atol=1e-13, dense_output=True)
    out = []
    for x in np.linspace(0, lna_span, 7):
        lt, u = sol.sol(x)
        g = TH.gas(np.exp(lt), u, -u, 0.0, with_s=with_s)
        out.append(dict(lna=float(x), t=float(np.exp(lt)), u=float(u), s_a3=float(g['s'] * np.exp(3 * x)),
                        nQ_a3=float(g['n_Q'] * np.exp(3 * x))))
    s0 = out[0]['s_a3']
    return dict(rows=out, max_rel_entropy_drift=max(abs(r['s_a3'] / s0 - 1) for r in out),
                max_rel_charge_drift=max(abs(r['nQ_a3'] / out[0]['nQ_a3'] - 1) for r in out), nfev=int(sol.nfev))


def check_parent_ledger(key='e1e-7_band'):
    """Integrate the parent decay in ln a from a = 1 to 30 a_unit: comoving daughters and conformal energy
    against the closed forms N_A and 2 N_A <q> (R = 40.5626 for the band)."""
    p = PARENTS[key]
    gam = p['Gamma_A_GeV'] / C.H_I
    au = gam**-0.5

    def rhs(lna, y):
        a = np.exp(lna)
        NA = y[0]
        d = gam * a * a * NA
        return [-d, d, d * a * p['mA_GeV'] / 2 / C.H_I * 2]      # N_A, N_+ (= N_-), conformal energy both signs /H_i^4
    sol = solve_ivp(rhs, (0, np.log(30 * au)), [p['N_A'], 0.0, 0.0], method='DOP853', rtol=1e-12, atol=1e-16)
    NA, Npm, Econf = sol.y[:, -1]
    R = Econf / C.E_X
    q0 = np.sqrt(p['mA_GeV']**2 / 4 - MX**2) * au
    R_closed = 2 * p['N_A'] * q0 * np.sqrt(np.pi / 2) / (C.H_I * C.E_X)
    return dict(key=key, N_A_left=float(NA), daughters_per_sign=float(Npm), closed_form=p['N_A'], R_ode=float(R),
                R_closed=float(R_closed), rel_err_number=float(Npm / p['N_A'] - 1), rel_err_R=float(R / R_closed - 1))


def check_object_work():
    d = json.loads((RES / 'objects_EQ.json').read_text())
    rows = sorted([r for r in d['rows'] if r['dQ_dOmega'] < 0], key=lambda r: r['Q'])
    Q = np.array([r['Q'] for r in rows]); E = np.array([r['E_GeV'] for r in rows]); Om = np.array([r['Omega'] for r in rows])
    # a population changing from (nL1, Q1) to (nL2, Q2): exact d(nL E) vs omega d(nL Q) + (E - Q omega) dnL
    # integrated along a straight path with the table's omega(Q)
    res = []
    for i in range(0, len(rows) - 6, 6):
        Q1, Q2 = Q[i], Q[i + 6]
        nL1, nL2 = 1.0, 1.3
        exact = nL2 * E[i + 6] - nL1 * E[i]
        xs = np.linspace(0, 1, 2001)
        Qp = Q1 + (Q2 - Q1) * xs; nLp = nL1 + (nL2 - nL1) * xs
        om = np.interp(Qp, Q, Om) * MX; Ep = np.interp(Qp, Q, E)
        integrand = om * (nLp * (Q2 - Q1) + Qp * (nL2 - nL1)) + (Ep - Qp * om) * (nL2 - nL1)
        approx = np.trapezoid(integrand, xs)
        linear = C.OMEGA_REF * MX * (nL2 * Q2 - nL1 * Q1)
        res.append(dict(Q1=float(Q1), Q2=float(Q2), exact_dE=float(exact), work_identity=float(approx),
                        rel_err=float(approx / exact - 1), linear_comparator=float(linear), linear_rel_err=float(linear / exact - 1)))
    return res


if __name__ == '__main__':
    print(check_isentropic())
    print(check_parent_ledger())
    print(check_parent_ledger('e1e-6_band'))
    for r in check_object_work():
        print(r)
    print(parent_source(3.14e10), hubble_of_a(A_REF), C.hubble(45.67), C.hubble(45.67, table=True))
