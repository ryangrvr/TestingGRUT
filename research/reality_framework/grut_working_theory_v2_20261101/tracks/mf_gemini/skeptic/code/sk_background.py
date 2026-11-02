"""mf_gemini SKEPTIC: independent background for the proposal's chi roll and M1r's oscillator, slice conventions, K.

Own integrator (scipy LSODA in cosmic time, with a radiation fluid fed by a constant width Gamma when asked).
Code units: M_P = 1, H0 = sqrt(V0/3) = 1 for the proposal (V = lam/4 (chi^2 - v^2)^2, v^2 = 8/3, lam = 27/16).
M1r regression: quadratic oscillator m = 3, at rest at phi = sqrt6/3 with H = 1, Gamma = 0.5 from t = 0;
dq_heat's engine gives K = a^2 H = 1.98581 (a = 1 at the start) [Q dq_heat NOTES].

Outputs (results/sk_background.json):
  - end of inflation (eps_H = 1): H_end/H0, chi_end/v, m/H_end;
  - S1 (one e-fold before eps_H = 1, the inflaton pass's slice): H/H0 and m/H;
  - first rest point (chi_dot = 0, all energy potential, the analogue of M1r's oscillator start): chi/v, H, m/H;
  - the first crossing of chi = v vs M1r's first zero crossing: m/H at each (the phase at which a holder on
    (chi - v) or on phi_I switches off);
  - K = a^2 H/H_slice at late times for a sudden decay at t and for constant widths switched on at eps_H = 1,
    in three slice normalisations (eps_H = 1, S1, first rest point).
Run: python3 -B sk_background.py
"""
import json
import os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'sk_background.json')
V2 = 8.0 / 3.0
LAM = 27.0 / 16.0
M = np.sqrt(2 * LAM * V2)          # = 3


def Vp(c):
    return 0.25 * LAM * (c * c - V2)**2


def dVp(c):
    return LAM * (c * c - V2) * c


def rhs_factory(V, dV, Gam_of_t):
    # y = (lna, phi, pi, rho_r); pi = phi_dot
    def f(t, y):
        lna, p, q, rr = y
        H = np.sqrt((0.5 * q * q + V(p) + rr) / 3.0)
        G = Gam_of_t(t)
        return [H, q, -(3 * H + G) * q - dV(p), -4 * H * rr + G * q * q]
    return f


def run(V, dV, y0, t1, Gam_of_t=lambda t: 0.0, events=None, max_step=0.01, dense=True):
    return solve_ivp(rhs_factory(V, dV, Gam_of_t), (0, t1), y0, method='LSODA', rtol=1e-10, atol=1e-14,
                     events=events, max_step=max_step, dense_output=dense)


def H_of(y, V):
    return np.sqrt((0.5 * y[2]**2 + V(y[1]) + y[3]) / 3.0)


def m1r_regression():
    V = lambda p: 0.5 * 9.0 * p * p
    dV = lambda p: 9.0 * p
    y0 = [0.0, np.sqrt(6.0) / 3.0, 0.0, 0.0]
    s = run(V, dV, y0, 120.0, Gam_of_t=lambda t: 0.5, max_step=0.02)
    yf = s.y[:, -1]
    K = np.exp(2 * yf[0]) * H_of(yf, V)
    # first zero crossing and eps = 1 time without decay (pure oscillator)
    s0 = run(V, dV, y0, 5.0, max_step=0.002)
    t = np.linspace(0, 5.0, 50001)
    Y = s0.sol(t)
    Hs = np.sqrt((0.5 * Y[2]**2 + V(Y[1])) / 3.0)
    eps = 1.5 * Y[2]**2 / (0.5 * Y[2]**2 + V(Y[1]))
    i0 = np.where(np.sign(Y[1][1:]) != np.sign(Y[1][:-1]))[0][0]
    ie = np.where(eps >= 1.0)[0][0]
    return dict(K_late_Gamma0p5=float(K), t_final=120.0, rho_phi_frac_final=float(1 - yf[3] / (3 * H_of(yf, V)**2)),
                m_over_H_first_zero_crossing=float(3.0 / Hs[i0]), t_first_zero=float(t[i0]),
                m_over_H_first_eps1=float(3.0 / Hs[ie]), t_first_eps1=float(t[ie]))


def main():
    out = dict(units='M_P = H0 = 1 (proposal); M1r oscillator: H_i = 1 at rest')
    out['M1r_regression'] = m1r_regression()
    # ------- proposal: hilltop roll from chi = 1e-6 with the growing-mode velocity -------
    s_rate = (-3 + np.sqrt(9 + 4 * LAM * V2)) / 2
    c0 = 1e-6
    y0 = [0.0, c0, s_rate * c0, 0.0]

    def ev_end(t, y):
        H = H_of(y, Vp)
        return 0.5 * y[2]**2 / H**2 - 1.0
    ev_end.terminal = True
    ev_end.direction = 1
    s1 = run(Vp, dVp, y0, 400.0, events=ev_end, max_step=0.05)
    te = s1.t_events[0][0]
    ye = s1.y_events[0][0]
    He = H_of(ye, Vp)
    # S1: one e-fold before the end
    Ns1 = ye[0] - 1.0
    ts1 = brentq(lambda t: s1.sol(t)[0] - Ns1, 0, te)
    ys1 = s1.sol(ts1)
    Hs1 = H_of(ys1, Vp)
    # after the end (no decay): first crossing of v, first rest point
    s2 = solve_ivp(rhs_factory(Vp, dVp, lambda t: 0.0), (te, te + 60.0), ye, method='LSODA', rtol=1e-11, atol=1e-14,
                   max_step=0.005, dense_output=True)
    t = np.linspace(te, te + 60.0, 240001)
    Y = s2.sol(t)
    H = np.sqrt((0.5 * Y[2]**2 + Vp(Y[1])) / 3.0)
    icv = np.where(np.sign(Y[1][1:] - np.sqrt(V2)) != np.sign(Y[1][:-1] - np.sqrt(V2)))[0]
    irest = np.where(np.sign(Y[2][1:]) != np.sign(Y[2][:-1]))[0]
    tr = t[irest[0]]
    yr = s2.sol(tr)
    Hr = H_of(yr, Vp)
    # oscillation-averaged w over the first few periods, and late
    rho = 0.5 * Y[2]**2 + Vp(Y[1])
    p = 0.5 * Y[2]**2 - Vp(Y[1])
    w = []
    for a, b in zip(irest[0:-2:2][:8], irest[2::2][:8]):
        w.append(float(np.trapezoid(p[a:b], t[a:b]) / np.trapezoid(rho[a:b], t[a:b])))
    # sudden-decay K table (a = 1 at eps_H = 1): K(t) = a^2 H/H_end
    a = np.exp(Y[0] - ye[0])
    Ktab = {f'{x:g}': float(np.interp(te + x / He, t, a**2 * H / He)) for x in [0.5, 1, 2, 2.4, 4, 8, 16]}
    out['proposal'] = dict(
        N_from_chi1em6_to_end=float(ye[0]), H_end_over_H0=float(He), chi_end_over_v=float(ye[1] / np.sqrt(V2)),
        m_over_H_end=float(M / He), H_S1_over_H0=float(Hs1), m_over_H_S1=float(M / Hs1),
        first_v_crossing=dict(t_after_end_Hend=float((t[icv[0]] - te) * He), H_over_H0=float(H[icv[0]]),
                              m_over_H=float(M / H[icv[0]])),
        first_rest_point=dict(t_after_end_Hend=float((tr - te) * He), chi_over_v=float(yr[1] / np.sqrt(V2)),
                              H_over_H0=float(Hr), m_over_H=float(M / Hr),
                              G1_identity_H=float(abs(yr[1] - np.sqrt(V2)) * M / np.sqrt(6.0)),
                              a_over_a_end=float(np.exp(yr[0] - ye[0]))),
        second_rest_point_chi_over_v=float(s2.sol(t[irest[1]])[1] / np.sqrt(V2)),
        w_first_periods=w, K_sudden_decay_a1_at_end=Ktab,
        chi_min_after_end_over_v=float(Y[1].min() / np.sqrt(V2)))
    # ------- constant widths switched on at eps_H = 1: late K in three normalisations -------
    Kw = {}
    for lab, G in [('0.5H0', 0.5), ('1H0', 1.0), ('0.5Hend', 0.5 * He), ('1Hend', He), ('2Hend', 2 * He),
                   ('5Hend', 5 * He), ('0.1Hend', 0.1 * He)]:
        s3 = solve_ivp(rhs_factory(Vp, dVp, lambda tt, G=G: G), (te, te + max(60.0, 40.0 / G)), ye, method='LSODA',
                       rtol=1e-10, atol=1e-14, max_step=0.01)
        yf = s3.y[:, -1]
        Hf = H_of(yf, Vp)
        frac_left = float(1 - yf[3] / (3 * Hf**2))
        K_end = float(np.exp(2 * (yf[0] - ye[0])) * Hf / He)
        Kw[lab] = dict(Gamma_over_Hend=float(G / He), K_a1_at_eps1=K_end,
                       K_a1_at_S1=float(K_end * np.exp(2.0) * He / Hs1),
                       K_a1_at_first_rest=float(np.exp(2 * (yf[0] - yr[0])) * Hf / Hr),
                       field_energy_fraction_left=frac_left)
    out['proposal']['K_constant_width_from_eps1'] = Kw
    # Gamma (in H_end) for K = 1.38 / 1.405 on the eps_H = 1 normalisation
    def Kof(G):
        s3 = solve_ivp(rhs_factory(Vp, dVp, lambda tt: G), (te, te + max(60.0, 40.0 / G)), ye, method='LSODA',
                       rtol=1e-9, atol=1e-13, max_step=0.01)
        yf = s3.y[:, -1]
        return float(np.exp(2 * (yf[0] - ye[0])) * H_of(yf, Vp) / He)
    gs = {}
    for Kc in [1.3815, 1.4051]:
        try:
            gs[f'{Kc}'] = float(brentq(lambda lg: Kof(np.exp(lg)) - Kc, np.log(0.05), np.log(0.7), xtol=1e-3))
            gs[f'{Kc}'] = float(np.exp(gs[f'{Kc}']) / He)
        except ValueError:
            gs[f'{Kc}'] = None
    out['proposal']['Gamma_over_Hend_needed_for_Kcrit'] = gs
    with open(OUT, 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
