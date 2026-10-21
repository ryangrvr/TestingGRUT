"""Toy-scaled flat-FRW integration of q-theory with dust + radiation + a matter-sector
vacuum term V_m that drops by Delta at a transition (latent heat -> radiation), mu conserved.

Units: K(q0_post) = 1 (toy Planck units). Toy functions (supplied, illustrative only):
  eps(q) = 1/2 + q^2/2         -> eps'' = 1,   q0_post = 1, mu0_post = 1  (for V_m,after = 0)
  K(q)   = 1 + kappa (q - 1)   -> K' = kappa,  m_q^2 = (2/3) K eps''/K'^2
The hierarchy m_q/H is O(10-10^3) here instead of ~10^60; NOTES section 5.4 explains what
carries over (attractor, remnant formula) and what does not (oscillation amplitude).

First-order system (NOTES section 3.3), y = (q, v=qdot, H, rho_d, rho_r, N=ln a):
  qdot = v
  K' vdot = -K'' v^2 - 3 H K' v - (T - 4 rhoV - K R(q))/3 ... (see rhs) with R(q) = 2(eps'-mu)/K'
  Hdot = R(q)/6 - 2 H^2
  rho_d' = -3 H rho_d ;  rho_r' = -4 H rho_r - dV_m/dt ;  N' = H
Constraint C = 3H^2 K + 3 H K' v - (rho_d + rho_r + V_m) - (eps - mu q);  Cdot = -4 H C (derived).
"""
import json, os, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

KAPPA = 0.2
LAMB = 0.5


def eps(q): return LAMB + 0.5 * q * q
def eps1(q): return q
def Kf(q): return 1.0 + KAPPA * (q - 1.0)
def K1(q): return KAPPA
def K2(q): return 0.0


def Vm_of_N(N, Vb, Va, Nc, w):
    x = (N - Nc) / w
    return Va + (Vb - Va) * 0.5 * (1.0 - np.tanh(x))


def dVm_dN(N, Vb, Va, Nc, w):
    x = (N - Nc) / w
    th = np.tanh(x)
    return -(Vb - Va) * 0.5 * (1.0 - th * th) / w


def make_rhs(mu, Vb, Va, Nc, w):
    def rhs(t, y):
        q, v, H, rd, rr, N = y
        K, k1, k2 = Kf(q), K1(q), K2(q)
        Vm = Vm_of_N(N, Vb, Va, Nc, w)
        Vdot = dVm_dN(N, Vb, Va, Nc, w) * H
        rhoV = eps(q) - mu * q
        rho_m = rd + rr + Vm
        p_m = rr / 3.0 - Vm
        Tm = -rho_m + 3.0 * p_m
        Rq = 2.0 * (eps1(q) - mu) / k1
        Kdd = -3.0 * H * k1 * v + (rho_m - 3.0 * p_m + 4.0 * rhoV - K * Rq) / 3.0
        vdot = (Kdd - k2 * v * v) / k1
        Hdot = Rq / 6.0 - 2.0 * H * H
        return [v, vdot, Hdot, -3.0 * H * rd, -4.0 * H * rr - Vdot, H]
    return rhs


def constraint(y, mu, Vb, Va, Nc, w):
    q, v, H, rd, rr, N = y
    Vm = Vm_of_N(N, Vb, Va, Nc, w)
    C = 3 * H**2 * Kf(q) + 3 * H * K1(q) * v - (rd + rr + Vm) - (eps(q) - mu * q)
    scale = abs(rd) + abs(rr) + abs(Vm) + abs(eps(q) - mu * q) + 3 * H**2 * Kf(q) + 1e-300
    return C, C / scale


def quasi_static_q(mu, rd, rr, Vm, qguess):
    # trace equation with Box K -> 0:  T_m - 4 rhoV + K R(q) = 0,  T_m = -rd - 4 Vm
    f = lambda q: (-rd - 4 * Vm) - 4 * (eps(q) - mu * q) + Kf(q) * 2 * (eps1(q) - mu) / K1(q)
    return brentq(f, qguess - 0.5, qguess + 0.5, xtol=1e-15)


def static_root(mu, Vinf, qguess):
    # de Sitter / AdS static point: eps' - mu = 2 K' rhoV_tot / K, rhoV_tot = eps + Vinf - mu q
    f = lambda q: eps1(q) - mu - 2 * K1(q) * (eps(q) + Vinf - mu * q) / Kf(q)
    qs = brentq(f, qguess - 0.3, qguess + 0.3, xtol=1e-15)
    return qs, eps(qs) + Vinf - mu * qs, Kf(qs)


def run_case(name, mu, Vb, Va, Nc=0.8, w=0.1, rr0=0.3, rd0=0.01, t_end=3000.0, max_step=0.05):
    Vm0 = Vm_of_N(0.0, Vb, Va, Nc, w)
    q_i = quasi_static_q(mu, rd0, rr0, Vm0, 1.0)
    q_nodust = quasi_static_q(mu, 0.0, rr0, Vm0, 1.0)
    rhoV = eps(q_i) - mu * q_i
    H_i = np.sqrt((rd0 + rr0 + Vm0 + rhoV) / (3 * Kf(q_i)))
    v_i = -3.0 * H_i * (q_i - q_nodust)          # quasi-static drift of the dust-sourced shift (prop. to rho_d ~ a^-3)
    # re-solve the constraint for H with v_i:  3 K H^2 + 3 K' v H - rho = 0
    Kq, k1 = Kf(q_i), K1(q_i)
    rho_tot = rd0 + rr0 + Vm0 + rhoV
    H_i = (-3 * k1 * v_i + np.sqrt(9 * k1**2 * v_i**2 + 12 * Kq * rho_tot)) / (6 * Kq)
    y0 = [q_i, v_i, H_i, rd0, rr0, 0.0]
    rhs = make_rhs(mu, Vb, Va, Nc, w)

    def ev_crunch(t, y):  # stop when contracting strongly
        return y[2] + 2.0
    ev_crunch.terminal = True

    def ev_turn(t, y):
        return y[2]
    ev_turn.terminal = False
    sol = solve_ivp(rhs, (0.0, t_end), y0, method='DOP853', rtol=1e-11, atol=1e-14,
                    dense_output=False, max_step=max_step, events=[ev_crunch, ev_turn])
    ys = sol.y
    ts = sol.t
    Cs = np.array([constraint(ys[:, i], mu, Vb, Va, Nc, w)[1] for i in range(ys.shape[1])])
    res = {'name': name, 'mu': mu, 'V_before': Vb, 'V_after': Va, 'Nc': Nc, 'w': w,
           'status': int(sol.status), 'n_steps': int(len(ts)), 't_final': float(ts[-1]),
           'max_rel_constraint': float(np.max(np.abs(Cs))),
           'H_turnaround_times': [float(x) for x in sol.t_events[1]]}
    return res, ts, ys


def osc_frequency(ts, qs, q0, tmin):
    sel = ts > tmin
    tt, dq = ts[sel], qs[sel] - q0
    # remove slow drift by subtracting running mean is unnecessary at late times (drift ~ rho_d/m^2)
    idx = np.where(np.sign(dq[:-1]) != np.sign(dq[1:]))[0]
    if len(idx) < 4:
        return None
    tz = tt[idx] - dq[idx] * (tt[idx + 1] - tt[idx]) / (dq[idx + 1] - dq[idx])
    period = 2 * np.mean(np.diff(tz))
    return 2 * np.pi / period


def osc_energy_a3(ts, ys, mu, Vb, Va, Nc, w, m2, frac=0.5):
    """Einstein-frame oscillation energy (3/4) m^2 K'^2 A_q^2 times a^3, A_q from the late-time envelope of
    q - q_qs, where q_qs is the instantaneous quasi-static solution (removes the dust-sourced drift)."""
    q, v, H, rd, rr, N = ys
    sel = np.where(ts > frac * ts[-1])[0][::5]
    dq = np.array([q[i] - quasi_static_q(mu, rd[i], rr[i], Vm_of_N(N[i], Vb, Va, Nc, w), 1.0) for i in sel])
    a3 = np.exp(3 * N[sel])
    env = np.abs(dq) * np.sqrt(a3)               # A_q a^{3/2} should be constant for dust-like decay
    nblk = 20
    blocks = np.array_split(np.arange(len(sel)), nblk)
    amp = np.array([env[b].max() for b in blocks if len(b)])
    E = 0.75 * m2 * K1(1.0)**2 * amp**2
    return float(np.mean(E)), float(np.std(E) / np.mean(E))


def main():
    Delta = 0.02
    m2_post = (2.0 / 3.0) * Kf(1.0) * 1.0 / KAPPA**2
    out = {'model': {'eps': '1/2 + q^2/2', 'K': '1 + %.3g (q-1)' % KAPPA, 'Delta': Delta,
                     'm_q^2(q0_post)': m2_post, 'units': 'K(q0_post)=1'}}
    # equilibrium values
    q0_post, mu0_post = 1.0, 1.0
    q0_pre = np.sqrt(2 * (LAMB + Delta)); mu0_pre = q0_pre
    out['model'].update({'q0_pre': q0_pre, 'mu0_pre': mu0_pre, 'q0_post': q0_post, 'mu0_post': mu0_post,
                         'dmu0_dc_check (finite diff vs 1/q0)': [(mu0_pre - mu0_post) / Delta, 1 / q0_post]})
    cases = []
    # A: mu = mu0[post-transition = asymptotic eps], adiabatic-ish and sudden transitions
    for w, label, D in [(0.1, 'A0 control: no transition (Delta=0)', 0.0), (0.1, 'A1 post-transition mu, w=0.1', Delta),
                        (0.003, 'A2 post-transition mu, sudden w=0.003', Delta),
                        (0.5, 'A3 post-transition mu, slow w=0.5', Delta)]:
        res, ts, ys = run_case(label, mu0_post, D, 0.0, w=w, t_end=3000.0)
        q, v, H, rd, rr, N = ys
        a = np.exp(N)
        late = ts > 0.5 * ts[-1]
        rhoV_tot = eps(q) - mu0_post * q + Vm_of_N(N, D, 0.0, 0.8, w)
        Eosc, Eosc_relspread = osc_energy_a3(ts, ys, mu0_post, D, 0.0, 0.8, w, m2_post)
        res.update({'H_final': float(H[-1]), 'tH_final': float(ts[-1] * H[-1]),
                    'rhoV_tot_final': float(rhoV_tot[-1]),
                    'rhoV_tot_late_mean': float(np.mean(rhoV_tot[late])),
                    'rho_d_final': float(rd[-1]),
                    'osc_energy_a3_late (Einstein frame, from q envelope)': Eosc,
                    'osc_energy_a3_rel_spread_over_late_blocks': Eosc_relspread,
                    'pred_sudden_rho_osc_a3 = (4/3)Delta^2/(K m^2) a_c^3': float((4 / 3) * D**2 / m2_post * np.exp(3 * 0.8)),
                    'quasi_static_dust_piece_a3 ~ rho_d^2/(12 K m^2) a^3': float(np.mean((rd[late]**2 / (12 * m2_post)) * a[late]**3)),
                    'omega_measured': osc_frequency(ts, q, 1.0, 0.5 * ts[-1]),
                    'H_at_Nc': float(np.interp(0.8, N, H)),
                    'adiabatic_form_factor [x/sinh x]^2, x = pi m w/(2 H_c)': float(
                        (lambda x: (x / np.sinh(x))**2)(np.pi * np.sqrt(m2_post) * w / (2 * np.interp(0.8, N, H)))),
                    'measured_suppression_vs_sudden_formula': float(Eosc / ((4 / 3) * D**2 / m2_post * np.exp(3 * 0.8))) if D else None,
                    'm_q_predicted': float(np.sqrt(m2_post))})
        cases.append(res)
        print(res)
    # B: mu = mu0[pre-transition eps] -> remnant -Delta + ...: AdS-type, flat FRW recollapses
    res, ts, ys = run_case('B pre-transition mu (drop Delta)', mu0_pre, Delta, 0.0, w=0.1, t_end=3000.0)
    q, v, H, rd, rr, N = ys
    rhoV_tot = eps(q) - mu0_pre * q
    i_turn = np.argmin(np.abs(H)) if res['H_turnaround_times'] else -1
    qs_, rV_, Ks_ = static_root(mu0_pre, 0.0, q0_pre)
    K_pre = Kf(q0_pre); m2_pre = (2 / 3) * K_pre / KAPPA**2
    res.update({'static_AdS_root_rhoV_exact': rV_, 'static_root_q': qs_,
                'series_-Delta+(4/3)Delta^2/(K m^2)': -Delta + (4 / 3) * Delta**2 / (K_pre * m2_pre),
                'rhoV_tot_at_turnaround': float(rhoV_tot[i_turn]),
                'final_H': float(H[-1])})
    cases.append(res); print(res)
    # C: vacuum energy RISES by Delta with pre-transition mu -> de Sitter attractor; compare to static root
    res, ts, ys = run_case('C pre-transition mu, rise by Delta (dS)', 1.0, 0.0, Delta, w=0.1, t_end=1500.0, max_step=0.05)
    q, v, H, rd, rr, N = ys
    qs_, rV_, Ks_ = static_root(1.0, Delta, 1.0)
    late = ts > 0.8 * ts[-1]
    rhoV_tot = eps(q) - 1.0 * q + Delta
    res.update({'static_dS_root_rhoV_exact': rV_, 'static_root_q': qs_, 'static_H2': rV_ / (3 * Ks_),
                'series_Delta+(4/3)Delta^2/(K m^2)': Delta + (4 / 3) * Delta**2 / (1.0 * m2_post),
                'series_H2 = Delta/(3K0) (2nd order cancels)': Delta / 3.0,
                'rhoV_tot_late': float(np.mean(rhoV_tot[late])), 'q_late': float(np.mean(q[late])),
                'H2_late': float(np.mean(H[late]**2)), 'N_final': float(N[-1])})
    cases.append(res); print(res)
    out['cases'] = cases
    with open(os.path.join(os.path.dirname(__file__), '..', 'results', 'frw_transition_toy.json'), 'w') as f:
        json.dump(out, f, indent=1, default=float)


if __name__ == '__main__':
    main()
