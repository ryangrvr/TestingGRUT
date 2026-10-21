"""v_minority: independent Nambu-Goto collapse of a super-horizon spherical wall in a ∝ tau^alpha.

Written from scratch by the verifier (does not import or copy modules/minority_relic).

Action  S = -4 pi sigma ∫ dtau a^3 r^2 sqrt(1 - r'^2)   (comoving radius r, conformal time tau)
EOM     d/dtau(a^3 r^2 gamma r') = -2 a^3 r / gamma
  with  u = gamma r':  du/dtau = -2 sqrt(1+u^2)/r - 3 (a'/a) u
Energy  E_phys = 4 pi sigma a^2 r^2 gamma   (sigma x proper area x Lorentz factor)
Plunge  (independent variable r):  eps = r^2 gamma,  d eps/d tau = -3 (a'/a) eps r'^2  (from dH_c/dtau = dH_c/dtau|explicit,
        H_c = sqrt((4 pi sigma a^3 r^2)^2 + p^2) = 4 pi sigma a^3 r^2 gamma);  r' = -sqrt(1 - r^4/eps^2).
Re-entry convention: physical radius of the frozen sphere equals the Hubble radius, a(tau_re) r0 = 1/H(tau_re)
  => tau_re = alpha r0 for a ∝ tau^alpha.
kappa_E = E_collapse / (sigma R_re^2), R_re = a(tau_re) r0.

Two independent integrations are compared: (A) (r,u) phase + (tau,eps)(r) plunge; (B) brute (r, r') second-order ODE
with eps extrapolated from r -> 0.
"""
import json, math, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gamma as Gfun

OUT = os.path.dirname(os.path.abspath(__file__))


def run(alpha, x=1e3, r0=1.0, r_sw=0.25):
    tau_i = r0 / x if alpha > 0 else 0.0
    aoa = (lambda t: alpha / t) if alpha > 0 else (lambda t: 0.0)

    def rhs1(t, z):
        r, u = z
        return [u / math.sqrt(1 + u * u), -2 * math.sqrt(1 + u * u) / r - 3 * aoa(t) * u]

    ev = lambda t, z: z[0] - r_sw
    ev.terminal = True
    ev.direction = -1
    s1 = solve_ivp(rhs1, [tau_i, 50 * r0], [r0, 0.0], method='DOP853', rtol=1e-12, atol=1e-14,
                   events=ev, dense_output=True)
    t_sw = s1.t_events[0][0]
    r_s, u_s = s1.y_events[0][0]
    eps_s = r_s ** 2 * math.sqrt(1 + u_s ** 2)

    def rhs2(r, z):
        t, eps = z
        rd = -math.sqrt(max(0.0, 1.0 - r ** 4 / eps ** 2))
        return [1.0 / rd, -3 * aoa(t) * eps * rd]

    s2 = solve_ivp(rhs2, [r_s, 0.0], [t_sw, eps_s], method='DOP853', rtol=1e-12, atol=1e-15, dense_output=True)
    tau_c, eps_c = s2.y[0, -1], s2.y[1, -1]

    out = dict(alpha=alpha, x=x, tau_c_over_r0=tau_c / r0, eps_c_over_r0sq=eps_c / r0 ** 2)
    a = (lambda t: t ** alpha) if alpha > 0 else (lambda t: 1.0)
    if alpha > 0:
        tau_re = alpha * r0
        kappa = 4 * math.pi * (a(tau_c) / a(tau_re)) ** 2 * eps_c / r0 ** 2
        out['tau_re_over_r0'] = tau_re / r0
        out['kappa_E'] = kappa
        out['kappa_E_over_4pi'] = kappa / (4 * math.pi)
        out['t_c_over_t_re_cosmic'] = (tau_c / tau_re) ** (alpha + 1)
        # trajectory diagnostics (phase 1 only covers r > r_sw; energy maximum is well before)
        tt = np.linspace(tau_i, t_sw, 200001)
        rr, uu = s1.sol(tt)
        E = (np.array([a(t) for t in tt]) / a(tau_re)) ** 2 * rr ** 2 * np.sqrt(1 + uu ** 2) / r0 ** 2  # /(R_re^2)
        k = int(np.argmax(E))
        out['E_max_over_4pi_sigma_Rre2'] = float(E[k])
        out['tau_at_E_max_over_r0'] = float(tt[k])
        if tau_re < t_sw:
            r_re, u_re = s1.sol(tau_re)
            out['r_at_tau_re_over_r0'] = float(r_re)
            out['E_at_tau_re_over_4pi_sigma_Rre2'] = float(r_re ** 2 * math.sqrt(1 + u_re ** 2))
            out['rdot_at_tau_re'] = float(u_re / math.sqrt(1 + u_re ** 2))
            out['kappa_E_relative_to_actual_radius_at_re'] = kappa / r_re ** 2
        # time when wall starts to move appreciably (r' = -0.1)
        evm = lambda t, z: z[1] / math.sqrt(1 + z[1] ** 2) + 0.1
        evm.terminal = True
        s3 = solve_ivp(rhs1, [tau_i, 50 * r0], [r0, 0.0], method='DOP853', rtol=1e-11, atol=1e-14, events=evm)
        out['tau_when_rdot_eq_-0.1'] = float(s3.t_events[0][0])
        # energy at the moment the wall is at r = 0.5 r0 (a "collapse onset" variant)
        evh = lambda t, z: z[0] - 0.5 * r0
        evh.terminal = True
        s4 = solve_ivp(rhs1, [tau_i, 50 * r0], [r0, 0.0], method='DOP853', rtol=1e-11, atol=1e-14, events=evh)
        th, (rh, uh) = s4.t_events[0][0], s4.y_events[0][0]
        out['E_at_r_half_over_4pi_sigma_Rre2'] = float((a(th) / a(tau_re)) ** 2 * rh ** 2 * math.sqrt(1 + uh ** 2))
        out['rest_part_at_r_half'] = float((a(th) / a(tau_re)) ** 2 * rh ** 2)
    return out


def run_brute(alpha, x=1e3, r0=1.0):
    """(B) second-order ODE in (r, v=r'), stop at r = 1e-4 r0; extrapolate eps=r^2 gamma (finite limit)."""
    tau_i = r0 / x
    def rhs(t, z):
        r, v = z
        return [v, -2 * (1 - v * v) / r - 3 * (alpha / t) * v * (1 - v * v)]
    res = []
    for rstop in [1e-2, 3e-3, 1e-3]:
        ev = lambda t, z: z[0] - rstop
        ev.terminal = True
        s = solve_ivp(rhs, [tau_i, 50], [r0, 0.0], method='Radau', rtol=1e-12, atol=1e-15, events=ev)
        t, (r, v) = s.t_events[0][0], s.y_events[0][0]
        res.append((rstop, t, r * r / math.sqrt(max(1e-300, 1 - v * v))))
    # final approach r' -> -1: tau_c ≈ t + r ; eps nearly constant
    rstop, t, eps = res[-1]
    return dict(alpha=alpha, tau_c_est=t + rstop, eps_at_rstop=[e for (_, _, e) in res],
                kappa_E=4 * math.pi * ((t + rstop) / (alpha * r0)) ** (2 * alpha) * eps)


if __name__ == '__main__':
    R = {}
    flat = run(0)
    flat['exact_tau_c'] = Gfun(0.25) ** 2 / (4 * math.sqrt(2 * math.pi))
    flat['rel_err'] = flat['tau_c_over_r0'] / flat['exact_tau_c'] - 1
    flat['eps_conservation_err'] = flat['eps_c_over_r0sq'] - 1.0
    R['flat'] = flat
    R['RD'] = {f'x={x:g}': run(1, x=x) for x in [1e2, 1e3, 1e4]}
    R['MD'] = {f'x={x:g}': run(2, x=x) for x in [1e2, 1e3, 1e4]}
    R['RD_switch_test'] = run(1, x=1e3, r_sw=0.1)
    R['brute_RD'] = run_brute(1)
    R['brute_MD'] = run_brute(2)
    k = R['RD']['x=1000']['kappa_E']
    R['C_E_sphere_RD'] = 3 * k / (4 * math.pi)
    R['C_E_sphere_MD'] = 3 * R['MD']['x=1000']['kappa_E'] / (4 * math.pi)
    with open(os.path.join(OUT, 'results_collapse.json'), 'w') as f:
        json.dump(R, f, indent=1)
    print(json.dumps(R, indent=1))
