"""Energy released by the collapse of a closed super-horizon minority domain after horizon re-entry.

Thin-wall (Nambu-Goto) spherical wall in a flat FRW background a ∝ tau^alpha (alpha = 1 radiation, 2 matter),
conformal time tau, comoving radius r(tau). Action S = -sigma ∫ d^3xi sqrt(-gamma) = -4 pi sigma ∫ dtau a^3 r^2 sqrt(1-r'^2).
Equation of motion [derived, NOTES §2]:
    r'' = -2 (1 - r'^2)/r - 3 (a'/a) r' (1 - r'^2)
With u = gamma r' (gamma = 1/sqrt(1-r'^2)):
    du/dtau = -2 sqrt(1+u^2)/r - 3 (a'/a) u,   dr/dtau = u/sqrt(1+u^2).
The conformal Hamiltonian is H_c = 4 pi sigma a^3 r^2 gamma = a E_phys, with dH_c/dtau = -dL/dtau|explicit, which gives
for eps = r^2 gamma (E_phys = 4 pi sigma a^2 eps):
    deps/dtau = -3 (a'/a) eps r'^2,   r' = -sqrt(1 - r^4/eps^2)      (regular at r -> 0).
The domain starts at rest at tau_i << r0 (super-horizon: frozen in comoving size). Horizon re-entry is defined by
a r0 = 1/H, i.e. tau_re = alpha r0 (H = alpha/(a tau)). The released energy is E_c = 4 pi sigma a_c^2 eps_c at r = 0, and

    E_c = kappa_E sigma R_re^2,   R_re = a(tau_re) r0 = 1/H_re,   kappa_E = 4 pi (tau_c/tau_re)^(2 alpha) eps_c/r0^2.

Also reported: kappa_E relative to the actual radius at tau_re (the wall has already started to move), and the
check tau_c(RD) = 1.867 r0 against network_sim (bubbles.json) and v_network (1.8670).
"""
import numpy as np
from scipy.integrate import solve_ivp
from mr_common import dump


def collapse(alpha=1.0, x=1e3, r_switch=0.02):
    r0 = 1.0
    ti = r0 / x
    hub = lambda t: alpha / t

    def f1(t, y):
        r, u = y
        g = np.sqrt(1.0 + u * u)
        return [u / g, -2.0 * g / r - 3.0 * hub(t) * u]

    ev1 = lambda t, y: y[0] - r_switch
    ev1.terminal = True
    ev1.direction = -1
    s1 = solve_ivp(f1, (ti, 50.0 * r0), [r0, 0.0], method="DOP853", rtol=1e-12, atol=1e-14, events=ev1,
                   dense_output=True)
    t1 = s1.t_events[0][0]
    r1, u1 = s1.y_events[0][0]
    eps1 = r1 * r1 * np.sqrt(1.0 + u1 * u1)

    def f2(t, y):
        r, e = y
        q = max(0.0, 1.0 - r ** 4 / e ** 2)
        return [-np.sqrt(q), -3.0 * hub(t) * e * q]

    ev2 = lambda t, y: y[0] - 1e-12
    ev2.terminal = True
    ev2.direction = -1
    s2 = solve_ivp(f2, (t1, t1 + 10.0), [r1, eps1], method="DOP853", rtol=1e-12, atol=1e-16, events=ev2)
    tc = s2.t_events[0][0]
    epsc = s2.y_events[0][0][1]
    tre = alpha * r0
    r_at_re = float(s1.sol(tre)[0]) if tre < t1 else None
    kE = 4 * np.pi * (tc / tre) ** (2 * alpha) * epsc / r0 ** 2
    # energy history: E(tau)/(4 pi sigma R_re^2) on the phase-1 solution
    ts = np.linspace(ti, t1, 400)
    rr, uu = s1.sol(ts)
    E = (ts / tre) ** (2 * alpha) * rr ** 2 * np.sqrt(1 + uu ** 2) / r0 ** 2
    return {"alpha": alpha, "x_r0_over_tau_i": x, "tau_c_over_r0": tc / r0, "tau_c_over_tau_re": tc / tre,
            "t_c_over_t_re (cosmic)": (tc / tre) ** (1 + alpha), "a_c_over_a_re": (tc / tre) ** alpha,
            "eps_c_over_r0sq": epsc / r0 ** 2, "kappa_E": kE, "kappa_E_over_4pi": kE / (4 * np.pi),
            "r(tau_re)/r0": r_at_re,
            "kappa_E_rel_actual_radius_at_re": (kE * (1.0 / r_at_re) ** 2 if r_at_re else None),
            "E_max_over_4pi_sigma_Rre2_before_switch": float(E.max()), "tau_at_Emax_over_r0": float(ts[np.argmax(E)])}


def main():
    out = {"module": "minority_relic", "part": "collapse energy factor (thin-wall NG sphere)",
           "labels": "derived equation, computed ODE (DOP853 rtol 1e-12)"}
    rows = []
    for alpha in (1.0, 2.0):
        for x in (1e2, 1e3, 1e4):
            rows.append(collapse(alpha, x))
    out["rows"] = rows
    # flat-space check: alpha = 0 -> tau_c = Gamma(1/4)^2/(4 sqrt(2 pi)) r0 = 1.3110 r0, eps conserved
    from math import gamma, pi, sqrt
    flat = collapse_flat()
    out["flat_check"] = {"tau_c_numeric": flat[0], "tau_c_exact": gamma(0.25) ** 2 / (4 * sqrt(2 * pi)),
                         "eps_c_over_r0sq": flat[1], "exact": 1.0}
    rd = [r for r in rows if r["alpha"] == 1.0 and r["x_r0_over_tau_i"] == 1e3][0]
    md = [r for r in rows if r["alpha"] == 2.0 and r["x_r0_over_tau_i"] == 1e3][0]
    out["adopted"] = {"kappa_E_RD": rd["kappa_E"], "kappa_E_MD": md["kappa_E"],
                      "C_E_sphere_RD = 3 kappa_E/(4 pi)": 3 * rd["kappa_E"] / (4 * np.pi),
                      "C_E_sphere_MD": 3 * md["kappa_E"] / (4 * np.pi),
                      "tau_c_RD_over_r0 (network_sim 1.867)": rd["tau_c_over_r0"]}
    dump("collapse.json", out)
    return out


def collapse_flat(r0=1.0):
    def f(t, y):
        r, u = y
        g = np.sqrt(1 + u * u)
        return [u / g, -2 * g / r]
    ev = lambda t, y: y[0] - 0.02
    ev.terminal = True
    s = solve_ivp(f, (0, 5), [r0, 0.0], method="DOP853", rtol=1e-12, atol=1e-14, events=ev)
    t1 = s.t_events[0][0]
    r1, u1 = s.y_events[0][0]
    e1 = r1 * r1 * np.sqrt(1 + u1 * u1)

    def f2(t, y):
        r, e = y
        return [-np.sqrt(max(0.0, 1 - r ** 4 / e ** 2)), 0.0]
    ev2 = lambda t, y: y[0] - 1e-12
    ev2.terminal = True
    s2 = solve_ivp(f2, (t1, t1 + 2), [r1, e1], method="DOP853", rtol=1e-12, atol=1e-16, events=ev2)
    return s2.t_events[0][0], s2.y_events[0][0][1]


if __name__ == "__main__":
    import json
    print(json.dumps(main(), indent=1))
