"""Independent re-integration of the Kaloper-Padilla linear-trigger history (review_physics track).

Different method from tracks/sequestering/code/seq_history.py: one integration in cosmic time t
with the second-order (Raychaudhuri) equation for H, so turnaround needs no branch switch; the
Friedmann constraint is monitored, not imposed.  Units M = 1, H0 = 1 (rho_crit0 = 3).

GR with source = matter + radiation + scalar phi, V = s*phi, no separate constant.
  tau = -rho_m + 0 (radiation) + phidot^2 - 4 s phi
  Sequestering consistency in this gauge: <tau> = 0  <=>  <R> = 0, with M^2 R = -tau.
  Lambda_res = <tau>/4.
The crunch is followed until rho = rho_s, then the w = 1 tail is added analytically:
  int a^3 tau dt (tail) = 2 int a^3 rho dt = 6 int a^3 H^2 dt = |a^3 H|_s * ln(rho_cut/rho_s).
"""
import json
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

h = 0.674
OM = 0.315
OR = 4.18e-5 / h**2
M_EV = 2.435e27
H0_EV = 2.1332e-33 * h
RHO_PL = (M_EV / H0_EV) ** 2          # M^4 in units M^2 H0^2
RHO_TEV = (1e12) ** 4 / (M_EV**2 * H0_EV**2)


def run(s, ok, rho_s=1e10, a_i=1e-6):
    """Phase 1: N = ln a with H from the Friedmann constraint (expansion) up to a = 1 and on to
    H = 0.5 (still expanding).  Phase 2: cosmic time t with Raychaudhuri for H, through
    turnaround to rho = rho_s.  Constraint violation is monitored in phase 2."""
    k = -ok
    oph = 1.0 - OM - OR - ok

    def parts(a, phi, psi):
        return 3 * OM / a**3, 3 * OR / a**4, 0.5 * psi**2 + s * phi, 0.5 * psi**2 - s * phi

    def Hc(a, phi, psi):
        rm, rr, rphi, _ = parts(a, phi, psi)
        return math.sqrt(max((rm + rr + rphi) / 3 - k / a**2, 1e-300))

    def rhsN(N, y):                      # y = [phi, psi, t, I1, Itau]
        a = math.exp(N)
        phi, psi = y[0], y[1]
        H = Hc(a, phi, psi)
        rm = 3 * OM / a**3
        tau = -rm + psi**2 - 4 * s * phi
        return [psi / H, (-3 * H * psi - s) / H, 1 / H, a**3 / H, a**3 * tau / H]

    def rhsT(t, y):                      # y = [a, H, phi, psi, I1, Itau]
        a, H, phi, psi = y[:4]
        rm, rr, rphi, pphi = parts(a, phi, psi)
        dH = -(rm + rr + rphi + rr / 3 + pphi) / 2 + k / a**2
        tau = -rm + psi**2 - 4 * s * phi
        return [a * H, dH, psi, -3 * H * psi - s, a**3, a**3 * tau]

    def y0(phi_in):
        t0 = a_i**2 / (2 * math.sqrt(OR))
        psi0 = -s * t0 / 2.5
        I1 = a_i**3 * t0 / 2.5
        Itau = (-3 * OM) * t0 + (-4 * s * phi_in) * I1
        return [phi_in, psi0, t0, I1, Itau]

    def to_today(phi_in):
        sol = solve_ivp(rhsN, (math.log(a_i), 0.0), y0(phi_in), method="DOP853", rtol=1e-11, atol=1e-16)
        if not sol.success:
            return None
        return sol.y[:, -1]

    def g(p):
        y = to_today(p)
        if y is None:
            return -1e9
        return 0.5 * y[1]**2 + s * y[0] - 3 * oph

    p0 = 3 * oph / s
    lo, hi = p0 - 0.3, p0 + 0.3
    while g(lo) > 0:
        lo -= 0.3
    while g(hi) < 0:
        hi += 0.3
    phi_in = brentq(g, lo, hi, xtol=1e-12, rtol=1e-12)
    yN = to_today(phi_in)
    t_now = yN[2]
    # continue in N until H < 0.5 (or a = 1e3)
    ev = lambda N, y: Hc(math.exp(N), y[0], y[1]) - 0.5
    ev.terminal = True
    sol = solve_ivp(rhsN, (0.0, math.log(1e3)), yN, method="DOP853", rtol=1e-11, atol=1e-16, events=ev)
    Nsw = sol.t[-1]
    phi, psi, t_sw, I1, Itau = sol.y[:, -1]
    a = math.exp(Nsw)
    yT = [a, Hc(a, phi, psi), phi, psi, I1, Itau]

    def ev_c(t, y):
        rm, rr, rphi, _ = parts(y[0], y[2], y[3])
        return (rm + rr + rphi) - rho_s
    ev_c.terminal = True
    ev_c.direction = 1
    solT = solve_ivp(rhsT, (t_sw, t_sw + 1e3), yT, method="DOP853", rtol=1e-12, atol=1e-20, events=ev_c)
    a, H, phi, psi, I1, Itau = solT.y_events[0][0]
    t_s = solT.t_events[0][0]
    rm, rr, rphi, _ = parts(a, phi, psi)
    rho = rm + rr + rphi
    constraint = 3 * (H**2 + k / a**2) / rho - 1
    kin_frac = 0.5 * psi**2 / rho
    C = abs(a**3 * H)
    amax = float(np.max(solT.y[0]))
    res = {"s": s, "Omega_k": ok, "phi_in": phi_in, "t_now_H0inv": t_now,
           "t_crunch_from_now_H0inv": t_s - t_now, "t_crunch_from_now_Gyr": (t_s - t_now) * 9.7779 / h,
           "ln_amax": math.log(amax), "friedmann_violation_at_rho_s": constraint,
           "kinetic_fraction_at_rho_s": kin_frac, "rho_DE0": 3 * oph}
    for name, rc in (("planck", RHO_PL), ("tev", RHO_TEV), ("none", rho_s)):
        tail = C * math.log(rc / rho_s)
        tau_avg = (Itau + tail) / I1
        res[name] = {"Lambda_res_over_rhoDE0": tau_avg / 4 / (3 * oph),
                     "Lambda_res_over_s_rhoDE0": tau_avg / 4 / (3 * oph) / s,
                     "R_avg_over_H0sq": -tau_avg,
                     "tail_fraction": tail / (Itau + tail)}
    return res


if __name__ == "__main__":
    out = {"note": "independent cosmic-time integration; compare tracks/sequestering Part B"}
    rows = []
    for s in (0.3, 1.0, 2.0, 3.0):
        rows.append(run(s, 0.0))
    out["flat"] = rows
    out["curv"] = [run(1.0, ok) for ok in (-0.5, -0.05, 0.0023)]
    # fit Lambda_res = m^3 M [A ln(rho_cut/(m^3 M)) + B] at s = 1 from planck and tev cutoffs
    r = rows[1]
    s = 1.0
    # in units M^2 H0^2: m^3 M = s, rho_DE0 = 3*0.685
    L_pl = r["planck"]["Lambda_res_over_rhoDE0"] * r["rho_DE0"] / s
    L_te = r["tev"]["Lambda_res_over_rhoDE0"] * r["rho_DE0"] / s
    x_pl, x_te = math.log(RHO_PL / s), math.log(RHO_TEV / s)
    A = (L_pl - L_te) / (x_pl - x_te)
    B = L_pl - A * x_pl
    out["fit_s1"] = {"A": A, "B": B, "L_planck_over_m3M": L_pl, "L_tev_over_m3M": L_te}
    path = sys.argv[1] if len(sys.argv) > 1 else "seq_trigger_check.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
