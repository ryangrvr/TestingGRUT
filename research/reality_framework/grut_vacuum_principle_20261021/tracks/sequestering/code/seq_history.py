"""Historic averages for global vacuum-energy sequestering on FRW histories.

Track: sequestering (grut_vacuum_principle_20261021).  Computed [C] numbers only.

Conventions (match Kaloper-Padilla arXiv:1309.6562, 1406.0711):
  signature (-,+,+,+);  reduced Planck mass M = (8 pi G)^(-1/2);
  global field equation   M^2 G^mu_nu = tau^mu_nu - (1/4) delta^mu_nu <tau>,
  <Q> = int d^4x sqrt(g) Q / int d^4x sqrt(g),   tau = tau^a_a = -rho + 3p.
  The residual enters Friedmann as an energy density Lambda_eff = <tau>/4:
      3 M^2 (H^2 + k/a^2) = rho + Lambda_eff .
  Tracing and averaging gives <R> = 0 (any matter content, any loop order of the
  protected sector).  For GR with a source tau_GR (no separate constant), M^2 R = -tau_GR,
  so <tau_GR>/4 = -M^2 <R>/4; the sequestering condition is <R> = 0.

Units: M = 1, H0 = 1 (H at the epoch labelled a = 1); densities in M^2 H0^2 (rho_crit0 = 3).

Part A  (toy, Lombriser 1805.05918 sec. 2 check): LambdaCDM expansion to a_c followed by a
        time-mirrored collapse; tau = -rho_m (radiation traceless); the observed Lambda is
        taken to BE Lambda_eff, so only matter enters <tau>.
Part B  (Kaloper-Padilla 'The End' arXiv:1409.7073 trigger): protected sector = matter +
        radiation + scalar phi with V = s*phi, s = m^3/(M H0^2); no separate constant
        (any constant is a shift of phi).  GR integration bang -> crunch, crunch cut at
        rho_cut; compute <R>, decompose <tau>/4, and solve <R> = 0 for Omega_k.
"""
import json
import math
import sys

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

# ---------------------------------------------------------------- inputs (stated)
H_LITTLE = 0.674                     # Planck 2018-like h (input, not predicted)
OM = 0.315                           # Omega_m today (input)
OR = 4.18e-5 / H_LITTLE**2           # photons + 3 massless nu (input)
M_PL_RED_EV = 2.435e27               # reduced Planck mass [eV]
H0_EV = 2.1332e-33 * H_LITTLE        # H0 [eV]
HUBBLE_TIME_GYR = 9.7779 / H_LITTLE  # 1/H0 [Gyr]
RHO_CUT_PLANCK = (M_PL_RED_EV / H0_EV) ** 2   # rho = M^4 in units M^2 H0^2
RHO_CUT_TEV = (1e12) ** 4 / (M_PL_RED_EV**2 * H0_EV**2)


# ================================================================= Part A
def part_a(om=OM, orad=OR):
    ol = 1.0 - om - orad
    E = lambda a: math.sqrt(orad / a**4 + om / a**3 + ol)

    def ratio(ac):
        t = quad(lambda a: 1.0 / (a * E(a)), 0.0, ac, limit=400)[0]
        v = quad(lambda a: a**2 / E(a), 0.0, ac, limit=400)[0]
        # <tau>/4 = -(1/4) rho_m0 * t / int a^3 dt  (mirror doubles both)
        return -0.25 * (3 * om) * t / v / (3 * ol)

    r1 = ratio(1.0)
    ac = brentq(lambda a: abs(ratio(a)) - 1.0, 0.2, 1.0)
    return {"Omega_m": om, "Lambda_eff_over_rhoLambda_if_mirror_collapse_at_a=1": r1,
            "a_c_where_|Lambda_eff|=rho_Lambda": ac,
            "sign_note": "negative: dust-dominated history gives Lambda_eff=-<rho_m>/4 (KP 1406.0711 eq.45)"}


# ================================================================= Part B
class History:
    def __init__(self, s, ok, om=OM, orad=OR, a_i=1e-7, h_sw=0.3, h_sw2=3.0, rho_s=1e12):
        self.s, self.ok, self.om, self.orad = s, ok, om, orad
        self.ophi = 1.0 - om - orad - ok
        self.a_i, self.h_sw, self.h_sw2, self.rho_s = a_i, h_sw, h_sw2, rho_s

    # densities
    def rho(self, a, phi, psi):
        return 3 * self.om / a**3 + 3 * self.orad / a**4 + 0.5 * psi**2 + self.s * phi

    def H2(self, a, phi, psi):
        return self.rho(a, phi, psi) / 3.0 + self.ok / a**2

    # ---- phase A: expansion, variable N = ln a; y = [phi, psi, t, I4, Jk, JV, Ia, IH2]
    def rhs_A(self, N, y):
        a = math.exp(N)
        phi, psi = y[0], y[1]
        h2 = self.H2(a, phi, psi)
        H = math.sqrt(max(h2, 1e-300))
        a3 = a**3
        return [psi / H, (-3 * H * psi - self.s) / H, 1 / H, a3 / H, a3 * psi**2 / H,
                -4 * self.s * phi * a3 / H, a / H, a3 * H]

    def initial(self, phi_in):
        a = self.a_i
        t = a**2 / (2 * math.sqrt(self.orad))
        psi = -self.s * t / 2.5                      # radiation-era attractor (KP 1409.7073 eq.11)
        I4 = a**3 * t / 2.5
        return [phi_in, psi, t, I4, 0.0, -4 * self.s * phi_in * I4, a * t / 1.5, 0.0]

    def rho_phi_today(self, phi_in):
        y0 = self.initial(phi_in)
        sol = solve_ivp(self.rhs_A, (math.log(self.a_i), 0.0), y0, method="DOP853",
                        rtol=1e-11, atol=1e-14)
        if not sol.success:
            return -1e9
        phi, psi = sol.y[0, -1], sol.y[1, -1]
        return 0.5 * psi**2 + self.s * phi

    def shoot(self):
        target = 3 * self.ophi
        g = lambda p: self.rho_phi_today(p) - target
        p0 = target / self.s
        lo, hi = p0 - 2.0, p0 + 2.0
        while g(lo) > 0:
            lo -= 2.0
        while g(hi) < 0:
            hi += 2.0
        return brentq(g, lo, hi, xtol=1e-13, rtol=1e-13)

    def run(self, rho_cuts=(RHO_CUT_PLANCK,)):
        phi_in = self.shoot()
        y0 = self.initial(phi_in)
        out = {"s": self.s, "Omega_k": self.ok, "Omega_phi": self.ophi, "phi_in": phi_in}

        # today: w0 and CPL-fitted wa over 0.5<a<1 (diagnostic only)
        Ngrid = np.linspace(math.log(0.5), 0.0, 41)
        solT = solve_ivp(self.rhs_A, (math.log(self.a_i), 0.0), y0, method="DOP853",
                         rtol=1e-11, atol=1e-14, dense_output=True)
        ys = solT.sol(Ngrid)
        kin = 0.5 * ys[1] ** 2
        V = self.s * ys[0]
        w = (kin - V) / (kin + V)
        ag = np.exp(Ngrid)
        A = np.vstack([np.ones_like(ag), 1 - ag]).T
        w0, wa = np.linalg.lstsq(A, w, rcond=None)[0]
        out.update({"w_today": float(w[-1]), "cpl_w0_fit": float(w0), "cpl_wa_fit": float(wa),
                    "t0_over_tH": float(solT.y[2, -1])})

        # phase A to H = h_sw
        evA = lambda N, y: math.sqrt(max(self.H2(math.exp(N), y[0], y[1]), 0)) - self.h_sw
        evA.terminal, evA.direction = True, -1
        solA = solve_ivp(self.rhs_A, (math.log(self.a_i), 50.0), y0, method="DOP853",
                         rtol=1e-11, atol=1e-14, events=evA)
        NA = solA.t[-1]
        a = math.exp(NA)
        phi, psi, t, I4, Jk, JV, Ia, IH2 = solA.y[:, -1]
        H = math.sqrt(self.H2(a, phi, psi))

        # phase B: time variable through turnaround
        s = self.s

        def rhs_B(tt, y):
            a, ad, phi, psi = y[0], y[1], y[2], y[3]
            Hh = ad / a
            a3 = a**3
            acc = -a / 6.0 * (3 * self.om / a**3 + 2 * 3 * self.orad / a**4 + 2 * psi**2 - 2 * s * phi)
            return [ad, acc, psi, -3 * Hh * psi - s, a3, a3 * psi**2, -4 * s * phi * a3, a, a3 * Hh**2]

        evB = lambda tt, y: y[1] / y[0] + self.h_sw2
        evB.terminal, evB.direction = True, -1
        evT = lambda tt, y: y[1]
        evT.terminal, evT.direction = False, -1
        yB0 = [a, a * H, phi, psi, I4, Jk, JV, Ia, IH2]
        solB = solve_ivp(rhs_B, (t, t + 200.0), yB0, method="DOP853", rtol=1e-11, atol=1e-14,
                         events=[evB, evT])
        if solB.status != 1:
            raise RuntimeError("no collapse within 200 Hubble times")
        t_turn = solB.t_events[1][0]
        a_max = solB.y_events[1][0][0]
        tB = solB.t[-1]
        a, ad, phi, psi, I4, Jk, JV, Ia, IH2 = solB.y[:, -1]
        fr_err = (ad / a) ** 2 - self.H2(a, phi, psi)

        # phase C: collapse, variable u = -ln a increasing
        def rhs_C(u, y):
            a = math.exp(-u)
            phi, psi = y[0], y[1]
            Hh = -math.sqrt(max(self.H2(a, phi, psi), 1e-300))
            a3 = a**3
            return [-psi / Hh, (3 * Hh * psi + s) / Hh, -1 / Hh, -a3 / Hh, -a3 * psi**2 / Hh,
                    4 * s * phi * a3 / Hh, -a / Hh, -a3 * Hh]

        evC = lambda u, y: self.rho(math.exp(-u), y[0], y[1]) - self.rho_s
        evC.terminal, evC.direction = True, 1
        yC0 = [phi, psi, tB, I4, Jk, JV, Ia, IH2]
        solC = solve_ivp(rhs_C, (-math.log(a), 200.0), yC0, method="DOP853", rtol=1e-11,
                         atol=1e-14, events=evC)
        uC = solC.t[-1]
        a_s = math.exp(-uC)
        phi_s, psi_s, t_s, I4, Jk, JV, Ia, IH2 = solC.y[:, -1]
        H_s = -math.sqrt(self.H2(a_s, phi_s, psi_s))
        rho_s = self.rho(a_s, phi_s, psi_s)
        kin_frac = 0.5 * psi_s**2 / rho_s

        # analytic stiff tail from rho_s to rho_cut (w = 1: a^3 H = const)
        C = a_s**3 * abs(H_s)
        tau_s = 1.0 / (3 * abs(H_s))
        cslope = psi_s / (3 * H_s)
        res_cuts = {}
        for rc in rho_cuts:
            L = math.log(rc / rho_s)
            I4t = I4 + a_s**3 * tau_s / 2
            Jkt = Jk + C * L
            JVt = JV - 4 * s * a_s**3 * tau_s * (phi_s / 2 - cslope / 4)
            ttot = t_s + tau_s
            Jm = -3 * self.om * ttot
            IH2t = IH2 + C * L / 6
            Iat = Ia + 0.75 * a_s * tau_s
            tau_avg = (Jm + Jkt + JVt) / I4t
            R_avg = -tau_avg
            # identity check: int a^3 R dt = 6[a^3H] + 6 int (k a - a^3 H^2) dt, k = -Omega_k
            a3H_f = -C
            a3H_i = self.a_i**3 * (math.sqrt(self.orad) / self.a_i**2)
            R_identity = (6 * (a3H_f - a3H_i) + 6 * (-self.ok * Iat - IH2t)) / I4t
            rhoL = 3 * self.ophi
            res_cuts[f"{rc:.3e}"] = {
                "ln(rho_cut/rho_s)": L,
                "Lambda_res_over_rho_DE0": tau_avg / 4 / rhoL,
                "parts_over_rho_DE0": {"matter": Jm / I4t / 4 / rhoL,
                                       "phi_kinetic_pre_tail": Jk / I4t / 4 / rhoL,
                                       "phi_kinetic_stiff_tail": C * L / I4t / 4 / rhoL,
                                       "phi_potential": JVt / I4t / 4 / rhoL},
                "R_avg_over_H0sq": R_avg,
                "R_avg_identity_check": R_identity,
                "four_volume_per_comoving_vol": I4t,
                "t_crunch_over_tH": ttot,
            }
        out.update({"t_turn_over_tH": t_turn, "a_max": a_max,
                    "t_turn_from_now_Gyr": (t_turn - out["t0_over_tH"]) * HUBBLE_TIME_GYR,
                    "rho_turn_scale_3H_age^2_proxy": 3.0 / t_turn**2,
                    "friedmann_drift_phaseB": fr_err, "kinetic_fraction_at_tail_start": kin_frac,
                    "cuts": res_cuts})
        return out


def R_avg(s, ok, rc=RHO_CUT_PLANCK):
    r = History(s, ok).run(rho_cuts=(rc,))
    return r["cuts"][f"{rc:.3e}"]["R_avg_over_H0sq"], r


def solve_ok(s, rc=RHO_CUT_PLANCK, lo=-0.95, hi=0.0):
    f = lambda ok: R_avg(s, ok, rc)[0]
    grid = np.linspace(lo, hi, 20)
    vals = []
    for g in grid:
        try:
            vals.append(f(g))
        except Exception as e:  # noqa: BLE001
            vals.append(np.nan)
    root = None
    for i in range(len(grid) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and vals[i] * vals[i + 1] < 0:
            root = brentq(f, grid[i], grid[i + 1], xtol=1e-6)
            break
    return root, list(zip([float(g) for g in grid], [float(v) for v in vals]))


def part_c(I4_total, label):
    """Phase-transition weight in the global version (KP 1406.0711 eqs. 49-52) [E inputs].

    After a sudden transition at t*, <V_after - V> = -DeltaV * int_0^t* a^3 dt / int a^3 dt.
    Radiation era, g_* changes ignored except in a* (order-of-magnitude only).
    """
    T0_eV = 2.348e-4
    gs0 = 3.91
    rhoDE0_eV4 = 3 * (1 - OM - OR) * (M_PL_RED_EV * H0_EV) ** 2
    out = {}
    for name, Tstar_eV, gs, dV_eV4 in [
        ("electroweak", 160e9, 106.75, (125.1e9) ** 2 * (246.2e9) ** 2 / 8.0),
        ("QCD", 0.15e9, 61.75, (0.2e9) ** 4),
    ]:
        a_s = T0_eV / Tstar_eV * (gs0 / gs) ** (1.0 / 3.0)
        t_s = a_s**2 / (2 * math.sqrt(OR))
        w = a_s**3 * t_s / 2.5 / I4_total
        out[name] = {"a_star": a_s, "DeltaV_over_rhoDE0": dV_eV4 / rhoDE0_eV4,
                     "volume_weight": w,
                     "residual_shift_over_rhoDE0": -dV_eV4 / rhoDE0_eV4 * w}
    out["history_used"] = label
    return out


if __name__ == "__main__":
    res = {"inputs": {"h": H_LITTLE, "Omega_m": OM, "Omega_r": OR,
                      "rho_cut_planck_units_M2H0^2": RHO_CUT_PLANCK,
                      "rho_cut_TeV_units_M2H0^2": RHO_CUT_TEV,
                      "rho_Lambda_obs_over_M^4": 3 * (1 - OM - OR) * (H0_EV / M_PL_RED_EV) ** 2}}
    res["partA"] = {"Om=0.315": part_a(0.315), "Om=0.308": part_a(0.308)}
    rows = []
    for s in [0.3, 0.6, 1.0, 1.5, 2.0, 3.0]:
        for ok in [0.0, -0.002, 0.0023]:
            r = History(s, ok).run(rho_cuts=(RHO_CUT_PLANCK, RHO_CUT_TEV, 1e12))
            rows.append(r)
            c = r["cuts"][f"{RHO_CUT_PLANCK:.3e}"]
            print(f"s={s:4.2f} Ok={ok:+.4f} w0={r['w_today']:+.3f} cpl(w0,wa)=({r['cpl_w0_fit']:+.3f},{r['cpl_wa_fit']:+.3f}) "
                  f"turn in {r['t_turn_from_now_Gyr']:.1f} Gyr  Lres/rhoDE={c['Lambda_res_over_rho_DE0']:+.3f} "
                  f"<R>/H0^2={c['R_avg_over_H0sq']:+.4f} (id {c['R_avg_identity_check']:+.4f})", flush=True)
    res["partB_grid"] = rows
    roots = {}
    for s in [0.6, 1.0, 1.5, 2.0]:
        root, scan = solve_ok(s)
        roots[str(s)] = {"Omega_k_root_Planck_cut": root, "scan": scan}
        print("s", s, "Omega_k root (<R>=0, Planck cut):", root, flush=True)
    res["partB_Ok_roots"] = roots
    for row in rows:
        if row["s"] == 1.0 and row["Omega_k"] == 0.0:
            I4 = row["cuts"][f"{RHO_CUT_PLANCK:.3e}"]["four_volume_per_comoving_vol"]
            res["partC_phase_transitions"] = part_c(I4, "s=1, Omega_k=0, Planck cut")
            print("partC", res["partC_phase_transitions"], flush=True)
    # slope-scaling summary: Lambda_res is linear in s (late dynamics scale-free under t->t/sqrt(s))
    summ = []
    for row in rows:
        if row["Omega_k"] != 0.0:
            continue
        cs = row["cuts"]
        summ.append({"s": row["s"], "w_today": row["w_today"], "cpl_w0": row["cpl_w0_fit"],
                     "cpl_wa": row["cpl_wa_fit"], "a_max": row["a_max"],
                     "turn_from_now_Gyr": row["t_turn_from_now_Gyr"],
                     "Lres_over_rhoDE0_planck": cs[f"{RHO_CUT_PLANCK:.3e}"]["Lambda_res_over_rho_DE0"],
                     "Lres_over_rhoDE0_TeV": cs[f"{RHO_CUT_TEV:.3e}"]["Lambda_res_over_rho_DE0"],
                     "Lres_over_rhoDE0_notail": cs[f"{1e12:.3e}"]["Lambda_res_over_rho_DE0"],
                     "Lres_over_3H_turn^2_planck": cs[f"{RHO_CUT_PLANCK:.3e}"]["Lambda_res_over_rho_DE0"]
                     * 3 * (1 - OM - OR) / (3.0 / row["t_turn_over_tH"] ** 2),
                     "R_avg_over_H0sq_planck": cs[f"{RHO_CUT_PLANCK:.3e}"]["R_avg_over_H0sq"]})
    res["summary_flat"] = summ
    path = sys.argv[1] if len(sys.argv) > 1 else "seq_history_results.json"
    with open(path, "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    print("wrote", path)
