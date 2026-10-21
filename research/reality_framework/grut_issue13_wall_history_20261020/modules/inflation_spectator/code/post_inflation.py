"""Part 4: post-inflation evolution of a large phi condensate (|phi_end| >> v), estimates with citations.

(a) Quartic oscillations in an (outside-action) radiation era that starts instantly at the end of inflation.
    Units H_inf = 1, a(tau) = tau (a = 1 at tau = 1), chi = a phi:
        chi'' + lam chi^3 - lam v^2 tau^2 chi = 0,   chi(1) = y_end, chi'(1) = y_end (phi frozen, phi' = 0),
    exact for RD (a'' = 0; conformal). Count zero crossings of phi until the physical energy density
    E = phi_dot^2/2 + V(phi) falls below the barrier lam v^4/4 (afterwards Hubble friction forbids crossings).
    Derived estimate (NOTES section 4): N_cross ~= (2 sqrt(lam)/7.4163) (y_end^2/(v/H)) sqrt(1 + 2/(lam y_end^2))
    -> (2 sqrt2/7.4163) phi_end/v = 0.3814 phi_end/v  for H^2 >> lam phi_end^2.
(b) Sign sensitivity: phi_end -> final sign map, flip spacing Delta y_end, compared with the isocurvature
    difference between neighbouring Hubble patches, delta y ~ sqrt(2 n) / (2 pi).
(c) Lame-equation Floquet chart for x'' + (kappa^2 + q cn^2(T, m=1/2)) x = 0 (Greene-Kofman-Linde-Starobinsky
    hep-ph/9705347): q = g_X^2/lam_phi = 1 (X production) and q = 3 (phi self-resonance).
(d) Back-reaction time from conformal vacuum fluctuations and the nonthermal-restoration criterion.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipj, ellipk
from common import dump

LAM = 0.01
G2 = 0.01


def oscillate(y_end, vH, lam=LAM, tau_max_factor=4.0):
    f0 = y_end * (1.0 + 2.0 / (lam * y_end ** 2)) ** 0.25
    tau_trap_est = f0 / vH
    tmax = max(50.0, tau_max_factor * tau_trap_est)

    def rhs(t, u):
        return [u[1], -lam * u[0] ** 3 + lam * vH ** 2 * t * t * u[0]]

    def zero(t, u):
        return u[0]

    def trapped(t, u):
        phi = u[0] / t
        phidot = (u[1] - u[0] / t) / t ** 2
        E = 0.5 * phidot ** 2 + 0.25 * lam * (phi ** 2 - vH ** 2) ** 2
        return E - 0.25 * lam * vH ** 4 * (1 - 1e-12)
    trapped.terminal = True
    trapped.direction = -1
    period = 7.4163 / (np.sqrt(lam) * f0)
    sol = solve_ivp(rhs, (1.0, tmax), [y_end, y_end], method="DOP853", rtol=1e-10, atol=1e-12,
                    events=(zero, trapped), max_step=period / 20)
    n_cross = len(sol.t_events[0])
    final_sign = np.sign(sol.y[0, -1])
    t_trap = sol.t_events[1][0] if len(sol.t_events[1]) else np.nan
    return n_cross, final_sign, t_trap


def ncross_formula(y_end, vH, lam=LAM):
    return 2 * np.sqrt(lam) / 7.4163 * (y_end ** 2 / vH) * np.sqrt(1 + 2 / (lam * y_end ** 2))


def floquet_chart(q, kap2, m=0.5):
    K = ellipk(m)
    Tp = 2 * K
    n = len(kap2)

    def rhs(t, u):
        sn, cn, dn, ph = ellipj(t, m)
        w2 = kap2 + q * cn ** 2
        x1, p1, x2, p2 = u.reshape(4, n)
        return np.concatenate([p1, -w2 * x1, p2, -w2 * x2])
    u0 = np.concatenate([np.ones(n), np.zeros(n), np.zeros(n), np.ones(n)])
    sol = solve_ivp(rhs, (0, Tp), u0, method="DOP853", rtol=1e-11, atol=1e-13)
    x1, p1, x2, p2 = sol.y[:, -1].reshape(4, n)
    tr = x1 + p2
    mu = np.where(np.abs(tr) > 2, np.log(np.abs(tr) / 2 + np.sqrt(np.maximum(tr ** 2 / 4 - 1, 0))) / Tp, 0.0)
    return mu, Tp


def bands(kap2, mu):
    on = mu > 1e-7
    edges = []
    i = 0
    while i < len(on):
        if on[i]:
            j = i
            while j + 1 < len(on) and on[j + 1]:
                j += 1
            k = i + np.argmax(mu[i:j + 1])
            edges.append({"kappa2_lo": float(kap2[i]), "kappa2_hi": float(kap2[j]), "mu_max": float(mu[k]),
                          "kappa2_at_max": float(kap2[k])})
            i = j + 1
        else:
            i += 1
    return edges


def backreaction_T(kap2, mu, qcoef, eps):
    """Smallest T at which the fluctuation-induced frequency^2 of the condensate reaches eps * lam f0^2.
    Conformal-vacuum seeds |x_k|^2 = 1/(2k) per real component, k = kappa sqrt(lam) f0:
        <x^2> = (lam f0^2/(4 pi^2)) I(T),  I(T) = int kappa e^{2 mu(kappa) T} dkappa.
    X (q = 1):   2 g^2 <|X|^2> = 2 g^2 <x^2>  >= eps lam f0^2  <=>  (2 g^2/(4 pi^2)) I >= eps   (qcoef = 2 g^2)
    phi (q = 3): 3 lam <dphi^2>               >= eps lam f0^2  <=>  (3 lam/(4 pi^2)) I >= eps   (qcoef = 3 lam)
    Superhorizon modes amplified during inflation start far above the conformal vacuum, so this T is an
    upper estimate."""
    kap = np.sqrt(kap2)
    for T in np.arange(1.0, 3000.0, 0.5):
        val = qcoef / (4 * np.pi ** 2) * np.trapezoid(kap * np.exp(2 * mu * T), kap)
        if val >= eps:
            return float(T)
    return None


def main():
    out = {"module": "inflation_spectator", "part": "4 post-inflation (estimates)",
           "labels": "computed ODE/Floquet numbers; everything about the inflaton, reheating and the RD start is "
                     "OUTSIDE-ACTION preparation data; back-reaction/restoration statements are estimates"}
    # (a) crossings
    rows = []
    for y_end in (1.0, 1.5, 2.5):
        for HV in (1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0):
            vH = 1.0 / HV
            n, sg, tt = oscillate(y_end, vH)
            rows.append({"y_end": y_end, "H_over_v": HV, "phi_end_over_v": y_end * HV, "N_cross": n,
                         "N_cross_formula": ncross_formula(y_end, vH), "final_sign": int(sg), "tau_trap": tt})
    out["zero_crossings"] = rows
    # fit N ~ c phi_end/v for H >> sqrt(lam) phi
    big = [r for r in rows if r["phi_end_over_v"] >= 30]
    c = np.array([r["N_cross"] / r["phi_end_over_v"] for r in big])
    out["N_cross_per_phi_end_over_v"] = {"computed_mean": float(c.mean()), "computed_min": float(c.min()),
                                         "computed_max": float(c.max()), "derived_limit": 2 * np.sqrt(2) / 7.4163}
    # (b) sign sensitivity: dN_cross/dy_end from two nearby end values; flip spacing = 1/(dN/dy)
    sens = []
    coef = float(c.mean())
    for HV in (10.0, 30.0, 100.0, 300.0):
        vH = 1.0 / HV
        n1 = oscillate(1.5, vH)[0]
        n2 = oscillate(1.7, vH)[0]
        dNdy = (n2 - n1) / 0.2
        dy_iso = np.sqrt(2.0) / (2 * np.pi)   # rms difference of phi_end/H between patches one e-fold apart
        sens.append({"H_over_v": HV, "N_cross_at_y1.5": n1, "N_cross_at_y1.7": n2, "dN_dy_end": dNdy,
                     "flip_spacing_dy": (1.0 / dNdy) if dNdy > 0 else None,
                     "isocurvature_dy_one_efold": dy_iso,
                     "expected_flips_between_neighbour_patches": dy_iso * dNdy})
    out["sign_sensitivity"] = sens
    out["sign_randomisation_criterion"] = {
        "statement": "neighbouring Hubble patches (one e-fold apart) end in uncorrelated vacua when "
                     "dN_cross/dy_end * sqrt(2)/(2 pi) >= 1 with dN_cross/dy_end = coef * H/v",
        "coef_computed": coef,
        "H_over_v_threshold": float(2 * np.pi / (coef * np.sqrt(2.0))),
        "note": "purely adiabatic curvature perturbations only shift the local clock; on uniform-density slices "
                "they do not change the number of crossings at trapping. Only the isocurvature delta phi does."}
    out["N_cross_estimates"] = {"trapping_at_amplitude_v": 2 * np.sqrt(2) / 7.4163,
                                "trapping_at_amplitude_sqrt2_v": 2 / 7.4163,
                                "computed": coef}
    # (c) Floquet
    kap2 = np.linspace(1e-4, 4.0, 4001)
    fl = {}
    for q in (1.0, 3.0):
        mu, Tp = floquet_chart(q, kap2)
        fl["q=%g" % q] = {"bands": bands(kap2, mu), "Floquet_period_T": Tp}
        fl["q=%g" % q]["_mu"] = mu
    out["floquet"] = {k: {"bands": v["bands"], "Floquet_period_T": v["Floquet_period_T"]} for k, v in fl.items()}
    out["floquet"]["literature_recall_unverified"] = ("GKLS hep-ph/9705347 (abstract read): analytic solutions exist "
                                                      "for g^2/lam = 1 and 3; band edges/mu_max NOT read from the paper")
    # (d) back-reaction (in zero crossings = T / (2K))
    K = ellipk(0.5)
    br = {}
    for q, qc, label in ((1.0, 2 * G2, "X (q=1): 2 g^2 <|aX|^2> vs lam f0^2"),
                         (3.0, 3 * LAM, "phi self (q=3): 3 lam <(a dphi)^2> vs lam f0^2")):
        mu = fl["q=%g" % q]["_mu"]
        res = {}
        for eps in (0.1, 1.0):
            T = backreaction_T(kap2, mu, qc, eps)
            res["eps=%g" % eps] = {"T_br": T, "N_cross_br": (T / (2 * K)) if T else None}
        br[label] = res
    out["backreaction"] = br
    out["restoration_criterion"] = {
        "Hartree": "m_eff^2(phi=0) = -lam v^2 + 3 lam <dphi^2> + 2 g_X^2 <|X|^2>; restoration iff "
                   "3 <dphi^2> + 2 (g_X^2/lam) <|X|^2> > v^2, i.e. <dphi^2> > v^2/3 or g_X^2 <|X|^2> > lam v^2/2",
        "scaling": "after back-reaction both variances and phi_amp^2 redshift as a^-2 (conformal quartic regime), "
                   "so with eps_fl = (3 lam <dphi^2> + 2 g^2 <|X|^2>)/(lam phi_amp^2) roughly constant, restoration "
                   "persists until phi_amp^2 ~ v^2/eps_fl; the final Z2 breaking then starts from a "
                   "fluctuation-dominated state (estimate)",
        "hidden_thermal_alternative": {
            "T_h_over_phi_amp_if_thermalised": float((30 * LAM / (4 * np.pi ** 2 * 4)) ** 0.25),
            "phi_amp_over_v_for_T_h_eq_Tc_1.56v": float(1.5571 / (30 * LAM / (4 * np.pi ** 2 * 4)) ** 0.25),
            "note": "g_h = 4 real hidden d.o.f.; T_c/v = 1.5571-1.5641 from modules/thermal_potential (quoted)"},
        "gravitational_production_bound": {
            "rho_grav_upper_dimensional": "<~ H_inf^4",
            "T_equiv_over_H": float((30 / (np.pi ** 2 * 4)) ** 0.25),
            "H_over_v_below_which_no_thermal_restoration": float(1.5571 / (30 / (np.pi ** 2 * 4)) ** 0.25)},
    }
    dump("post_inflation.json", out)
    return out


if __name__ == "__main__":
    o = main()
    import json
    print(json.dumps({k: v for k, v in o.items() if k != "zero_crossings"}, indent=1, default=float)[:4000])
    for r in o["zero_crossings"]:
        print(r)
