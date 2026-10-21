"""Item 1 and 3: zero-chemical-potential hidden-sector transition of the unchanged action.
For every scheme (oneloop, AE, Parwani, DR) and renormalisation scale (mu0/2, mu0, 2 mu0, v), with couplings either
RG-evolved from mu0 ('running') or held at their mu0 values ('fixed'), find the curvature zero T0, any coexisting broken
minimum, T_c, phi_c/T_c, spinodals, latent heat, the effective quartic and the perturbative expansion parameters.
Writes results/transition.json."""
import json
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from common import (BENCH, MU0, SCHEMES, SCHEME_DOC, Veff, JBSeries, JBTable, debye_coeffs, run_couplings, to_builtin,
                    PI, V_GEV)

warnings.filterwarnings("ignore")
OUT = Path(__file__).resolve().parent.parent / "results"
PHI = np.geomspace(1e-7, 2.0, 6000)
S = PHI**2


def broken_minima(V, T):
    """Return list of (phi_min, Delta W = W(phi_min)-W(0)) for local minima with phi > 0."""
    g = V.dW(S, T)
    idx = np.where((g[:-1] < 0) & (g[1:] >= 0))[0]
    out = []
    for i in idx:
        s_r = brentq(lambda s: V.dW(s, T)[0], S[i], S[i + 1], xtol=1e-18, rtol=1e-15, maxiter=500)
        out.append((np.sqrt(s_r), float(V.W(s_r, T)[0] - V.W(0.0, T)[0])))
    return out


def curv0(V, T):
    return 2.0 * V.dW(0.0, T)[0]


def m2_at(V, s, T, h=1e-9):
    """phi-direction curvature d^2V/dphi^2 = 2 W'(s) + 4 s W''(s)."""
    if s == 0.0:
        return curv0(V, T)
    w2 = (V.dW(s * (1 + 1e-6), T)[0] - V.dW(s * (1 - 1e-6), T)[0]) / (2e-6 * s)
    return 2 * V.dW(s, T)[0] + 4 * s * w2


def lambda_eff(V, T, h=1e-6):
    """Effective phi^4 coefficient lambda_eff (V ~ ... + lambda_eff/4 phi^4) at phi = 0: lambda_eff = 2 W''(0)."""
    return 2.0 * (V.dW(h, T)[0] - V.dW(0.0, T)[0]) / h


def analyse(p, mu, scheme, jb):
    V = Veff(p, mu, scheme, jb)
    T0 = brentq(lambda T: curv0(V, T), 1.2, 2.0, xtol=1e-15, rtol=1e-15)
    offs = np.concatenate([[0.0], np.geomspace(1e-10, 2e-2, 60)])
    exists = [len(broken_minima(V, T0 + o)) > 0 for o in offs[1:]]
    r = {"scheme": scheme, "mu_over_v": mu, "T0_curvature_zero_over_v": T0}
    pic = debye_coeffs(p)
    r["Pi_over_T2"] = {"phi": pic[0], "X_per_real_mode": pic[1], "D": pic[2]}
    if any(exists):
        # first order (perturbatively): upper spinodal T1 = largest T with a broken minimum
        last = max(i for i, e in enumerate(exists) if e)
        lo, hi = T0 + offs[1:][last], T0 + (offs[1:][last + 1] if last + 1 < len(offs) - 1 else 0.05)
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if len(broken_minima(V, mid)) > 0:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-14:
                break
        T1 = lo

        def dW(T):
            b = broken_minima(V, T)
            return min(x[1] for x in b)
        Tc = brentq(dW, T0 + 1e-15, T1, xtol=1e-15, rtol=1e-15)
        bm = min(broken_minima(V, Tc), key=lambda x: x[1])
        phic = bm[0]
        sb = phic**2
        # latent heat via envelope theorem: L = Tc * d/dT [W(sb,T) - W(0,T)] at fixed sb
        h = 1e-7
        dd = ((V.W(sb, Tc + h)[0] - V.W(0.0, Tc + h)[0]) - (V.W(sb, Tc - h)[0] - V.W(0.0, Tc - h)[0])) / (2 * h)
        L = Tc * dd
        rho_h = PI**2 / 30 * 4 * Tc**4
        # barrier: maximum of W between 0 and sb at Tc
        res_b = minimize_scalar(lambda x: -(V.W(x * x, Tc)[0] - V.W(0.0, Tc)[0]), bounds=(0, phic), method="bounded",
                                options={"xatol": 1e-12})
        r.update(order_perturbative="first", Tc_over_v=Tc, T1_upper_spinodal_over_v=T1, T0_lower_spinodal_over_v=T0,
                 phic_over_Tc=phic / Tc, phic_over_v=phic, latent_heat_over_v4=L, latent_heat_GeV4=L * V_GEV**4,
                 latent_heat_over_rho_hidden=L / rho_h, barrier_height_over_Tc4=-res_b.fun / Tc**4,
                 barrier_location_over_Tc=res_b.x / Tc,
                 m_phi_broken_at_Tc_over_Tc=np.sqrt(max(m2_at(V, sb, Tc), 0)) / Tc,
                 m_phi_symmetric_at_Tc_over_Tc=np.sqrt(max(curv0(V, Tc), 0)) / Tc)
    else:
        Tc = T0
        below = []
        for o in [1e-8, 1e-6, 1e-4, 1e-3, 1e-2]:
            b = broken_minima(V, T0 * (1 - o))
            below.append({"eps": -o, "n_minima": len(b), "phi_min_over_Tc": (b[0][0] / (T0 * (1 - o))) if b else None})
        r.update(order_perturbative="second (no coexisting minimum found for T0 < T <= T0+0.02 v on a phi grid from 1e-7 v)",
                 Tc_over_v=Tc, T1_upper_spinodal_over_v=T0, T0_lower_spinodal_over_v=T0, phic_over_Tc=0.0,
                 latent_heat_over_v4=0.0, latent_heat_over_rho_hidden=0.0, broken_minimum_below_Tc=below)
    # derivative of the phi=0 curvature with respect to T^2 (kappa), effective quartic, expansion parameters
    dT = 1e-6
    kappa = (curv0(V, Tc + dT) - curv0(V, Tc - dT)) / (2 * dT) / (2 * Tc)
    lam = lambda_eff(V, Tc)
    mX = np.sqrt(p["mx"] + pic[1] * Tc**2)
    mD = np.sqrt(p["md"] + pic[2] * Tc**2)
    if scheme == "oneloop":
        lam = None     # X is massless at phi = 0 without resummation: V has a |phi|^3 term, W''(0) is not defined
    r.update(kappa_dm2_dT2_at_Tc=kappa, kappa_over_leading_Pi_phi=kappa / pic[0],
             lambda_eff_at_phi0_Tc=lam, lambda_eff_over_lambda_phi=(lam / p["lp"]) if lam is not None else None,
             tree_mX2_at_mu_over_v2=p["mx"],
             m_X_debye_over_Tc=mX / Tc, m_D_thermal_over_Tc=mD / Tc,
             eps_X=max(p["c"], 4 * p["lx"]) * Tc / (4 * PI * mX),
             eps_D=max(p["d"], 3 * p["ld"]) * Tc / (4 * PI * mD))
    if scheme in ("AE", "Parwani"):
        M2 = p["mp"] + pic[0] * Tc**2
        r["phi_debye_mass_M_at_Tc_over_Tc"] = np.sqrt(max(M2, 0)) / Tc
        r["eps_phi_soft_at_Tc_3lpT_over_4piM"] = 3 * p["lp"] * Tc / (4 * PI * np.sqrt(M2)) if M2 > 0 else np.inf
    return r, V


def main():
    t0 = time.time()
    jb = JBSeries()
    out = {"conventions": {"units": "v = 1 (v = 1000 GeV); mu0 = g_X v = 0.1 v",
                           "schemes": SCHEME_DOC,
                           "rg": "running: all ten coefficients evolved from mu0 with the one-loop betas of "
                                 "grut_candidate_1_quantum_20261018; fixed: mu0 values kept, mu changed only in the CW logs"}}
    pic = debye_coeffs(BENCH)
    out["leading_high_T"] = {"Pi_phi_over_T2": pic[0], "check_lp/4+gX^2/6+hD^2/12": 0.01 / 4 + 0.01 / 6 + 1e-4 / 12,
                             "Pi_X_over_T2_per_real_mode": pic[1], "check_gX^2/12+lambda_X/3": 0.01 / 12 + 0.001 / 3,
                             "Pi_D_over_T2": pic[2], "check_hD^2/12+lambda_D/4": 1e-4 / 12 + 0.001 / 4,
                             "T_curvature_zero_over_v": np.sqrt(0.01 / pic[0])}
    # zero temperature check
    V0 = Veff(BENCH, MU0, "oneloop", jb)
    sv = brentq(lambda s: V0.dW(s, 0.0)[0], 0.8, 1.2, xtol=1e-16)
    out["zero_T_check"] = {"phi_min_over_v_numerical_V0+V1": np.sqrt(sv),
                           "strict_one_loop_archived": 1.00012193740,
                           "curvature_over_v2_numerical": m2_at(V0, sv, 0.0),
                           "strict_one_loop_curvature_archived": 0.020012778428,
                           "note": "numerical minimisation of truncated V0+V1 differs from the strict expansion at O(loop^2)"}
    runs = []
    mus = {"mu0/2": MU0 / 2, "mu0": MU0, "2mu0": 2 * MU0, "v": 1.0}
    for rg in ("running", "fixed"):
        for lab, mu in mus.items():
            p = run_couplings(mu) if rg == "running" else dict(BENCH)
            for sch in SCHEMES:
                r, _ = analyse(p, mu, sch, jb)
                r["scale_label"], r["rg"] = lab, rg
                runs.append(r)
                print("%-8s %-6s %-8s Tc/v=%.7f order=%s phic/Tc=%.4g L/rho=%.3g lam_eff/lp=%s" % (
                    rg, lab, sch, r["Tc_over_v"], r["order_perturbative"][:6], r["phic_over_Tc"],
                    r["latent_heat_over_rho_hidden"], r["lambda_eff_over_lambda_phi"]))
    out["runs"] = runs
    # evaluator cross-check with the quadrature table
    tab_path = OUT / "jb_table.npz"
    tab = JBTable() if not tab_path.exists() else None
    if tab is None:
        tab = JBTable.__new__(JBTable)
        d = np.load(tab_path)
        from common import jb_nonanalytic
        from scipy.interpolate import CubicSpline
        tab.grid, tab.jq = d["grid"], d["jq"]
        tab.R = tab.jq - jb_nonanalytic(tab.grid)
        tab.spl = CubicSpline(tab.grid, tab.R)
        tab.dspl = tab.spl.derivative()
        tab.ymin, tab.ymax = tab.grid[0], tab.grid[-1]
    xc = {}
    for sch in SCHEMES:
        a, _ = analyse(BENCH, MU0, sch, jb)
        b, _ = analyse(BENCH, MU0, sch, tab)
        xc[sch] = {"Tc_series": a["Tc_over_v"], "Tc_table": b["Tc_over_v"], "diff": a["Tc_over_v"] - b["Tc_over_v"],
                   "phic_over_Tc_series": a["phic_over_Tc"], "phic_over_Tc_table": b["phic_over_Tc"]}
    out["evaluator_crosscheck_mu0"] = xc
    # adopted range: resummed schemes (AE, Parwani) and DR, running couplings, four scales
    sel = [r for r in runs if r["scheme"] in ("AE", "Parwani", "DR") and r["rg"] == "running"]
    tcs = [r["Tc_over_v"] for r in sel]
    allres = [r["Tc_over_v"] for r in runs if r["scheme"] in ("AE", "Parwani", "DR")]
    ol = [r for r in runs if r["scheme"] == "oneloop"]
    dr0 = [r for r in runs if r["scheme"] == "DR" and r["rg"] == "running" and r["scale_label"] == "mu0"][0]
    kappa = dr0["kappa_dm2_dT2_at_Tc"]
    Tc = dr0["Tc_over_v"]
    lp = BENCH["lp"]
    tG = (3 * lp / (4 * PI)) ** 2 / (2 * kappa)
    # two-loop soft sunset estimate of the T_c shift: delta m^2 ~ 6 lp^2 T^2 L/(16 pi^2), L in [1, ln(Tc/m_G)]
    mG = 3 * lp * Tc / (4 * PI)
    Lmax = np.log(Tc / mG)
    dTc = [6 * lp**2 * Lx / (16 * PI**2) / (2 * kappa) for Lx in (1.0, Lmax)]
    out["adopted"] = {
        "Tc_over_v_range_resummed_running": [min(tcs), max(tcs)],
        "Tc_GeV_range_resummed_running": [min(tcs) * V_GEV, max(tcs) * V_GEV],
        "Tc_over_v_range_resummed_running_and_fixed": [min(allres), max(allres)],
        "Tc_over_v_central_AE_mu0": [r for r in runs if r["scheme"] == "AE" and r["rg"] == "running" and r["scale_label"] == "mu0"][0]["Tc_over_v"],
        "Tc_over_v_DR_mu0": Tc,
        "order_resummed_schemes": sorted(set(r["order_perturbative"].split(" ")[0] for r in sel)),
        "order_oneloop_unresummed": sorted(set(r["order_perturbative"].split(" ")[0] for r in ol)),
        "phic_over_Tc_oneloop_range": [min(r["phic_over_Tc"] for r in ol), max(r["phic_over_Tc"] for r in ol)],
        "latent_heat_over_rho_hidden_oneloop_range": [min(r["latent_heat_over_rho_hidden"] for r in ol),
                                                       max(r["latent_heat_over_rho_hidden"] for r in ol)],
        "shift_from_leading_estimate_percent": [100 * (min(tcs) / out["leading_high_T"]["T_curvature_zero_over_v"] - 1),
                                                100 * (max(tcs) / out["leading_high_T"]["T_curvature_zero_over_v"] - 1)],
        "ginzburg": {"definition": "eps_phi = 3 lambda_phi T/(4 pi m_phi(T)) = 1 with mean-field m^2 = 2 kappa Tc^2 |t|",
                     "kappa_DR": kappa, "t_G": tG, "t_G_band_O1_criteria": [tG / 36, 4 * tG],
                     "m_phi_at_Ginzburg_over_Tc": mG / Tc},
        "two_loop_soft_sunset_Tc_shift_estimate_fraction": dTc,
        "nonperturbative_order_statement": (
            "DR effective quartic at Tc is positive (lambda_eff/lambda_phi = %.4f); X and D are Debye-screened and "
            "perturbative (eps_X = %.4f, eps_D = %.4f).  The soft phi zero mode is a 3D real scalar with Z2 and a "
            "positive quartic: 3D Ising universality, continuous transition.  Any perturbative first-order signal "
            "below lies inside the Ginzburg region." % (dr0["lambda_eff_over_lambda_phi"], dr0["eps_X"], dr0["eps_D"])),
    }
    # oneloop first-order signal inside Ginzburg region?
    for r in ol:
        if r["order_perturbative"].startswith("first"):
            r["eps_phi_at_phic_broken"] = 3 * lp * r["Tc_over_v"] / (4 * PI * r["m_phi_broken_at_Tc_over_Tc"] * r["Tc_over_v"])
            r["(T1-T0)/Tc"] = (r["T1_upper_spinodal_over_v"] - r["T0_lower_spinodal_over_v"]) / r["Tc_over_v"]
            r["window_inside_ginzburg_t_G"] = r["(T1-T0)/Tc"] < tG
    out["adopted"]["oneloop_ir_sensitivity"] = [
        {"rg": r["rg"], "scale": r["scale_label"], "tree_mX2_over_v2": r["tree_mX2_at_mu_over_v2"],
         "order": r["order_perturbative"].split(" ")[0], "phic_over_Tc": r["phic_over_Tc"]} for r in ol]
    out["adopted"]["lambda_eff_over_lambda_phi_range"] = {
        sch: [min(r["lambda_eff_over_lambda_phi"] for r in runs if r["scheme"] == sch),
              max(r["lambda_eff_over_lambda_phi"] for r in runs if r["scheme"] == sch)] for sch in ("AE", "Parwani", "DR")}
    out["adopted"]["eps_phi_soft_at_Tc_range_AE_Parwani"] = [
        min(r["eps_phi_soft_at_Tc_3lpT_over_4piM"] for r in runs if r["scheme"] in ("AE", "Parwani")),
        max(r["eps_phi_soft_at_Tc_3lpT_over_4piM"] for r in runs if r["scheme"] in ("AE", "Parwani"))]
    out["seconds"] = time.time() - t0
    (OUT / "transition.json").write_text(json.dumps(to_builtin(out), indent=1))
    print(json.dumps(to_builtin(out["adopted"]), indent=1))
    print(json.dumps(to_builtin(out["evaluator_crosscheck_mu0"]), indent=1))
    print("time %.1fs" % out["seconds"])


if __name__ == "__main__":
    main()
