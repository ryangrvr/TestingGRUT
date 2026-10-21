"""Cross-checks of minority_relic against verified packet numbers and closed forms. Writes results/checks.json."""
import math
import json
import os
import numpy as np
import mr_common as C
import production as P
import relic as R

R.load_collapse()


def rel(a, b):
    return abs(a - b) / abs(b)


def main():
    checks = []

    def add(name, value, ref, tol, src):
        ok = rel(value, ref) <= tol
        checks.append({"check": name, "value": value, "reference": ref, "rel_dev": rel(value, ref), "tol": tol,
                       "pass": bool(ok), "reference_source": src})

    # 1. fresh-start p_- against the verifier's independent backward-Kolmogorov values (v_inflation VERIFY §2)
    ver = {250.0: (2.412e-9, 5.278e-9), 335.0: (1.820e-4, 3.402e-4), 1000.0: (0.18381, 0.20755),
           2000.0: (0.33011, 0.34626)}
    L1ref = {1000.0: 6.32e-3, 501.0: 2.13e-3}
    for H, (a, b) in ver.items():
        r = P.production(H)
        add("p_minus fresh N50 at %g GeV" % H, r["p_minus_fresh_N50"], a, 3e-3, "v_inflation VERIFY.md §2 table")
        add("p_minus fresh N60 at %g GeV" % H, r["p_minus_fresh_N60"], b, 3e-3, "v_inflation VERIFY.md §2 table")
        if H in L1ref:
            add("Lambda_1 (eigenvalue) at %g GeV" % H, r["Lambda1_eig"], L1ref[H], 5e-3, "v_inflation corrected column")
    r = P.production(501.0)
    add("Lambda_1 (eigenvalue) at 501 GeV", r["Lambda1_eig"], 2.13e-3, 5e-3, "v_inflation corrected column")
    # 2. collapse
    d = json.load(open(os.path.join(C.RESULTS, "collapse.json")))
    add("NG sphere collapse tau_c/r0 (RD, x>>1)", d["adopted"]["tau_c_RD_over_r0 (network_sim 1.867)"], 1.867, 1e-3,
        "network_sim bubbles.json; v_network 1.8670")
    add("flat-space sphere collapse tau_c/r0", d["flat_check"]["tau_c_numeric"], d["flat_check"]["tau_c_exact"], 1e-9,
        "Gamma(1/4)^2/(4 sqrt(2 pi))")
    add("flat-space energy conservation eps_c/r0^2", d["flat_check"]["eps_c_over_r0sq"], 1.0, 1e-9, "exact")
    # 3. re-entry bookkeeping vs inflation_spectator (step g_*, h = 0.674) - expect small differences
    add("n'(1 MeV)", C.nprime_of_T(1e-3), 17.766, 2e-3, "inflation_spectator reentry_efolds 17.766 (step g_*)")
    add("n'(8.0 keV)", C.nprime_of_T(7.984e-6), 12.692, 2e-3, "inflation_spectator 12.692")
    for H, ref in ((150.0, 48.115), (300.0, 48.462), (1000.0, 49.064)):
        add("N_0 instant at %g GeV" % H, C.N0_instant(H), ref, 2e-3, "inflation_spectator Nstar_instant_reheating")
    add("T_RH instant at 300 GeV", C.T_RH_instant(300.0), 1.46e10, 1e-2, "v_inflation §4.5 (g = 106.75)")
    add("T_dom,iso sigma/(M_P^2 H) = 1", C.T_dom_iso(), 7.97e-6, 2e-3, "critic I2 (7.97 keV)")
    add("Y_N,max at m_X from rho/s max", C.RHO_OVER_S_MAX / C.M_X, 4.373e-12, 1e-3, "v_relic strict 4.373e-12")
    # 4. cosmo_observable CHC-3: constant per-ln-t minority fraction p, d rho/d ln t = c_b p sigma/t, f c_b p = 1:
    #    rho/s = 2 f c_b p sigma/(t s) at T_end with t = 1/(2H)  ->  T_end = 60.3 MeV for rho/s = 4.349e-10 GeV.
    from scipy.optimize import brentq
    for fcp, ref, lo, hi in ((1e-6, 126e-9, 1e-8, 1e-5), (1e-3, 0.107e-3, 1e-6, 1e-2), (1.0, 0.0603, 1e-4, 1.0)):
        g = lambda lT: math.log(4 * fcp * C.SIGMA * C.H_of_T(math.exp(lT)) / C.s_of_T(math.exp(lT))) - math.log(4.349e-10)
        Tend = math.exp(brentq(g, math.log(lo), math.log(hi)))
        add("CHC-3 T_end for f c_b p_- = %g (t = 1/(2H) at T_end, my background)" % fcp, Tend, ref,
            0.03 if fcp < 1 else 0.2, "cosmo_observable NOTES §7 table; at 60 MeV their t(T) includes the g_* history "
            "(QCD/muon steps), so t != 1/(2H) there and a ~14% offset is expected")
    # identity: my rho/s = f C_E sigma int (dp/dn) H/s dn'; with dp/dn = 2 p (d ln t = 2 dn) and C_E = 2 c_b this is
    # 4 f c_b p sigma (H/s)(T_end) = 2 f c_b p sigma/(t_end s_end) exactly in RD (derived; see NOTES §3).
    # 5. closed form of the re-entry map for instant reheating vs the inflation_spectator classification table (H = 250 GeV)
    H = 250.0
    N0 = C.N0_instant(H)
    TRH = C.T_RH_instant(H)
    for m, ref in ((10.0, 6.1e5), (30.0, 1.8e-3), (40.0, 1.06e-7)):
        Tnum = float(C.T_of_nprime(N0 - m))
        f_g = lambda T: C.g_rho(T) ** 0.5 * C.g_s(T) ** (-1.0 / 3.0)
        Tcf = TRH * math.exp(-m)
        for _ in range(5):
            Tcf = TRH * math.exp(-m) * f_g(TRH) / f_g(Tcf)
        checks.append({"check": "T_re for a scale exiting %g e-folds before the end, H = 250 GeV" % m,
                       "numeric_GeV": Tnum, "closed_form_GeV": Tcf, "module_classification_table_GeV": ref,
                       "rel_dev_numeric_vs_closed_form": rel(Tnum, Tcf), "pass": bool(rel(Tnum, Tcf) < 0.05),
                       "note": "module table uses step g_* and instant RD (estimate); <= 2x agreement expected below 1 GeV"})
    # 6. relaxed-start closed form for the relic: rho/s ~ f C_E sigma (Lambda_1/2) (H/s)(T_cut) (RD, flat dp/dn)
    prod = json.load(open(os.path.join(C.RESULTS, "production.json")))
    row = min(prod["rows"], key=lambda rr: abs(rr["H_GeV"] - 200.0))
    ns = np.round(np.arange(0.0, 62.0 + 1e-9, 0.05), 10)
    q = R.quantities(ns, np.array(row["log10_dpdn_relaxed"]) * math.log(10), 50.0, 0.0, T_cuts=R.T_CUTS)
    num = 10 ** q["log10_rho_over_s_unit_Tdom_iso"]
    Tc = C.T_DOM_ISO
    cf = C.SIGMA * row["Lambda1_adopted"] / 2 * C.H_of_T(Tc) / C.s_of_T(Tc)
    add("relaxed relic integral vs closed form sigma (Lambda_1/2) H/s(T_cut), H = %.1f GeV" % row["H_GeV"], num, cf,
        0.10, "derived closed form (flat dp/dn = Lambda_1/2, RD, g_* const near 8 keV); the relaxed start's slow "
        "transient makes dp/dn exceed Lambda_1/2 by ~5% at B ~ 40 and ~20% at B ~ 11")
    # 7. regime B/C boundary p_-(end) = 1e-3 from the exact-start runs
    Hs = [rr["H_GeV"] for rr in prod["rows"]]
    for N in (50, 60):
        t = R.threshold(Hs, [math.log10(rr["p_minus_fresh_N%d" % N]) for rr in prod["rows"]], 1e-3)
        add("H at p_-(fresh, N%d) = 1e-3" % N, t["H_GeV"], {50: 367.8, 60: 355.07}[N], 5e-3, "v_inflation thresholds")
        t = R.threshold(Hs, [math.log10(rr["p_minus_fresh_N%d" % N]) for rr in prod["rows"]], 1e-10)
        add("H at p_-(fresh, N%d) = 1e-10" % N, t["H_GeV"], {50: 239.1, 60: 236.6}[N], 5e-3,
            "inflation_spectator / v_inflation")
    out = {"module": "minority_relic", "part": "checks", "n_checks": len(checks),
           "n_pass": sum(1 for c in checks if c["pass"]), "checks": checks}
    C.dump("checks.json", out)
    return out


if __name__ == "__main__":
    o = main()
    for c in o["checks"]:
        print(("PASS " if c["pass"] else "FAIL ") + c["check"], {k: v for k, v in c.items() if k not in ("check", "pass")})
    print(o["n_pass"], "/", o["n_checks"])
