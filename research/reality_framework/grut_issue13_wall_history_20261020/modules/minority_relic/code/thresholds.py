"""Thresholds in H_inf for every criterion, start condition, N_* and e-fold convention; systematics; spectra.

Reads results/production.json and results/collapse.json (this module) and the R_X entropic factors from
inflation_spectator/results/stochastic_1d.json (read-only, [Q]). Writes results/thresholds.json,
results/spectra.json, results/headline.json.
"""
import math
import json
import os
import numpy as np
import mr_common as C
import relic as R

R.load_collapse()
LN10 = math.log(10.0)
LOG_RHOS_MAX = math.log10(C.RHO_OVER_S_MAX)


def load():
    d = json.load(open(os.path.join(C.RESULTS, "production.json")))
    rows = sorted(d["rows"], key=lambda r: r["H_GeV"])
    ns = np.round(np.arange(0.0, rows[0]["n_grid"]["nmax"] + 1e-9, rows[0]["n_grid"]["dn"]), 10)
    return rows, ns


def RX_interp():
    p = os.path.join(C.PACKET, "modules", "inflation_spectator", "results", "stochastic_1d.json")
    d = json.load(open(p))
    H = np.array([r["H_GeV"] for r in d["rows"] if r["R_X_entropic"] is not None])
    RX = np.array([r["R_X_entropic"] * r["R_D_entropic"] for r in d["rows"] if r["R_X_entropic"] is not None])
    return lambda h: float(np.exp(np.interp(math.log(h), np.log(H), np.log(RX), left=np.nan, right=np.nan)))


RXf = RX_interp()
# 2D (phi, |X|) / 1D odd-eigenvalue ratio vs B_HM, v_inflation VERIFY.md §3 (corrected small-B values) [Q]
_B2D = np.array([0.11, 1.16, 2.71, 5.73, 10.9, 17.1])
_R2D = np.array([1.86, 2.49, 3.55, 6.40, 11.7, 16.7])


def FX_factor(kind, H, B):
    """Multiplicative X-fluctuation correction to dp/dn (1D X = 0 path -> X fluctuating).
    '1': none. '2': fresh-start 2D/1D p_- ratio measured at 330 GeV (1.91 N50, 2.17 N60; v_inflation).
    '2D': 2D/1D Lambda_1 ratio interpolated in ln B (v_inflation), above B = 17.1 continued as 0.70 R_X (Langer).
    'RX': Langer entropic R_X R_D where defined (B >= 3; an over-correction = upper bound), else the 2D ratio."""
    if kind == "1":
        return 1.0
    if kind == "2":
        return 2.0
    rx = RXf(H)
    if B <= _B2D[-1]:
        r2d = float(np.exp(np.interp(math.log(max(B, _B2D[0])), np.log(_B2D), np.log(_R2D))))
    else:
        r2d = 0.70 * rx if np.isfinite(rx) else _R2D[-1]
    if kind == "2D":
        return r2d
    if kind == "RX":
        return rx if (np.isfinite(rx) and B >= 3.0) else r2d
    raise ValueError(kind)
VOL_RATIO = None   # set in main (observable-universe / Hubble volume)

# criteria: name -> (quantity key, log10 offset added to the quantity, target)
F_LIST = (1.0, 0.1, 0.01)


def criteria_values(q):
    """Return dict criterion -> log10(value/target) (crossing zero = threshold)."""
    ce = R.C_E_SPHERE
    out = {}
    for f in F_LIST:
        out["relic_Omega_le_0.12_f%g_CEsphere" % f] = q["log10_rho_over_s_unit_Tdom_iso"] + math.log10(f * ce) - LOG_RHOS_MAX
    out["relic_Omega_le_0.12_f1_CE3x"] = q["log10_rho_over_s_unit_Tdom_iso"] + math.log10(3 * ce) - LOG_RHOS_MAX
    out["relic_f1_cut_Tdom_net"] = q["log10_rho_over_s_unit_Tdom_net"] + math.log10(ce) - LOG_RHOS_MAX
    out["relic_f1_only_collapses_before_1MeV"] = q["log10_rho_over_s_unit_BBN_1MeV"] + math.log10(ce) - LOG_RHOS_MAX
    out["relic_f1_injection_to_0.03eV_(critic_style)"] = q["log10_rho_over_s_unit_late_0.03eV"] + math.log10(ce) - LOG_RHOS_MAX
    out["N_reenter_after_Tdom_iso_lt_1"] = q["log10_N_after_Tdom_iso"]
    out["N_reenter_after_Tdom_iso_lt_1e-3"] = q["log10_N_after_Tdom_iso"] + 3
    out["N_reenter_after_Tdom_net_lt_1"] = q["log10_N_after_Tdom_net"]
    out["Omega_late_enclosed_max_le_Omega_DM"] = q["log10_rho_late_over_s_max"] - LOG_RHOS_MAX
    out["Omega_late_enclosed_max_le_1e-3_Omega_DM"] = q["log10_rho_late_over_s_max"] - LOG_RHOS_MAX + 3
    out["Omega_late_enclosed_max_le_1e-6_Omega_DM"] = q["log10_rho_late_over_s_max"] - LOG_RHOS_MAX + 6
    out["N_reenter_after_Teq_lt_1"] = q["log10_N_after_Teq"]
    out["N_reenter_after_Trec_lt_1"] = q["log10_N_after_Trec"]
    out["N_reenter_after_BBN_1MeV_lt_1"] = q["log10_N_after_BBN_1MeV"]
    out["N_reenter_after_BBN_1MeV_lt_1e-3"] = q["log10_N_after_BBN_1MeV"] + 3
    out["N_reenter_after_BBN_70keV_lt_1"] = q["log10_N_after_BBN_70keV"]
    VR = math.log10(VOL_RATIO)
    out["N_obs_universe_reenter_after_Tdom_iso_lt_1"] = q["log10_N_after_Tdom_iso"] + VR
    out["N_obs_universe_reenter_after_BBN_1MeV_lt_1"] = q["log10_N_after_BBN_1MeV"] + VR
    out["N_obs_universe_reenter_after_Trec_lt_1"] = q["log10_N_after_Trec"] + VR
    out["N_superhorizon_wall_in_Hubble_volume_today_lt_1"] = q["log10_N_wall_superhorizon_today"]
    out["N_all_sizes_in_Hubble_volume_lt_1e-3_(regime_A)"] = q["log10_N_all_sizes_in_Hubble_volume"] + 3
    return out


def scenario_DeltaN(conv, Nstar, H, DeltaN=None):
    if conv == "instant":
        return Nstar - C.N0_instant(H)
    if conv == "exit":
        return 0.0
    return DeltaN


def run_scenario(rows, ns, start, Nstar, conv, DeltaN=None, s_off=0.0, FX="1"):
    vals = {}
    per_H = []
    for r in rows:
        H = r["H_GeV"]
        ld = np.array(r["log10_dpdn_" + start]) * LN10
        fx = FX_factor(FX, H, r["B_HM"])
        dN = scenario_DeltaN(conv, Nstar, H, DeltaN)
        q = R.quantities(ns, ld, Nstar, dN, s_off=s_off, FX=fx, T_cuts=R.T_CUTS)
        cv = criteria_values(q)
        for k, v in cv.items():
            vals.setdefault(k, []).append(v)
        per_H.append({"H_GeV": H, "DeltaN": dN, **{k: (None if not np.isfinite(v) else v) for k, v in q.items()}})
    Hs = [r["H_GeV"] for r in rows]
    thr = {k: R.threshold(Hs, v, 1.0) for k, v in vals.items()}
    return thr, per_H


def spectra(rows, ns, Hsel=(200.0, 230.0, 240.0, 260.0, 300.0, 400.0, 600.0, 1000.0)):
    out = []
    c_over_H0_Mpc = 299792.458 / (100 * cosmo_h())
    MSUN_GEV = 1.98847e30 * 5.60958865e26
    for Ht in Hsel:
        r = min(rows, key=lambda rr: abs(math.log(rr["H_GeV"] / Ht)))
        H = r["H_GeV"]
        for start in ("fresh", "relaxed"):
            for Nstar in (50.0, 60.0):
                dN = Nstar - C.N0_instant(H)
                ld = np.array(r["log10_dpdn_" + start])
                tab = []
                for nn in np.arange(0.0, Nstar + 1e-9, 2.0):
                    i = int(round(nn / 0.05))
                    npr = nn - dN
                    T = float(C.T_of_nprime(npr)) if npr >= 0 else None
                    row = {"n_from_start": nn, "efolds_before_end": Nstar - nn, "nprime": npr,
                           "log10_dpdn": ld[i],
                           "comoving_radius_today_Mpc": c_over_H0_Mpc * math.exp(-npr)}
                    if T is not None and np.isfinite(T):
                        Hre = C.H_of_T(T)
                        E = R.KAPPA_E_RD * C.SIGMA / Hre ** 2
                        row.update({"T_re_GeV": T, "E_released_per_domain_GeV": E,
                                    "supercritical_NG_energy_invalid": bool(C.SIGMA / (C.M_P ** 2 * Hre) > 1.0),
                                    "E_released_per_domain_Msun": E / MSUN_GEV,
                                    "wall_energy_fraction_sigma_over_MP2H": C.SIGMA / (C.M_P ** 2 * Hre),
                                    "r_s_over_wall_width": E / (4 * math.pi * C.M_P ** 2) / C.DELTA_WALL,
                                    "horizon_mass_Msun": 4 * math.pi * C.M_P ** 2 / Hre / MSUN_GEV,
                                    "log10_domains_per_Hubble_volume_per_efold": ld[i] + 3 * npr / LN10,
                                    "log10_rho_over_s_per_efold_GeV_fCE1": ld[i] + math.log10(C.SIGMA * Hre / C.s_of_T(T))})
                    tab.append(row)
                out.append({"H_GeV": H, "start": start, "N_star": Nstar, "DeltaN_instant": dN,
                            "T_RH_instant_GeV": C.T_RH_instant(H), "N0_instant": C.N0_instant(H), "table": tab})
    return out


def cosmo_h():
    return C.H0 / C.cosmo_core.H100_GEV


def Nstar_reach():
    """N_0 range for H_inf in the regime-B/C band under w <= 1 post-inflation histories [P]."""
    res = []
    for H in (150.0, 250.0, 400.0, 1000.0, 1600.0):
        rho_inf = 3 * C.M_P ** 2 * H ** 2
        N0 = C.N0_instant(H)
        rho4 = C.cosmo_core.rho_rad(4e-3)
        rho10 = C.cosmo_core.rho_rad(1e-2)
        res.append({"H_GeV": H, "N0_instant": N0,
                    "N0_max_kination_until_4MeV": N0 + math.log(rho_inf / rho4) / 12.0,
                    "N0_min_matter_dominated_reheating_T_RH_10MeV": N0 - math.log(rho_inf / rho10) / 12.0})
    return res


def particle_horizon_volume_ratio():
    """(chi_p a0 H0)^3: observable-universe comoving volume over today's Hubble volume (LambdaCDM calibration)."""
    from scipy.integrate import quad
    Om = C.cosmo_core.RHO_M0 / (3 * C.M_P ** 2 * C.H0 ** 2)
    Or = C.cosmo_core.RHO_R0 / (3 * C.M_P ** 2 * C.H0 ** 2)
    OL = 1 - Om - Or
    f = lambda a: 1.0 / math.sqrt(Or + Om * a + OL * a ** 4)
    chi = quad(f, 0, 1, limit=200)[0]
    return chi, chi ** 3


def main():
    rows, ns = load()
    Hs = [r["H_GeV"] for r in rows]
    out = {"module": "minority_relic", "part": "thresholds",
           "labels": "computed; inflation/reheating/N_*/start conditions are outside-action preparation data [P]",
           "inputs": {"sigma_GeV3": C.SIGMA, "kappa_E_RD": R.KAPPA_E_RD, "C_E_sphere": R.C_E_SPHERE,
                      "rho_over_s_max_GeV": C.RHO_OVER_S_MAX, "T_dom_iso_GeV": C.T_DOM_ISO,
                      "T_dom_iso_recomputed_GeV": C.T_dom_iso(), "T_dom_net_GeV": C.T_DOM_NET,
                      "T_eq_GeV": R.T_CUTS["counts"]["Teq"], "T_rec_GeV": C.T_REC,
                      "nprime": {k: C.nprime_of_T(v) for k, v in R.T_CUTS["counts"].items()},
                      "H_grid_GeV": [Hs[0], Hs[-1], len(Hs)]}}
    global VOL_RATIO
    chi, vr = particle_horizon_volume_ratio()
    VOL_RATIO = vr
    out["inputs"]["particle_horizon_chi_a0H0"] = chi
    out["inputs"]["observable_over_Hubble_volume"] = vr
    out["Nstar_reach"] = Nstar_reach()
    main_thr = {}
    per_H_store = {}
    for start in ("fresh", "relaxed"):
        for Nstar in (50.0, 60.0):
            for conv in ("instant", "exit"):
                key = "%s_N%d_%s" % (start, Nstar, conv)
                thr, per_H = run_scenario(rows, ns, start, Nstar, conv)
                main_thr[key] = thr
                per_H_store[key] = per_H
    out["thresholds_GeV"] = main_thr
    # systematics on the instant-reheating scenarios
    sysd = {}
    for start in ("fresh", "relaxed"):
        for Nstar in (50.0, 60.0):
            for label, kw in (("FX2", {"FX": "2"}), ("FX_2Drate", {"FX": "2D"}), ("FX_RX_bound", {"FX": "RX"}),
                              ("s_off_+3", {"s_off": 3.0}), ("s_off_-2", {"s_off": -2.0})):
                key = "%s_N%d_instant_%s" % (start, Nstar, label)
                thr, _ = run_scenario(rows, ns, start, Nstar, "instant", **kw)
                sysd[key] = thr
    out["systematics_thresholds_GeV"] = sysd
    # DeltaN scan (fresh start; DeltaN - s_off is the only combination that enters)
    scan = {}
    for dN in (-8.0, -4.0, 0.0, 2.0, 4.0, 8.0, 12.0, 16.0, 21.0):
        thr, _ = run_scenario(rows, ns, "fresh", 60.0, "fixed", DeltaN=dN)
        scan["DeltaN=%g" % dN] = {k: thr[k] for k in ("relic_Omega_le_0.12_f1_CEsphere", "relic_Omega_le_0.12_f0.01_CEsphere",
                                                       "N_reenter_after_Tdom_iso_lt_1", "N_reenter_after_BBN_1MeV_lt_1",
                                                       "N_reenter_after_Trec_lt_1", "Omega_late_enclosed_max_le_Omega_DM")}
    out["fresh_DeltaN_scan_N60window"] = scan
    # Lambda_1 table and critic comparison
    L1 = [r["Lambda1_adopted"] for r in rows]
    out["critic_comparison"] = {
        "critic_estimate": "relaxed start, f c_b ~ 1, injection to late times needs Lambda_1/2 <~ 1e-11, i.e. H_inf <~ 240 GeV [E]",
        "H_at_Lambda1_over_2_eq_1e-11": R.threshold(Hs, [math.log10(l / 2) for l in L1], 1e-11),
        "H_at_p_minus_fresh_N60_eq_1e-10": R.threshold(Hs, [math.log10(r["p_minus_fresh_N60"]) for r in rows], 1e-10),
        "H_at_p_minus_fresh_N50_eq_1e-10": R.threshold(Hs, [math.log10(r["p_minus_fresh_N50"]) for r in rows], 1e-10)}
    # closed forms for the relaxed start (dp/dn ~ Lambda_1/2, RD, n'-integrals done analytically) [derived]
    sHs_dom = C.SIGMA * C.H_of_T(C.T_DOM_ISO) / C.s_of_T(C.T_DOM_ISO)
    lam = [math.log10(l) for l in L1]
    nB, nD, nR = (C.nprime_of_T(C.T_BBN), C.nprime_of_T(C.T_DOM_ISO), C.nprime_of_T(C.T_REC))
    cf = {"Lambda1_max_N_after_BBN_lt_1": 6 * math.exp(-3 * nB), "Lambda1_max_N_after_Tdom_lt_1": 6 * math.exp(-3 * nD),
          "Lambda1_max_N_after_Trec_lt_1": 6 * math.exp(-3 * nR),
          "Lambda1_max_relic_f1_CEsphere": 2 * C.RHO_OVER_S_MAX / (R.C_E_SPHERE * sHs_dom),
          "sigma_H_over_s_at_Tdom_iso_GeV": sHs_dom}
    cf["H_thresholds_from_closed_forms_GeV"] = {k.replace("Lambda1_max_", ""): R.threshold(Hs, lam, v)
                                                for k, v in cf.items() if k.startswith("Lambda1_max")}
    out["relaxed_closed_forms"] = cf
    out["p_minus_check"] = [{"H_GeV": r["H_GeV"], "B_HM": r["B_HM"], "Lambda1_adopted": r["Lambda1_adopted"],
                             "Lambda1_eig": r["Lambda1_eig"], "Lambda1_mfpt": r["Lambda1_mfpt"],
                             "p_fresh_N50": r["p_minus_fresh_N50"], "p_fresh_N60": r["p_minus_fresh_N60"],
                             "p_relaxed_N50": r["p_minus_relaxed_N50"], "p_relaxed_N60": r["p_minus_relaxed_N60"],
                             "N0_instant": C.N0_instant(r["H_GeV"]), "T_RH_instant": C.T_RH_instant(r["H_GeV"]),
                             "R_X_times_R_D": RXf(r["H_GeV"])} for r in rows]
    C.dump("thresholds.json", out)
    C.dump("per_H_quantities.json", per_H_store)
    C.dump("spectra.json", {"module": "minority_relic", "part": "minority-domain spectra vs comoving size",
                            "labels": "computed; instant reheating [P]; c_R = 1 (s_off = 0); kappa_E (RD sphere)",
                            "rows": spectra(rows, ns)})
    return out


if __name__ == "__main__":
    o = main()
    for k, v in o["thresholds_GeV"].items():
        print(k)
        for c, t in v.items():
            print("   %-60s %s" % (c, t))
