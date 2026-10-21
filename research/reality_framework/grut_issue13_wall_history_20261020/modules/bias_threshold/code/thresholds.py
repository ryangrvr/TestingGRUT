"""thresholds.py -- critical systematic bias that separates (i) persistence to wall-radiation parity, (ii) decay between
the relic-safe temperature and parity (stable-relic overproduction), (iii) harmless early decay.

Scaling law (Bray, Adv. Phys. 43 (1994) 357, §4.2 eqs. (91)-(95) and §5.2 eq. (175), with t^(d/4) -> L^(d - lambda_bar)):
    m(K) = erf( y_f K^gamma ),   y_f = erfinv(m_f),   gamma = d - lambda_bar,
m_f = systematic population bias of the domains at formation (m = 2p - 1), K = comoving growth factor of the network
length from formation.  The network falls below 10% of the unbiased area (tau_10) when y_f K^gamma = x10.
Calibration [C]: 3D PRS collapse of network_sim's biased runs, gamma = 1.406 +- 0.053, x10 = 1.54 (fit_K10.json);
3D lattice autocorrelation lambda_bar = 1.61-1.69 (lattice_response.json).  Band gamma in [1.30, 1.50] (1.50 is the
Fisher-Huse bound lambda_bar >= d/2 for a scaling state).
    m_f < eps_persist = erf(x10 K_par^-gamma)            -> network survives to parity            (i)
    eps_persist < m_f < eps_safe = erf(x10 K_safe^-gamma) -> decays after T_safe: relic too large   (ii)
    m_f > eps_safe                                          -> decays before T_safe                   (iii)
The extrapolation is from lattice K <= 2.7-5 (3D) to K ~ 1e11-1e18: it rests on the scaling theory, not on runs.
Writes results/thresholds.json and results/headline.json.
"""
from __future__ import annotations

import itertools
import json
import math
import os

import numpy as np
from scipy.special import erf, erfinv

import common as C

OUT = os.path.join(C.HERE, "..", "results")

# ------------------------------------------------------------------ calibration (from this module's fits)
fitK = json.load(open(os.path.join(OUT, "fit_K10.json")))
latt = json.load(open(os.path.join(OUT, "lattice_response.json")))
FF3 = fitK["prs3d_bias_N128_L128"]["formation_fit"]
GAMMA_C = FF3["gamma"]                      # 1.406
GAMMA_BAND = (1.30, 1.50)
GAMMA_LIST = [1.30, GAMMA_C - FF3["sig_gamma"], GAMMA_C, GAMMA_C + FF3["sig_gamma"], 1.50]
X10_TABLE = {float(k): v for k, v in FF3["x10_at_fixed_gamma"].items()}


def x10_for(gamma):
    """x10 consistent with gamma: refit of ln y_f + gamma ln K10 (median over the 9 3D runs)."""
    pts = FF3["points"]
    return float(np.exp(np.median([math.log(p["y_f"]) + gamma * math.log(p["K10"]) for p in pts])))


def eps_of_K(K, gamma, x10):
    return float(erf(x10 * K ** (-gamma)))


def K10_of_m(m, gamma, x10):
    return (x10 / float(erfinv(min(abs(m), 0.999999)))) ** (1.0 / gamma)


# ------------------------------------------------------------------ inputs
A_LIST = [0.7, 0.8, 0.9]              # [Q] HKS arXiv:1309.5001 A = 0.8 +- 0.1; network_sim 0.78-0.88
F_LIST = [0.01, 1.0]                  # [P] rest-mass fraction of the wall energy ending in stable h, X, D
KAPPA_LIST = [1.0, 2.0]               # [E] earlier scaling emission adds ~ rho_w(T_d) once more (geometric sum)
GGEO_LIST = [1.0, 3.0]                # [E] network length V/Area at formation = g_geo x (KZ correlation length)
K0_OVER_K10 = [1.0, 2.3]              # [C network_sim] tau_0 (no wall left) vs tau_10 in K: 1.6-2.3
T_PAR_BACKREACTION_KEV = 6.37         # [Q cosmo_observable] parity with backreaction (A = 0.8)


def kz_row(xi_target=5e-4):
    kz = C.load_json("thermal_potential/results/kibble_zurek.json")
    rows = kz["rows"]
    lx = np.array([math.log(r["xi"]) for r in rows])
    i = int(np.argmin(np.abs(lx - math.log(xi_target))))
    r = rows[i]
    out = {"xi_row": r["xi"], "T_SM_GeV": r["T_SM_GeV"], "H_GeV": r["H_GeV"], "variants": {}}
    for k, v in r.items():
        if isinstance(v, dict) and v.get("xi_KZ_GeV^-1"):
            out["variants"][k] = {"xi_KZ_GeV^-1": v["xi_KZ_GeV^-1"], "log10_N_H": v["log10_N_H"]}
    out["z1_log10_N_H_tau0_x1e-2_x1e2"] = r.get("z1_log10_N_H_tau0_x1e-2_x1e2")
    out["z2_log10_N_H_tau0_x1e-2_x1e2"] = r.get("z2_log10_N_H_tau0_x1e-2_x1e2")
    return out


def history_cases():
    """Formation epoch and formation length scale for the three reference histories (all [P]/[Q]/[E] labelled)."""
    cases = {}
    kz = kz_row(5e-4)
    # H-star: thermal exit at T_h = T_c, xi = T_h/T_SM = 5e-4 [P]; KZ correlation length [Q thermal_potential]
    main = {k: v for k, v in kz["variants"].items() if k in ("z1_ising_ginzburg_matched", "z2_ising_ginzburg_matched")}
    ext = dict(kz["variants"])
    cases["H_star"] = {"label": "H-star: thermal restoration, exit at T_SM = 3.1e6 GeV (xi = 5e-4 [P]); KZ domains",
                       "T_f": kz["T_SM_GeV"], "xi0_main": [v["xi_KZ_GeV^-1"] for v in main.values()],
                       "xi0_ext": [v["xi_KZ_GeV^-1"] for v in ext.values()], "sQ": [1.0],
                       "provenance": "T_f, xi_KZ from modules/thermal_potential/results/kibble_zurek.json row xi = %.4g "
                                     "(main: z1/z2 Ginzburg-matched; ext: all six variants)" % kz["xi_row"],
                       "kz": kz}
    # H-star-Q: charged exit at T_h = 0.29 v -> T_SM = 0.29 v / xi = 580 TeV [Q REPORT/charge_phase];
    # formation scale 21-29 ell simulated (fast quench), KZ estimate ~960 ell at 200 TeV [E charged_exit_sim]
    ell = 1.0 / C.M_X
    T_fQ = 0.29 * C.V / 5e-4
    xi_KZQ = 960.0 * (200e3 / T_fQ) ** (2.0 / 3.0) * ell
    cases["H_star_Q"] = {"label": "H-star-Q: charged exit at T_SM ~ 580 TeV; domain scale 21 ell - 960 ell",
                         "T_f": T_fQ, "xi0_main": [21 * ell, 960 * ell], "xi0_ext": [21 * ell, 29 * ell, xi_KZQ, 960 * ell],
                         "sQ": [1.0, 1.3, 2.0],
                         "provenance": "T_h = 0.29 v (REPORT §2); xi0 = 21-29 ell measured at simulated H/m ~ 1e7 x "
                                       "physical, ~960 ell KZ estimate at 200 TeV (charged_exit_sim NOTES §6.1, §7; "
                                       f"rescaled to {T_fQ:.3g} GeV: {xi_KZQ / ell:.0f} ell); charged slowdown "
                                       "sQ = 1.3-2 (charged_exit_sim), applied as A_eff = sQ A"}
    # inflationary randomization, H_inf >> v [P]: spectator oscillation starts at H ~ sqrt(lambda) phi_end and is trapped
    # when its amplitude reaches ~v: H_trap ~ sqrt(lambda) v^2/phi_end, phi_end ~ 1.5 H_inf [E inflation_spectator §4];
    # instant visible reheating [P]; restoration by preheating may delay formation to T_trap sqrt(eps_fl),
    # eps_fl in [0.1, 1] [E].  Domain scale: Hubble-patch randomization 1/sqrt(H_inf H_trap) or 1-10/m_phi [E].
    for Hinf in (2.4e4, 1e6, 1e8, 1e10, 1e13):
        phi_end = 1.5 * Hinf
        H_trap = math.sqrt(C.LAMBDA_PHI) * C.V ** 2 / phi_end
        T_trap = C.T_of_H(H_trap)
        xi_list = [1.0 / math.sqrt(Hinf * H_trap), 1.0 / C.M_PHI, 10.0 / C.M_PHI]
        cases[f"infl_Hinf_{Hinf:.0e}"] = {
            "label": f"inflationary randomization, H_inf = {Hinf:.1e} GeV [P], instant reheating [P]",
            "T_f": T_trap, "T_f_range": [T_trap * math.sqrt(0.1), T_trap], "H_trap_GeV": H_trap,
            "xi0_main": xi_list, "xi0_ext": xi_list, "sQ": [1.0],
            "provenance": "H_trap = sqrt(lambda_phi) v^2/phi_end, phi_end = 1.5 H_inf [E, inflation_spectator §4.1-4.4]"}
    return cases


def K_targets(T_f, L_f, A, sQ, f, kappa):
    Tpar = C.T_parity(A, sQ)
    Tsafe = C.T_relic_safe(A, f, kappa, sQ)
    return {"T_par": Tpar, "T_safe": Tsafe,
            "K_par": C.K_of_T(Tpar, T_f, L_f, A, sQ),
            "K_par_backreaction_T": C.K_of_T(T_PAR_BACKREACTION_KEV * 1e-6, T_f, L_f, A, sQ),
            "K_safe": C.K_of_T(Tsafe, T_f, L_f, A, sQ),
            "K_horizon": C.K_of_T(T_f, T_f, L_f, A, sQ)}   # network length t/A already at T_f: sub-horizon phase factor


def band_for_case(case, which="xi0_main"):
    """Scan all nuisance inputs; return central values and min/max of log10 thresholds."""
    T_list = case.get("T_f_range", [case["T_f"]])
    rows = []
    for T_f, xi0, ggeo, A, sQ, f, kappa in itertools.product(T_list, case[which], GGEO_LIST, A_LIST, case["sQ"],
                                                             F_LIST, KAPPA_LIST):
        L_f = ggeo * xi0
        Kt = K_targets(T_f, L_f, A, sQ, f, kappa)
        for gamma in GAMMA_LIST:
            x10 = x10_for(gamma)
            for k0 in K0_OVER_K10:
                Kp = max(Kt["K_par"], Kt["K_par_backreaction_T"]) if k0 > 1 else min(Kt["K_par"], Kt["K_par_backreaction_T"])
                # persistence uses the 'gone' criterion (tau_0 ~ k0 x tau_10 in K): survive iff K0 > K_par
                eps_p = eps_of_K(Kp / k0, gamma, x10)
                eps_s = eps_of_K(Kt["K_safe"], gamma, x10)   # relic-safe uses tau_10 (energy released by then)
                rows.append({"T_f": T_f, "xi0": xi0, "ggeo": ggeo, "A": A, "sQ": sQ, "f": f, "kappa": kappa,
                             "gamma": gamma, "x10": x10, "k0": k0, "K_par": Kp, "K_safe": Kt["K_safe"],
                             "K_horizon": Kt["K_horizon"], "T_safe": Kt["T_safe"], "T_par": Kt["T_par"],
                             "eps_persist": eps_p, "eps_safe": eps_s})
    def rng(key, sel=None):
        vals = [r[key] for r in rows if (sel is None or sel(r))]
        return [float(min(vals)), float(max(vals))]
    # central: gamma_c, x10(gamma_c), A = 0.8, ggeo = 1, k0 = 1, f = 1, kappa = 1, sQ = first, middle xi0
    xi_c = sorted(case[which])[len(case[which]) // 2] if len(case[which]) % 2 == 1 else math.sqrt(
        sorted(case[which])[0] * sorted(case[which])[-1])
    T_fc = case["T_f"]
    Kc = K_targets(T_fc, xi_c, 0.8, case["sQ"][0], 1.0, 1.0)
    x10c = x10_for(GAMMA_C)
    central = {"T_f_GeV": T_fc, "xi0_GeV^-1": xi_c, "A": 0.8, "g_geo": 1.0, "f": 1.0, "kappa": 1.0,
               "sQ": case["sQ"][0], "gamma": GAMMA_C, "x10": x10c,
               "T_par_keV": Kc["T_par"] * 1e6, "T_safe_MeV": Kc["T_safe"] * 1e3,
               "log10_K_par": math.log10(Kc["K_par"]), "log10_K_safe": math.log10(Kc["K_safe"]),
               "log10_K_horizon_subhorizon_phase": math.log10(Kc["K_horizon"]),
               "eps_persist": eps_of_K(Kc["K_par"], GAMMA_C, x10c), "eps_safe": eps_of_K(Kc["K_safe"], GAMMA_C, x10c)}
    lg = lambda v: [math.log10(v[0]), math.log10(v[1])]
    band = {"log10_K_par": lg(rng("K_par")), "log10_K_safe": lg(rng("K_safe")),
            "log10_eps_persist": lg(rng("eps_persist")), "log10_eps_safe": lg(rng("eps_safe")),
            "log10_eps_persist_gamma_central_only": lg(rng("eps_persist", lambda r: abs(r["gamma"] - GAMMA_C) < 1e-9)),
            "log10_eps_safe_gamma_central_only": lg(rng("eps_safe", lambda r: abs(r["gamma"] - GAMMA_C) < 1e-9)),
            "log10_eps_persist_inputs_at_central_gamma_x_f_kappa_fixed": lg(rng(
                "eps_persist", lambda r: abs(r["gamma"] - GAMMA_C) < 1e-9 and r["k0"] == 1.0)),
            "log10_eps_safe_f1": lg(rng("eps_safe", lambda r: r["f"] == 1.0)),
            "log10_eps_safe_f0.01": lg(rng("eps_safe", lambda r: r["f"] == 0.01)),
            "T_safe_MeV": [1e3 * x for x in rng("T_safe")], "n_combinations": len(rows)}
    return central, band


def classify(m, case_central):
    """Fate of a given systematic formation bias m for a history (central inputs)."""
    g, x10 = case_central["gamma"], case_central["x10"]
    K10 = K10_of_m(m, g, x10)
    lK = math.log10(K10)
    if lK < case_central["log10_K_safe"]:
        fate = "(iii) decays before the relic-safe temperature (harmless for relic/parity)"
        if lK < case_central["log10_K_horizon_subhorizon_phase"]:
            fate += "; already before the network reaches horizon size"
    elif lK < case_central["log10_K_par"]:
        fate = "(ii) decays between T_safe and parity: stable-relic overproduction"
    else:
        fate = "(i) survives to parity"
    return {"m_f": m, "log10_K10": lK, "fate": fate}


def main():
    cases = history_cases()
    out = {"calibration": {
        "law": "m(K) = erf(y_f K^gamma), decay (tau_10) at y_f K^gamma = x10; eps(K) = erf(x10 K^-gamma)",
        "gamma_central_3D_PRS_collapse": GAMMA_C, "sig_gamma": FF3["sig_gamma"], "x10_central": x10_for(GAMMA_C),
        "gamma_band": GAMMA_BAND, "gamma_list_scanned": GAMMA_LIST,
        "x10_for_gamma": {f"{g:.3f}": x10_for(g) for g in GAMMA_LIST},
        "3D_lattice_lambda_bar_safe": latt["prs3d_N160"].get("lambda_bar_safe"),
        "3D_lattice_lambda_bar_ext": latt["prs3d_N160"].get("lambda_bar_ext"),
        "3D_lattice_direct_response_gamma_safe": latt["prs3d_N160"].get("gamma_Ry_safe_m0=0.002"),
        "quench_amplification_3D_lattice_R_y_form": latt["prs3d_N160"].get("R_y_at_formation_m0=0.002"),
        "inputs": {"A [Q]": A_LIST, "f [P]": F_LIST, "kappa_emit [E]": KAPPA_LIST, "g_geo [E]": GGEO_LIST,
                   "K0/K10 [C network_sim]": K0_OVER_K10, "rho/s max GeV [Q]": C.RHO_OVER_S_MAX,
                   "T_par backreaction keV [Q cosmo_observable]": T_PAR_BACKREACTION_KEV}}}
    hist = {}
    for name, case in cases.items():
        central, band = band_for_case(case, "xi0_main")
        _, band_ext = band_for_case(case, "xi0_ext")
        hist[name] = {"label": case["label"], "provenance": case["provenance"], "central": central, "band": band,
                      "band_extended_formation_scale": band_ext}
        if "kz" in case:
            hist[name]["kz_inputs"] = case["kz"]
        if "H_trap_GeV" in case:
            hist[name]["H_trap_GeV"] = case["H_trap_GeV"]
    out["histories"] = hist

    # extrapolation span (decades in L/xi0) relative to the largest lattice K
    K_latt = {"network_sim 3D K10 max": 2.69, "this module 3D lattice K (t <= L/2)": 5.4,
              "network_sim 2D K10 max": 10.8}
    hs = hist["H_star"]["band"]
    out["extrapolation_decades"] = {
        "H_star_parity": [hs["log10_K_par"][0] - math.log10(5.4), hs["log10_K_par"][1] - math.log10(2.69)],
        "H_star_relic_safe": [hs["log10_K_safe"][0] - math.log10(5.4), hs["log10_K_safe"][1] - math.log10(2.69)],
        "lattice_K_reached": K_latt}

    # statistical (exact Z2) vs systematic
    hc = hist["H_star"]["central"]
    Kp = 10 ** hc["log10_K_par"]
    out["statistical_vs_systematic"] = {
        "statistical_bias_on_scale_K: (xi0/L)^(d/2)": "K^-1.5",
        "H_star_K_par_central": Kp, "statistical_floor_at_parity_K_par^-1.5": Kp ** -1.5,
        "eps_persist_central": hc["eps_persist"],
        "amplified_statistical_bias_at_scale_L: K^(gamma - d/2)": {f"{g:.2f}": Kp ** (g - 1.5) for g in GAMMA_LIST},
        "note": "exact Z2: the random magnetization of a region of the current domain size ~ K^(gamma - d/2) <= O(1) "
                "never grows (Fisher-Huse lambda_bar >= d/2), so it cannot remove the network; a systematic bias is "
                "the same on every scale and grows as K^gamma."}

    # candidate residual biases [Q review finding / modules], classified for H-star (and inflationary case)
    cands = {
        "thermal remnant 1e-3 of soft-mode level, undamped [E, review finding]": float(erf(1e-3 / math.sqrt(2))),
        "post-inflation trapping p_- = 0.4999 [Q inflation_spectator]": 2e-4,
        "post-inflation trapping p_- = 0.4997 [Q inflation_spectator]": 6e-4,
        "forced charge 2 n_c, p(+) = 0.71 [Q charged_exit_sim]": 0.42,
        "relaxed start p_- = 0.23 at end of inflation [Q]": 0.54,
    }
    out["candidates_H_star"] = {k: classify(v, hc) for k, v in cands.items()}
    ci = hist["infl_Hinf_1e+08"]["central"]
    out["candidates_infl_Hinf_1e8"] = {k: classify(v, ci) for k, v in cands.items() if "trapping" in k}
    # damping needed for an undamped 1e-3 thermal remnant to reach the persistence regime
    m_rem = cands["thermal remnant 1e-3 of soft-mode level, undamped [E, review finding]"]
    out["damping_required_H_star"] = {
        "m_remnant_undamped": m_rem,
        "amplitude_efolds_to_reach_eps_safe": math.log(m_rem / hc["eps_safe"]),
        "amplitude_efolds_to_reach_eps_persist": math.log(m_rem / hc["eps_persist"]),
        "band_efolds_to_eps_persist": [math.log(m_rem / 10 ** hs["log10_eps_persist"][1]),
                                       math.log(m_rem / 10 ** hs["log10_eps_persist"][0])],
        "band_efolds_to_eps_safe": [math.log(m_rem / 10 ** hs["log10_eps_safe"][1]),
                                    math.log(m_rem / 10 ** hs["log10_eps_safe"][0])]}

    # one-at-a-time sensitivity for H-star (other inputs central)
    Kp_c, Ks_c = 10 ** hc["log10_K_par"], 10 ** hc["log10_K_safe"]
    sens = {"gamma": {f"{g:.3f}": {"eps_persist": eps_of_K(Kp_c, g, x10_for(g)), "eps_safe": eps_of_K(Ks_c, g, x10_for(g))}
                      for g in GAMMA_LIST + [1.327, 1.39]}}
    hstar = cases["H_star"]
    sens["f_kappa"] = {}
    for f in (1.0, 0.1, 0.01):
        for kap in (1.0, 2.0):
            Kt = K_targets(hstar["T_f"], hc["xi0_GeV^-1"], 0.8, 1.0, f, kap)
            sens["f_kappa"][f"f={f:g},kappa={kap:g}"] = {"T_safe_MeV": Kt["T_safe"] * 1e3,
                                                        "log10_K_safe": math.log10(Kt["K_safe"]),
                                                        "eps_safe": eps_of_K(Kt["K_safe"], GAMMA_C, x10_for(GAMMA_C))}
    sens["A"] = {}
    for A in A_LIST:
        Kt = K_targets(hstar["T_f"], hc["xi0_GeV^-1"], A, 1.0, 1.0, 1.0)
        sens["A"][f"{A}"] = {"T_par_keV": Kt["T_par"] * 1e6, "log10_K_par": math.log10(Kt["K_par"]),
                             "eps_persist": eps_of_K(Kt["K_par"], GAMMA_C, x10_for(GAMMA_C)),
                             "eps_safe": eps_of_K(Kt["K_safe"], GAMMA_C, x10_for(GAMMA_C))}
    sens["xi0_variants"] = {}
    for k, v in hstar["kz"]["variants"].items():
        Kt = K_targets(hstar["T_f"], v["xi_KZ_GeV^-1"], 0.8, 1.0, 1.0, 1.0)
        sens["xi0_variants"][k] = {"log10_N_H": v["log10_N_H"], "log10_K_par": math.log10(Kt["K_par"]),
                                   "eps_persist": eps_of_K(Kt["K_par"], GAMMA_C, x10_for(GAMMA_C)),
                                   "eps_safe": eps_of_K(Kt["K_safe"], GAMMA_C, x10_for(GAMMA_C))}
    sens["Ymax_4.35e-10_instead_of_4.373e-10_T_safe_ratio"] = (
        C.T_relic_safe(0.8, 1.0, 1.0, 1.0, 4.35e-10) / C.T_relic_safe(0.8, 1.0, 1.0, 1.0))
    out["sensitivity_H_star"] = sens
    # bias expressed as a coherent zero-mode offset over the fluctuation level at formation: phibar/sigma = sqrt2 erfinv(m)
    out["zero_mode_equivalent_H_star"] = {"phibar_over_sigma_at_eps_persist": math.sqrt(2) * float(erfinv(hc["eps_persist"])),
                                          "phibar_over_sigma_at_eps_safe": math.sqrt(2) * float(erfinv(hc["eps_safe"]))}

    # comparison: reviewer's naive 2D-fit extrapolation and gamma = 0.84 (2D relativistic) stress values
    out["comparison"] = {
        "naive_2D_fit_(1/2-p)_at_K=1e17": (1e17 / 0.262) ** (-1 / 1.255),
        "naive_2D_fit_m0_at_K_par_H_star_central": 2 * (Kp / 0.262) ** (-1 / 1.255),
        "network_sim_3D_fit_initial_m0_at_K_par_H_star_central (c=0.49, alpha=0.722)":
            (Kp / fitK["prs3d_bias_N128_L128"]["fit_all"]["c"]) ** (-1 / fitK["prs3d_bias_N128_L128"]["fit_all"]["alpha"]),
        "stress_gamma_0.84_eps_persist_H_star": eps_of_K(Kp, 0.84, 1.55),
        "stress_gamma_0.84_eps_safe_H_star": eps_of_K(10 ** hc["log10_K_safe"], 0.84, 1.55)}

    with open(os.path.join(OUT, "thresholds.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    return out


if __name__ == "__main__":
    o = main()
    print("calibration gamma %.3f +- %.3f, x10 %.3f" % (o["calibration"]["gamma_central_3D_PRS_collapse"],
                                                      o["calibration"]["sig_gamma"], o["calibration"]["x10_central"]))
    for n, h in o["histories"].items():
        c, b = h["central"], h["band"]
        print(f"{n}: T_f {c['T_f_GeV']:.3g} GeV, T_par {c['T_par_keV']:.2f} keV, T_safe {c['T_safe_MeV']:.1f} MeV, "
              f"log10 K_par {c['log10_K_par']:.2f} [{b['log10_K_par'][0]:.2f}, {b['log10_K_par'][1]:.2f}], "
              f"log10 K_safe {c['log10_K_safe']:.2f} [{b['log10_K_safe'][0]:.2f}, {b['log10_K_safe'][1]:.2f}], "
              f"log10 K_hor {c['log10_K_horizon_subhorizon_phase']:.2f}")
        print(f"     eps_persist {c['eps_persist']:.2e}  band log10 [{b['log10_eps_persist'][0]:.2f}, {b['log10_eps_persist'][1]:.2f}]"
              f" (gamma central only [{b['log10_eps_persist_gamma_central_only'][0]:.2f}, {b['log10_eps_persist_gamma_central_only'][1]:.2f}]);"
              f" eps_safe {c['eps_safe']:.2e} band [{b['log10_eps_safe'][0]:.2f}, {b['log10_eps_safe'][1]:.2f}]"
              f" (f=1 [{b['log10_eps_safe_f1'][0]:.2f}, {b['log10_eps_safe_f1'][1]:.2f}], f=0.01 [{b['log10_eps_safe_f0.01'][0]:.2f}, {b['log10_eps_safe_f0.01'][1]:.2f}])")
        be = h["band_extended_formation_scale"]
        print(f"     extended formation scale: eps_persist [{be['log10_eps_persist'][0]:.2f}, {be['log10_eps_persist'][1]:.2f}], eps_safe [{be['log10_eps_safe'][0]:.2f}, {be['log10_eps_safe'][1]:.2f}]")
    print("extrapolation decades:", o["extrapolation_decades"])
    print("stat vs syst:", {k: v for k, v in o["statistical_vs_systematic"].items() if k != "note"})
    for k, v in o["candidates_H_star"].items():
        print("  H*:", k, "log10 K10 %.2f" % v["log10_K10"], v["fate"])
    for k, v in o["candidates_infl_Hinf_1e8"].items():
        print("  infl 1e8:", k, "log10 K10 %.2f" % v["log10_K10"], v["fate"])
    print("damping:", o["damping_required_H_star"])
    print("comparison:", o["comparison"])
    print("sensitivity gamma:", {k: ("%.2e" % v["eps_persist"], "%.2e" % v["eps_safe"]) for k, v in o["sensitivity_H_star"]["gamma"].items()})
    print("sensitivity f,kappa:", {k: (round(v["T_safe_MeV"], 3), round(v["log10_K_safe"], 2), "%.2e" % v["eps_safe"]) for k, v in o["sensitivity_H_star"]["f_kappa"].items()})
    print("sensitivity A:", {k: (round(v["T_par_keV"], 2), round(v["log10_K_par"], 2), "%.2e" % v["eps_persist"], "%.2e" % v["eps_safe"]) for k, v in o["sensitivity_H_star"]["A"].items()})
    print("sensitivity xi0:", {k: (v["log10_N_H"], round(v["log10_K_par"], 2), "%.2e" % v["eps_persist"], "%.2e" % v["eps_safe"]) for k, v in o["sensitivity_H_star"]["xi0_variants"].items()})
    print("Ymax 4.35e-10:", o["sensitivity_H_star"]["Ymax_4.35e-10_instead_of_4.373e-10_T_safe_ratio"])
    print("zero-mode equivalent:", o["zero_mode_equivalent_H_star"])
