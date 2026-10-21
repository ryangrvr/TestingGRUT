"""Part 3: classification of the inflation-prepared sign field and the minority-domain size distribution.

Inputs: results/stochastic_1d.json (and stochastic_2d.json when present). Re-computes, for a few H_inf,
the delta-start crossing-rate profile dp_-/dn and maps comoving scales to re-entry temperatures (instant
reheating into the visible SM bath: OUTSIDE-ACTION preparation data). Also evaluates the Gaussian-field
approximation requested by the brief and shows where it fails (it misses barrier trapping).
"""
import json
import os
import numpy as np
from scipy.special import erfc
from common import Spectator, V, MP, SIGMA_TREE, PC_SITE_SC, RESULTS, BETA_FP, dump
import fp1d
import cosmo


def load(name):
    p = os.path.join(RESULTS, name)
    return json.load(open(p)) if os.path.exists(p) else None


def rate_profile(H, Nstar=50.0):
    sp = Spectator(H)
    yy, h = fp1d.make_grid(sp)
    m0 = int(np.argmin(np.abs(yy - V / H)))
    ns, logp, floor = fp1d.p_minus_delta_curve(sp, yy, h, m0, Nmax=Nstar)
    L1, _ = fp1d.lambda1_mfpt(sp, yy, h)
    p = np.exp(logp)
    rate = np.gradient(p, ns)
    return sp, ns, p, rate, L1


def gaussian_fraction(sp, n_before_end, Nstar):
    """Gaussian (linearised about +v) volume fraction with phi < 0 for the field smoothed on the scale that
    exited n e-folds before the end: sigma^2 = (3/(8 pi^2 mu2)) (e^{-nu n} - e^{-nu N*}), nu = 2 mu2/3,
    mu2 = w''(y_min) = 2 lam y_min^2 (units H = 1)."""
    ym = sp.ymin()
    mu2 = 2 * sp.lp * ym ** 2
    nu = 2 * mu2 / 3
    var = 3 / (8 * np.pi ** 2 * mu2) * (np.exp(-nu * n_before_end) - np.exp(-nu * Nstar))
    return 0.5 * erfc(ym / np.sqrt(2 * np.maximum(var, 1e-300))), ym ** 2 / (2 * 3 / (8 * np.pi ** 2 * mu2))


def main():
    s1 = load("stochastic_1d.json")
    s2 = load("stochastic_2d.json")
    out = {"module": "inflation_spectator", "part": "3 classification",
           "labels": "computed from the 1D/2D Fokker-Planck results; percolation threshold quoted (Wang et al. "
                     "arXiv:1302.0421, abstract); re-entry mapping is an estimate with instant reheating "
                     "(outside-action preparation data)"}
    # regimes for N* = 50 and 60, delta start (phi = +v)
    reg = {}
    for N in ("N50", "N60"):
        t = s1["thresholds_GeV"][N]["delta_start"]
        tr = s1["thresholds_GeV"][N]["relaxed_start"]
        reg[N] = {
            "A_no_opposite_region_expected": "H_inf < %.4g GeV (all-size count < 1e-3; relaxed start %.4g)"
                                             % (t["H_star_Nreg_all_1e-3"], tr["H_star_Nreg_all_1e-3"]),
            "B_rare_closed_minority_regions": "%.4g < H_inf < %.4g GeV (p_- < 1e-3)" % (
                t["H_star_Nreg_all_1e-3"], t["H_star_p_minus_1e-3"]),
            "C_closed_minority_domains_below_pc": "%.4g < H_inf < %.4g GeV (1e-3 < p_- < p_c)" % (
                t["H_star_p_minus_1e-3"], t["H_perc_p_minus_eq_pc"]),
            "D_both_vacua_percolate_network": "H_inf > %.4g GeV (p_- > p_c = %.4f); p_- -> 1/2 for H_inf >> v" % (
                t["H_perc_p_minus_eq_pc"], PC_SITE_SC),
            "late_reentry_counts": "N(minority domains re-entering after 1 MeV) < 1e-3 needs H_inf < %.4g GeV "
                                   "(delta start) / %.4g GeV (relaxed start)" % (
                                       t["H_star_Nreg_after_1MeV_1e-3"], tr["H_star_Nreg_after_1MeV_1e-3"]),
        }
    out["regimes"] = reg
    # high-H summary
    hi = [r for r in s1["rows"] if r["H_GeV"] >= 1e4]
    out["high_H_summary_1D"] = [{"H_GeV": r["H_GeV"], "p_minus_fresh_N50": r["N50"]["delta_start"]["p_minus"],
                                 "p_minus_fresh_N60": r["N60"]["delta_start"]["p_minus"],
                                 "Lambda1": r["Lambda1_mfpt"],
                                 "sign_field_spectral_tilt_2Lambda1": 2 * r["Lambda1_mfpt"],
                                 "eq_draws_N50": r["N50"].get("equilibrium_draws")} for r in hi[::3]]
    if s2:
        out["high_H_summary_2D"] = s2.get("high_H")
    # size distribution for small p_-: rate profile and re-entry temperatures
    prof = []
    for H in (200.0, 250.0, 300.0):
        sp, ns, p, rate, L1 = rate_profile(H, 50.0)
        Tr = cosmo.T_reheat_instant(H)
        entries = []
        for nb in (0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 45.0, 50.0):     # e-folds before the end of inflation
            n_after = 50.0 - nb                                          # e-folds after the patch start
            i = int(round(n_after / 0.25))
            i = min(max(i, 1), len(ns) - 2)
            k_rel = np.exp(-nb)
            Tre = cosmo.T_reentry(k_rel, H)
            entries.append({"efolds_before_end": nb, "comoving_size_over_end_Hubble": float(np.exp(nb)),
                            "T_reentry_GeV": Tre, "RD_mapping_valid": bool(Tre > 1e-8),
                            "dp_dn_delta_start": float(rate[i]), "dp_dn_relaxed": L1 / 2,
                            "domains_per_Hubble_volume_at_reentry_per_efold_delta": float(rate[i]),
                            "wall_energy_fraction_one_wall_per_Hubble_patch": SIGMA_TREE / (MP ** 2 * cosmo.H_of_T(Tre))})
        gf = [gaussian_fraction(sp, nb, 50.0)[0] for nb in (0.0, 5.0, 10.0)]
        prof.append({"H_GeV": H, "T_reheat_instant_GeV": Tr, "Lambda1": L1, "p_minus_N50_delta": float(p[-1]),
                     "profile": entries,
                     "gaussian_approx_fraction_smoothed_at_n_0_5_10": gf,
                     "gaussian_exponent_over_HM_exponent": gaussian_fraction(sp, 0, 50)[1] / sp.barrier_B()})
    out["minority_domain_profiles"] = prof
    out["notes"] = [
        "Stationary barrier crossing creates minority domains at a constant rate per Hubble volume per e-fold, "
        "so the minority volume fraction per unit ln k is flat (= Lambda_1/2) and the number per Hubble volume at "
        "re-entry is ~ Lambda_1/2 per e-fold of scale. A fresh start at +v suppresses the largest scales (the "
        "first 20-45 e-folds after the patch start).",
        "The Gaussian approximation (linearised about +v) puts all opposite-sign volume on scales that exited in "
        "the last ~3H^2/(2 m^2) e-folds and has an exponent 4x the Hawking-Moss/stochastic barrier exponent: it "
        "misses trapping in the other well and underestimates p_- by many orders of magnitude where p_- is small.",
        "Re-entry: in RD with constant g_*, T_re = T_reheat e^{-n} for a scale that exited n e-folds before the end.",
        "Site percolation p_c = 0.3116 is for uncorrelated sites; the inflationary sign field is correlated "
        "(power-law, tilt 2 Lambda_1), so p_c is only a proxy for the network/no-network boundary.",
    ]
    dump("classification.json", out)
    return out


if __name__ == "__main__":
    o = main()
    print(json.dumps(o["regimes"], indent=1))
    for pr in o["minority_domain_profiles"]:
        print(pr["H_GeV"], pr["Lambda1"], pr["p_minus_N50_delta"], pr["gaussian_approx_fraction_smoothed_at_n_0_5_10"],
              pr["gaussian_exponent_over_HM_exponent"])
        for e in pr["profile"]:
            print("   ", e["efolds_before_end"], "%.3e" % e["T_reentry_GeV"], "%.3e" % e["dp_dn_delta_start"],
                  "%.3e" % e["dp_dn_relaxed"], "%.3e" % e["wall_energy_fraction_one_wall_per_Hubble_patch"])
