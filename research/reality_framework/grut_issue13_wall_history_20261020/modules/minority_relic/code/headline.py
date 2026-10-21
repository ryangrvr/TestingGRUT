"""Headline table: thresholds per criterion/scenario with systematic ranges, counts and relic at fixed H_inf.
Reads results/thresholds.json and results/per_H_quantities.json; writes results/headline.json."""
import json
import math
import os
import numpy as np
import mr_common as C

LOG_RHOS_MAX = math.log10(C.RHO_OVER_S_MAX)


def main():
    T = json.load(open(os.path.join(C.RESULTS, "thresholds.json")))
    Q = json.load(open(os.path.join(C.RESULTS, "per_H_quantities.json")))
    thr, sysd = T["thresholds_GeV"], T["systematics_thresholds_GeV"]
    crit = ["relic_Omega_le_0.12_f1_CEsphere", "relic_Omega_le_0.12_f0.1_CEsphere", "relic_Omega_le_0.12_f0.01_CEsphere",
            "N_reenter_after_Tdom_iso_lt_1", "Omega_late_enclosed_max_le_Omega_DM",
            "Omega_late_enclosed_max_le_1e-3_Omega_DM", "Omega_late_enclosed_max_le_1e-6_Omega_DM", "N_reenter_after_Trec_lt_1",
            "N_reenter_after_BBN_1MeV_lt_1", "N_all_sizes_in_Hubble_volume_lt_1e-3_(regime_A)",
            "relic_f1_injection_to_0.03eV_(critic_style)"]
    g = lambda t: (round(t["H_GeV"], 1) if (t and t.get("H_GeV")) else None)
    table = {}
    for c in crit:
        row = {}
        for sc in ("fresh_N50", "fresh_N60", "relaxed_N50", "relaxed_N60"):
            base = g(thr[sc + "_instant"][c])
            ex = g(thr[sc + "_exit"][c])
            allv = [base, ex] + [g(sysd[k][c]) for k in sysd if k.startswith(sc + "_") and c in sysd[k]]
            calib = "FX2" if sc.startswith("fresh") else "FX_2Drate"
            row[sc] = {"instant_baseline": base, "exit_convention": ex,
                       "with_X_calibrated (%s)" % calib: g(sysd["%s_instant_%s" % (sc, calib)][c]),
                       "range_over_systematics": [min([v for v in allv if v is not None], default=None),
                                                  max([v for v in allv if v is not None], default=None)],
                       "n_variants_not_reached": sum(1 for v in allv if v is None)}
        table[c] = row
    # quantities at fixed H (interpolated in ln H)
    def at(key, H, qname):
        rows = Q[key]
        Hs = np.array([r["H_GeV"] for r in rows])
        v = np.array([r[qname] if r[qname] is not None else -np.inf for r in rows], float)
        ok = np.isfinite(v)
        return float(np.interp(math.log(H), np.log(Hs[ok]), v[ok])) if ok.sum() > 1 else None
    fixed = {}
    for H in (180.0, 200.0, 220.0, 240.0, 260.0, 300.0, 360.0):
        d = {}
        for key in Q:
            relic = at(key, H, "log10_rho_over_s_unit_Tdom_iso")
            d[key] = {"log10_N_after_Tdom_iso": at(key, H, "log10_N_after_Tdom_iso"),
                      "log10_N_after_BBN_1MeV": at(key, H, "log10_N_after_BBN_1MeV"),
                      "log10_N_after_Trec": at(key, H, "log10_N_after_Trec"),
                      "log10_Omega_h2_collapse_relic_f1_CEsphere": (relic + math.log10(C_E()) + math.log10(0.12)
                                                                    - LOG_RHOS_MAX if relic is not None else None),
                      "log10_Omega_late_max_over_Omega_DM": (at(key, H, "log10_rho_late_over_s_max") - LOG_RHOS_MAX)}
        fixed["H=%g" % H] = d
    out = {"module": "minority_relic", "part": "headline",
           "statement": "Conditional history constraints on H_inf (outside-action preparation data) for the class: "
                        "inflation precedes the hidden history, hidden sector never reheated above T_c, no nonthermal "
                        "restoration, observable patch at phi = +v at the start of the last N_* e-folds (fresh) or in "
                        "the restricted + equilibrium (relaxed); instant reheating unless 'exit convention'. "
                        "Not predictions; not an exclusion or validation of Candidate 1.",
           "collapse_energy": {"kappa_E_RD": T["inputs"]["kappa_E_RD"], "C_E_sphere": T["inputs"]["C_E_sphere"]},
           "thresholds_GeV": table, "quantities_at_fixed_H": fixed,
           "critic_comparison": T["critic_comparison"], "Nstar_reach": T["Nstar_reach"],
           "upper_end_regime_C (percolation, v_inflation)": {"end_of_inflation": C.H_PERC, "after_trapping": C.H_PERC_TRAP}}
    C.dump("headline.json", out)
    return out


def C_E():
    d = json.load(open(os.path.join(C.RESULTS, "collapse.json")))
    return 3 * d["adopted"]["kappa_E_RD"] / (4 * math.pi)


if __name__ == "__main__":
    o = main()
    print(json.dumps(o["thresholds_GeV"], indent=1))
    print(json.dumps(o["quantities_at_fixed_H"]["H=240"], indent=1))
