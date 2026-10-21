"""(1) T_par with Claude's own cosmo_observable.core.history (imported read-only from a scratch extraction of c7b8480)
for area-count A and energy-weighted A = R A, including R computed from a downward-shifted s_E (estimator read-high
case).  Checks T_par proportional to R^(1/2).
(2) Codex charged-run velocity proxies straight from the raw CSVs in raw_diagnostics.tar.gz (own parser; no reuse of
audit/estimators/code/proxy_series.py).  Formation = first sample where both raw vacuum fractions (|S/a| > 0.5)
reach 0.1, and alternatively Codex's RESULTS.json formation a_first_sample.

Usage: python3 -I -B verify_tpar_proxies.py <claude_packet> <codex_packet> <raw_dir> <out_json>
"""
import csv
import json
import math
import sys


def tpar(packet):
    sys.path.insert(0, packet + "/modules/cosmo_observable/code")
    import core
    s = 0.456 ** 2
    Rs = {"area": 1.0, "E_lo": (1 - s) ** -0.5, "E_hi": 1 / (1 - s),
          "E_lo_sE/1.10": (1 - s / 1.10) ** -0.5, "E_hi_sE/1.10": 1 / (1 - s / 1.10)}
    rows = []
    for A in (0.7, 0.8, 0.9):
        for conv in ("t", "H"):
            for br in (False, True):
                row = {"A": A, "conv": conv, "br": br}
                for lab, R in Rs.items():
                    row[lab] = core.history(core.SIGMA_BENCH, A=A * R, backreaction=br, convention=conv)["T_switch_GeV"] * 1e6
                row["ratio_E_lo"] = row["E_lo"] / row["area"]
                row["ratio_E_hi"] = row["E_hi"] / row["area"]
                rows.append(row)
    rng = {lab: [min(r[lab] for r in rows), max(r[lab] for r in rows)] for lab in Rs}
    ratios = [min(min(r["ratio_E_lo"], r["ratio_E_hi"]) for r in rows), max(max(r["ratio_E_lo"], r["ratio_E_hi"]) for r in rows)]
    return {"R": Rs, "sqrtR": {k: math.sqrt(v) for k, v in Rs.items()}, "rows": rows, "range_keV": rng,
            "T_ratio_range": ratios}


def proxies(codex, raw):
    res_codex = json.load(open(codex + "/RESULTS.json"))["runs"]
    names = ["main_s101", "main_s202", "main_s303", "cold_s101", "kappa4_s101", "dt_half_s101", "grid256_s101",
             "grid384_s101", "large_s101", "gaussian_s101", "gaussian_s202", "gaussian_large_s101"]
    out = {}
    for nm in names:
        rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(open(f"{raw}/{nm}.csv"))]
        a = [r["a"] for r in rows]
        vc = [r["v2_core_proxy"] for r in rows]
        vs = [r["v2_smooth_proxy"] for r in rows]
        f_own = next((i for i, r in enumerate(rows) if r["fplus_05"] >= 0.1 and r["fminus_05"] >= 0.1), None)
        fc = res_codex.get(nm, {}).get("formation")
        a_form_codex = fc.get("a_first_sample") if isinstance(fc, dict) else None
        i_form = f_own if f_own is not None else len(rows)
        pre = [i for i in range(i_form) if vc[i] > 1]
        post = list(range(f_own, len(rows))) if f_own is not None else []
        out[nm] = {
            "final_v2_core": vc[-1], "final_v2_smooth": vs[-1],
            "codex_RESULTS_final_v2_core": res_codex.get(nm, {}).get("final", {}).get("v2_core_proxy"),
            "a_formation_own": a[f_own] if f_own is not None else None, "a_formation_codex": a_form_codex,
            "n_preformation_core_gt1": len(pre), "first_a_core_gt1": a[pre[0]] if pre else None,
            "max_core_pre": max(vc[:i_form]) if i_form > 0 else None,
            "a_max_core_pre": a[vc.index(max(vc[:i_form]))] if i_form > 0 else None,
            "post_core_values": [round(vc[i], 3) for i in post],
            "first_post_core": vc[post[0]] if post else None, "first_post_a": a[post[0]] if post else None,
            "q_fraction_core_final": rows[-1]["q_fraction_core"],
            "max_abs_S_preformation_first_gt1": None,
        }
        if pre:
            r = rows[pre[0]]
            out[nm]["max_abs_S_preformation_first_gt1"] = max(abs(r["min_S"]), abs(r["max_S"]))
            out[nm]["KG_phi_at_first_gt1"] = r["K_phi"] / r["G_phi"] if r["G_phi"] else None
    for nm in ("v_kink_N128", "v_kink_N256"):
        r = [{k: float(v) for k, v in x.items()} for x in csv.DictReader(open(f"{raw}/{nm}.csv"))][0]
        out[nm] = {"v2_core": r["v2_core_proxy"], "v2_smooth": r["v2_smooth_proxy"],
                   "rel_core": r["v2_core_proxy"] / 0.456 ** 2 - 1, "rel_smooth": r["v2_smooth_proxy"] / 0.456 ** 2 - 1}
    return out


def main(packet, codex, raw, out_path):
    res = {"tpar": tpar(packet), "proxies": proxies(codex, raw)}
    json.dump(res, open(out_path, "w"), indent=1)
    t = res["tpar"]
    print(json.dumps({"R": t["R"], "range_keV": t["range_keV"], "T_ratio_range": t["T_ratio_range"]}, indent=1))
    for r in t["rows"]:
        if r["A"] == 0.8 and r["conv"] == "t":
            print(r)
    for k, v in res["proxies"].items():
        print(k, v)


if __name__ == "__main__":
    main(*sys.argv[1:5])
