"""Read Codex's raw diagnostic CSVs (extracted read-only from a69bd86 raw_diagnostics.tar.gz) and tabulate the
v2_core_proxy / v2_smooth_proxy time series of every production run.

Usage: python3 -I proxy_series.py <raw_dir> <out_json>
"""
import csv, json, sys, math, pathlib

RUNS = ["main_s101", "main_s202", "main_s303", "cold_s101", "kappa4_s101", "dt_half_s101", "grid256_s101",
        "grid384_s101", "large_s101", "gaussian_s101", "gaussian_s202", "gaussian_large_s101"]


def load(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return [{k: float(v) for k, v in r.items()} for r in rows]


def main(raw, out):
    raw = pathlib.Path(raw)
    res = {}
    for name in RUNS:
        p = raw / (name + ".csv")
        if not p.exists():
            res[name] = "missing"
            continue
        rows = load(p)
        a = [r["a"] for r in rows]
        vc = [r["v2_core_proxy"] for r in rows]
        vs = [r["v2_smooth_proxy"] for r in rows]
        area = [r["area_core_com_over_V"] for r in rows]
        gp = [r["genuine_plus"] for r in rows]
        gm = [r["genuine_minus"] for r in rows]
        qc = [r["q_fraction_core"] for r in rows]
        Kphi = [r["K_phi"] for r in rows]
        Gphi = [r["G_phi"] for r in rows]
        n = len(rows)
        over_c = [i for i in range(n) if vc[i] > 1]
        over_s = [i for i in range(n) if vs[i] > 1]
        # formation (both raw fractions >= 0.1 at threshold 0.5), as in Codex's criterion
        form = next((i for i in range(n) if rows[i]["fplus_05"] >= 0.1 and rows[i]["fminus_05"] >= 0.1), None)
        post = range(form, n) if form is not None else range(0)
        def stats(v, idx):
            v = [v[i] for i in idx if math.isfinite(v[i])]
            if not v:
                return None
            return {"min": min(v), "max": max(v), "mean": sum(v) / len(v), "n": len(v)}
        res[name] = {
            "n_samples": n,
            "a_range": [a[0], a[-1]],
            "a_formation_first_sample": a[form] if form is not None else None,
            "v2_core_final": vc[-1], "v2_smooth_final": vs[-1],
            "v2_core_max": max(vc), "a_at_v2_core_max": a[vc.index(max(vc))],
            "v2_smooth_max": max(vs), "a_at_v2_smooth_max": a[vs.index(max(vs))],
            "n_core_gt1": len(over_c), "first_a_core_gt1": a[over_c[0]] if over_c else None,
            "n_smooth_gt1": len(over_s), "first_a_smooth_gt1": a[over_s[0]] if over_s else None,
            "v2_core_postformation": stats(vc, post),
            "v2_smooth_postformation": stats(vs, post),
            "area_core_final": area[-1], "genuine_plus_final": gp[-1], "genuine_minus_final": gm[-1],
            "q_fraction_core_final": qc[-1],
            "Kphi_over_Gphi_final_whole_box": Kphi[-1] / Gphi[-1] if Gphi[-1] else None,
            "series_every_10": [[a[i], vc[i], vs[i], area[i]] for i in range(0, n, max(1, n // 25))],
        }
    # planar controls
    for name in ["v_kink_N128", "v_kink_N256"]:
        p = raw / (name + ".csv")
        if p.exists():
            r = load(p)[0]
            res[name] = {"v2_core_proxy": r["v2_core_proxy"], "v2_smooth_proxy": r["v2_smooth_proxy"],
                         "target_0.456^2": 0.456 ** 2,
                         "rel_err_core": r["v2_core_proxy"] / 0.456 ** 2 - 1,
                         "rel_err_smooth": r["v2_smooth_proxy"] / 0.456 ** 2 - 1}
    json.dump(res, open(out, "w"), indent=1)
    for k, v in res.items():
        if isinstance(v, dict) and "v2_core_final" in v:
            print(f"{k:22s} form_a={v['a_formation_first_sample']} core_final={v['v2_core_final']:.3f} "
                  f"core_max={v['v2_core_max']:.3f}@a={v['a_at_v2_core_max']:.2f} n>1={v['n_core_gt1']}/{v['n_samples']} "
                  f"first>1 a={v['first_a_core_gt1']} smooth_max={v['v2_smooth_max']:.3f} n>1={v['n_smooth_gt1']} "
                  f"post_core={v['v2_core_postformation']} area_core={v['area_core_final']:.4f} KG={v['Kphi_over_Gphi_final_whole_box']:.3f}")
        else:
            print(k, v)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
