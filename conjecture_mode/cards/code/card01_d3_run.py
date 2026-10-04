#!/usr/bin/env python3
"""CARD-01 v1 D3 single minimization job (one start).  Uses ONLY the frozen run configuration.

usage: card01_d3_run.py COMBO MODEL EPS START OUTDIR
  COMBO  : key of card01_run_config.COMBINATIONS
  MODEL  : card  (eps fixed at EPS; eps = 0 is the nested LambdaCDM null)
           free  (eps free on the primary branch [-1, 0]; EPS = starting value)
           cpl   (w0waCDM comparator; EPS ignored)
  START  : start index 0..2 (seeds the random starting point; frozen protocol = 3 starts)

Statistic (owner ruling S2): chi2 = -2 ln L_total + nuisance Gaussian-prior terms, i.e.
  chi2_eff = 2*minuslogpost - 2*sum_{uniform-prior params} ln(width)
so flat cosmological priors contribute nothing; Gaussian likelihood-nuisance priors (Planck calibration /
foreground priors shipped with the likelihoods) are retained, identically for every model in a combination.
"""
import json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import card01_run_config as RC  # noqa: E402

PACKAGES = os.environ["CARD01_PACKAGES"]

# CR-4: start-point distributions (DESI DR2 input-yaml ref widths) so that 3 starts are distinct
REF = {"logA": (3.036, 0.001), "ns": (0.9649, 0.004), "H0": (67.5, 0.5), "ombh2": (0.02237, 0.0001),
       "omch2": (0.12, 0.001), "tau": (0.0544, 0.006), "w": (-1.0, 0.02), "wa": (0.0, 0.05)}


def build_info(combo, model, eps, start):
    params = {}
    for k, v in RC.BASE_PARAMS.items():
        if isinstance(v, dict) and "prior" in v:
            v = dict(v)
            if k in REF:
                v["ref"] = {"dist": "norm", "loc": REF[k][0], "scale": REF[k][1]}
        params[k] = v
    info = {"likelihood": {k: (dict(v) if v else None) for k, v in RC.COMBINATIONS[combo].items()},
            "params": params, "packages_path": PACKAGES, "debug": False}
    if model == "card":
        info["theory"] = RC.theory_card01()
        params["card01_eps"] = float(eps)
    elif model == "free":
        info["theory"] = RC.theory_card01()
        lo, hi = RC.EPS_PRIMARY
        params["card01_eps"] = {"prior": {"min": lo, "max": hi},
                                "ref": {"dist": "norm", "loc": float(eps), "scale": 0.01}, "proposal": 0.01}
    elif model == "cpl":
        info["theory"] = RC.theory_camb()
        for k, v in RC.CPL_PARAMS.items():
            v = dict(v); v["ref"] = {"dist": "norm", "loc": REF[k][0], "scale": REF[k][1]}
            params[k] = v
        info["prior"] = dict(RC.CPL_PRIOR_CONSTRAINT)
    else:
        raise ValueError(model)
    info["sampler"] = {"minimize": {"method": "bobyqa", "best_of": 1, "ignore_prior": False,
                                    "seed": 1000 * start + 17}}
    return info


def uniform_log_widths(model):
    tot = 0.0
    for name, pdf in zip(model.prior.params, model.prior.pdf):
        if getattr(getattr(pdf, "dist", None), "name", "") == "uniform":
            lo, hi = pdf.interval(1.0)
            tot += np.log(hi - lo)
    return tot


def main():
    combo, model_kind, eps, start, outdir = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    from cobaya.run import run
    tag = f"{combo}__{model_kind}__{eps:+.3f}__s{start}"
    out_json = os.path.join(outdir, tag + ".json")
    if os.path.exists(out_json):
        return
    info = build_info(combo, model_kind, eps, start)
    t0 = time.time()
    upd, sampler = run(info, stop_at_error=True)
    mn = sampler.products()["minimum"]
    lw = uniform_log_widths(sampler.model)
    names = list(mn.keys()) if hasattr(mn, "keys") else list(mn.data.columns)
    row = {c: float(mn[c]) for c in names if not c.startswith("weight")}
    chi2_lik = {c: v for c, v in row.items() if c.startswith("chi2__")}
    res = {"tag": tag, "combo": combo, "model": model_kind, "eps_input": eps, "start": start,
           "minuslogpost": row["minuslogpost"], "minuslogprior": row.get("minuslogprior"),
           "uniform_log_width_sum": lw, "chi2_eff": 2 * row["minuslogpost"] - 2 * lw,
           "chi2_breakdown": chi2_lik, "point": row, "seconds": time.time() - t0}
    os.makedirs(outdir, exist_ok=True)
    with open(out_json, "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps({k: res[k] for k in ("tag", "chi2_eff", "seconds")}))


if __name__ == "__main__":
    main()
