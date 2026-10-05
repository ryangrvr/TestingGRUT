#!/usr/bin/env python3
"""CARD #1 v1R — POST-DATA CONVERGENCE VERIFICATION: one target, one arm (CR-5 refinement, owner ruling v1R).

usage: card01_v1r_refine.py TARGET ARM OUTDIR
  TARGET : card:<eps>  (fixed eps; eps = 0 is the LambdaCDM null) | free | cpl
  ARM    : A  start = best valid v1 solution for the target (v1 outputs, d3_out)
           B  distinct admissible start for the independent convergence check:
              card/free -> fresh seeded draw from the v1 ref distributions (free: eps starts at FREE_B_EPS_START) (seed V1R_SEED_B, distinct from the v1
                           seeds 17/1017/2017); cpl -> highest-likelihood official DESI w0wa sample (v1r_cpl_startB.json)
Each cycle: Cobaya BOBYQA, rhoend = 0.005 (v1: 0.05), proposal covariance from the official DESI chains, started
exactly at the current best point.  Cycles repeat until one full cycle improves chi2 by < 0.1, or MAX_CYCLES is hit
(=> target NOT CONVERGED; no state is forced).  chi2_eff is the v1 definition (CR-3), unchanged.
"""
import json, os, sys, time, copy, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import card01_run_config as RC        # noqa: E402  (frozen v1 configuration, unchanged)
import card01_d3_run as V1            # noqa: E402  (frozen v1 info builder and chi2_eff definition)
import card01_v1r_config as C         # noqa: E402

PACKAGES = os.environ["CARD01_PACKAGES"]
NUIS = json.load(open(os.path.join(HERE, "v1r_inputs", "v1r_nuisance_param_defaults.json")))


def v1_best(model, eps):
    pat = os.path.join(C.V1_OUT, f"{RC.PRIMARY}__{model}__{eps:+.3f}__s*.json") if model == "card" else \
        os.path.join(C.V1_OUT, f"{RC.PRIMARY}__{model}__*.json")
    rows = [json.load(open(f)) for f in glob.glob(pat)]
    return min(rows, key=lambda r: r["chi2_eff"])


def start_point(model, eps, arm, sampled):
    if arm == "A":
        p = v1_best(model, eps)["point"]
        return {k: p[k] for k in sampled}
    if model == "cpl":
        st = json.load(open(os.path.join(HERE, "v1r_inputs", "v1r_cpl_startB.json")))["start"]
        return {k: st[k] for k in sampled}
    rng = np.random.default_rng(C.V1R_SEED_B + int(round(1000 * abs(eps))) + (500000 if model == "free" else 0))
    out = {}
    for k in sampled:
        if k in V1.REF:
            out[k] = rng.normal(*V1.REF[k])
        elif k == "card01_eps":
            out[k] = C.FREE_B_EPS_START     # fixed distinct admissible start for the free-eps arm B
        else:
            r = NUIS[k]["ref"]
            out[k] = rng.normal(r["loc"], r["scale"]) if isinstance(r, dict) else float(r)
    return out


def build(model, eps, point):
    info = V1.build_info(RC.PRIMARY, model, eps if model != "free" else point.get("card01_eps", eps), 0)
    for k, v in point.items():                       # exact deterministic start at `point`
        if k in NUIS:
            d = copy.deepcopy(NUIS[k]); d["ref"] = float(v); info["params"][k] = d
        elif isinstance(info["params"].get(k), dict):
            info["params"][k] = dict(info["params"][k]); info["params"][k]["ref"] = float(v)
    cov = C.COVMAT_W0WA if model == "cpl" else C.COVMAT_LCDM
    info["sampler"] = {"minimize": {"method": "bobyqa", "best_of": 1, "ignore_prior": False,
                                    "override_bobyqa": {"rhoend": C.RHOEND},
                                    "covmat": os.path.join(HERE, "v1r_inputs", cov)}}
    return info


def one_cycle(model, eps, point):
    from cobaya.run import run
    info = build(model, eps, point)
    _, sampler = run(info, stop_at_error=True)
    mn = sampler.products()["minimum"]
    names = list(mn.keys()) if hasattr(mn, "keys") else list(mn.data.columns)
    row = {c: float(mn[c]) for c in names if not c.startswith("weight")}
    chi2 = 2 * row["minuslogpost"] - 2 * V1.uniform_log_widths(sampler.model)
    sampled = list(sampler.model.parameterization.sampled_params())
    return chi2, {k: row[k] for k in sampled}, row


def main():
    target, arm, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    model, eps = (target.split(":")[0], float(target.split(":")[1])) if target.startswith("card:") else (target, 0.0)
    tag = f"V1R__{model}__{eps:+.3f}__{arm}"
    out = os.path.join(outdir, tag + ".json")
    if os.path.exists(out):
        return
    from cobaya.model import get_model
    probe = build(model, eps, {}); probe.pop("sampler")
    sampled = list(get_model(probe).parameterization.sampled_params())
    point = start_point(model, eps, arm, sampled)
    if model == "free" and arm == "A":
        eps = point["card01_eps"]
    hist, t0, prev, conv = [], time.time(), None, False
    for cyc in range(1, C.MAX_CYCLES + 1):
        chi2, point, row = one_cycle(model, eps, point)
        hist.append(chi2)
        print(f"{tag} cycle {cyc}: chi2_eff = {chi2:.4f}", flush=True)
        if prev is not None and prev - chi2 < C.CYCLE_IMPROVEMENT_STOP:
            conv = True
            break
        prev = chi2
    res = {"tag": tag, "model": model, "eps_target": eps if model == "card" else None, "arm": arm,
           "cycle_chi2": hist, "chi2_eff": min(hist), "converged_cycles": conv, "n_cycles": len(hist),
           "point": point, "final_row": row, "seconds": time.time() - t0,
           "start": {"A": "best valid v1 solution", "B": "distinct admissible start"}[arm]}
    os.makedirs(outdir, exist_ok=True)
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps({k: res[k] for k in ("tag", "chi2_eff", "converged_cycles", "n_cycles", "seconds")}), flush=True)


if __name__ == "__main__":
    main()
