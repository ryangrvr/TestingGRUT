#!/usr/bin/env python3
"""CARD #1 v1R2: refine ONE start with the unchanged v1R local procedure (BOBYQA rhoend 0.005, official-chain
proposal covariance, cycles until a full cycle improves chi2 < 0.1, cap 6; one subprocess per model build).
usage: card01_v1r2_run.py TAG MODEL EPS START_JSON OUTDIR
  MODEL card (EPS fixed) | free (EPS = start eps; point carries card01_eps) | cpl
"""
import json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import card01_v1r_refine as F    # noqa: E402  (unchanged v1R machinery: build(), one_cycle(), _sub())
import card01_v1r_config as C    # noqa: E402


def main():
    tag, model, eps, start_json, outdir = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4], sys.argv[5]
    out = os.path.join(outdir, tag + ".json")
    if os.path.exists(out):
        return
    os.makedirs(outdir, exist_ok=True)
    point = json.load(open(start_json))
    if model == "free":
        eps = point["card01_eps"]
    hist, t0, prev, conv, row = [], time.time(), None, False, None
    pj = os.path.join(outdir, tag + ".point.tmp.json")
    for cyc in range(1, C.MAX_CYCLES + 1):
        json.dump(point, open(pj, "w"))
        r = F._sub(["--cycle", model, repr(eps), pj])
        chi2, point, row = r["chi2"], r["point"], r["row"]
        hist.append(chi2)
        print(f"{tag} cycle {cyc}: chi2_eff = {chi2:.4f}", flush=True)
        if prev is not None and prev - chi2 < C.CYCLE_IMPROVEMENT_STOP:
            conv = True
            break
        prev = chi2
    os.path.exists(pj) and os.unlink(pj)
    res = {"tag": tag, "model": model, "eps_target": eps if model == "card" else None, "cycle_chi2": hist,
           "chi2_eff": min(hist), "converged_cycles": conv, "n_cycles": len(hist), "point": point,
           "final_row": row, "seconds": time.time() - t0, "start_file": os.path.basename(start_json)}
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps({k: res[k] for k in ("tag", "chi2_eff", "converged_cycles", "n_cycles", "seconds")}), flush=True)


if __name__ == "__main__":
    main()
