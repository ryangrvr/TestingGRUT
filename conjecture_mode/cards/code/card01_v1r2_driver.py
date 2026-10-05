#!/usr/bin/env python3
"""CARD #1 v1R2 — LOWER-MINIMUM REPRODUCIBILITY AUDIT driver (PRIMARY only; final local-optimizer campaign).

Stage A (four core minima: LCDM, CPL, fixed eps=-0.08, free eps):
  round 1: routes 2 and 3 for each target; round 2: route 4 only for targets not yet reproduced.
  If any core target is not reproduced after its 4 genuinely distinct starts (route 1 = retained v1R lowest):
  TERMINAL "CARD-01-v1R2 - COMPUTATIONALLY UNRESOLVED AT AVAILABLE MINIMIZATION METHOD"; stop (no Stage B).
Stage B (only if Stage A passes): continuation sweeps from the branch-minimum grid point toward eps=0 and toward
  eps=-1 (two sweeps in parallel).  At each grid point the start is the converged neighbour solution with eps changed;
  the best existing solution at that eps is retained in the pool.  Then load-bearing verification (branch minimum,
  crossing brackets, LCDM, CPL, free) with the same criterion; extra routes for a failing grid point: reverse-direction
  continuation, then the free-eps low solution projected to that eps; for a failing free fit: the lowest grid solution
  with eps released.  Load-bearing identification is repeated at most twice.  Then STOP FOR OWNER REVIEW.
usage: card01_v1r2_driver.py V1R_OUT V1R2_OUT LOGDIR NWORKERS
"""
import json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import card01_run_config as RC                       # noqa: E402
import card01_v1r2_lib as L                          # noqa: E402
from card01_d3_analyze import profile_set            # noqa: E402

V1R, OUT, LOGS, NW = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
STARTS = json.load(open(os.path.join(HERE, "v1r2_inputs", "stageA_starts.json")))
SDIR = os.path.join(OUT, "starts"); os.makedirs(SDIR, exist_ok=True); os.makedirs(LOGS, exist_ok=True)


def log(m):
    print(m, flush=True)


def run(job):
    tag, model, eps, point = job
    sj = os.path.join(SDIR, tag + ".start.json")
    json.dump(point, open(sj, "w"))
    env = dict(os.environ, OMP_NUM_THREADS=str(max(1, (os.cpu_count() or 4) // NW)))
    with open(os.path.join(LOGS, tag + ".log"), "a") as lf:
        rc = subprocess.call([sys.executable, os.path.join(HERE, "card01_v1r2_run.py"), tag, model, repr(eps), sj, OUT],
                             stdout=lf, stderr=subprocess.STDOUT, env=env)
    log(f"done {tag} rc={rc}")


def run_all(jobs):
    with ThreadPoolExecutor(NW) as ex:
        list(ex.map(run, jobs))


def model_eps(target):
    if target == "cpl":
        return "cpl", 0.0
    if target == "free":
        return "free", 0.0
    return "card", float(target.split(":")[1])


def stageA_job(target, route):
    m, e = model_eps(target)
    s = STARTS[target][route - 1]
    return (f"V1R2A__{target.replace(':', '_')}__r{route}", m, e, s["point"])


def pools():
    return L.pool(V1R, OUT)


def stage_a():
    run_all([stageA_job(t, r) for t in STARTS for r in (2, 3)])
    P = pools()
    todo = [t for t in STARTS if not L.reproduced(P.get(t, []))[0]]
    log(f"STAGE A round 1: not yet reproduced {todo}")
    if todo:
        run_all([stageA_job(t, 4) for t in todo])
    P = pools()
    status = {t: L.reproduced(P.get(t, [])) for t in STARTS}
    for t, (ok, lo, why) in status.items():
        log(f"STAGE A {t}: {'REPRODUCED' if ok else 'NOT REPRODUCED'} - {why}")
    return all(v[0] for v in status.values())


def pt_card(row, eps):
    return {k: v for k, v in row["point"].items() if k != "card01_eps"}


def stage_b():
    P = pools()
    grid = sorted(RC.EPS_GRID_PRIMARY)                      # ascending: -1 ... 0
    prof = {e: L.best(P[f"card:{round(e, 3)}"])["chi2_eff"] for e in grid}
    m = min(grid, key=lambda e: prof[e])
    up = [e for e in grid if e > m]                         # toward 0
    down = [e for e in reversed(grid) if e < m]             # toward -1
    log(f"STAGE B: branch-minimum grid point {m}; sweeps toward 0 {up} and toward -1 {down}")

    def sweep(seq, name):
        prev = m
        for e in seq:
            Pn = pools()
            src = L.best(Pn[f"card:{round(prev, 3)}"])
            run((f"V1R2B__card_{e}__cont_{name}", "card", e, pt_card(src, e)))
            prev = e
    with ThreadPoolExecutor(2) as ex:
        list(ex.map(lambda a: sweep(*a), [(up, "up"), (down, "down")]))
    log("STAGE B sweeps complete")

    extra_done = set()
    for it in range(2):
        P = pools()
        prof = {e: L.best(P[f"card:{round(e, 3)}"])["chi2_eff"] for e in grid}
        eps = np.array(grid); chi = np.array([prof[e] for e in grid])
        i = int(np.argmin(chi)); fb = L.best(P["free"])["chi2_eff"]
        cmin = min(chi[i], fb)
        lo, hi, lo_e, hi_e = profile_set(eps, chi, cmin, RC.PROFILE_SET_DCHI2)
        lb = {f"card:{round(float(eps[i]), 3)}", "card:0.0", "cpl", "free"}
        for val, hit in ((lo, lo_e), (hi, hi_e)):
            if not hit:
                lb |= {f"card:{round(float(e), 3)}" for e in (max([x for x in grid if x <= val], default=None),
                                                              min([x for x in grid if x >= val], default=None)) if e is not None}
        fails = [t for t in sorted(lb) if not L.reproduced(P.get(t, []))[0]]
        free_low = (chi[i] < fb - L.TOL)                    # a fixed-eps point beats the free fit: free not at its minimum
        if free_low and "free" not in fails:
            fails.append("free")
        log(f"LB iteration {it + 1}: load-bearing {sorted(lb)}; failing {fails}")
        jobs = []
        for t in fails:
            if t == "cpl":
                continue                                    # CPL is unaffected by Stage B (verified in Stage A)
            if t == "free":
                if ("free", "rel") not in extra_done:
                    g = L.best(P[f"card:{round(float(eps[i]), 3)}"]); pt = pt_card(g, 0); pt["card01_eps"] = float(eps[i])
                    jobs.append((f"V1R2B__free__release_eps{eps[i]}", "free", 0.0, pt)); extra_done.add(("free", "rel"))
                continue
            e = float(t.split(":")[1]); k = grid.index(e)
            # reverse-direction continuation: the neighbour the sweep did NOT come from (both, at the branch minimum)
            idx = (k + 1,) if e > m else ((k - 1,) if e < m else (k - 1, k + 1))
            nb = [grid[j] for j in idx if 0 <= j < len(grid)]
            for n in nb:
                if (t, n) not in extra_done:
                    jobs.append((f"V1R2B__card_{e}__from_{n}", "card", e, pt_card(L.best(P[f"card:{round(n, 3)}"]), e)))
                    extra_done.add((t, n)); break
            else:
                if (t, "free") not in extra_done:
                    jobs.append((f"V1R2B__card_{e}__from_free", "card", e, pt_card(L.best(P["free"]), e)))
                    extra_done.add((t, "free"))
        if not jobs:
            break
        run_all(jobs)
    log("STOP FOR OWNER REVIEW (v1R2 primary complete)")


def main():
    if not stage_a():
        log("TERMINAL: CARD-01-v1R2 - COMPUTATIONALLY UNRESOLVED AT AVAILABLE MINIMIZATION METHOD (Stage A failed; no Stage B)")
        log("STOP FOR OWNER REVIEW")
        return
    log("STAGE A PASSED - all four core minima reproduced")
    stage_b()


if __name__ == "__main__":
    main()
