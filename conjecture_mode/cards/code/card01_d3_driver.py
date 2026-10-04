#!/usr/bin/env python3
"""CARD-01 v1 D3 driver: runs the frozen job list with N parallel single-thread workers; resumable
(skips jobs whose JSON exists).  Phases, in order of verdict priority:
  P1  PRIMARY: eps grid x 3 starts, CPL comparator x 3
  P2  PRIMARY: free-eps minimization x 3 seeded at the best grid eps   (needs P1)
  P3  extensions (Pantheon+, Union3, DESY5): grid x 3, CPL x 3, then free x 3
  P4  CONTROL (PRIMARY combination, eps > 0 grid x 3) — CONTROL — NOT A GRUT CLAIM
usage: card01_d3_driver.py OUTDIR LOGDIR NWORKERS
"""
import glob, json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import card01_run_config as RC  # noqa: E402

OUT, LOGS, NW = sys.argv[1], sys.argv[2], int(sys.argv[3])
PY = sys.executable
EXTS = ["EXT_PANTHEONPLUS", "EXT_UNION3", "EXT_DESY5"]


def tag(combo, model, eps, s):
    return f"{combo}__{model}__{eps:+.3f}__s{s}"


def run_job(job):
    combo, model, eps, s = job
    t = tag(combo, model, eps, s)
    if os.path.exists(os.path.join(OUT, t + ".json")):
        return t, 0
    env = dict(os.environ, OMP_NUM_THREADS="1")
    with open(os.path.join(LOGS, t + ".log"), "w") as lf:
        rc = subprocess.call([PY, os.path.join(HERE, "card01_d3_run.py"), combo, model, repr(eps), str(s), OUT],
                             stdout=lf, stderr=subprocess.STDOUT, env=env)
    print(f"done {t} rc={rc}", flush=True)
    return t, rc


def run_all(jobs):
    with ThreadPoolExecutor(NW) as ex:
        return list(ex.map(run_job, jobs))


def best_grid_eps(combo):
    best = None
    for f in glob.glob(os.path.join(OUT, f"{combo}__card__*.json")):
        r = json.load(open(f))
        if r["eps_input"] <= 0 and (best is None or r["chi2_eff"] < best[0]):
            best = (r["chi2_eff"], r["eps_input"])
    return best[1] if best else 0.0


def grid_jobs(combo, grid):
    return [(combo, "card", e, s) for e in grid for s in range(RC.N_STARTS)]


def main():
    os.makedirs(OUT, exist_ok=True); os.makedirs(LOGS, exist_ok=True)
    P = RC.PRIMARY
    # P1
    run_all([(P, "cpl", 0.0, s) for s in range(RC.N_STARTS)] + grid_jobs(P, RC.EPS_GRID_PRIMARY))
    # P2
    e0 = best_grid_eps(P)
    run_all([(P, "free", e0, s) for s in range(RC.N_STARTS)])
    print(f"PHASE P1+P2 COMPLETE (primary); best grid eps = {e0}", flush=True)
    # P3
    jobs = []
    for c in EXTS:
        jobs += [(c, "cpl", 0.0, s) for s in range(RC.N_STARTS)] + grid_jobs(c, RC.EPS_GRID_PRIMARY)
    run_all(jobs)
    run_all([(c, "free", best_grid_eps(c), s) for c in EXTS for s in range(RC.N_STARTS)])
    print("PHASE P3 COMPLETE (extensions)", flush=True)
    # P4 control
    run_all(grid_jobs(P, [e for e in RC.EPS_GRID_CONTROL if e > 0]))
    print("PHASE P4 COMPLETE (control)", flush=True)
    print("ALL PHASES COMPLETE", flush=True)


if __name__ == "__main__":
    main()
