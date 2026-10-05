#!/usr/bin/env python3
"""CARD #1 v1R driver (PRIMARY only; resumable).  Stage A: all 22 targets.  Stage B: the load-bearing targets,
identified mechanically from the best-of-arms profile; re-identified up to MAX_LB_ITERATIONS times so that any newly
load-bearing target also gets its independent check.  Then STOP FOR OWNER REVIEW (no extensions, no control, no MCMC).
usage: card01_v1r_driver.py V1R_OUT LOGDIR NWORKERS
"""
import os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import card01_v1r_config as C   # noqa: E402

OUT, LOGS, NW = sys.argv[1], sys.argv[2], int(sys.argv[3])


def job(target_arm):
    target, arm = target_arm
    lf = os.path.join(LOGS, f"v1r_{target.replace(':', '_')}_{arm}.log")
    env = dict(os.environ, OMP_NUM_THREADS=str(max(1, (os.cpu_count() or 4) // NW)))
    with open(lf, "a") as f:
        rc = subprocess.call([sys.executable, os.path.join(HERE, "card01_v1r_refine.py"), target, arm, OUT],
                             stdout=f, stderr=subprocess.STDOUT, env=env)
    print(f"done {target} {arm} rc={rc}", flush=True)


def run_all(jobs):
    with ThreadPoolExecutor(NW) as ex:
        list(ex.map(job, jobs))


def main():
    os.makedirs(OUT, exist_ok=True); os.makedirs(LOGS, exist_ok=True)
    run_all([(t, "A") for t in C.TARGETS])
    def tagname(t, arm):
        m, e = (t.split(":")[0], float(t.split(":")[1])) if t.startswith("card:") else (t, 0.0)
        return os.path.join(OUT, f"V1R__{m}__{e:+.3f}__{arm}.json")
    missing = [t for t in C.TARGETS if not os.path.exists(tagname(t, "A"))]
    if missing:
        print(f"STAGE A INCOMPLETE {missing} - stopping (no stage B on an incomplete profile)", flush=True)
        return
    print("STAGE A COMPLETE", flush=True)
    from card01_v1r_analyze import load, load_bearing
    done_b = set()
    for it in range(C.MAX_LB_ITERATIONS):
        lb, _ = load_bearing(load(OUT))
        todo = sorted(lb - done_b)
        print(f"STAGE B iteration {it + 1}: load-bearing {sorted(lb)}; new {todo}", flush=True)
        if not todo:
            break
        run_all([(t, "B") for t in todo])
        done_b |= set(todo)
    print("STOP FOR OWNER REVIEW (v1R primary complete; no extensions/control/MCMC)", flush=True)


if __name__ == "__main__":
    main()
