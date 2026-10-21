"""Run the whole minority_relic module (about 6 min wall with 2 single-threaded workers).

    cd research/reality_framework/grut_issue13_wall_history_20261020/modules/minority_relic/code
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 run_all.py
"""
import os
import sys
import time
import json

sys.dont_write_bytecode = True
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

if __name__ == "__main__":
    t = {}
    import collapse
    t0 = time.time(); collapse.main(); t["collapse"] = time.time() - t0
    import production
    t0 = time.time(); production.main(); t["production"] = time.time() - t0
    import thresholds
    t0 = time.time(); thresholds.main(); t["thresholds"] = time.time() - t0
    import headline
    t0 = time.time(); headline.main(); t["headline"] = time.time() - t0
    import checks_mr
    t0 = time.time(); c = checks_mr.main(); t["checks"] = time.time() - t0
    import mr_common
    mr_common.dump("runtimes.json", t)
    print(json.dumps(t, indent=1), c["n_pass"], "/", c["n_checks"], "checks pass")
