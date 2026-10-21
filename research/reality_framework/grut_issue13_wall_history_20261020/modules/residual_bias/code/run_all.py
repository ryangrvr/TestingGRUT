"""Run every part of module residual_bias.  numpy/scipy only, 2 threads, deterministic.

    cd research/reality_framework/grut_issue13_wall_history_20261020/modules/residual_bias/code
    OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1 python3 run_all.py

Wall times on the shared 4-CPU machine: a_damping about 75 s, b_inflation about 11 min, d_charged about 12 min
(two single-threaded workers; existing run files in results/d_runs are reused), e_synthesis 1 s.
"""
import hashlib
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(os.path.dirname(HERE), "results")
env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
times = {}


def run(args):
    t0 = time.time()
    subprocess.run([sys.executable] + args, cwd=HERE, env=env, check=True)
    times[" ".join(args)] = round(time.time() - t0, 1)


run(["a_damping.py"])
run(["b_inflation.py"])
v = subprocess.run([sys.executable, "d_charged.py", "validate"], cwd=HERE, env=env, check=True,
                   capture_output=True, text=True).stdout
with open(os.path.join(RES, "d_validation.json"), "w") as fh:
    fh.write(v)
t0 = time.time()
procs = [subprocess.Popen([sys.executable, "d_charged.py", "worker", str(i), "2"], cwd=HERE, env=env)
         for i in range(2)]
for p in procs:
    assert p.wait() == 0
times["d_charged workers"] = round(time.time() - t0, 1)
run(["e_synthesis.py"])
sha = {}
for fn in sorted(os.listdir(HERE)):
    if fn.endswith(".py"):
        with open(os.path.join(HERE, fn), "rb") as fh:
            sha[fn] = hashlib.sha256(fh.read()).hexdigest()
with open(os.path.join(RES, "runtimes.json"), "w") as fh:
    json.dump(dict(wall_s=times, sha256=sha), fh, indent=1)
print(times)
