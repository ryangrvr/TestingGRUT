"""Run every charge_phase calculation in order (about 6-8 minutes on one core).
Usage:  cd modules/charge_phase/code && OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python3 run_all.py
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
for script in ["t0_phase.py", "thermal_phase.py", "z2_memory.py", "cosmo_map.py", "codex_consistency.py"]:
    t = time.time()
    r = subprocess.run([sys.executable, os.path.join(HERE, script)], cwd=HERE, capture_output=True, text=True)
    print("%-22s exit=%d  %.1fs" % (script, r.returncode, time.time() - t), flush=True)
    if r.returncode != 0:
        print(r.stderr[-3000:])
        sys.exit(r.returncode)
print("all done; results in ../results/")
