"""Runs the full pre-registered grid and writes results/*.json.  See ../README.md."""
import json, os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "..", "results")
os.makedirs(RES, exist_ok=True)
BETAS = "0.6,0.9,1.0,1.1"; NTRAJ4 = {"3": "48", "4": "24", "5": "16"}; NTRAJ2 = "64"; TMAX = "100000"
jobs = [
    [sys.executable, "toric4d.py", "3", BETAS, NTRAJ4["3"], TMAX, "1", f"{RES}/toric4d_L3.json"],
    [sys.executable, "toric4d.py", "4", BETAS, NTRAJ4["4"], TMAX, "1", f"{RES}/toric4d_L4.json"],
    [sys.executable, "toric4d.py", "5", BETAS, NTRAJ4["5"], TMAX, "1", f"{RES}/toric4d_L5.json"],
    [sys.executable, "toric2d.py", "4", BETAS, NTRAJ2, TMAX, "1", f"{RES}/toric2d_L4.json"],
    [sys.executable, "toric2d.py", "6", BETAS, NTRAJ2, TMAX, "1", f"{RES}/toric2d_L6.json"],
    [sys.executable, "toric2d.py", "8", BETAS, NTRAJ2, TMAX, "1", f"{RES}/toric2d_L8.json"],
    [sys.executable, "toric2d.py", "12", BETAS, NTRAJ2, TMAX, "1", f"{RES}/toric2d_L12.json"],
    [sys.executable, "toric2d_quantum_otoc.py", "0.3", f"{RES}/toric2d_quantum_otoc.json"],
]
if __name__ == "__main__":
    t0 = time.time()
    for j in jobs:
        print(">>", " ".join(j[1:]), flush=True)
        subprocess.run(j, cwd=HERE, check=True)
    print(f"all done in {time.time()-t0:.0f}s")
