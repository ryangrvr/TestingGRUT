"""Network job list and runner. Usage:
    python3 run_networks.py list            # print job names
    python3 run_networks.py run [pattern]   # run all (or matching) jobs, 2 worker processes, skip finished
Each job writes results/runs/<name>.json. Each job is capped at 24 min wall (max_seconds)."""
import os, sys, json, fnmatch
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results", "runs")
PS = [0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2]
SEEDS = [1, 2, 3]


def jobs():
    J = []
    # ---- (i) physical equation, 2D (3+1 walls independent of z): alpha = 2, beta = 2, a = tau
    phys2d = dict(d=2, alpha=2, beta=2, tau_a=1.0, tau_i=2.0, s_c=0.25, amp=0.1, courant=1 / 3, n_out=32)
    for s in SEEDS:
        J.append(dict(phys2d, name=f"phys2d_N2048_L64_s{s}", N=2048, L=64.0, tau_end=32.0, N_base=1024, seed=s))
    J.append(dict(phys2d, name="phys2d_N1024_L64_s1", N=1024, L=64.0, tau_end=32.0, N_base=1024, seed=1))
    J.append(dict(phys2d, name="phys2d_N4096_L64_s1", N=4096, L=64.0, tau_end=32.0, N_base=1024, seed=1, n_out=20))
    for p in PS:
        for s in SEEDS:
            J.append(dict(phys2d, name=f"phys2d_bias_N1024_L32_p{p:.2f}_s{s}", N=1024, L=32.0, tau_end=16.0,
                          N_base=1024, seed=s, p=p, stop_when_empty=True, n_out=28))
    # ---- (i') physical equation, 3D
    phys3d = dict(d=3, alpha=2, beta=2, tau_a=1.0, tau_i=2.0, s_c=0.3, amp=0.1, courant=1 / 3, n_out=20)
    for s in SEEDS:
        J.append(dict(phys3d, name=f"phys3d_N256_L24_s{s}", N=256, L=24.0, tau_end=12.0, N_base=128, seed=s))
    J.append(dict(phys3d, name="phys3d_N128_L24_s1", N=128, L=24.0, tau_end=12.0, N_base=128, seed=1))
    # ---- (ii) PRS 2D, beta = 0, alpha = 2 (standard 2+1 PRS) and alpha = 3 (damping of z-invariant 3+1 walls)
    prs2d = dict(d=2, beta=0, tau_i=2.0, s_c=1.5, amp=0.1, courant=1 / 3, n_out=32)
    for al in (2, 3):
        for s in SEEDS:
            J.append(dict(prs2d, name=f"prs2d_a{al}_N2048_L2048_s{s}", alpha=al, N=2048, L=2048.0, tau_end=1024.0,
                          N_base=2048, seed=s))
    J.append(dict(prs2d, name="prs2d_a2_N4096_L2048_s1", alpha=2, N=4096, L=2048.0, tau_end=1024.0, N_base=2048,
                  seed=1, n_out=24))
    J.append(dict(prs2d, name="prs2d_a2_N4096_L4096_s1", alpha=2, N=4096, L=4096.0, tau_end=2048.0, seed=1,
                  n_out=28))
    J.append(dict(prs2d, name="prs2d_a2_N2048_L2048_late_s1", alpha=2, N=2048, L=2048.0, tau_i=20.0,
                  tau_end=1024.0, N_base=2048, seed=1))
    for p in PS:
        for s in SEEDS:
            J.append(dict(prs2d, name=f"prs2d_bias_a2_N1024_L1024_p{p:.2f}_s{s}", alpha=2, N=1024, L=1024.0,
                          tau_end=512.0, seed=s, p=p, stop_when_empty=True, n_out=32))
    # late formation (tau_i = 100 >> s_c): domains form deep inside the horizon, closer to the cosmological
    # hierarchy xi_formation << H^-1; used for the biased-decay law
    for p in PS:
        for s in SEEDS:
            J.append(dict(prs2d, name=f"prs2d_biaslate_a2_N1024_L1024_p{p:.2f}_s{s}", alpha=2, N=1024, L=1024.0,
                          tau_i=100.0, tau_end=512.0, seed=s, p=p, stop_when_empty=True, n_out=32))
    # same, with a damped preparation phase (extra friction 0.5 phi' for tau < 120) that removes the post-quench
    # oscillations which otherwise produce non-wall sign changes in the weakly damped late-start runs
    for p in PS:
        for s in SEEDS:
            J.append(dict(prs2d, name=f"prs2d_biasdamp_a2_N1024_L1024_p{p:.2f}_s{s}", alpha=2, N=1024, L=1024.0,
                          tau_i=100.0, tau_end=512.0, seed=s, p=p, stop_when_empty=True, n_out=32,
                          damp_until=120.0, damp_gamma=0.5))
    # ---- (ii') PRS 3D, alpha = 3, beta = 0
    prs3d = dict(d=3, alpha=3, beta=0, tau_i=2.0, s_c=1.5, amp=0.1, courant=1 / 3, n_out=24)
    for s in SEEDS:
        J.append(dict(prs3d, name=f"prs3d_N256_L256_s{s}", N=256, L=256.0, tau_end=128.0, seed=s))
    J.append(dict(prs3d, name="prs3d_N256_L128_s1", N=256, L=128.0, tau_end=64.0, N_base=128, seed=1))
    J.append(dict(prs3d, name="prs3d_N384_L384_s1", N=384, L=384.0, tau_end=192.0, seed=1, n_out=20))
    for p in PS:
        for s in SEEDS:
            J.append(dict(prs3d, name=f"prs3d_bias_N128_L128_p{p:.2f}_s{s}", N=128, L=128.0, tau_end=64.0,
                          seed=s, p=p, stop_when_empty=True, n_out=24))
    return J


# rough cost ordering (largest first) so that the 2-worker pool finishes evenly
def cost(j):
    import math
    dx = j["L"] / j["N"]
    steps = (j["tau_end"] - j["tau_i"]) / (j["courant"] * dx)
    return steps * j["N"] ** j["d"]


def run_one(j):
    os.environ["OMP_NUM_THREADS"] = "1"
    sys.path.insert(0, HERE)
    import wallsim as ws
    path = os.path.join(OUT, j["name"] + ".json")
    if os.path.exists(path):
        return j["name"], "exists"
    r = ws.evolve(j, outdir=OUT, verbose=True)
    return j["name"], f"done {r['wall_seconds']:.0f}s stopped_early={r['stopped_early']}"


if __name__ == "__main__":
    J = jobs()
    if sys.argv[1] == "list":
        for j in sorted(J, key=cost, reverse=True):
            print(j["name"], f"{cost(j):.2e}")
        sys.exit()
    pat = sys.argv[2] if len(sys.argv) > 2 else "*"
    todo = [j for j in J if fnmatch.fnmatch(j["name"], pat)]
    todo.sort(key=cost, reverse=True)
    os.makedirs(OUT, exist_ok=True)
    with Pool(int(os.environ.get("NWORKERS", "2")), maxtasksperchild=1) as pool:
        for name, status in pool.imap_unordered(run_one, todo):
            print("JOB", name, status, flush=True)
