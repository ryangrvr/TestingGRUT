"""Run list and worker for charged_exit_sim.

usage:  python3 runs.py list
        python3 runs.py worker <worker_id> <n_workers>     (single-threaded; runs its share of the queue)
Each job is capped at MAX_WALL seconds; a run that is not finished is checkpointed (scratch directory) and resumed
by a later job.  Results: results/runs/<name>.json (compact time series, floats rounded to 6 significant digits).
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cesim as C  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results", "runs")
CKPT = os.environ.get("CES_CKPT", "/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/"
                      "scratchpad/ces_ckpt")
MAX_WALL = 22 * 60

BASE2 = dict(d=2, N=2048, dx=0.25, dt=0.05, s_c=1.0, eps_S=1e-2, eps_X=1e-3, block=16)


def cfg(name, **kw):
    c = dict(BASE2)
    c.update(kw)
    c["name"] = name
    return c


RUNS = [
    # (A) charged exit from the homogeneous phi = 0 condensate, n_i = kappa n_c at a = 1
    cfg("A_k2_t40_s1", kind="charged_exit", kappa=2.0, tau_i=40.0, a_final=3.0, seed=1, dtau_out=1.0,
        qlist_at_a=[1.6, 2.0, 2.5]),
    cfg("B_t40_s1", kind="neutral", kappa=0.0, tau_i=40.0, a_final=3.0, seed=1, dtau_out=1.0),
    cfg("A_k2_t150_s1", kind="charged_exit", kappa=2.0, tau_i=150.0, a_final=2.5, seed=1, dtau_out=2.5,
        qlist_at_a=[1.6, 2.0]),
    cfg("B_t150_s1", kind="neutral", kappa=0.0, tau_i=150.0, a_final=2.5, seed=1, dtau_out=2.5),
    cfg("C2_forced_k4_t150", kind="memory_forced", kappa=4.0, tau_i=150.0, a_final=2.5, seed=1, dtau_out=2.5),
    cfg("C1_sub_k0.9_t150", kind="memory_sub", kappa=0.9, tau_i=150.0, a_final=2.5, seed=1, dtau_out=2.5),
    cfg("A_k2_t40_s2", kind="charged_exit", kappa=2.0, tau_i=40.0, a_final=3.0, seed=2, dtau_out=1.0),
    cfg("A_k2_t40_s3", kind="charged_exit", kappa=2.0, tau_i=40.0, a_final=3.0, seed=3, dtau_out=1.0),
    cfg("A_k4_t40_s1", kind="charged_exit", kappa=4.0, tau_i=40.0, a_final=3.0, seed=1, dtau_out=1.0),
    cfg("A_k4_t40_s2", kind="charged_exit", kappa=4.0, tau_i=40.0, a_final=3.0, seed=2, dtau_out=1.0),
    cfg("B_t40_s2", kind="neutral", kappa=0.0, tau_i=40.0, a_final=3.0, seed=2, dtau_out=1.0),
    cfg("A_k2_t150_s2", kind="charged_exit", kappa=2.0, tau_i=150.0, a_final=2.5, seed=2, dtau_out=2.5),
    cfg("A_k4_t150_s1", kind="charged_exit", kappa=4.0, tau_i=150.0, a_final=2.8, seed=1, dtau_out=2.5),
    cfg("C2_forced_k2_t150", kind="memory_forced", kappa=2.0, tau_i=150.0, a_final=2.5, seed=1, dtau_out=2.5),
    # H -> 0 continuation: expansion stopped at a = 1.6 (after the exit), then Minkowski evolution
    cfg("F_k2_t150_freeze1.6", kind="charged_exit", kappa=2.0, tau_i=150.0, a_final=None, a_freeze=1.6,
        tau_final=500.0, seed=1, dtau_out=2.5, init_from=os.path.join(CKPT, "state_k2_t150_s1_a1.6.npz"),
        snap_dir=CKPT, snap_at_tau=[240.0, 370.0, 500.0]),
    cfg("F_neutral_t150_freeze1.6", kind="neutral", kappa=0.0, tau_i=150.0, a_final=None, a_freeze=1.6,
        tau_final=500.0, seed=1, dtau_out=2.5),
    # resolution pair: same continuum initial field (N_base = 512), L = 128
    cfg("R_lo_k2_t40", kind="charged_exit", kappa=2.0, tau_i=40.0, a_final=3.0, seed=1, dtau_out=1.0,
        N=512, dx=0.25, dt=0.05, N_base=512),
    cfg("R_hi_k2_t40", kind="charged_exit", kappa=2.0, tau_i=40.0, a_final=3.0, seed=1, dtau_out=1.0,
        N=1024, dx=0.125, dt=0.025, N_base=512),
    # small 3D runs (true Q-balls; 3D percolation)
    cfg("D3_k2_t40", kind="charged_exit", kappa=2.0, tau_i=40.0, a_final=2.2, seed=1, dtau_out=2.0,
        d=3, N=256, dx=0.3, dt=0.06),
    cfg("D3_neutral_t40", kind="neutral", kappa=0.0, tau_i=40.0, a_final=2.2, seed=1, dtau_out=2.0,
        d=3, N=256, dx=0.3, dt=0.06),
    # extra seeds (2D torus topology lottery, sign statistics, coarsening-rate comparison)
    cfg("A_k2_t40_s4", kind="charged_exit", kappa=2.0, tau_i=40.0, a_final=3.0, seed=4, dtau_out=1.0),
    cfg("A_k2_t40_s5", kind="charged_exit", kappa=2.0, tau_i=40.0, a_final=3.0, seed=5, dtau_out=1.0),
    cfg("A_k2_t40_s6", kind="charged_exit", kappa=2.0, tau_i=40.0, a_final=3.0, seed=6, dtau_out=1.0),
    cfg("B_t40_s3", kind="neutral", kappa=0.0, tau_i=40.0, a_final=3.0, seed=3, dtau_out=1.0),
    cfg("B_t40_s4", kind="neutral", kappa=0.0, tau_i=40.0, a_final=3.0, seed=4, dtau_out=1.0),
    cfg("B_t40_s5", kind="neutral", kappa=0.0, tau_i=40.0, a_final=3.0, seed=5, dtau_out=1.0),
    cfg("B_t40_s6", kind="neutral", kappa=0.0, tau_i=40.0, a_final=3.0, seed=6, dtau_out=1.0),
    cfg("B_t150_s2", kind="neutral", kappa=0.0, tau_i=150.0, a_final=2.5, seed=2, dtau_out=2.5),
    cfg("C2_forced_k4_t150_s2", kind="memory_forced", kappa=4.0, tau_i=150.0, a_final=2.5, seed=2, dtau_out=2.5),
    # 'decharged' controls: identical phi configuration after the charged exit, X and its momentum set to zero
    cfg("S_k2_t150_s1_to1.6", kind="charged_exit", kappa=2.0, tau_i=150.0, a_final=1.6, seed=1, dtau_out=2.5,
        save_state=os.path.join(CKPT, "state_k2_t150_s1_a1.6.npz")),
    cfg("DC_k2_t150_s1_from1.6", kind="charged_exit", kappa=2.0, tau_i=150.0, a_final=2.5, seed=1, dtau_out=2.5,
        init_from=os.path.join(CKPT, "state_k2_t150_s1_a1.6.npz"), zero_X=True),
    cfg("F_dc_k2_t150_freeze1.6", kind="charged_exit", kappa=2.0, tau_i=150.0, a_final=None, a_freeze=1.6,
        tau_final=500.0, seed=1, dtau_out=2.5, init_from=os.path.join(CKPT, "state_k2_t150_s1_a1.6.npz"),
        zero_X=True, snap_dir=CKPT, snap_at_tau=[370.0, 500.0]),
    cfg("S_k2_t40_s1_to1.8", kind="charged_exit", kappa=2.0, tau_i=40.0, a_final=1.8, seed=1, dtau_out=1.0,
        save_state=os.path.join(CKPT, "state_k2_t40_s1_a1.8.npz")),
    cfg("DC_k2_t40_s1_from1.8", kind="charged_exit", kappa=2.0, tau_i=40.0, a_final=3.0, seed=1, dtau_out=1.0,
        init_from=os.path.join(CKPT, "state_k2_t40_s1_a1.8.npz"), zero_X=True),
    # charged continuations from the same saved states (identical to the A runs, now with the L_vor diagnostic)
    cfg("CC_k2_t150_s1_from1.6", kind="charged_exit", kappa=2.0, tau_i=150.0, a_final=2.5, seed=1, dtau_out=2.5,
        init_from=os.path.join(CKPT, "state_k2_t150_s1_a1.6.npz")),
    cfg("CC_k2_t40_s1_from1.8", kind="charged_exit", kappa=2.0, tau_i=40.0, a_final=3.0, seed=1, dtau_out=1.0,
        init_from=os.path.join(CKPT, "state_k2_t40_s1_a1.8.npz")),
]


def rnd(x):
    if isinstance(x, float):
        if x == 0 or x != x:
            return x
        return float(f"{x:.6g}")
    if isinstance(x, dict):
        return {k: rnd(v) for k, v in x.items()}
    if isinstance(x, list):
        return [rnd(v) for v in x]
    return x


def do_run(c):
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(CKPT, exist_ok=True)
    path = os.path.join(OUT, c["name"] + ".json")
    if os.path.exists(path):
        return "exists"
    if c.get("a_final") is None:
        c = dict(c)
        c["a_final"] = c["tau_final"] / c["tau_i"]
    tmp = path + ".tmp"
    st = C.run(c, tmp, ckpt_dir=CKPT, max_wall=MAX_WALL, verbose=True, snap_dir=c.get("snap_dir"))
    if st == "done":
        with open(tmp) as f:
            res = json.load(f)
        with open(path, "w") as f:
            json.dump(rnd(res), f, separators=(",", ":"))
        os.remove(tmp)
    return st


if __name__ == "__main__":
    if sys.argv[1] == "list":
        for c in RUNS:
            print(c["name"], os.path.exists(os.path.join(OUT, c["name"] + ".json")))
    elif sys.argv[1] == "worker":
        wid, nw = int(sys.argv[2]), int(sys.argv[3])
        mine = [c for i, c in enumerate(RUNS) if i % nw == wid]
        if len(sys.argv) > 4:
            names = sys.argv[4].split(",")
            byname = {c["name"]: c for c in RUNS}
            mine = [byname[x] for x in names]          # in the order given (state-producing runs first)
        pending = list(mine)
        while pending:
            c = pending.pop(0)
            t = time.time()
            st = do_run(c)
            print(f"JOB {c['name']} -> {st} in {time.time() - t:.0f} s", flush=True)
            if st == "incomplete":
                pending.append(c)
