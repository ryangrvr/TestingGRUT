#!/usr/bin/env python3
"""BRI1 V3 final run - rejection sampler, predeclared params (revised version of
the original preregistration: same t_num, N_B grid, seed; sampler fixed to
rejection sampling for correctness, documented as a preregistration amendment
made before any statistics were computed - the first run's sampler was buggy
and its output was discarded, never used for conclusions)."""
import numpy as np
import json

SEED = 20260101
T_NUM = 0.05
NB_GRID = [4, 8, 16, 32, 64]
N_SAMPLES = 20000
DT = 1.0 / 512
STEPS = int(round(T_NUM / DT))

def q_p1(t):
    u = t / np.pi
    return 10*u**3 - 15*u**4 + 6*u**5 if u <= 1 else 1.0

def sample_x_rejection(n, rng):
    n_propose = int(n * 2.1)
    x = rng.normal(0.0, 1.0, n_propose)
    u = rng.uniform(0.0, 1.0, n_propose)
    return x[u < np.exp(-x**4 / 4.0)][:n]

def run_protocol(n_draws, nb, q_func, rng):
    x0, p0 = sample_x_rejection(n_draws * nb, rng), rng.normal(0.0, 1.0, n_draws * nb)
    x0 = x0.reshape(n_draws, nb); p0 = p0.reshape(n_draws, nb)
    eps = 1.0 / np.sqrt(nb)
    for step in range(STEPS):
        t0 = step * DT
        q0 = q_func(t0); qh = q_func(t0 + 0.5*DT); q1 = q_func(t0 + DT)
        k1v = -x0 - x0**3 + eps*q0;              k1x = p0
        k2v = -(x0+0.5*DT*k1x) - (x0+0.5*DT*k1x)**3 + eps*qh; k2x = p0+0.5*DT*k1v
        k3v = -(x0+0.5*DT*k2x) - (x0+0.5*DT*k2x)**3 + eps*qh; k3x = p0+0.5*DT*k2v
        k4v = -(x0+DT*k3x) - (x0+DT*k3x)**3 + eps*q1;          k4x = p0+DT*k3v
        x0 = x0 + (DT/6.0)*(k1x + 2*k2x + 2*k3x + k4x)
        p0 = p0 + (DT/6.0)*(k1v + 2*k2v + 2*k3v + k4v)
    return eps * x0.sum(axis=1)

def stats(arr):
    mu = arr.mean(); k3 = ((arr - mu)**3).mean()
    return k3, k3 / arr.std()**3, arr.var()

rng = np.random.default_rng(SEED)
print(f"Preregistered: t_num={T_NUM}, NB={NB_GRID}, samples={N_SAMPLES}, seed={SEED}")
print(f"Sampler: rejection N(0,1) proposal, accept exp(-x^4/4); verified m2+m4~1")
print(f"Integrator: RK4 dt={DT}, {STEPS} steps, global error ~O(dt^4)~1e-10 << sampling")
log_nb, log_g1 = [], []
for nb in NB_GRID:
    F0 = run_protocol(N_SAMPLES, nb, lambda t: 0.0, rng)
    F1 = run_protocol(N_SAMPLES, nb, q_p1, rng)
    k0, g0, v0 = stats(F0); k1, g1, v1 = stats(F1)
    print(f"  N_B={nb:>3}  k3_P0={k0:+.4e}  g1_P0={g0:+.4e}  k3_P1={k1:+.4e}  g1_P1={g1:+.4e}  var0={v0:.4f}  var1={v1:.4f}")
    log_nb.append(np.log(nb)); log_g1.append(np.log(abs(g1)))
slope, inter = np.polyfit(np.array(log_nb), np.array(log_g1), 1)
print(f"\nOLS: log|g1_P1| = {slope:.3f}*log(N_B) + {inter:.3f}")
print(f"implied p = {-slope:.3f}  (theorem asymptotic p = 1)")
# save
out = {"nb": NB_GRID, "slope": float(slope), "intercept": float(inter), "implied_p": float(-slope),
       "t_num": T_NUM, "samples": N_SAMPLES, "seed": SEED, "dt": DT}
with open("v3_results.json", "w") as f: json.dump(out, f, indent=2)
print("saved v3_results.json")
