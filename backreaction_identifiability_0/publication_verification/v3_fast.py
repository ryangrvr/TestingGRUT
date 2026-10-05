#!/usr/bin/env python3
"""V3 fast version: reduced samples for feasibility. Same physics, fewer samples.
Preregistered change documented in report: N_SAMPLES = 20000 (10x reduction)."""
import numpy as np

SEED = 20260101
T_NUM = 0.05
NB_GRID = [4, 8, 16, 32, 64]
N_SAMPLES = 20000   # reduced from 200k for compute feasibility; documented
DT = 1.0 / 512
STEPS = int(round(T_NUM / DT))
INTEGRATOR_TOL_NOTE = "RK4 fixed-step, dt=1/512; local error O(dt^5), global O(dt^4) ~ 1e-10; integration error negligible vs sampling"

def s(u): return 10*u**3 - 15*u**4 + 6*u**5

def q_p1(t): return s(t/np.pi) if t <= np.pi else 1.0

def sample_x_p(n, rng):
    """Vectorized rejection sampler for x with density ~ exp(-x^2/2 - x^4/4)."""
    # Use proposal: N(0, sigma) with sigma chosen so the tail is covered.
    # exp(-x^2/2 - x^4/4) <= exp(-x^2/2), so x's density <= N(0,1) up to constant.
    # Rejection: propose x ~ N(0, 1), accept with prob exp(-x^4/4).
    out = np.empty(n)
    filled = 0
    while filled < n:
        m = max(int((n - filled) * 1.3), 100)
        prop = rng.normal(0.0, 1.0, m)
        u = rng.uniform(0.0, 1.0, m)
        acc = prop[u < np.exp(-prop**4 / 4.0)]
        take = min(len(acc), n - filled)
        out[filled:filled+take] = acc[:take]
        filled += take
    p = rng.normal(0.0, 1.0, n)
    return out, p

def run_protocol(nb, q_func, rng):
    n = N_SAMPLES
    eps = 1.0 / np.sqrt(nb)
    x, p = sample_x_p(n * nb, rng)
    x = x.reshape(n, nb); p = p.reshape(n, nb)
    # Integrate all n*nb oscillators simultaneously to T_NUM with RK4 vectorized
    for step in range(STEPS):
        t0 = step * DT; q0 = q_func(t0); qh = q_func(t0 + 0.5*DT); q1 = q_func(t0 + DT)
        def dv(x, v, q): return -x - x**3 + eps*q
        k1v = dv(x, p, q0);          k1x = p
        k2v = dv(x + 0.5*DT*k1x, p + 0.5*DT*k1v, qh); k2x = p + 0.5*DT*k1v
        k3v = dv(x + 0.5*DT*k2x, p + 0.5*DT*k2v, qh); k3x = p + 0.5*DT*k2v
        k4v = dv(x + DT*k3x, p + DT*k3v, q1);          k4x = p + DT*k3v
        x = x + (DT/6.0)*(k1x + 2*k2x + 2*k3x + k4x)
        p = p + (DT/6.0)*(k1v + 2*k2v + 2*k3v + k4v)
    F = eps * x.sum(axis=1)
    return F

def k3(a):
    return ((a - a.mean())**3).mean()

def skew(a):
    return k3(a) / a.std()**3

def main():
    rng = np.random.default_rng(SEED)
    print(f"Preregistered: t_num={T_NUM}, NB={NB_GRID}, samples={N_SAMPLES}, seed={SEED}")
    print(f"Integrator: RK4 fixed dt={DT}, steps={STEPS}. {INTEGRATOR_TOL_NOTE}")
    print(f"{'N_B':>5} {'kappa3_P0':>13} {'gamma1_P0':>13} {'kappa3_P1':>13} {'gamma1_P1':>13} {'var_P0':>11} {'var_P1':>11}")
    log_nb, log_g = [], []
    for nb in NB_GRID:
        F0 = run_protocol(nb, lambda t: 0.0, rng)
        F1 = run_protocol(nb, q_p1, rng)
        print(f"{nb:>5} {k3(F0):>13.5e} {skew(F0):>13.5e} {k3(F1):>13.5e} {skew(F1):>13.5e} {F0.var():>11.4e} {F1.var():>11.4e}")
        log_nb.append(np.log(nb)); log_g.append(np.log(abs(skew(F1))))
    ln, lg = np.array(log_nb), np.array(log_g)
    slope, intercept = np.polyfit(ln, lg, 1)
    print(f"\nOLS: log|gamma1_P1| = {slope:.3f}*log(N_B) + {intercept:.3f}")
    print(f"implied p = {-slope:.3f}  (theorem asymptotic p = 1)")
    np.savez("v3_results.npz", log_nb=ln, log_g=lg, slope=slope, intercept=intercept,
             nb=np.array(NB_GRID), samples=N_SAMPLES, t_num=T_NUM, seed=SEED)

if __name__ == "__main__":
    main()
