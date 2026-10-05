#!/usr/bin/env python3
"""
BRI1 V3 - Illustrative numerical reproduction (Path A).
Preregistered before execution. Single file, deterministic seed.

System (frozen X1): x_j'' + x_j + x_j^3 = q(t)/sqrt(N_B), force F = (1/sqrt(N_B)) sum x_j.
Gibbs prep: (x,p) iid, rho ~ exp(-H0), H0 = p^2/2 + x^2/2 + x^4/4.
Protocols: P0 (q=0), P1 (q = s(t/pi) on [0,pi], then 1), s(u)=10u^3-15u^4+6u^5.

Preregistered choices (fixed BEFORE any run):
  - Path A small-time witness: t_num = 0.05 (inside small-t regime; justification in report)
  - N_B grid: {4, 8, 16, 32, 64}
  - samples per (N_B, protocol): 200,000 Gibbs draws, integrated with RK4, dt = 1/512
  - seed: 20260101 (single fixed seed; numpy Generator PCG64)
  - statistic: third central moment of F_q(t_num) and standardized skewness gamma1
  - slope fit: OLS on log|gamma1(P1)| vs log(N_B), weights none
  - exclusion criteria: none (all results reported)
"""
import numpy as np

SEED = 20260101
T_NUM = 0.05
NB_GRID = [4, 8, 16, 32, 64]
N_SAMPLES = 200_000
DT = 1.0 / 512
STEPS = int(round(T_NUM / DT))

def s(u):
    return 10*u**3 - 15*u**4 + 6*u**5

def q_p1(t):
    return s(t/np.pi) if t <= np.pi else 1.0

def sample_gibbs(n, rng, M=12):
    """Gibbs sampler for rho(x,p) ~ exp(-H0). Returns (x, p) arrays of shape (n,)."""
    x = rng.normal(0.0, 1.0, n)
    p = rng.normal(0.0, 1.0, n)
    for _ in range(M):
        # p | x is N(0,1): p ~ N(0,1) independent of x given H0 split. Actually H0 = p^2/2 + V(x);
        # p and x are coupled only through normalization, but exp(-p^2/2) factors out:
        # rho(x,p) = Z^{-1} exp(-p^2/2) exp(-V(x)). p ~ N(0,1) independent of x. But we
        # need consistency for both coordinates. Use 1D slice sampling for x (non-standard 1D density),
        # and exact normal for p.
        # Slice sample x from exp(-V(x)) with V = x^2/2 + x^4/4
        for i in range(n):
            y = -np.log(rng.uniform()) - 0.5*x[i]**2  # auxiliary
            L = -1.0
            while 0.5*L*L + 0.25*L**4 > y:
                L -= 0.5
            R = 1.0
            while 0.5*R*R + 0.25*R**4 > y:
                R += 0.5
            while True:
                xc = rng.uniform(L, R)
                if 0.5*xc*xc + 0.25*xc**4 <= y:
                    break
                if xc < x[i]:
                    L = xc
                else:
                    R = xc
            x[i] = xc
        p = rng.normal(0.0, 1.0, n)
    return x, p

def run_protocol(n, nb, q_func, rng):
    """Simulate nb oscillators for each of n Gibbs draws to time T_NUM, return F(t_num) array."""
    # We simulate n independent "bath instances", each with nb oscillators.
    # For each instance j-th oscillator has own initial (x,p).
    # Initials: shape (n, nb)
    F = np.zeros(n)
    # To avoid huge memory, process in chunks over draws
    chunk = max(1, 2000 // nb)
    for start in range(0, n, chunk):
        m = min(chunk, n - start)
        x, p = sample_gibbs(m * nb, rng)
        x = x.reshape(m, nb)
        p = p.reshape(m, nb)
        eps = 1.0 / np.sqrt(nb)
        # RK4 integrate x'' = -x - x^3 + eps*q(t) as first-order system
        for step in range(STEPS):
            t0 = step * DT
            q0 = q_func(t0)
            qh = q_func(t0 + 0.5*DT)
            q1 = q_func(t0 + DT)
            # state (x, v)
            def deriv(xx, vv, qq):
                return vv, -xx - xx**3 + eps*qq
            k1x, k1v = deriv(x, p, q0)
            k2x, k2v = deriv(x + 0.5*DT*k1x, p + 0.5*DT*k1v, qh)
            k3x, k3v = deriv(x + 0.5*DT*k2x, p + 0.5*DT*k2v, qh)
            k4x, k4v = deriv(x + DT*k3x, p + DT*k3v, q1)
            x = x + (DT/6.0)*(k1x + 2*k2x + 2*k3x + k4x)
            p = p + (DT/6.0)*(k1v + 2*k2v + 2*k3v + k4v)
        F[start:start+m] = eps * x.sum(axis=1)
    return F

def third_central_moment(arr):
    mu = arr.mean()
    return ((arr - mu)**3).mean()

def skewness(arr):
    return third_central_moment(arr) / arr.std()**3

def main():
    rng = np.random.default_rng(SEED)
    print(f"Preregistered: t_num={T_NUM}, NB={NB_GRID}, samples={N_SAMPLES}, seed={SEED}, dt={DT}, steps={STEPS}")
    print(f"{'N_B':>5} {'kappa3_P0':>14} {'gamma1_P0':>14} {'kappa3_P1':>14} {'gamma1_P1':>14} {'var_P0':>12} {'var_P1':>12}")
    log_nb, log_g = [], []
    for nb in NB_GRID:
        F0 = run_protocol(N_SAMPLES, nb, lambda t: 0.0, rng)
        F1 = run_protocol(N_SAMPLES, nb, q_p1, rng)
        k0, k1 = third_central_moment(F0), third_central_moment(F1)
        g0, g1 = skewness(F0), skewness(F1)
        v0, v1 = F0.var(), F1.var()
        print(f"{nb:>5} {k0:>14.6e} {g0:>14.6e} {k1:>14.6e} {g1:>14.6e} {v0:>12.4e} {v1:>12.4e}")
        log_nb.append(np.log(nb))
        log_g.append(np.log(abs(g1)) if g1 != 0 else -np.inf)
    # slope fit
    log_nb = np.array(log_nb)
    log_g = np.array(log_g)
    mask = np.isfinite(log_g)
    if mask.sum() >= 2:
        slope, intercept = np.polyfit(log_nb[mask], log_g[mask], 1)
        print(f"\nOLS fit: log|gamma1_P1| = {slope:.3f} * log(N_B) + {intercept:.3f}")
        print(f"Implied power p = {-slope:.3f} (theorem asymptotic p = 1)")
    # save results
    np.save("v3_results_F0.npy", F0)
    np.save("v3_results_F1.npy", F1)

if __name__ == "__main__":
    main()
