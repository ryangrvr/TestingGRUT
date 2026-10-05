#!/usr/bin/env python3
"""BRI1 V3 redo — the analytical signal at t=0.05 is ~1e-13, far below the
Monte Carlo noise (~1e-3). t_num must be increased so the signal dominates.

Preregistered for THIS run (after the diagnostic, before any statistics):
  t_num = 1.0   (in the small-t regime where the t^7 term dominates; the
                t^7 coefficient dominates the t^8 correction for t < |C7|/(2CR),
                which is existential; the explicit safe choice here is
                justified separately below)
  N_B grid: {4, 8, 16, 32, 64}
  samples: 40000
  seed: 20260101
  integrator: RK4, dt = 1/512
  statistic: kappa3 and gamma1 of F(t_num)
  slope: OLS log|g1| vs log N_B

Justification for t_num=1.0:
  At t_num=1.0 the leading-order prediction is kappa3(F_P1) ~ K(1)/N_B with
  K(1) = 3*c(1).  c(t) ~ C7 t^7 with C7 = -0.313/(28*pi^3) = -3.58e-4.
  So c(1) ~ -3.58e-4 and K(1) = 3*c(1) ~ -1.07e-3.
  kappa3(F_P1,N(1)) ~ -1.07e-3 / N_B.
  For N_B=4: kappa3 ~ -2.7e-4.  For N_B=64: ~ -1.7e-5.
  The Monte Carlo noise (std of kappa3 estimator) at 40000 samples is ~
  sqrt(15/40000)*var^{3/2} ~ 0.019*0.32 ~ 6.2e-3, which is LARGER than
  the signal at N_B=64. So we still won't see a clean signal.
  -> We need many more samples, OR a larger t_num where c(t) is larger.

  Better choice: t_num = 2.0.  c(2) ~ C7 * 2^7 = -3.58e-4 * 128 = -4.58e-2
  (rough, using only leading term — but the t^8 correction and higher terms
  may be significant at t=2).  Still, the sign should be negative.
  kappa3(F_P1,2) ~ 3*c(2)/N_B ~ -0.137/N_B.
  For N_B=4: ~ -3.4e-2.  For N_B=64: ~ -2.1e-3.
  Monte Carlo noise at 40000 samples: ~ 6.2e-3.  Signal-to-noise at N_B=4:
  ~5.5. At N_B=64: ~0.34.  Borderline.

  Better: t_num=3.0?  But 3.0 > pi, outside the small-t regime where the
  t^7 expansion is valid.  Need to stay within the small-t regime.
  Actually the V1 proof says c(t) = C7 t^7 + R(t) with |R| <= C_R t^8,
  and delta = min(1, |C7|/(2 C_R)). Since C_R is not computed explicitly,
  the V1 proof doesn't give a numerical delta. But t=1.0 is a "safe" choice
  for the leading term (t^7 vs t^8 correction is a 1/t ratio = 1, so the
  correction is comparable to the leading term at t=1). The *sign* of c(t)
  for small enough t is guaranteed negative by the V1 proof; for t=1.0
  it's not guaranteed by the proof alone.
  
  Path A requires the leading t^7 term to dominate. At t_num, the ratio of
  the t^8 correction to the t^7 term is ~ (C_R/C_7) * t. For this to be
  small, need t << |C7|/C_R. Since neither is computed numerically, we
  adopt Path B: illustrative coefficient-level check.
  
REVISED PREREGISTRATION (amendment 2, before any statistics):
  t_num = 1.0, Path B (illustrative consistency check, not claiming membership
  in the proven delta interval).
  N_B grid {4, 8, 16, 32, 64, 128}
  samples: 200,000 (increased for signal-to-noise)
  seed: 20260101
"""
import numpy as np, json, sys

SEED = 20260101
T_NUM = 1.0
NB_GRID = [4, 8, 16, 32, 64, 128]
N_SAMPLES = 200000
DT = 1.0 / 256
STEPS = int(round(T_NUM / DT))

def q_p1(t):
    u = t / np.pi
    return 10*u**3 - 15*u**4 + 6*u**5 if u <= 1 else 1.0

def sample_x_rejection(n, rng):
    n_propose = int(n * 2.2)
    x = rng.normal(0.0, 1.0, n_propose)
    u = rng.uniform(0.0, 1.0, n_propose)
    return x[u < np.exp(-x**4 / 4.0)][:n]

def run_protocol(n_draws, nb, q_func, rng):
    x0 = sample_x_rejection(n_draws * nb, rng)
    p0 = rng.normal(0.0, 1.0, n_draws * nb)
    x0 = x0.reshape(n_draws, nb); p0 = p0.reshape(n_draws, nb)
    eps = 1.0 / np.sqrt(nb)
    for step in range(STEPS):
        t0 = step * DT
        q0 = q_func(t0); qh = q_func(t0 + 0.5*DT); q1 = q_func(t0 + DT)
        xh = x0 + 0.5*DT*p0
        v1 = -x0 - x0**3 + eps*q0
        x2 = x0 + 0.5*DT*(p0 + 0.5*DT*v1)
        v2 = -x2 - x2**3 + eps*qh
        x3 = x0 + 0.5*DT*(p0 + 0.5*DT*v2)
        v3 = -x3 - x3**3 + eps*qh
        x4 = x0 + DT*(p0 + DT*v3)
        v4 = -x4 - x4**3 + eps*q1
        x0 = x0 + (DT/6.0)*(p0 + 2*(p0 + 0.5*DT*v1) + 2*(p0 + 0.5*DT*v2) + (p0 + DT*v3))
        p0 = p0 + (DT/6.0)*(v1 + 2*v2 + 2*v3 + v4)
    return eps * x0.sum(axis=1)

def stats(a):
    mu = a.mean(); k3 = ((a-mu)**3).mean()
    return k3, k3/a.std()**3, a.var()

rng = np.random.default_rng(SEED)
print(f"Amended preregistration: t_num={T_NUM}, NB={NB_GRID}, samples={N_SAMPLES}, seed={SEED}, dt={DT}, steps={STEPS}")
print(f"{'N_B':>5} {'k3_P0':>13} {'g1_P0':>13} {'k3_P1':>13} {'g1_P1':>13} {'var0':>9} {'var1':>9} {'k3_P1*N_B':>13}")
log_nb, log_g = [], []
for nb in NB_GRID:
    F0 = run_protocol(N_SAMPLES, nb, lambda t: 0.0, rng)
    F1 = run_protocol(N_SAMPLES, nb, q_p1, rng)
    k0, g0, v0 = stats(F0); k1, g1, v1 = stats(F1)
    print(f"{nb:>5} {k0:>13.4e} {g0:>13.4e} {k1:>13.4e} {g1:>13.4e} {v0:>9.4f} {v1:>9.4f} {k1*nb:>13.4e}")
    log_nb.append(np.log(nb)); log_g.append(np.log(abs(g1)))
slope, inter = np.polyfit(np.array(log_nb), np.array(log_g), 1)
print(f"\nOLS log|g1_P1| = {slope:.3f}*log(N_B) + {inter:.3f};  implied p = {-slope:.3f}")
json.dump({"nb": NB_GRID, "slope": float(slope), "implied_p": float(-slope), "t_num": T_NUM,
           "samples": N_SAMPLES, "seed": SEED}, open("v3_results.json","w"), indent=2)
print("saved v3_results.json")
