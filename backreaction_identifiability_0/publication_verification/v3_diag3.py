import numpy as np
rng = np.random.default_rng(20260101)
# Rejection sampler: target rho(x) ~ exp(-x^2/2 - x^4/4). Use proposal g = N(0, sigma^2) with
# sigma chosen so that ratio <= 1 everywhere. exp(-x^4/4) <= 1 so rho <= N(0,1) density shape.
# Actually: rho(x) = exp(-x^2/2) * exp(-x^4/4). Let g(x) = N(0,1) density. Then
# rho(x)/[C*g(x)] with proper C: rho(x) = exp(-x^2/2) * exp(-x^4/4), and g(x) proportional to exp(-x^2/2).
# So rho(x)/g(x) ∝ exp(-x^4/4) <= 1. So: sample from N(0,1), accept with prob exp(-x^4/4).
n_target = 400_000
n_propose = int(n_target * 2.0)
x = rng.normal(0.0, 1.0, n_propose)
u = rng.uniform(0.0, 1.0, n_propose)
acc = x[u < np.exp(-x**4/4.0)][:n_target]
print(f"accepted {len(acc)} of {n_propose}, target {n_target}")
m2 = (acc**2).mean(); m4 = (acc**4).mean()
print(f"m2 = {m2:.6f}, m4 = {m4:.6f}, m2+m4 = {m2+m4:.6f}")
print(f"Var(x^2) = {m4 - m2**2:.6f}  (should be > 0)")
