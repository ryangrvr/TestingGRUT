# P-15 control D: the boundaries are exactly unattainable (Beta(2,2) entrance); show the MC "absorption" is Euler leakage -> 0 with dt.
import numpy as np
rng = np.random.default_rng(16)
for dt in (2e-4, 5e-5, 1.25e-5):
    n, x, t, dead = 2000, np.full(2000, 0.25), 0.0, np.zeros(2000, bool)
    while t < 5.0:
        a = ~dead
        x[a] = x[a] + 4*(0.5-x[a])*dt + np.sqrt(2*np.clip(x[a]*(1-x[a]), 0, None))*np.sqrt(dt)*rng.standard_normal(a.sum())
        dead |= (x <= 0) | (x >= 1); t += dt
    print(f"dt={dt:.2e}: fraction crossing a boundary by t=5: {dead.mean():.4f}")
