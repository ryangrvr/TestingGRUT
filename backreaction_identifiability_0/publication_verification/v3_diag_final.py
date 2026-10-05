import numpy as np
rng = np.random.default_rng(20260101)
def sample_x_rejection(n, rng):
    n_propose = int(n * 2.1)
    x = rng.normal(0.0, 1.0, n_propose)
    u = rng.uniform(0.0, 1.0, n_propose)
    return x[u < np.exp(-x**4 / 4.0)][:n]
def run_protocol(n_draws, nb, q_func, rng):
    x0, p0 = sample_x_rejection(n_draws * nb, rng), rng.normal(0.0, 1.0, n_draws * nb)
    x0 = x0.reshape(n_draws, nb); p0 = p0.reshape(n_draws, nb)
    eps = 1.0 / np.sqrt(nb)
    DT = 1.0/512; STEPS = int(round(0.05/DT))
    for step in range(STEPS):
        t0 = step * DT
        q0 = q_func(t0); qh = q_func(t0 + 0.5*DT); q1 = q_func(t0 + DT)
        k1v = -x0 - x0**3 + eps*q0;              k1x = p0
        k2v = -(x0+0.5*DT*k1x) - (x0+0.5*DT*k1x)**3 + eps*qh; k2x = p0+0.5*DT*k1v
        k3v = -(x0+0.5*DT*k2x) - (x0+0.5*DT*k2x)**3 + eps*qh; k3x = p0+0.5*DT*k2v
        k4v = -(x0+DT*k3x) - (x0+DT*k3x)**3 + eps*q1;          k4x = p0+DT*k3v
        x0 = x0 + (DT/6.0)*(k1x + 2*k2x + 2*k3x + k4x)
        p0 = p0 + (DT/6.0)*(k1v + 2*k2v + 2*k3v + k4v)
    return eps * x0.sum(axis=1), x0, p0

def q_p1(t):
    u = t / np.pi
    return 10*u**3 - 15*u**4 + 6*u**5 if u <= 1 else 1.0

def stats(arr):
    mu = arr.mean(); k3 = ((arr - mu)**3).mean()
    return k3, k3 / arr.std()**3, arr.var()

# ANALYTIC VALIDATION at t_num=0.05, N_B large (epsilon ~ 0): check kappa3(F_P1) vs kappa3(x0) should match at eps=0.
# At eps=0, F_P1 = eps*sum(x0_j) = sqrt(N_B)*eps*x0 -> wait, x_j(0) are drawn iid but F = eps*sum x_j;
# at eps=0, F_P1(t) = eps * sum x0_j(t) (no driving). But k3(F) = eps^3 * N_B * k3(x0(t)) = k3(x0(t)) / N_B^(1/2)?? No:
# F = eps * sum x0_j, eps = 1/sqrt(nb). k3(F) = eps^3 * nb * k3(x0) = eps^3 * nb * 0 = 0 by symmetry.
# So k3(F_P0) should be exactly 0. The measured k3_P0 values are small and random-sign - consistent with 0.
# For P1: the expected k3 at small t should be negative with magnitude ~ t^7 * K ~ 0.05^7 * 0.5 / 28/pi^3 ~ 3e-11, tiny!
# So the k3_P1 measured at t=0.05 should be dominated by noise. Let me check magnitude.
print("Expected analytic |kappa3(F_P1)| ~ Var(x0^2)*t^7/(28*pi^3) at eps -> 0")
t = 0.05
expected = 0.313 * t**7 / (28 * np.pi**3)
print(f"  = {expected:.3e}")
print(f"  With N_B = {4}: expected kappa3(F_P1) ~ kappa3(X^eps) * N_B^(1-3/2) = (eps*K/N_B...) need care")
# From V1: kappa3(F_P1,N(t)) = K(t)/N_B + O(1/N_B^2) where K(t) = 3*c(t) and c(t) = C7*t^7
K = 3 * 0.313 * t**7 / (28 * np.pi**3)
print(f"  K(t) = 3*c(t) = {K:.3e}")
for nb in [4, 64]:
    print(f"  N_B={nb}: predicted kappa3(F_P1) = K/N_B = {K/nb:.3e}")
print(f"\nSampling noise at N_SAMPLES=20000: std(kappa3) ~ N^{-1} * const")
print(f"  For var=0.47, k3 sampling std ~ sqrt(15/N)*var^1.5 ~ {np.sqrt(15/20000)*0.47**1.5:.3e}")
