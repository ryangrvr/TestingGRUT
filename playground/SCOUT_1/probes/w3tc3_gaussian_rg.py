"""SCOUT-1 W3-TC3: what does Gaussian (quadratic) coarse-graining do to the quotient data?
 1. exact decimation (Schur complement) preserves the retained resolvent G_rr(z) exactly -> mu_r invariant.
 2. IR forgetting inside the Gaussian class: NNN couplings with a quadratic band minimum leave gamma = -1/2 (bulk):
    higher-order dispersion data are irrelevant (genuine IR reduction).
 3. but the Gaussian class has a CONTINUUM of fixed points: long-range hopping J(r) ~ r^{-(1+s)} (0<s<2) gives
    dispersion ~ |k|^s and bulk gamma = 1/s - 1, continuous in s; decimation does not move s.
"""
import numpy as np

rng = np.random.default_rng(3)

print("=== 1. Schur-complement decimation preserves G_rr(z) ===")
n = 200
a = 2.3 + 0.2 * rng.standard_normal(n); b = -(1 + 0.2 * rng.random(n - 1))
K = np.diag(a) + np.diag(b, 1) + np.diag(b, -1)
keep = np.arange(0, n, 2)              # decimate every other site (retained site 0 kept)
drop = np.setdiff1d(np.arange(n), keep)
for z in (-0.5 + 0.1j, 1.0 + 0.05j, 3.0 + 0.01j):
    G_full = np.linalg.inv(z * np.eye(n) - K)[0, 0]
    A = K[np.ix_(keep, keep)]; B = K[np.ix_(keep, drop)]; D = K[np.ix_(drop, drop)]
    Keff = A + B @ np.linalg.solve(z * np.eye(len(drop)) - D, B.T)   # energy-dependent effective coupling
    G_dec = np.linalg.inv(z * np.eye(len(keep)) - Keff)[0, 0]
    print(f"  z = {z}: |G_full - G_decimated| = {abs(G_full - G_dec):.2e}")
print("  decimation = exact reduction: the retained quotient mu_r is invariant (no flow of quotient data)")


def bulk_gamma_from_dispersion(eps_k, k):
    """bulk LDOS rho(E) = (1/2pi) int dk delta(E - eps(k)); fit near the band bottom"""
    e = eps_k - eps_k.min()
    W = e.max()
    hist, edges = np.histogram(e, bins=np.logspace(np.log10(W) - 4.5, np.log10(W) - 1.5, 25))
    c = np.sqrt(edges[1:] * edges[:-1]); dens = hist / np.diff(edges)
    m = hist > 50
    return np.polyfit(np.log(c[m]), np.log(dens[m]), 1)[0]


print("\n=== 2. IR forgetting: nearest + next-nearest hopping, bulk gamma (quadratic minimum => -1/2) ===")
k = np.linspace(-np.pi, np.pi, 4_000_001)
for t2 in (0.0, 0.1, 0.2, -0.2):
    eps_k = 2 * (1 - np.cos(k)) + 2 * t2 * (1 - np.cos(2 * k))
    print(f"  t2 = {t2:+.1f}: gamma = {bulk_gamma_from_dispersion(eps_k, k):+.3f}")
eps_k = 2 * (1 - np.cos(k)) - 0.5 * (1 - np.cos(2 * k))   # t2 = -1/4: quadratic term cancels -> quartic minimum
print(f"  t2 = -0.25 (fine-tuned: k^2 term cancels, k^4 minimum): gamma = {bulk_gamma_from_dispersion(eps_k, k):+.3f} (predicted -3/4)")

print("\n=== 3. continuum of Gaussian fixed points: long-range hopping J(r) = r^{-(1+s)} ===")
rmax = 20000
r = np.arange(1, rmax + 1)
kk = np.linspace(1e-5, np.pi, 200001)
for s in (0.5, 1.0, 1.5, 1.8):
    J = r ** (-(1.0 + s))
    eps_k = np.array([2 * np.sum(J * (1 - np.cos(q * r))) for q in kk[::400]])
    kq = kk[::400]
    sl = np.polyfit(np.log(kq[1:40]), np.log(eps_k[1:40]), 1)[0]
    print(f"  s = {s}: dispersion exponent ~ {sl:.3f} (predicted s) -> bulk gamma = 1/s - 1 = {1 / s - 1:+.3f}")
print("  decimating every other site of a power-law chain leaves J_eff(r) ~ r^{-(1+s)} (same tail): s is a")
print("  fixed-point LABEL, not a flowing coupling. The Gaussian class is a continuous family of fixed points.")
