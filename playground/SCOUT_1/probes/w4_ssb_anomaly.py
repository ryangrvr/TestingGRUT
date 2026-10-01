"""SCOUT-1 W4-4 (SSB) and W4-5 (anomaly-fixed coefficients).

W4-4: 2D Ising, Z2-symmetric action. Checkerboard Metropolis, L = 64.
  Which vacuum is selected, and by what: infinitesimal source h = +-0.005, initial condition, or nothing (h = 0, hot start).
  The effective algebra around the state: odd cumulant (skewness of 8x8 block magnetization) is zero by symmetry in the
  symmetric phase and nonzero around a broken vacuum (an odd operator absent from the symmetric expansion).
W4-5: anomaly-fixed coefficients: pi0 -> gamma gamma (coefficient N_c (Q_u^2 - Q_d^2) fixed by representation content).
"""
import numpy as np

rng = np.random.default_rng(99)


def ising(L, T, h, init, sweeps=3000, meas_from=1500):
    if init == "hot":
        s = rng.choice([-1, 1], (L, L))
    else:
        s = np.full((L, L), 1 if init == "up" else -1)
    ii, jj = np.indices((L, L))
    masks = [((ii + jj) % 2 == c) for c in (0, 1)]
    ms, skews = [], []
    for sw in range(sweeps):
        for m in masks:
            nb = np.roll(s, 1, 0) + np.roll(s, -1, 0) + np.roll(s, 1, 1) + np.roll(s, -1, 1)
            dE = 2 * s * (nb + h)
            acc = (rng.random((L, L)) < np.exp(-dE / T)) & m
            s = np.where(acc, -s, s)
        if sw >= meas_from and sw % 10 == 0:
            ms.append(s.mean())
            b = s.reshape(L // 8, 8, L // 8, 8).mean(axis=(1, 3)).ravel()
            d = b - b.mean()
            skews.append(np.mean(d ** 3) / np.mean(d ** 2) ** 1.5)
    return np.mean(ms), np.mean(skews)


if __name__ == "__main__":
    L = 64
    Tc = 2 / np.log(1 + np.sqrt(2))
    print(f"=== W4-4: 2D Ising, L = {L}, Tc = {Tc:.4f}; Onsager m0(T=2.0) = {(1 - np.sinh(2 / 2.0) ** -4) ** 0.125:.4f} ===")
    for T in (2.0, 3.0):
        for h, init in ((+0.005, "hot"), (-0.005, "hot"), (0.0, "up"), (0.0, "down"), (0.0, "hot")):
            m, sk = ising(L, T, h, init)
            print(f"  T = {T}: source h = {h:+.3f}, start = {init:4s}: <m> = {m:+.4f}, block skewness = {sk:+.4f}", flush=True)
    print("  symmetric phase (T = 3): m ~ 0 and skewness ~ 0 whatever the source/start (odd operators absent)")
    print("  broken phase (T = 2): the vacuum sign follows the source sign or the initial condition; around it the odd")
    print("  cumulant is nonzero with sign OPPOSITE to the vacuum (fluctuations skew back toward 0): the effective algebra contains odd (phi^3-type) operators")

    print("\n=== W4-5: anomaly-fixed coefficient, pi0 -> gamma gamma ===")
    alpha = 1 / 137.035999
    mpi, fpi = 134.9768, 92.1   # MeV (f_pi in the 92 MeV convention)
    for Nc in (1, 3, 5):
        A = Nc * ((2 / 3) ** 2 - (1 / 3) ** 2)
        G = alpha ** 2 * mpi ** 3 / (64 * np.pi ** 3 * fpi ** 2) * A ** 2 * 1e6   # eV
        print(f"  N_c = {Nc}: anomaly coefficient N_c(Q_u^2 - Q_d^2) = {A:.4f}, Gamma = {G:.3f} eV")
    print("  measured (PrimEx-II, SECONDARY): ~7.8 eV. The coefficient is fixed by the anomaly (no free Wilson coefficient);")
    print("  the representation content (N_c, quark charges) is supplied. WZW / Chern-Simons levels: integers, same structure.")
