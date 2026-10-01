"""SCOUT-1 W4-OG sharper diagnostics (the width exponents of w4og_generation.py were too noisy to decide).

D-1  coarse current-density relation: stationary mean flux j(u0) = <J_LF + diffusive part> measured on a grid of
     conserved sectors u0; fit j = c0 + c1 u0 + c2 u0^2 + c3 u0^3. c2 = lambda_2,eff (the GENERATED quadratic current).
     Cases: (lam3, a) = (0.5, 0) Z2 control; (0.5, 0.8) Z2 broken only by noise; (0, 0.8) no nonlinear current.
D-2  KPZ signature: skewness of the time-integrated current Q_i(t) through each bond (EW: 0 by symmetry;
     KPZ: nonzero with sign set by lambda_2). W4-2 sectors u0 = 0, 0.3, 0.6 and the W4-3 cases.
"""
import sys
import numpy as np

rng = np.random.default_rng(5)


def step_flux(u, lam3, a, T, dt):
    up = np.roll(u, -1)
    J = lam3 * u ** 3; Jn = lam3 * up ** 3
    amax = np.maximum(3 * lam3 * u ** 2, 3 * lam3 * up ** 2)
    det = -(up - u) + 0.5 * (J + Jn) - 0.5 * amax * (up - u)
    Tl = T * np.clip(1 + a * 0.5 * (u + up), 0.05, None)
    F = det + np.sqrt(2 * Tl / dt) * rng.standard_normal(len(u))
    return F, det


def mean_current(lam3, a, u0, N=1024, dt=0.02, burn=50.0, tmeas=400.0, T=1.0):
    u = np.full(N, float(u0))
    acc, n = 0.0, 0
    for s in range(int((burn + tmeas) / dt)):
        F, det = step_flux(u, lam3, a, T, dt)
        u = u - dt * (F - np.roll(F, 1))
        if s * dt > burn:
            acc += det.mean(); n += 1
    return acc / n


def q_skew(lam3, a, u0, N=4096, dt=0.02, t=600.0, T=1.0, reps=3):
    sk = []
    for _ in range(reps):
        u = np.full(N, float(u0)); Q = np.zeros(N)
        for s in range(int(t / dt)):
            F, _ = step_flux(u, lam3, a, T, dt)
            u = u - dt * (F - np.roll(F, 1))
            Q += F * dt
        d = Q - Q.mean()
        sk.append(np.mean(d ** 3) / np.mean(d ** 2) ** 1.5)
    return np.mean(sk), np.std(sk) / np.sqrt(reps)


if __name__ == "__main__":
    part = sys.argv[1]
    if part == "D1":
        grid = np.linspace(-0.4, 0.4, 9)
        print("=== D-1 coarse current-density relation j(u0); fit c0 + c1 u + c2 u^2 + c3 u^3 ===")
        for lam3, a, lab in ((0.5, 0.0, "Z2 control"), (0.5, 0.8, "noise breaks Z2"), (0.0, 0.8, "no nonlinear current")):
            j = np.array([mean_current(lam3, a, u0) for u0 in grid])
            c3, c2, c1, c0 = np.polyfit(grid, j, 3)
            print(f"  lam3 = {lam3}, a = {a} ({lab:20s}): j(u0) = {np.round(j, 4)}", flush=True)
            print(f"      fit: c1 = {c1:+.4f}, c2 = lambda_2,eff = {c2:+.4f}, c3 = {c3:+.4f}", flush=True)
    elif part == "D2":
        print("=== D-2 skewness of the integrated current Q(t = 600), N = 4096, 3 runs (mean +- s.e.) ===")
        for lam3, a, u0, lab in ((0.5, 0.0, 0.0, "W4-2 u0 = 0"), (0.5, 0.0, 0.3, "W4-2 u0 = 0.3"),
                                 (0.5, 0.0, 0.6, "W4-2 u0 = 0.6"), (0.5, 0.0, -0.6, "W4-2 u0 = -0.6"),
                                 (0.5, 0.8, 0.0, "W4-3 a = 0.8"), (0.0, 0.8, 0.0, "W4-3 lam3 = 0 control"),
                                 (0.0, 0.0, 0.0, "linear EW reference")):
            m, e = q_skew(lam3, a, u0)
            print(f"  {lab:24s}: skew(Q) = {m:+.4f} +- {e:.4f}", flush=True)
