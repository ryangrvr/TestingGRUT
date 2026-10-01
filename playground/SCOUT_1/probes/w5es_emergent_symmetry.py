"""SCOUT-1 W5-ES: emergent symmetry (C3 hostile).

Part 1 (W5-2 / W5-4): O(N) model with cubic anisotropy, one-loop in d = 4 - eps (standard normalization,
  H = u/4! (sum phi^2)^2 + v/4! sum phi_i^4):
    beta_u = -eps u + (N+8)/6 u^2 + u v
    beta_v = -eps v + 3/2 v^2 + 2 u v
  RG eigenvalue of v at the O(N) point: y_v = eps - 2 u* = eps (N - 4)/(N + 8).
  Flows integrated toward the IR (dl > 0 means coarse-graining: dg/dl = -beta_g).
Part 2 (W5-2 lattice): 2D Ising on the square lattice at T_c (Wolff clusters): the C4-only microscopic symmetry
  vs an emergent SO(2) at long distance: C(r) along the axis vs along the diagonal at equal |r|.
  Hostile: uniaxial anisotropy Jy/Jx = 0.5 at its own T_c: the anisotropy is REDUNDANT (absorbed into a coordinate
  rescaling), so the ratio of correlation lengths survives as geometry.
"""
import numpy as np
from collections import deque
from scipy.integrate import solve_ivp


def part1():
    eps = 1.0

    def flow(l, g, N):
        u, v = g
        bu = -eps * u + (N + 8) / 6 * u ** 2 + u * v
        bv = -eps * v + 1.5 * v ** 2 + 2 * u * v
        return [-bu, -bv]

    print("=== Part 1: O(N) + cubic anisotropy, one loop, eps = 1 ===")
    for N in (2, 3, 3.5, 4.5, 6):
        y = eps * (N - 4) / (N + 8)
        ustar = 6 * eps / (N + 8)
        print(f"  N = {N}: y_v at the O(N) point = {y:+.4f} ({'IRRELEVANT -> O(N) emerges' if y < 0 else 'RELEVANT -> no enhancement'})")
        for u0, v0 in ((0.2, 0.15), (0.2, -0.1), (0.05, 0.3)):
            sol = solve_ivp(flow, (0, 200), [u0, v0], args=(N,), rtol=1e-9, atol=1e-12, dense_output=True,
                            events=lambda l, g, N: abs(g[0]) + abs(g[1]) - 50)
            l_end = sol.t[-1]
            u, v = sol.y[:, -1]
            tag = "RUNAWAY (fluct.-induced 1st order)" if l_end < 199 else f"-> (u, v) = ({u:.4f}, {v:.4f})"
            print(f"     UV (u, v) = ({u0}, {v0:+.2f}): {tag}" + ("" if l_end < 199 else f"   [O(N): u* = {ustar:.4f}, v = 0]"))
    print("  one-loop N_c = 4; higher-loop / bootstrap estimates N_c ~ 2.9 (SECONDARY): for N = 3 the cubic")
    print("  perturbation is in fact weakly RELEVANT (y_v ~ +0.01), i.e. N = 3 is NOT robustly enhanced (contested/marginal)")


def wolff_ising(L, Jx, Jy, T, nclusters, rng):
    s = rng.choice([-1, 1], (L, L))
    px, py = 1 - np.exp(-2 * Jx / T), 1 - np.exp(-2 * Jy / T)
    acc = np.zeros((L, L))
    nmeas = 0
    for c in range(nclusters):
        i0, j0 = rng.integers(0, L, 2)
        sp = s[i0, j0]
        q = deque([(i0, j0)]); s[i0, j0] = -sp
        while q:
            i, j = q.popleft()
            for di, dj, p in ((1, 0, px), (-1, 0, px), (0, 1, py), (0, -1, py)):
                a, b = (i + di) % L, (j + dj) % L
                if s[a, b] == sp and rng.random() < p:
                    s[a, b] = -sp
                    q.append((a, b))
        if c > nclusters // 5 and c % 5 == 0:
            f = np.fft.fft2(s)
            acc += np.real(np.fft.ifft2(f * np.conj(f))) / L ** 2
            nmeas += 1
    return acc / nmeas


def part2():
    rng = np.random.default_rng(3)
    L = 128
    print("\n=== Part 2: 2D Ising at T_c, emergent rotation symmetry (L = 128, Wolff) ===")
    Tc = 2 / np.log(1 + np.sqrt(2))
    C = wolff_ising(L, 1.0, 1.0, Tc, 6000, rng)
    print("  isotropic J: ratio C_axis(r) / C_diag(r) at equal |r| (diag points (k,k), |r| = k sqrt 2; axis interpolated)")
    ax = np.array([C[r, 0] for r in range(L // 2)])
    for k in (1, 2, 4, 8, 16):
        rd = k * np.sqrt(2)
        cax = np.interp(rd, np.arange(L // 2), ax)
        print(f"    |r| = {rd:6.2f}: C_axis = {cax:.4f}, C_diag = {C[k, k]:.4f}, ratio = {cax / C[k, k]:.4f}")
    # anisotropic: Jy = 0.5 Jx; critical point sinh(2Jx/T) sinh(2Jy/T) = 1
    from scipy.optimize import brentq
    Jx, Jy = 1.0, 0.5
    Tca = brentq(lambda T: np.sinh(2 * Jx / T) * np.sinh(2 * Jy / T) - 1, 0.5, 5)
    Ca = wolff_ising(L, Jx, Jy, Tca, 6000, rng)
    print(f"  uniaxial anisotropy Jy/Jx = 0.5 at its T_c = {Tca:.4f}: C along x vs along y at the same r")
    for r in (2, 4, 8, 16):
        print(f"    r = {r:2d}: C_x = {Ca[r, 0]:.4f}, C_y = {Ca[0, r]:.4f}, ratio = {Ca[r, 0] / Ca[0, r]:.4f}")
    print("  -> the anisotropy does NOT wash out: it is redundant (a coordinate rescaling x -> x / xi_x), so the IR is")
    print("     rotation-invariant only in rescaled coordinates; the aspect ratio is kept as geometry, not erased.")


if __name__ == "__main__":
    part1()
    part2()
