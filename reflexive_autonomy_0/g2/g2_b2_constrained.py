"""RA0 / G2 B2 repair [RA1-03] -- exact constrained mean-density sup defect.

The original g2_asymptotic.py `mean_defect` maximized each source class independently (an UPPER BOUND). Here the
supremum over pairs x, x' in the same block (same occupancy n_i of every class) is computed exactly under the global
particle-number constraint  sum_i n_i = N = L/2,  0 <= n_i <= |class i|:

    D_j = max_{n: sum n_i = N}  sum_i  g_i(n_i),   g_i(k) = top_k(w|class i) - bottom_k(w|class i),

where w_y = P(particle starting at y lies in target class j at time tau) (exact one-particle heat kernel, SSEP duality).
The joint maximization is an exact dynamic program over classes. Both numbers are reported; the old log is preserved.
Independent code path, not independent reviewer. Illustration.
"""
import numpy as np


def heat_kernel(L, tau):
    q = np.arange(L)
    lam = 1 - (2.0 / L) * (1 - np.cos(2 * np.pi * q / L))
    return np.real(np.fft.ifft(np.exp(tau * np.log(lam.astype(complex)))))


def defects(L, classes, tau):
    k0 = heat_kernel(L, tau); m = classes.max() + 1; N = L // 2
    ub_worst = ex_worst = 0.0
    for j in range(m):
        tgt = np.where(classes == j)[0]
        w = np.array([k0[(tgt - y) % L].sum() for y in range(L)])
        gs = []
        for i in range(m):
            v = np.sort(w[classes == i])
            g = np.concatenate([[0.0], np.cumsum(v[::-1]) - np.cumsum(v)])     # g(k), k = 0..|class|
            gs.append(g)
        ub = sum(g.max() for g in gs)
        dp = np.full(N + 1, -np.inf); dp[0] = 0.0
        for g in gs:                                                           # exact knapsack-style DP
            new = np.full(N + 1, -np.inf)
            for k in range(len(g)):
                if k > N:
                    break
                cand = dp[: N + 1 - k] + g[k]
                new[k:] = np.maximum(new[k:], cand)
            dp = new
        ub_worst = max(ub_worst, ub / len(tgt)); ex_worst = max(ex_worst, dp[N] / len(tgt))
    return ub_worst, ex_worst


if __name__ == "__main__":
    print("B2 repaired: (old unconstrained upper bound) / (exact constrained sup), tau = 0.05 L^3, N = L/2")
    for L in [64, 256, 1024, 4096]:
        tau = 0.05 * L ** 3
        c4 = (np.arange(L) * 4) // L
        ell = int(round(np.sqrt(L))); cs = np.arange(L) // ell
        i4 = np.arange(L) % 4
        r = [defects(L, c, tau) for c in (c4, cs, i4)]
        print(f"  L={L:5d}: 4 contiguous {r[0][0]:.4f} / {r[0][1]:.4f}   sqrt(L) cells {r[1][0]:.4f} / {r[1][1]:.4f}"
              f"   4 interleaved {r[2][0]:.4f} / {r[2][1]:.4f}")
