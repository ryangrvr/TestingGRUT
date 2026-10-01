"""SCOUT-1 W2-FP: criticality without tuning? BTW sandpile (2D, open boundaries).
Slow drive: one grain at a random site, relax completely (parallel toppling; BTW is abelian, so the avalanche
size s = number of topplings is order-independent).
 1. conservative bulk, infinitely slow drive: P(s) power law with a cutoff s_c growing with L (no tuned parameter).
 2. hostile: bulk dissipation eps (each transferred grain lost with prob eps): cutoff set by eps, not L.
 3. hostile: finite drive rate r (grains keep arriving during avalanches): avalanches merge.
Measured: <s>, s_c = <s^2>/<s>, local slope of the log-binned P(s).
"""
import sys
import numpy as np

rng = np.random.default_rng(5)


def relax(z, eps=0.0, r=0.0, max_iter=200000):
    s = 0
    it = 0
    L = z.shape[0]
    while True:
        t = z >= 4
        n = int(t.sum())
        if r > 0:
            k = rng.poisson(r * L * L)
            if k:
                np.add.at(z, (rng.integers(0, L, k), rng.integers(0, L, k)), 1)
                t = z >= 4
                n = int(t.sum())
        if n == 0:
            return s, False
        s += n
        it += 1
        z -= 4 * t
        ti = t.astype(np.int64)
        if eps == 0.0:
            z[1:, :] += ti[:-1, :]; z[:-1, :] += ti[1:, :]
            z[:, 1:] += ti[:, :-1]; z[:, :-1] += ti[:, 1:]
        else:
            for sl_to, sl_from in (((slice(1, None), slice(None)), (slice(None, -1), slice(None))),
                                   ((slice(None, -1), slice(None)), (slice(1, None), slice(None))),
                                   ((slice(None), slice(1, None)), (slice(None), slice(None, -1))),
                                   ((slice(None), slice(None, -1)), (slice(None), slice(1, None)))):
                src = ti[sl_from]
                keep = rng.random(src.shape) >= eps
                z[sl_to] += src * keep
        if it > max_iter:
            return s, True


def run(L, n_aval, eps=0.0, r=0.0, warm=None):
    z = rng.integers(0, 4, (L, L))
    warm = warm if warm is not None else 3 * L * L
    sizes = []
    for k in range(warm + n_aval):
        i, j = rng.integers(0, L, 2)
        z[i, j] += 1
        s, runaway = relax(z, eps, r)
        if runaway:
            return np.array(sizes), True
        if k >= warm:
            sizes.append(s)
    return np.array(sizes), False


def summarize(label, s):
    s = s[s > 0]
    m1, m2 = s.mean(), (s.astype(float) ** 2).mean()
    bins = np.logspace(0, np.log10(s.max() + 1), 25)
    h, e = np.histogram(s, bins=bins)
    c = np.sqrt(e[1:] * e[:-1]); p = h / np.diff(e) / len(s)
    ok = h > 20
    # slope in the scaling window: between 3 and s_c/10
    sc = m2 / m1
    win = ok & (c > 3) & (c < sc / 10)
    slope = np.polyfit(np.log(c[win]), np.log(p[win]), 1)[0] if win.sum() >= 3 else np.nan
    print(f"  {label:40s} n={len(s):6d}  <s>={m1:9.2f}  s_c=<s^2>/<s>={sc:10.1f}  max={s.max():7d}  P(s) slope (3..s_c/10)={slope:+.3f}")
    return m1, sc, slope


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 30000
    print("=== 1. conservative BTW, infinitely slow drive (no tuned parameter) ===")
    res = {}
    for L in (16, 32, 64):
        s, _ = run(L, n)
        res[L] = summarize(f"L={L}", s)
    Ls = np.array(sorted(res)); sc = np.array([res[L][1] for L in Ls]); m1 = np.array([res[L][0] for L in Ls])
    print(f"  finite-size scaling: s_c ~ L^{np.polyfit(np.log(Ls), np.log(sc), 1)[0]:.3f};"
          f"  <s> ~ L^{np.polyfit(np.log(Ls), np.log(m1), 1)[0]:.3f} (exact BTW: <s> ~ L^2)")

    print("\n=== 2. hostile: bulk dissipation eps (L = 64) ===")
    for eps in (0.01, 0.03, 0.1):
        s, _ = run(64, n // 2, eps=eps)
        summarize(f"eps={eps}", s)
    print("  and at eps = 0.03 for L = 32 vs 64 (cutoff set by eps, not L):")
    for L in (32, 64):
        s, _ = run(L, n // 2, eps=0.03)
        summarize(f"eps=0.03, L={L}", s)

    print("\n=== 3. hostile: finite drive rate r (grains per site per parallel step), L = 32 ===")
    for r in (1e-4, 1e-3, 3e-3):
        s, runaway = run(32, n // 3, r=r)
        if runaway or len(s) == 0:
            print(f"  r={r}: activity never stops (merged / infinite avalanche)")
        else:
            summarize(f"r={r}", s)
