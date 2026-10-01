"""SCOUT-1 W3-NL: does a nonzero nonlinearity in a conserved stochastic field flow to a universal fixed point,
removing a free dimensionless datum (a classical analogue of dimensional transmutation)?

Exclusion process on a ring (L sites, half filling) <-> single-step interface h_i = sum_{j<=i}(1 - 2 eta_j).
Symmetric (SSEP, linear current)            -> Edwards-Wilkinson class:  W(t) ~ t^{1/4}
Asymmetric (ASEP, nonlinear current, any eps) -> KPZ class:              W(t) ~ t^{1/3}
Checkerboard bond updates; per half-sweep a (1,0) bond hops right w.p. (1+eps)/4, a (0,1) bond hops left w.p. (1-eps)/4.
Reports: local growth exponents d ln W / d ln t, and the crossover time t_x(eps) (KPZ theory: t_x ~ eps^{-4}).
"""
import numpy as np

rng = np.random.default_rng(31)


def run(L, eps, tmax, reps):
    times = np.unique(np.round(np.logspace(0, np.log10(tmax), 30)).astype(int))
    W2 = np.zeros(len(times))
    p, q = (1 + eps) / 4, (1 - eps) / 4
    for _ in range(reps):
        eta = np.zeros(L, dtype=np.int8); eta[::2] = 1
        k = 0
        for t in range(1, tmax + 1):
            for off in (0, 1):
                i = np.arange(off, L, 2); j = (i + 1) % L
                a, b = eta[i], eta[j]
                u = rng.random(len(i))
                right = (a == 1) & (b == 0) & (u < p)
                left = (a == 0) & (b == 1) & (u < q)
                mv = right | left
                eta[i[mv]], eta[j[mv]] = eta[j[mv]], eta[i[mv]]
            if t == times[k]:
                h = np.cumsum(1 - 2 * eta.astype(np.int64))
                W2[k] += h.var()
                k += 1
                if k == len(times):
                    break
    W = np.sqrt(W2 / reps)
    return times, W


if __name__ == "__main__":
    L, tmax, reps = 8192, 4000, 6
    print(f"L = {L}, t_max = {tmax} sweeps, {reps} runs each (flat initial condition)")
    res = {}
    for eps in (0.0, 0.25, 0.5, 1.0):
        t, W = run(L, eps, tmax, reps)
        res[eps] = (t, W)
        lt, lW = np.log(t), np.log(W)
        # local exponent over successive decades
        def beta(t1, t2):
            m = (t >= t1) & (t <= t2)
            return np.polyfit(lt[m], lW[m], 1)[0]
        print(f"  eps={eps:4.2f}: W(t) at t=10,100,1000,4000 = {np.round(np.interp([10, 100, 1000, 4000], t, W), 3)}; "
              f"beta[10-100] = {beta(10, 100):.3f}, beta[100-1000] = {beta(100, 1000):.3f}, beta[1000-4000] = {beta(1000, 4000):.3f}")
    print("  references: EW beta = 0.25, KPZ beta = 1/3 = 0.333")
    print("  KPZ scaling: W^3 / t -> const * eps^... ; nonlinearity strength sets only the amplitude / crossover scale")
    for eps in (0.25, 0.5, 1.0):
        t, W = res[eps]
        print(f"   eps={eps}: W(4000)^3/4000 = {np.interp(4000, t, W) ** 3 / 4000:.4f}")


def extended():
    """small asymmetry: the crossover to KPZ happens later (t_x ~ eps^-4 in KPZ theory)"""
    print("\nextended run: eps = 0.25 and eps = 0.5 to t = 40000 (L = 16384, 3 runs)")
    for eps in (0.25, 0.5):
        t, W = run(16384, eps, 40000, 3)
        lt, lW = np.log(t), np.log(W)
        out = []
        for t1, t2 in ((100, 1000), (1000, 5000), (5000, 15000), (15000, 40000)):
            m = (t >= t1) & (t <= t2)
            out.append(f"beta[{t1}-{t2}] = {np.polyfit(lt[m], lW[m], 1)[0]:.3f}")
        print(f"  eps={eps}: " + ", ".join(out))
