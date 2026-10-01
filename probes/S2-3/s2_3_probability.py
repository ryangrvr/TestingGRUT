"""SCOUT-2 S2-3: probability from deterministic microdynamics + coarse access?

Access: a fixed coarse partition (cell A = [0, 1/2)).  'Statistics' = frequency of A.
FIREWALL: an invariant measure existing is not a probability rule; we ask what FORCES the frequencies.

1. Doubling map x -> 2x mod 1 (mixing, positive entropy), exact in binary (bits shift).
   Preparations (ensembles of initial points):
     (a) Lebesgue (iid fair bits)                         (b) absolutely continuous density 2x (x = sqrt(u))
     (c) singular invariant Bernoulli(0.3) (iid bits, P(1) = 0.3)   (d) Dirac on the fixed point 0
   -> ensemble frequency of A vs time.
2. Single trajectories (time averages, no ensemble):
     doubling map from a Lebesgue-random point; from a constructed point whose digit blocks grow (0^1 1^2 0^4 1^8 ...)
     (Baire-generic-type: running frequency oscillates forever).
3. Irrational rotation x -> x + phi mod 1 (uniquely ergodic, not mixing):
     time averages from ARBITRARY initial points (no measure needed); ensemble from an a.c. preparation: does NOT relax.
4. Affinity: the ensemble frequency of a mixture of preparations equals the mixture of frequencies (pushforward is
   linear) -> affine structure is inherited from 'preparations are measures' (DEFINITIONAL).
"""
import numpy as np

rng = np.random.default_rng(13)
NB = 200          # bits per point (exact doubling for 200 steps)


def bits_iid(n, p1):
    return (rng.random((n, NB)) < p1).astype(np.int8)


def bits_from_density_2x(n):
    # x = sqrt(u) has density 2x; get 52-bit binary expansion, pad with fair bits (tail irrelevant at early times)
    x = np.sqrt(rng.random(n))
    out = np.zeros((n, NB), dtype=np.int8)
    y = x.copy()
    for k in range(52):
        y *= 2
        b = (y >= 1).astype(np.int8)
        out[:, k] = b
        y -= b
    out[:, 52:] = rng.integers(0, 2, (n, NB - 52))
    return out


def ensemble_freq(B, times):
    # after t doublings, x_t starts with bit t; x_t in A = [0,1/2) iff bit t == 0
    return [1 - B[:, t].mean() for t in times]


if __name__ == "__main__":
    N = 200000
    times = [0, 1, 2, 4, 8, 16, 32, 45]
    print("=== 1. doubling map: ensemble frequency of A = [0, 1/2) after t steps ===")
    preps = {"(a) Lebesgue": bits_iid(N, 0.5), "(b) a.c. density 2x": bits_from_density_2x(N),
             "(c) singular invariant Bernoulli(0.3)": bits_iid(N, 0.3), "(d) Dirac at fixed point 0": np.zeros((N, NB), np.int8)}
    for name, B in preps.items():
        print(f"  {name:40s}: " + " ".join(f"{f:.4f}" for f in ensemble_freq(B, times)))
    print("  -> a.c. preparations forget their shape and converge to 1/2 (Lebesgue/SRB value);")
    print("     singular invariant preparations stay at their own value forever (0.7, 1.0): the limit is preparation-class dependent")

    print("\n=== 2. single trajectories (time averages over 190 steps; exact) ===")
    b = rng.integers(0, 2, NB)
    run = np.cumsum(1 - b[:190]) / np.arange(1, 191)
    print(f"  Lebesgue-random point: running frequency at n = 10, 50, 100, 190: {np.round(run[[9, 49, 99, 189]], 3)}")
    blocks = []
    k, v = 0, 0
    while sum(len(x) for x in blocks) < NB:
        blocks.append([v] * (2 ** k)); v ^= 1; k += 1
    g = np.array([x for blk in blocks for x in blk][:NB])
    run_g = np.cumsum(1 - g[:190]) / np.arange(1, 191)
    ends = np.cumsum([2 ** j for j in range(8)]) - 1
    print("  constructed point (blocks 0^1 1^2 0^4 1^8 ...): running frequency at block ends:",
          np.round(run_g[ends[ends < 190]], 3))
    print("  -> no limiting frequency: such points form a residual (Baire-generic) set. 'Generic' in the topological")
    print("     sense picks NON-convergent points; 'generic' in the Lebesgue sense picks 1/2. The notion of generic is supplied.")

    print("\n=== 3. irrational rotation (uniquely ergodic, not mixing) ===")
    phi = (np.sqrt(5) - 1) / 2
    for x0 in (0.0, 0.123456, 0.5, 1 / 3):
        n = np.arange(200000)
        xs = (x0 + n * phi) % 1.0
        print(f"  time average from x0 = {x0:.6f}: {np.mean(xs < 0.5):.5f}   (every initial point; no measure used)")
    x0s = 0.1 * rng.random(N)                       # a.c. preparation concentrated on [0, 0.1)
    fr = [np.mean(((x0s + t * phi) % 1.0) < 0.5) for t in (0, 1, 2, 5, 10, 50, 100, 1000)]
    print(f"  ensemble from a.c. preparation on [0, 0.1): frequency at t = 0,1,2,5,10,50,100,1000: {np.round(fr, 3)}")
    print("  -> the measure is EARNED (unique ergodicity: time averages converge for every point), but ensembles never relax")

    print("\n=== 4. affinity ===")
    Ba, Bc = preps["(a) Lebesgue"], preps["(c) singular invariant Bernoulli(0.3)"]
    lam = 0.3
    mix = np.vstack([Ba[: int(lam * N)], Bc[: N - int(lam * N)]])
    for t in (0, 8, 45):
        f_mix = ensemble_freq(mix, [t])[0]
        f_aff = lam * ensemble_freq(Ba, [t])[0] + (1 - lam) * ensemble_freq(Bc, [t])[0]
        print(f"  t = {t:2d}: frequency of the mixture = {f_mix:.4f}; mixture of frequencies = {f_aff:.4f}")
    print("  -> affine automatically: preparations are measures and the dynamics acts by linear pushforward (DEFINITIONAL)")
