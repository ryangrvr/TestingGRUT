"""SCOUT-2 S2-6: can dynamics select basin / preparation / measure data without supplying a measure?

1. Coexisting attractors (deterministic gradient flow x' = x - x^3 + 0.3, attractors near -0.88 and +1.13):
   'probability' of each basin under different reference measures on initial data.
2. Maximum entropy: unconstrained MaxEnt on [0,1] in coordinate x vs y = x^2 -> different densities (Jaynes).
3. Physical measure: logistic map r = 4; two different a.c. preparations converge to the same statistics (arcsine);
   a Dirac preparation on the period-2 orbit never does.
4. Shared reference frames (link to S2-5): local-tomography signal for the J-correlation = 2 <f_A f_B>, where f = +-1 are
   the parties' frame orientations. Frames = spins of a 2D Ising 'reference field' (Metropolis, L = 64):
   distant-frame correlation below vs above T_c. Sign of the order (which vacuum) is irrelevant (J -> -J is complex
   conjugation); only EXISTENCE of long-range order matters.
"""
import numpy as np

rng = np.random.default_rng(21)

print("=== 1. basin probabilities depend on the reference measure ===")
x = rng.uniform(-2, 2, 100000)
x_exp = np.clip(rng.exponential(0.6, 100000) * np.sign(rng.random(100000) - 0.7) , -2, 2)   # an asymmetric reference
for name, x0 in (("Lebesgue on [-2,2]", x), ("asymmetric exponential reference", x_exp)):
    y = x0.copy()
    for _ in range(4000):
        y = y + 0.01 * (y - y ** 3 + 0.3)
    print(f"  {name:34s}: P(basin of +1.13) = {np.mean(y > 0):.4f}")
print("  -> the attractors are dynamical; their 'probabilities' are the reference-measure volumes of the basins (supplied)")

print("\n=== 2. MaxEnt is coordinate (reference-measure) dependent ===")
u = rng.random(1000000)
in_x = u                       # MaxEnt uniform in x
in_y = np.sqrt(u)              # MaxEnt uniform in y = x^2  -> x = sqrt(y)
for name, s in (("uniform in x", in_x), ("uniform in y = x^2", in_y)):
    print(f"  {name:20s}: P(x < 0.5) = {np.mean(s < 0.5):.4f}, <x> = {s.mean():.4f}")

print("\n=== 3. logistic map r = 4: physical measure forgets a.c. preparations; singular ones are not forgotten ===")
f = lambda z: 4 * z * (1 - z)
preps = {"a.c. uniform": rng.random(100000), "a.c. narrow Gaussian at 0.3": np.clip(rng.normal(0.3, 0.02, 100000), 1e-9, 1 - 1e-9),
         "Dirac on period-2 orbit": np.full(100000, (5 - np.sqrt(5)) / 8)}
for name, z in preps.items():
    zz = z.copy()
    for _ in range(30):
        zz = f(zz)
    print(f"  {name:30s}: P(z < 0.5) after 30 steps = {np.mean(zz < 0.5):.4f}  (arcsine value 0.5000)")
print("  (the period-2 value is exact in exact arithmetic; float roundoff eventually ejects it, a numerical artefact)")

print("\n=== 4. shared reference frames from long-range order (2D Ising frame field, L = 64) ===")


def ising_corr(L, T, sweeps=2500, meas_from=1000, r=32):
    s = np.ones((L, L), dtype=int)
    ii, jj = np.indices((L, L))
    masks = [((ii + jj) % 2 == c) for c in (0, 1)]
    acc, n = 0.0, 0
    for sw in range(sweeps):
        for m in masks:
            nb = np.roll(s, 1, 0) + np.roll(s, -1, 0) + np.roll(s, 1, 1) + np.roll(s, -1, 1)
            dE = 2 * s * nb
            flip = (rng.random((L, L)) < np.exp(-dE / T)) & m
            s = np.where(flip, -s, s)
        if sw >= meas_from and sw % 5 == 0:
            acc += np.mean(s * np.roll(s, r, 1)); n += 1
    return acc / n


for T in (1.8, 2.0, 2.2, 2.6, 3.2):
    c = ising_corr(64, T)
    print(f"  T = {T}: <f_A f_B> at distance 32 = {c:+.3f}  ->  local-tomography signal 2<f_A f_B> = {2 * c:+.3f}")
print("  -> in the ordered phase the parties share a frame and the J-correlation becomes locally accessible;")
print("     in the disordered phase it does not. EXISTENCE of the shared frame = being in the ordered basin")
print("     (T < T_c, d >= 2); its SIGN is irrelevant (complex conjugation).")
