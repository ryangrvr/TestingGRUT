import numpy as np
SEED=20260101
rng=np.random.default_rng(SEED)
# 1) Check the slice sampler correctness on the target density:
#    Direct rejection test: sample 200k x from density ∝ exp(-x^2/2-x^4/4), compare m2, m4 to
#    high-accuracy values: m2+m4=1 (exact). Also compute Var(x^2)=1-m2-m2^2 > 0.
#    If m2+m4 = 1 holds, the sampler is correct.
def sample_x(n, rng):
    x = rng.normal(0.0, 1.0, n)
    for _ in range(15):
        y = -np.log(rng.uniform(size=n)) - 0.5*x**2
        L = np.full(n, -1.0); R = np.full(n, 1.0)
        # expand brackets vectorized
        def V(z): return 0.5*z*z + 0.25*z**4
        cont = V(L) > y
        while cont.any():
            L[cont] -= 0.5
            cont = (V(L) > y) & (L > -50)
        cont = V(R) > y
        while cont.any():
            R[cont] += 0.5
            cont = (V(R) > y) & (R < 50)
        # sample within brackets
        xc = rng.uniform(L, R)
        acc = V(xc) <= y
        for _ in range(30):
            if acc.all(): break
            redo = ~acc
            xc[redo] = rng.uniform(L[redo], R[redo])
            acc[redo] = V(xc[redo]) <= y[redo]
        x = xc
    return x
x = sample_x(200_000, rng)
m2 = (x**2).mean(); m4 = (x**4).mean()
print(f"m2 = {m2:.6f}, m4 = {m4:.6f}, m2+m4 = {m2+m4:.6f} (should be 1.000000)")
print(f"Var(x^2) = {m4 - m2**2:.6f} (should be > 0)")
# True m2 for this density ~ 0.55 (approximately); check
