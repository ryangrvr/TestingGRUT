import numpy as np
rng = np.random.default_rng(20260101)
# --- Correct slice sampler for rho(x) ~ exp(-V(x)), V(x)=x^2/2+x^4/4 ---
def sample_x(n, rng, n_iter=20):
    x = np.zeros(n)
    def V(z): return 0.5*z*z + 0.25*z**4
    for _ in range(n_iter):
        # slice level
        y = -np.log(rng.uniform(size=n)) - V(x)   # log of uniform * exp(-V)
        L = x.copy(); R = x.copy()
        # expand left
        step = 1.0
        cont = V(L) > y
        while cont.any():
            L[cont] -= step
            cont = V(L) > y
        cont = V(R) > y
        while cont.any():
            R[cont] += step
            cont = V(R) > y
        # shrink
        for _ in range(50):
            xc = rng.uniform(L, R)
            ok = V(xc) <= y
            x[ok] = xc[ok]
            L[~ok & (xc < x)] = xc[~ok & (xc < x)]
            R[~ok & (xc >= x)] = xc[~ok & (xc >= x)]
            if (V(x) <= y).all(): break
    return x

x = sample_x(200000, rng, n_iter=15)
m2 = (x**2).mean(); m4 = (x**4).mean()
print(f"m2 = {m2:.6f}, m4 = {m4:.6f}, m2+m4 = {m2+m4:.6f}")
print(f"Var(x^2) = {m4 - m2**2:.6f}")
