# ZOOM_OUT_03 check: edge exponent of the retained local spectral measure -> kernel tail exponent (Karamata/Tauberian).
# Chain bath, pin = 0 (C1-a family, amenable): mu_r(lam) ~ lam^{1/2} near 0 => k(tau) ~ tau^{-3/2} (= S5-1's native t^{-3/2}).
import numpy as np
for n in (383, 1535):
    K = np.zeros((n, n))
    for i in range(n-1):
        K[i, i] += 1; K[i+1, i+1] += 1; K[i, i+1] -= 1; K[i+1, i] -= 1
    K[0, 0] += 1
    lam, U = np.linalg.eigh(K); w = U[0]**2
    taus = np.array([20.0, 40.0, 80.0, 160.0])
    k = np.array([np.sum(w*np.exp(-lam*t)) for t in taus])
    sl = np.diff(np.log(k))/np.diff(np.log(taus))
    edge = w[lam < 0.01].sum()/0.01**1.5, w[lam < 0.04].sum()/0.04**1.5
    print(f"n={n}: local slope d ln k/d ln tau on 20-160: {np.round(sl, 3)}  (-> -1.5);  mu_r([0,x])/x^1.5 at x=0.01,0.04: {edge[0]:.3f}, {edge[1]:.3f} (const => cumulative mass ~ x^3/2, density ~ lam^1/2)")
