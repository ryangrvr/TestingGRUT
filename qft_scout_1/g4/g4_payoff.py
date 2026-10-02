"""QFT-SCOUT-1 G4 -- payoff firewall numerics (illustrations; comparator = ordinary relativistic QFT / AQFT).
Q2: vary mass, field content (central charge), lattice regulator, buffer scale for the G1 entanglement candidates.
Q3: vary mass (and state) for the near-cut 2*pi modular slope."""
import numpy as np
import mpmath as mp
def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)
L = 6.0
def kmat(N, mu, stencil):
    if stencil == "2nd":
        return (mu ** 2 + 2) * np.eye(N) - np.eye(N, k=1) - np.eye(N, k=-1)
    # 4th-order improved Laplacian: -(1/12) f_{j-2} + (4/3) f_{j-1} - (5/2) f_j + (4/3) f_{j+1} - (1/12) f_{j+2}
    return (mu ** 2 + 2.5) * np.eye(N) - 4 / 3 * (np.eye(N, k=1) + np.eye(N, k=-1)) + 1 / 12 * (np.eye(N, k=2) + np.eye(N, k=-2))
def vac(K):
    lam, V = np.linalg.eigh(K); return 0.5 * (V * lam ** -0.5) @ V.T, 0.5 * (V * lam ** 0.5) @ V.T
def S(X, P, idx):
    nu = np.sqrt(np.clip(np.linalg.eigvals(X[np.ix_(idx, idx)] @ P[np.ix_(idx, idx)]).real, 0.25, None)); nu = nu[nu > 0.5 + 1e-12]
    return float(np.sum((nu + 0.5) * np.log(nu + 0.5) - (nu - 0.5) * np.log(nu - 0.5)))

hdr("Q2  entanglement candidates under variation of supplied theory data (mass, field content, regulator, buffer)")
Ns = [240, 480, 960, 1920]
for stencil in ("2nd", "4th"):
    for m in (0.5, 1.0, 2.0):
        Ss, dE, MI = [], [], []
        for N in Ns:
            a = 2 * L / N; X, P = vac(kmat(N, m * a, stencil)); h = N // 2
            Ss.append(S(X, P, list(range(h))))
            K = kmat(N, m * a, stencil); A, B = list(range(h)), list(range(h, N))
            dE.append(-np.sum(K[np.ix_(A, B)] * X[np.ix_(A, B)]) / a)          # energy cost of the product of vacuum marginals
            k = int(round(1.0 / (2 * a))); MI.append(S(X, P, list(range(h - k))) + S(X, P, list(range(h + k, N))) - S(X, P, list(range(h - k)) + list(range(h + k, N))))
        inc = np.diff(Ss) / np.log(2)
        fin = Ss[-1] - np.log(1 / (2 * L / Ns[-1])) / 6
        print(f"  stencil {stencil}, m={m}: dS per halving {np.round(inc, 4)} (-> c/6=0.1667); finite part S-(1/6)ln(1/a) = {fin:+.4f};"
              f" dE_prod(a->) {np.round(dE, 1)}; I(A:B) d=1: {MI[-1]:.5f}")
print("  field content: n decoupled fields are additive -> sharp-cut coefficient n/6 (c = n); e.g. 2 fields: increment = 2 x 0.1667 = 0.3333")
print("  -> universal: the log coefficient c/6 (fixed by the supplied CFT data c, i.e. field content); non-universal: finite parts,")
print("     buffered MI values (mass- and buffer-dependent), and the product-state energy coefficient (regulator-dependent).")

hdr("Q3  near-cut 2*pi modular slope under mass variation (80-digit precision; interval of 30 sites)")
mp.mp.dps = 80
def near_cut_ratio(Nn, Rn, mu, beta=None):
    K = mp.matrix(Nn, Nn)
    for i in range(Nn):
        K[i, i] = mu ** 2 + 2
        if i + 1 < Nn: K[i, i + 1] = K[i + 1, i] = -1
    lam, V = mp.eigsy(K)
    f = lambda l: (1 if beta is None else mp.coth(beta * mp.sqrt(l) / 2))
    R = list(range(Nn // 2, Nn // 2 + Rn)); X = mp.matrix(Rn, Rn); Pm = mp.matrix(Rn, Rn)
    for a_, i in enumerate(R):
        for b_, j in enumerate(R):
            X[a_, b_] = sum(V[i, k] * V[j, k] * lam[k] ** -0.5 * f(lam[k]) for k in range(Nn)) / 2
            Pm[a_, b_] = sum(V[i, k] * V[j, k] * lam[k] ** 0.5 * f(lam[k]) for k in range(Nn)) / 2
    w, O = mp.eigsy(X); Xh = O * mp.diag([mp.sqrt(x) for x in w]) * O.T
    nu2, Q = mp.eigsy(Xh * Pm * Xh); nu = [mp.sqrt(x) for x in nu2]
    eps = [mp.log((2 * v + 1) / (2 * v - 1)) for v in nu]
    Hp = Xh * Q * mp.diag([e / v for e, v in zip(eps, nu)]) * Q.T * Xh
    return [float(Hp[j, j] / (2 * mp.pi * (j + 0.5))) for j in (0, 1, 2)]
for mu in ("0.02", "0.2", "0.5"):
    print(f"  vacuum, mu = m a = {mu}: H_p[j,j]/(2 pi (j+1/2)) at j = 0, 1, 2: {np.round(near_cut_ratio(120, 30, mp.mpf(mu)), 3)}")
print(f"  thermal beta=4, mu=0.02 (state change): {np.round(near_cut_ratio(120, 30, mp.mpf('0.02'), mp.mpf(4)), 3)}")
print("  -> the 2*pi slope at the entangling point is mass-independent (vacuum, any mass) but STATE-dependent; it is the standard")
print("     Bisognano-Wichmann / Unruh normalization of relativistic QFT (comparator), not a GRUT-distinctive number.")
