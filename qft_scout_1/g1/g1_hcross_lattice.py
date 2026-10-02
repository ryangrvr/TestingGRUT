"""QFT-SCOUT-1 G1 -- lattice ILLUSTRATION (not proof) of H_cross = 0 admissibility: sharp vs split-buffer localization.
Free massive scalar in 1+1D, mass m = 1, Dirichlet box [-L, L], lattice spacing a -> 0 (quantum lift QP-5 supplied, hbar = 1).
H = (1/a) [ sum p_j^2 / 2 + q^T K q / 2 ],  K = (m a)^2 I + (2 I - T);  vacuum covariances X = K^{-1/2}/2, P = K^{1/2}/2.
Type-I control: at every finite a the algebras are type I (finite tensor factors) and product states are always admissible
and normal; the question is what happens to their cost as a -> 0."""
import numpy as np
m, L = 1.0, 6.0
def vacuum(N, a):
    K = (m * a) ** 2 * np.eye(N) + 2 * np.eye(N) - np.eye(N, k=1) - np.eye(N, k=-1)
    lam, V = np.linalg.eigh(K)
    X = 0.5 * (V * lam ** -0.5) @ V.T; P = 0.5 * (V * lam ** 0.5) @ V.T
    return X, P
def entropy(X, P, idx):
    nu2 = np.linalg.eigvals(X[np.ix_(idx, idx)] @ P[np.ix_(idx, idx)]).real
    nu = np.sqrt(np.clip(nu2, 0.25, None)); nu = nu[nu > 0.5 + 1e-12]
    return float(np.sum((nu + 0.5) * np.log(nu + 0.5) - (nu - 0.5) * np.log(nu - 0.5)))
print(f"m = {m}, box [-{L}, {L}], Dirichlet; sites N = 2L/a")
print(" a        N     S_sharp(A)   (1/6)ln(1/ma)   I_sharp=2S   dE_prod (phys)   I(A:B) d=0.5   I(A:B) d=1.0   I(A:B) d=2.0")
rows = []
for N in (240, 480, 960, 1920, 3840):
    a = 2 * L / N; X, P = vacuum(N, a); h = N // 2
    A = list(range(h)); B = list(range(h, N))
    S = entropy(X, P, A)
    dE = X[h - 1, h] / a                                    # energy of rho_A (x) rho_B minus vacuum energy (physical units)
    mi = []
    for d in (0.5, 1.0, 2.0):
        k = int(round(d / (2 * a)))
        Ab = list(range(0, h - k)); Bb = list(range(h + k, N))
        mi.append(entropy(X, P, Ab) + entropy(X, P, Bb) - entropy(X, P, Ab + Bb))
    rows.append((a, S, dE, mi))
    print(f" {a:.5f} {N:5d}   {S:9.5f}     {np.log(1/(m*a))/6:9.5f}     {2*S:9.5f}   {dE:13.3f}   {mi[0]:12.6f}   {mi[1]:12.6f}   {mi[2]:12.6f}")
sl = [(rows[i + 1][1] - rows[i][1]) / np.log(2) for i in range(len(rows) - 1)]
print("\nsharp entropy increment per halving of a (expected -> c/6 = 0.1667 for c = 1, one boundary point):", np.round(sl, 4))
print("energy cost of the sharp product state, ratio per halving of a:", np.round([rows[i + 1][2] / rows[i][2] for i in range(len(rows) - 1)], 3))
for j, d in enumerate((0.5, 1.0, 2.0)):
    v = [r[3][j] for r in rows]
    print(f"buffer d = {d}: I(A:B) changes per halving of a: {np.round(np.diff(v), 7)}  (converging -> finite continuum value ~ {v[-1]:.5f})")
print("\nreading: sharp cut -> the product of vacuum marginals costs unbounded relative entropy (I = 2S ~ (1/3) ln(1/a)) and unbounded")
print("energy (~ ln(1/a)/a) as a -> 0, consistent with no NORMAL product state across a sharp complementary cut in the continuum;")
print("with a buffer the vacuum correlations stay finite, consistent with the split property (normal product states exist).")
print("The vacuum itself is correlated across every tested split (I > 0 for all d): H_cross != 0 for this state, as Reeh-Schlieder suggests.")
