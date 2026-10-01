# P-02b: does interaction supply a density -> IR-velocity mechanism that restores the SF-1 sector split?
# Parent: H = -sum_j (b_j^dag b_{j+1} + h.c.) + (U/2) sum_j n_j(n_j - 1), periodic ring (U = NEW ASSUMPTION).
# U = 0 is the unmodified P-02 control (free bosons). Same readout rho_q, same support object as SF-1.
# Exact sector ground state (unique for every U by Perron-Frobenius: hopping off-diagonals <= 0, connected).
# Lowest supported excitation at q = lowest Ritz value of H on the Krylov space of rho_q|0> (full reorth).
# Requires numpy, scipy.
import itertools, math, numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

def basis(L, N, nmax):
    out = []
    def rec(pref, left, sites):
        if sites == 0:
            if left == 0: out.append(tuple(pref))
            return
        for x in range(min(left, nmax), -1, -1): rec(pref+[x], left-x, sites-1)
    rec([], N, L); return out

def build(L, N, U, hardcore=False):
    nmax = 1 if hardcore else N
    B = basis(L, N, nmax); idx = {c: i for i, c in enumerate(B)}
    rows, cols, vals = [], [], []
    for c, i in idx.items():
        if not hardcore and U:
            rows.append(i); cols.append(i); vals.append(U/2*sum(x*(x-1) for x in c))
        for a in range(L):
            for (s, t) in ((a, (a+1) % L), ((a+1) % L, a)):
                if c[t] == 0 or c[s] + 1 > nmax: continue
                d = list(c); amp = math.sqrt(d[t]); d[t] -= 1; amp *= math.sqrt(d[s]+1); d[s] += 1
                rows.append(idx[tuple(d)]); cols.append(i); vals.append(-amp)
    H = csr_matrix((vals, (rows, cols)), shape=(len(B), len(B)))
    occ = np.array(B, dtype=float)
    return H, occ

def ground(H):
    if H.shape[0] < 400:
        E, V = np.linalg.eigh(H.toarray()); return E[0], V[:, 0], E[1]-E[0]
    E, V = eigsh(H, k=2, which='SA', tol=1e-12); o = np.argsort(E)
    return E[o[0]], V[:, o[0]], E[o[1]]-E[o[0]]

def lowest_supported(H, v, E0, m=200):
    nv = np.linalg.norm(v)
    if nv < 1e-12: return None
    Q = [v/nv]; al, be = [], []
    w = H @ Q[0]
    for j in range(m):
        a = np.vdot(Q[j], w).real; al.append(a)
        w = w - a*Q[j] - (be[-1]*Q[j-1] if j else 0)
        for qv in Q: w -= np.vdot(qv, w)*qv          # full reorthogonalization
        b = np.linalg.norm(w)
        if b < 1e-10: break
        be.append(b); Q.append(w/b); w = H @ Q[-1]
    T = np.diag(al) + np.diag(be[:len(al)-1], 1) + np.diag(be[:len(al)-1], -1)
    th, S = np.linalg.eigh(T)
    for t, s0 in zip(th, S[0]):
        if s0**2 > 1e-10: return t - E0
    return None

def edges(L, N, U, hardcore=False, ms=None):
    H, occ = build(L, N, U, hardcore)
    E0, g, gap = ground(H)
    ms = ms or range(1, L)
    res = {}
    for m in ms:
        q = 2*np.pi*m/L
        rho = (occ*np.exp(-1j*q*np.arange(L))[None, :]).sum(1)
        res[m] = lowest_supported(H.astype(complex), rho*g, E0)
    return H.shape[0], gap, res

print("=== D-type sector (nu = 1/2, N = L/2 odd): omega_1(L) = lowest support at q = 2pi/L, and the q = pi edge ===")
DL = ((6, 3), (10, 5), (14, 7))
table = {}
for lab, U, hc in (("U=0 (control)", 0.0, False), ("U=1", 1.0, False), ("U=4", 4.0, False), ("U=inf (hard-core)", 0.0, True)):
    w1, wpi = [], []
    for L, N in DL:
        dim, gap, r = edges(L, N, U, hc, ms=[1, L//2])
        w1.append(r[1]); wpi.append(r[L//2])
    zs = [-math.log(w1[i+1]/w1[i])/math.log(DL[i+1][0]/DL[i][0]) for i in range(2)]
    table[lab] = (w1, wpi, zs)
    print(f"  {lab:18}: omega_1 = {[round(x, 5) for x in w1]}  local z_P' = {[round(z, 3) for z in zs]}  "
          f"omega(q=pi) = {[round(x, 5) for x in wpi]}")
print("  U=0: z=2, no soft point at pi (free-boson collapse).  U>0: z -> 1 and omega(q=pi) -> 0 (Lieb II / 2k_F soft point)")

print("\n=== E-type sector (fixed N = 3, L grows): interaction cannot generate an IR velocity at vanishing density ===")
EL = ((6, 3), (10, 3), (14, 3))
for lab, U, hc in (("U=0 (control)", 0.0, False), ("U=4", 4.0, False), ("U=inf (hard-core)", 0.0, True)):
    w1 = [edges(L, N, U, hc, ms=[1])[2][1] for L, N in EL]
    zs = [-math.log(w1[i+1]/w1[i])/math.log(EL[i+1][0]/EL[i][0]) for i in range(2)]
    print(f"  {lab:18}: omega_1 = {[round(x, 6) for x in w1]}  local z_P' = {[round(z, 3) for z in zs]}  (-> 2)")

print("\n=== Full support lower edge, D sector L=10 N=5 (U=4): soft points ===")
dim, gap, r = edges(10, 5, 4.0)
print("  m:", list(r.keys())); print("  omega^-(q):", [round(r[m], 4) for m in r])
print("  ground state unique (gap):", round(gap, 4))

print("\n=== Bogoliubov (mean-field, weak U; exact in d >= 2 condensates, NOT the 1D exact answer) ===")
for n in (0.5, 0.0):
    for qq in (0.2, 0.1, 0.05):
        eq = 4*math.sin(qq/2)**2; w = math.sqrt(eq*(eq + 2*1.0*n))
        print(f"  n={n}: q={qq}: omega_B = {w:.6f}", end="")
    print()
print("  finite n: omega_B ~ sqrt(2Un) q (z=1) with soft set {0} only -> class (1,1); n -> 0: omega_B -> eps_q (z=2), class (2,1)")
