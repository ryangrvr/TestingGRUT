# P-01: do the E-1 floor edges (gap -> P_memory, passivity -> P_positivity) port off the C1-a chain?
# Retained site r coupled by one unit spring to bath vertex b0; bath block K_b = pin*I + L_bath + e_b0 e_b0^T
# (the coupling spring's diagonal, as build_K does); kernel k(tau) = [exp(-K_b tau)]_{b0,b0} (k(0) = 1).
# Requires numpy. Deterministic (fixed seed). See P01_RESULT.md.
import numpy as np
rng = np.random.default_rng(1)
TAU = np.arange(1.0, 40.01, 0.5)

def kernel(Kb, b0=0, taus=TAU):
    lam, U = np.linalg.eigh(Kb)
    w = U[b0, :]**2
    return lam, w, np.array([np.sum(w*np.exp(-lam*t)) for t in taus])
def grade(k):   # L0-1a frozen comparator: max-abs residual of ln k vs tau (exp) and ln k vs ln tau (alg)
    y = np.log(np.maximum(k, 1e-300))
    r = []
    for x in (TAU, np.log(TAU)):
        A = np.vstack([x, np.ones_like(x)]).T; c, *_ = np.linalg.lstsq(A, y, rcond=None)
        r.append(np.max(np.abs(A @ c - y)))
    return ("EXPONENTIAL" if r[0] < r[1] else "ALGEBRAIC"), r
def chain_bath(n, pin):     # C1-a bath block: path on n sites, coupling spring at site 0
    K = pin*np.eye(n)
    for i in range(n-1):
        K[i, i] += 1; K[i+1, i+1] += 1; K[i, i+1] -= 1; K[i+1, i] -= 1
    K[0, 0] += 1
    return K
def random_regular(n, d):   # configuration model, reject loops / multi-edges
    while True:
        stubs = np.repeat(np.arange(n), d); rng.shuffle(stubs); E = stubs.reshape(-1, 2)
        if np.any(E[:, 0] == E[:, 1]): continue
        key = set()
        ok = True
        for a, b in E:
            t = (min(a, b), max(a, b))
            if t in key: ok = False; break
            key.add(t)
        if ok: return E
def graph_bath(n, E, pin):
    K = pin*np.eye(n)
    for a, b in E:
        K[a, a] += 1; K[b, b] += 1; K[a, b] -= 1; K[b, a] -= 1
    K[0, 0] += 1
    return K

print("=== C0 replication anchor: chain, pin=0.3, bath n=23 (C1-a) ===")
lam, w, k = kernel(chain_bath(23, 0.3))
print(f"  lam_min={lam.min():.4f} (>= pin)  k(40)={k[-1]:.3e} (<= e^-12 = 6.1e-6)  grade={grade(k)[0]}")

print("\n=== C1 GAP edge: delete the pin (pin = 0) on different bath graphs ===")
print("  chain (amenable, 1D):")
for n in (23, 95, 383):
    lam, w, k = kernel(chain_bath(n, 0.0))
    print(f"    n={n:4d}: lam_min={lam.min():.2e}  local mass below 0.15 = {w[lam < 0.15].sum():.3f}  k(40)={k[-1]:.3e}  grade={grade(k)[0]}")
for d in (3, 4):
    bot = d - 2*np.sqrt(d-1)
    print(f"  random {d}-regular bath (non-amenable limit; Kesten-McKay bottom d-2sqrt(d-1) = {bot:.4f}):")
    for n in (200, 800, 3200):
        E = random_regular(n, d)
        lam, w, k = kernel(graph_bath(n, E, 0.0))
        g, r = grade(k)
        soft = lam < 0.9*bot
        kreg = np.array([np.sum(w[~soft]*np.exp(-lam[~soft]*t)) for t in TAU])
        A = np.vstack([TAU, np.ones_like(TAU)]).T; rate = -np.linalg.lstsq(A, np.log(kreg), rcond=None)[0][0]
        print(f"    n={n:5d}: lam_min={lam.min():.2e}  soft mass (lam<{0.9*bot:.3f}) = {w[soft].sum():.2e}  n*soft = {n*w[soft].sum():.3f}  "
              f"raw window grade={g} | soft-removed kernel: grade={grade(kreg)[0]}, fitted rate={rate:.3f}, min visible lam={lam[~soft].min():.3f}")
print("  -> chain: soft mass is n-independent (a real algebraic tail). Expander: soft mass ~ 1/n -> 0, and the")
print("     remaining kernel is exponential with rate >= the Kesten-McKay bottom: pin-free memory is exponential in the limit.")
print("  -> the carrier is the bottom of the retained")
print("     site's local spectral measure, which the pin supplies only when the graph itself does not")

print("\n=== C2 PASSIVITY edge, symmetric families ===")
pin, wneg = 0.3, -1.5
# bath = {h, a, b}: h-a, h-b unit springs, a-b spring w < 0; retained site couples to h (b0 = 0)
K = pin*np.eye(3); K[0, 0] += 1
for (i, j, s) in ((0, 1, 1.0), (0, 2, 1.0), (1, 2, wneg)):
    K[i, i] += s; K[j, j] += s; K[i, j] -= s; K[j, i] -= s
lam, w, k = kernel(K)
print(f"  two-branch bath with a-b spring w={wneg}: eigenvalues {np.round(lam, 4)} (NOT passive)")
print(f"    retained weights {np.round(w, 6)}  -> unstable mode weight = {w[lam < 0].sum():.1e} (antisymmetric, invisible)")
mono = np.all(np.diff(k) <= 1e-12); print(f"    kernel strictly decreasing on window: {mono}; k(40) = {k[-1]:.3e} > 0; visible Bernstein weights all at lam >= 0: {np.all(lam[w > 1e-14] >= 0)}")
Kp = K.copy(); Kp[0, 1] -= 0.05; Kp[1, 0] -= 0.05; Kp[0, 0] += 0.05; Kp[1, 1] += 0.05   # break the symmetry slightly
lam2, w2, k2 = kernel(Kp)
print(f"  same graph, symmetry broken by 0.05 on h-a: unstable weight = {w2[lam2 < 0].sum():.2e}, k(40)/k(1) = {k2[-1]/k2[0]:.2e} (growth: positivity lost)")
print("  chain (Jacobi): every eigenvector has nonzero endpoint component -> min retained weight over modes:",
      f"{kernel(chain_bath(23, 0.3))[1].min():.2e} > 0, so on the chain an unstable mode is always visible")

print("\n=== C3 PASSIVITY edge, asymmetric (directed) Laplacian: accretive but not CM ===")
n, pin = 12, 0.05
P = np.roll(np.eye(n), 1, axis=1)
Kd = pin*np.eye(n) + (np.eye(n) - P); Kd[0, 0] += 1
sym_min = np.linalg.eigvalsh(Kd + Kd.T).min()
lamd = np.linalg.eigvals(Kd)
taus = np.linspace(0, 40, 801)
from scipy.linalg import expm
kd = np.array([expm(-Kd*t)[0, 0] for t in taus])
print(f"  min eig(K+K^T) = {sym_min:.3f} >= 0 (accretive/passive); eigenvalues complex: {np.any(np.abs(lamd.imag) > 1e-9)}")
print(f"  kernel min = {kd.min():+.3e} (sign change: {kd.min() < 0}); increasing somewhere: {np.any(np.diff(kd) > 1e-12)} -> not CM (record E-4: affinity -> not CM)")
