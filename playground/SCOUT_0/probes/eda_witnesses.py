# EDGE-DATA AUDIT witnesses (requires numpy). Same earned predicates, different supplied primitive,
# different edge datum, different effective law. See EDGE_DATA_AUDIT_01.md.
import numpy as np

def kernel_stats(K, b0, taus):
    lam, U = np.linalg.eigh(K); w = U[b0]**2
    k = np.array([np.sum(w*np.exp(-lam*t)) for t in taus])
    return lam, w, k

def chain(n, pin):
    K = pin*np.eye(n)
    for i in range(n-1):
        K[i, i] += 1; K[i+1, i+1] += 1; K[i, i+1] -= 1; K[i+1, i] -= 1
    K[0, 0] += 1; return K, 0

def half_plane(Lx, Ly, pin):        # square lattice, retained site couples to boundary-middle vertex (0, Ly//2)
    n = Lx*Ly; idx = lambda x, y: x*Ly + y
    K = pin*np.eye(n)
    for x in range(Lx):
        for y in range(Ly):
            for dx, dy in ((1, 0), (0, 1)):
                X, Y = x+dx, y+dy
                if X < Lx and Y < Ly:
                    a, b = idx(x, y), idx(X, Y)
                    K[a, a] += 1; K[b, b] += 1; K[a, b] -= 1; K[b, a] -= 1
    b0 = idx(0, Ly//2); K[b0, b0] += 1; return K, b0

taus = np.array([5.0, 10.0, 20.0, 40.0, 80.0])
print("=== W-A: same earned predicates (local NN coupling, passive, gapped, linear), different supplied graph ===")
for pin in (0.3, 0.0):
    for name, (K, b0) in (("chain n=400", chain(400, pin)), ("half-plane 46x46", half_plane(46, 46, pin))):
        lam, w, k = kernel_stats(K, b0, taus)
        k0 = k*np.exp(pin*taus)                    # exact factorization k_pin = e^{-pin tau} k_0
        sl = np.diff(np.log(k0))/np.diff(np.log(taus))
        preds = dict(local=True, passive=bool(lam.min() >= -1e-12), gap=bool(lam.min() >= 0.05),
                     CM=bool(np.all(np.diff(k) < 0)))
        print(f"  pin={pin}: {name:17} predicates {preds}  lam_min={lam.min():.4f}  "
              f"local slope of k_0 (tau 5-80): {np.round(sl, 3)}")
print("  -> chain: edge exponent gamma=1/2, k_0 ~ tau^-3/2.  half-plane boundary site with a point coupling defect:")
print("     2D-type edge (k_0 ~ tau^-1 up to a slowly drifting logarithmic point-defect correction; slope ~ -1.16 here).")
print("     At pin=0.3 both models satisfy every E-1/E-2/E-3 predicate; edge exponent and memory law differ.")

print("\n=== W-B: short-time coefficients read spectral MOMENTS, not edge data (equivalent invariant) ===")
for name, (K, b0) in (("chain", chain(400, 0.3)), ("half-plane", half_plane(46, 46, 0.3))):
    lam, w, _ = kernel_stats(K, b0, taus)
    print(f"  {name:10}: int lam dmu_r = {np.sum(w*lam):.6f}  vs  (K_b)_{{b0,b0}} = {K[b0, b0]:.6f}"
          f"   (C-B Delta c3 = 24 beta T1 a (44 beta a^2 + 5 K11): K11 is the first moment of mu_r)")

print("\n=== W-C: S5-1 Markov coefficient kappa = 1/(2 sqrt(K11)) moves with the supplied pin ===")
for pin in (0.3, 0.1, 0.5):
    K11 = pin + 2.0
    print(f"  pin={pin}: K11={K11:.1f}  kappa={1/(2*np.sqrt(K11)):.6f}" + ("   (record: 1/(2 sqrt 2.3))" if pin == 0.3 else ""))
