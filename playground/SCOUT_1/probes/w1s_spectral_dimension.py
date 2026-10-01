"""SCOUT-1 W1-S: does spectral dimension fix the retained-site edge exponent gamma?
mu_r(d lambda) ~ (lambda - lambda_0)^gamma d lambda near the edge  <=>  k_0(tau) ~ tau^{-(gamma+1)}.

 1. translation-invariant hypercubic lattices d = 1, 2, 3: bulk / Neumann face / Dirichlet face / Dirichlet
    corner sites, isotropic and anisotropic couplings (exact return probabilities via Bessel functions).
 2. a point defect V at a bulk site (threshold Green's function): relevant (d=1), marginal (d=2), irrelevant
    unless resonant (d=3).
 3. local, passive, bounded-coupling chains realizing a CONTINUUM of edge exponents (radial graphs of real
    dimension D): gamma = D/2 - 1 for any real D.
"""
import numpy as np
from scipy.special import ive
from scipy.integrate import quad
import scipy.sparse as sps

np.set_printoptions(precision=4, suppress=True)
taus = np.array([50.0, 200.0, 800.0, 3200.0])


def slope(f):
    v = np.array([f(t) for t in taus])
    return np.diff(np.log(v)) / np.diff(np.log(taus))


# 1D factors of the continuous-time return probability with coupling J (generator J*(2 - shift - shift^-1))
bulk = lambda t, J=1.0: ive(0, 2 * J * t)                        # e^{-2Jt} I_0(2Jt)
dirichlet = lambda t, J=1.0: ive(0, 2 * J * t) - ive(2, 2 * J * t)  # image method, site next to a removed site
neumann = lambda t, J=1.0: ive(0, 2 * J * t) + ive(1, 2 * J * t)    # reflecting boundary site

print("=== 1. hypercubic lattices: local slope of k_0(tau) at tau = 50..3200; predicted -(gamma+1) ===")
cases = [
    ("d=1 bulk", lambda t: bulk(t), -0.5),
    ("d=1 Dirichlet end", lambda t: dirichlet(t), 0.5),
    ("d=1 Neumann end", lambda t: neumann(t), -0.5),
    ("d=2 bulk", lambda t: bulk(t) ** 2, 0.0),
    ("d=2 bulk anisotropic Jx=1, Jy=0.2", lambda t: bulk(t) * bulk(t, 0.2), 0.0),
    ("d=2 Neumann edge", lambda t: neumann(t) * bulk(t), 0.0),
    ("d=2 Dirichlet edge", lambda t: dirichlet(t) * bulk(t), 1.0),
    ("d=2 Dirichlet corner", lambda t: dirichlet(t) ** 2, 2.0),
    ("d=3 bulk", lambda t: bulk(t) ** 3, 0.5),
    ("d=3 bulk anisotropic 1, 0.5, 0.1", lambda t: bulk(t) * bulk(t, 0.5) * bulk(t, 0.1), 0.5),
    ("d=3 Dirichlet face", lambda t: dirichlet(t) * bulk(t) ** 2, 1.5),
    ("d=3 Dirichlet corner", lambda t: dirichlet(t) ** 3, 3.5),
]
for name, f, g in cases:
    print(f"  {name:36s} slopes {slope(f)}   -> gamma ~ {-slope(f)[-1] - 1:+.3f}   predicted {g:+.1f} (= d/2 - 1 + #Dirichlet)")


print("\n=== 2. point defect V at a bulk site: rho_V(E) = rho_0(E) / |1 - V G_0(E)|^2 near the edge E -> 0 ===")
# G_0(-s) = -int_0^inf e^{-s t} P(t) dt (resolvent just below the band bottom, E measured from the edge)
for d in (1, 2, 3):
    P = lambda t, d=d: bulk(t) ** d
    vals = []
    for s in (1e-2, 1e-3, 1e-4, 1e-5):
        G, _ = quad(lambda t: np.exp(-s * t) * P(t), 0, np.inf, limit=800)
        vals.append(G)
    print(f"  d={d}: |G_0| at distance s = 1e-2..1e-5 below the edge: {np.round(vals, 4)}")
print("  d=1: |G_0| ~ s^{-1/2} diverges -> rho_V ~ E^{+1/2}: any V != 0 is RELEVANT (gamma -1/2 -> +1/2)")
print("  d=2: |G_0| ~ log(1/s) diverges slowly -> rho_V ~ 1/log^2 E: MARGINAL (gamma = 0 with log corrections)")
print("  d=3: |G_0| -> G_c finite (Watson) -> gamma = 1/2 unchanged unless V = -1/G_c (threshold resonance: gamma -> -1/2)")


print("\n=== 3. continuum of exponents: radial graphs of real dimension D (local, passive, bounded couplings) ===")


def radial_chain_sparse(n, D):
    c = (np.arange(n - 1) + 1.0) ** (D - 1)
    m = (np.arange(n) + 1.0) ** (D - 1)
    diag = np.zeros(n); diag[:-1] += c; diag[1:] += c
    off = -c / np.sqrt(m[:-1] * m[1:])
    return sps.diags([off, diag / m, off], [-1, 0, 1])


def radial_chain(n, D):
    """birth-death chain with conductances c_k=(k+1)^(D-1) between shells k,k+1 and masses m_k=(k+1)^(D-1);
    symmetrized generator K = M^{-1/2} L_c M^{-1/2}"""
    c = (np.arange(n - 1) + 1.0) ** (D - 1)
    m = (np.arange(n) + 1.0) ** (D - 1)
    diag = np.zeros(n); diag[:-1] += c; diag[1:] += c
    K = np.diag(diag / m)
    off = -c / np.sqrt(m[:-1] * m[1:])
    K += np.diag(off, 1) + np.diag(off, -1)
    return K


import scipy.sparse as sps
from scipy.sparse.linalg import expm_multiply
n = 20000
tt = np.array([100.0, 400.0, 1600.0, 6400.0])
for D in (1.0, 1.5, 2.0, 2.5, 3.0, 3.7):
    K = sps.csr_matrix(radial_chain_sparse(n, D))
    e0 = np.zeros(n); e0[0] = 1.0
    k = []
    for t in tt:
        k.append(expm_multiply(-t * K, e0)[0])
    k = np.array(k)
    sl = np.diff(np.log(k)) / np.diff(np.log(tt))
    offd = -K.diagonal(1)
    print(f"  D={D:3.1f}: slopes (tau 100..6400) {np.round(sl, 3)} -> gamma ~ {-sl[-1] - 1:+.3f} (predicted D/2-1 = {D / 2 - 1:+.2f});"
          f"  couplings in [{offd.min():.3f}, {offd.max():.3f}]")
print("  every member: nearest-neighbour (local), symmetric PSD (passive), bounded couplings; add a pin for a gap.")
print("  => without translation invariance + integer dimension, gamma ranges over a continuum.")
