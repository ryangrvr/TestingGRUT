"""SCOUT-1 W2-ETH: does a closed chain thermalize a subsystem to a single-temperature KMS state without a
supplied passive/Gibbs preparation?

Mixed-field Ising chain, PBC: H = sum Z_i Z_{i+1} + g sum X_i + h sum Z_i
  non-integrable: (g, h) = (0.9045, 0.8090)     integrable (free-fermion): (g, h) = (0.9045, 0)
Initial product states |theta> = prod (cos(theta/2)|0> + sin(theta/2)|1>) (rotated about y).
Infinite-time average = diagonal ensemble rho_DE = sum |c_n|^2 |n><n| (exact).
Compare the 2-site reduced density matrix of rho_DE with that of the Gibbs state at beta fixed by <H>.
"""
import numpy as np
from scipy.optimize import brentq

X = np.array([[0, 1], [1, 0]], float); Z = np.diag([1.0, -1.0]); I2 = np.eye(2)


def op(o, i, L):
    out = np.array([[1.0]])
    for k in range(L):
        out = np.kron(out, o if k == i else I2)
    return out


def ham(L, g, h):
    H = np.zeros((2 ** L, 2 ** L))
    Zs = [op(Z, i, L) for i in range(L)]
    for i in range(L):
        H += Zs[i] @ Zs[(i + 1) % L] + g * op(X, i, L) + h * Zs[i]
    return H


def product_state(L, theta):
    s = np.array([np.cos(theta / 2), np.sin(theta / 2)])
    v = np.array([1.0])
    for _ in range(L):
        v = np.kron(v, s)
    return v


def reduced_2site(rho, L, i=0):
    # keep sites i, i+1 (adjacent; i = 0 by translation invariance)
    R = rho.reshape([2] * (2 * L))
    keep = [i, i + 1]
    trace_out = [k for k in range(L) if k not in keep]
    for k in sorted(trace_out, reverse=True):
        R = np.trace(R, axis1=k, axis2=k + R.ndim // 2)
    return R.reshape(4, 4)


def tdist(a, b):
    return 0.5 * np.sum(np.abs(np.linalg.eigvalsh(a - b)))


def analyse(L, g, h, thetas):
    H = ham(L, g, h)
    E, V = np.linalg.eigh(H)
    res = []
    for th in thetas:
        psi = product_state(L, th)
        c = V.T @ psi
        e0 = float(c ** 2 @ E)
        # Gibbs beta from energy
        f = lambda b: (np.exp(-b * (E - E.min())) @ E) / np.exp(-b * (E - E.min())).sum() - e0
        beta = brentq(f, -5, 20)
        w = np.exp(-beta * (E - E.min())); w /= w.sum()
        rho_G = (V * w) @ V.T
        rho_DE = (V * c ** 2) @ V.T   # diagonal ensemble (degeneracies: off-diagonal terms within degenerate blocks dropped)
        rG, rD = reduced_2site(rho_G, L), reduced_2site(rho_DE, L)
        rho0 = np.outer(psi, psi)
        r0 = reduced_2site(rho0, L)
        res.append((th, e0 / L, beta, tdist(rD, rG), tdist(r0, rG)))
    return res


if __name__ == "__main__":
    thetas = (0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4)
    for label, g, h in (("NON-INTEGRABLE (g,h)=(0.9045,0.809)", 0.9045, 0.8090),
                        ("INTEGRABLE (g,h)=(0.9045,0)", 0.9045, 0.0)):
        print(f"\n=== {label} ===")
        print("   L  theta   E/L     beta    D(rho_DE,rho_Gibbs)_2site   D(initial,Gibbs)_2site")
        for L in (8, 10, 12):
            for th, el, b, d, d0 in analyse(L, g, h, thetas):
                print(f"  {L:2d}  {th:5.3f}  {el:+.4f}  {b:+.4f}   {d:.5f}                    {d0:.5f}")
    print("\nNote: the diagonal ensemble keeps only |c_n|^2; degenerate blocks (momentum/parity) would add"
          " off-diagonal terms for a true time average; product states are translation- and parity-invariant, so"
          " they couple only to the symmetric sector where degeneracies are rare.")
