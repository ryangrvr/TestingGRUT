"""Quantum illustration (NOT a test of the bridge claim): exact 2D toric code on a 3x3 torus
(18 qubits) with a uniform transverse+longitudinal field h.  Infinite-temperature OTOC
C(t) = 1 - Re <W(t) V W(t) V> for local Paulis, and local vs logical infinite-T autocorrelations,
by stochastic trace with Haar-random states.  2D is not self-correcting; this only shows that a
perturbed stabilizer code scrambles.
"""
import json, sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply

L = 3; n = 2 * L * L; D = 2 ** n
I2 = sp.identity(2, format="csr"); X = sp.csr_matrix([[0, 1], [1, 0]]); Z = sp.csr_matrix([[1, 0], [0, -1]])


def site_op(op, i):
    mats = [I2] * n; mats[i] = op
    out = mats[0]
    for m in mats[1:]:
        out = sp.kron(out, m, format="csr")
    return out


def edge(x, y, d):
    return ((x % L) * L + (y % L)) * 2 + d


def build(h):
    Xs = [site_op(X, i) for i in range(n)]; Zs = [site_op(Z, i) for i in range(n)]
    H = sp.csr_matrix((D, D))
    stars, plaqs = [], []
    for x in range(L):
        for y in range(L):
            es = [edge(x, y, 0), edge(x - 1, y, 0), edge(x, y, 1), edge(x, y - 1, 1)]
            A = Xs[es[0]] @ Xs[es[1]] @ Xs[es[2]] @ Xs[es[3]]
            ep = [edge(x, y, 0), edge(x, y, 1), edge(x + 1, y, 1), edge(x, y + 1, 0)]
            B = Zs[ep[0]] @ Zs[ep[1]] @ Zs[ep[2]] @ Zs[ep[3]]
            stars.append(A); plaqs.append(B); H = H - A - B
    for i in range(n):
        H = H - h * (Xs[i] + Zs[i])
    # logical Z-bar: Z on the horizontal edges of row y=0 (a primal non-contractible loop)
    Zbar = Zs[edge(0, 0, 0)] @ Zs[edge(1, 0, 0)] @ Zs[edge(2, 0, 0)]
    return H.tocsr(), Xs, Zs, Zbar, stars, plaqs


def main(h=0.3, times=None, n_states=2, seed=0):
    if times is None:
        times = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 14.0, 20.0]
    H, Xs, Zs, Zbar, stars, plaqs = build(h)
    rs = np.random.RandomState(seed)
    W = Zs[edge(0, 0, 0)]                # local probe
    V = Xs[edge(1, 2, 1)]                # a distant local perturbation (graph distance >= 2)
    Zloc = Zs[edge(0, 0, 0)]
    out = {"h": h, "times": times, "otoc": [], "acf_local": [], "acf_logical": [], "n_states": n_states}
    t0 = time.time()
    for t in times:
        o, al, ag = [], [], []
        for s in range(n_states):
            psi = rs.randn(D) + 1j * rs.randn(D); psi /= np.linalg.norm(psi)
            # a(t) = W(t) V psi ; b(t) = V W(t) psi ; OTOC = 1 - Re<b|a>
            def Ut(v):  # exp(-iHt) v
                return expm_multiply(-1j * t * H, v) if t > 0 else v
            def Udag(v):
                return expm_multiply(1j * t * H, v) if t > 0 else v
            a = Udag(W @ Ut(V @ psi)); b = V @ Udag(W @ Ut(psi))
            o.append(float(1.0 - np.real(np.vdot(b, a))))
            # infinite-T autocorrelations <O(t) O> = <psi| U^dag O U O |psi>
            al.append(float(np.real(np.vdot(psi, Udag(Zloc @ Ut(Zloc @ psi))))))
            ag.append(float(np.real(np.vdot(psi, Udag(Zbar @ Ut(Zbar @ psi))))))
        out["otoc"].append(float(np.mean(o))); out["acf_local"].append(float(np.mean(al))); out["acf_logical"].append(float(np.mean(ag)))
        print(f"t={t:5.1f} OTOC={out['otoc'][-1]:.4f} acf_local={out['acf_local'][-1]:.4f} acf_logical={out['acf_logical'][-1]:.4f} ({time.time()-t0:.0f}s)", flush=True)
    # sanity: unperturbed code (h=0) has OTOC = 0 for commuting W,V and conserved logical
    return out


if __name__ == "__main__":
    h = float(sys.argv[1]) if len(sys.argv) > 1 else 0.3
    outp = sys.argv[2] if len(sys.argv) > 2 else "toric2d_quantum_otoc.json"
    json.dump(main(h), open(outp, "w"), indent=1)
