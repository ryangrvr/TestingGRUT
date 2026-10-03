"""RA0 / G2 -- asymptotic autonomy.  NUMERICAL ILLUSTRATION; independent code path, not independent reviewer.

No epsilon is chosen anywhere: every defect is reported as a sequence in N (or L) and judged only by its limit trend.

A  independent two-state particles (product kernel)            -- symmetry / slow-mode control
B  symmetric simple exclusion (SSEP) on a ring, discrete time   -- interacting model with a known hydrodynamic limit
   B1 exact-block (sup) lumpability defect of the two-cell count partition at diffusive times (small L, exact)
   B2 mean-density sup defect of cell partitions via the exact one-particle heat kernel (SSEP duality; large L)
   B3 slow spectrum of the full many-particle kernel vs the one-particle kernel; slow eigenvectors linear in occupations?
C  two-lane ladder (one-particle generator): does the slow set depend on the supplied family scaling gamma_L = L^-beta?
D  twin test: identical SSEP dynamics, contiguous cells vs interleaved classes
"""
import itertools
import numpy as np
import scipy.sparse as sp

np.set_printoptions(precision=4, suppress=True)


# ---------------------------------------------------------------- A: independent particles
def control_A():
    print("== A  independent two-state particles, flip probs a=0.1 (0->1), b=0.25 (1->0) ==")
    a, b = 0.1, 0.25
    p = np.array([[1 - a, a], [b, 1 - b]])
    for N in [2, 4, 6, 8, 10]:
        P = p
        for _ in range(N - 1):
            P = np.kron(P, p)
        states = list(itertools.product([0, 1], repeat=N))
        occ = np.array([sum(s) for s in states]); first = np.array([s[0] for s in states])

        def defect(lab):
            k = lab.max() + 1
            S = np.column_stack([P[:, lab == j].sum(1) for j in range(k)])
            return max(np.abs(S[lab == i] - S[lab == i][0]).max() for i in range(k) if np.any(lab == i))
        rng = np.random.default_rng(N); rnd = rng.integers(0, 3, len(states))
        ev = np.sort(np.abs(np.linalg.eigvals(P)))[::-1]
        lam1 = 1 - a - b
        mult = int(np.sum(np.isclose(ev, abs(lam1), atol=1e-8)))
        print(f"  N={N:2d}: one-step defect  occupation-number {defect(occ):.1e}  first-particle {defect(first):.1e}"
              f"  random 3-block {defect(rnd):.3f}   slowest nontrivial |eig| {abs(lam1):.3f} (N-independent)"
              f" multiplicity {mult} (= N single-particle modes)")


# ---------------------------------------------------------------- B: SSEP on a ring
def ssep(L, N):
    confs = [c for c in itertools.combinations(range(L), N)]
    idx = {c: i for i, c in enumerate(confs)}
    rows, cols, vals = [], [], []
    for i, c in enumerate(confs):
        occ = set(c)
        for bnd in range(L):
            u, v = bnd, (bnd + 1) % L
            if (u in occ) != (v in occ):
                new = set(occ); new.symmetric_difference_update({u, v})
                rows.append(i); cols.append(idx[tuple(sorted(new))]); vals.append(1.0 / L)
            else:
                rows.append(i); cols.append(i); vals.append(1.0 / L)
    P = sp.csr_matrix((vals, (rows, cols)), shape=(len(confs), len(confs)))
    return confs, P


def control_B1():
    print("== B1 SSEP ring, N=L/2: exact-block sup defect of the left-half-count partition at time tau = c L^3 ==")
    for c in [0.02, 0.1]:
        out = []
        for L in [6, 8, 10, 12, 14]:
            N = L // 2; confs, P = ssep(L, N)
            lab = np.array([sum(1 for s in cf if s < L // 2) for cf in confs]); K = lab.max() + 1
            E = np.column_stack([(lab == j).astype(float) for j in range(K)])
            R = E.copy(); tau = int(round(c * L ** 3))
            for _ in range(tau):
                R = P @ R
            sup = max(np.abs(R[lab == i][:, None, :] - R[lab == i][None, :, :]).max() for i in range(K))
            avg = np.mean([np.abs(R[lab == i] - R[lab == i].mean(0)).max(1).mean() for i in range(K)])
            out.append(f"L={L}:sup {sup:.3f} / avg {avg:.3f}")
        print(f"  c={c}: " + "  ".join(out))


def heat_kernel(L, tau):
    """one-particle kernel of the same bond dynamics: stay 1-2/L, hop +-1 w.p. 1/L each; tau steps, exact via FFT"""
    q = np.arange(L)
    lam = 1 - (2.0 / L) * (1 - np.cos(2 * np.pi * q / L))
    return np.real(np.fft.ifft(np.exp(tau * np.log(lam.astype(complex)))))    # K^tau(0, d), d = 0..L-1


def mean_defect(L, classes, tau):
    """sup over configurations with equal class counts of |E[count in class j at tau | x] - same | x'| / |class|"""
    k0 = heat_kernel(L, tau)
    m = classes.max() + 1; worst = 0.0
    for j in range(m):
        tgt = np.where(classes == j)[0]
        w = np.array([k0[(tgt - y) % L].sum() for y in range(L)])      # P(particle from y ends in class j)
        tot = 0.0
        for i in range(m):
            v = np.sort(w[classes == i])
            top = np.cumsum(v[::-1]); bot = np.cumsum(v)
            tot += np.max(top - bot)
        worst = max(worst, tot / len(tgt))
    return worst


def control_B2():
    print("== B2 SSEP mean-density sup defect via exact heat kernel (duality), diffusive time tau = 0.05 L^3 ==")
    for L in [64, 256, 1024, 4096]:
        tau = 0.05 * L ** 3
        contig4 = (np.arange(L) * 4) // L
        ell = int(round(np.sqrt(L))); contig_sqrt = np.arange(L) // ell
        inter4 = np.arange(L) % 4
        print(f"  L={L:5d}: 4 contiguous cells {mean_defect(L, contig4, tau):.4f}   "
              f"sqrt(L)-size cells (m={L // ell}) {mean_defect(L, contig_sqrt, tau):.4f}   "
              f"4 interleaved classes {mean_defect(L, inter4, tau):.4f}")


def control_B3():
    print("== B3 slow spectrum: full SSEP kernel vs one-particle kernel; slow eigenvectors linear in occupations? ==")
    for L in [8, 10, 12]:
        N = L // 2; confs, P = ssep(L, N)
        Pd = P.toarray(); w, V = np.linalg.eigh(Pd)
        gaps = np.sort(1 - w)
        one = np.sort((2.0 / L) * (1 - np.cos(2 * np.pi * np.arange(L) / L)))
        occ = np.array([[1.0 if s in cf else 0.0 for s in range(L)] for cf in confs])
        lin = np.column_stack([np.ones(len(confs)), occ])
        Q, _ = np.linalg.qr(lin)
        nz = np.where(1 - w > 1e-9)[0]
        order = nz[np.argsort((1 - w)[nz])]
        res = []
        for idx in order[:6]:
            v = V[:, idx]; r = np.linalg.norm(v - Q @ (Q.T @ v)) / np.linalg.norm(v)
            res.append(r)
        print(f"  L={L:2d} ({len(confs)} states): lowest many-body gaps {gaps[1:6]}  one-particle {one[1:6]}")
        print(f"         residual of 6 slowest eigenvectors outside span(1, eta_x): {np.array(res)}")
    print("  one-particle ratio lambda_q/lambda_1 (q=1..5) as L grows:")
    for L in [10, 100, 1000, 10000]:
        lam = (2.0 / L) * (1 - np.cos(2 * np.pi * np.arange(1, 6) / L))
        print(f"    L={L:5d}: {lam / lam[0]}   fast-mode ratio lambda_(L/2)/lambda_1 = "
              f"{(2.0 / L) * 2 / lam[0]:.3e}")


# ---------------------------------------------------------------- C: two-lane ladder, family scaling
def control_C():
    print("== C  two-lane ladder (one particle): lane-imbalance mode vs along-lane slow modes, gamma_L = L^-beta ==")
    for beta in [0, 1, 2, 3]:
        row = []
        for L in [50, 200, 800, 3200]:
            g = L ** (-float(beta))
            lam1 = 2 * (1 - np.cos(2 * np.pi / L))          # slowest along-lane density mode
            lane = 2 * g                                     # lane-imbalance mode (q=0, antisymmetric)
            row.append(f"L={L}: {lane / lam1:.3g}")
        print(f"  beta={beta}: ratio (lane-imbalance rate)/(slowest density rate): " + "  ".join(row))


# ---------------------------------------------------------------- D: twin test
def control_D():
    print("== D  twin test on identical SSEP dynamics: class-difference relaxation rate / slowest rate ==")
    for L in [64, 256, 1024, 4096]:
        lam1 = (2.0 / L) * (1 - np.cos(2 * np.pi / L))
        # contiguous 4 cells: class-difference observables contain the q=1 Fourier mode
        contig = 1.0
        # interleaved 4 classes (site mod 4): class-difference observables live at q = L/4, L/2, 3L/4
        lam_inter = (2.0 / L) * (1 - np.cos(2 * np.pi * (L // 4) / L))
        print(f"  L={L:5d}: contiguous-cell differences slowest ratio {contig:.1f}; "
              f"interleaved-class differences ratio {lam_inter / lam1:.3e}")


if __name__ == "__main__":
    control_A(); control_B1(); control_B2(); control_B3(); control_C(); control_D()
