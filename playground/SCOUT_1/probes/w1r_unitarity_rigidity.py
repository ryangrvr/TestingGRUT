"""SCOUT-1 W1-R: does unitarity + conformal symmetry with c < 1 quantize the IR central charge, and does a
lattice parent land on the quantized value independent of microscopic couplings?

 A. Kac/FQS table: c(m) = 1 - 6/(m(m+1)), h_{r,s}; discreteness and accumulation at c = 1.
 B. free-fermion XY chains (BdG correlation matrices, N = 400 OBC): entanglement fit of c
    Ising line h = 1 for anisotropies g = 0.3, 0.6, 1.0 (-> 1/2), XX line (-> 1), off-critical (-> 0).
 C. interacting, non-integrable Ising-class point: TFIM + NNN coupling J2; locate h_c by the L*gap crossing,
    fit c from the PBC half-chain entanglement (ED).
 D. hostile at c = 1: XXZ chain, scaling dimension of the S^z=1 sector varies continuously with Delta
    (ED vs Bethe ansatz) -> no rigidity at c >= 1.
"""
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
from fractions import Fraction

np.set_printoptions(precision=5, suppress=True)

# ---------------- A. Kac table ----------------
print("=== A. FQS/Kac unitary minimal models (c < 1) ===")
for m in range(3, 9):
    c = Fraction(1) - Fraction(6, m * (m + 1))
    hs = sorted({Fraction(((m + 1) * r - m * s) ** 2 - 1, 4 * m * (m + 1)) for r in range(1, m) for s in range(1, r + 1)})
    print(f"  m={m}: c = {str(c):6s} = {float(c):.5f}; primaries h = {', '.join(str(h) for h in hs[:7])}{' ...' if len(hs) > 7 else ''}")
print("  c(m) -> 1 as m -> inf; unitary c < 1 is exactly this discrete set (FQS 1984, GKO coset construction)")


# ---------------- B. free-fermion XY chains ----------------
def xy_entropies(N, g, h):
    A = np.diag(np.full(N, -2.0 * h)) - (np.eye(N, k=1) + np.eye(N, k=-1))
    B = -g * (np.eye(N, k=1) - np.eye(N, k=-1))
    Hb = np.block([[A, B], [-B, -A]])
    E, W = np.linalg.eigh(Hb)
    Wp = W[:, E > 0]
    C = Wp @ Wp.conj().T          # <Psi Psi^dagger>, Psi = (c, c^dagger)
    S = []
    for l in range(1, N):
        idx = np.r_[0:l, N:N + l]
        nu = np.linalg.eigvalsh(C[np.ix_(idx, idx)])
        nu = np.clip(nu, 1e-15, 1 - 1e-15)
        S.append(-0.5 * np.sum(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)))
    return np.array(S)


def fit_c_obc(S, N, lo=0.1, hi=0.9):
    l = np.arange(1, N)
    x = np.log((2 * N / np.pi) * np.sin(np.pi * l / N))
    sel = (l > lo * N) & (l < hi * N)
    # average even/odd to suppress parity oscillations
    p = np.polyfit(x[sel], S[sel], 1)
    pe = np.polyfit(x[sel & (l % 2 == 0)], S[sel & (l % 2 == 0)], 1)
    po = np.polyfit(x[sel & (l % 2 == 1)], S[sel & (l % 2 == 1)], 1)
    return 6 * p[0], 6 * pe[0], 6 * po[0]


print("\n=== B. free-fermion XY chains, N = 400 OBC: c from S(l) = (c/6) log[(2N/pi) sin(pi l/N)] + s0 ===")
N = 400
for label, g, h in (("Ising line g=1.0, h=1", 1.0, 1.0), ("Ising line g=0.6, h=1", 0.6, 1.0),
                    ("Ising line g=0.3, h=1", 0.3, 1.0), ("XX line g=0, h=0", 0.0, 0.0),
                    ("XX line g=0, h=0.5", 0.0, 0.5), ("off-critical g=1, h=1.3", 1.0, 1.3),
                    ("off-critical g=0.6, h=0.7 (ordered)", 0.6, 0.7)):
    S = xy_entropies(N, g, h)
    c_all, c_e, c_o = fit_c_obc(S, N)
    print(f"  {label:38s}: c_fit = {c_all:.4f} (even l {c_e:.4f}, odd l {c_o:.4f});  S(N/2) = {S[N // 2 - 1]:.4f}")


# ---------------- C/D. ED tools (spin-1/2, PBC) ----------------
def bits(L):
    return np.arange(2 ** L, dtype=np.int64)


def tfim_nnn(L, J=1.0, J2=0.0, h=1.0):
    """H = -J sum ZZ_{i,i+1} - J2 sum ZZ_{i,i+2} - h sum X_i, PBC, full 2^L basis"""
    s = bits(L)
    z = 1 - 2 * ((s[:, None] >> np.arange(L)) & 1)
    diag = -J * np.sum(z * np.roll(z, -1, axis=1), axis=1) - J2 * np.sum(z * np.roll(z, -2, axis=1), axis=1)
    rows = [s]; cols = [s]; vals = [diag.astype(float)]
    for i in range(L):
        rows.append(s ^ (1 << i)); cols.append(s); vals.append(np.full(len(s), -h))
    return sps.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(2 ** L,) * 2)


def half_entropy(psi, L):
    M = psi.reshape(2 ** (L - L // 2), 2 ** (L // 2))  # bit order: low bits = first sites
    sv = np.linalg.svd(M, compute_uv=False) ** 2
    sv = sv[sv > 1e-16]
    return -np.sum(sv * np.log(sv))


def low(H, k=2):
    E, V = spla.eigsh(H, k=k, which="SA")
    o = np.argsort(E)
    return E[o], V[:, o]


print("\n=== C. interacting Ising-class point: TFIM + NNN J2 (non-integrable) ===")
for J2 in (0.0, 0.3):
    # locate h_c from the crossing of L*gap (Z2-even ground to first excited across sectors -> use 2 lowest)
    Ls = (10, 12, 14)
    def Lgap(L, h):
        E, _ = low(tfim_nnn(L, J2=J2, h=h), k=3)
        # at criticality the ground state is unique; the two lowest excitations are O(1/L)
        return L * (E[1] - E[0])
    hs = np.linspace(0.9, 2.2, 27) if J2 > 0 else np.linspace(0.85, 1.15, 13)
    curves = {L: np.array([Lgap(L, h) for h in hs]) for L in Ls}
    # crossing of L=12 and L=14
    d = curves[14] - curves[12]
    k = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
    hc = None
    if len(k):
        k = k[0]
        hc = hs[k] - d[k] * (hs[k + 1] - hs[k]) / (d[k + 1] - d[k])
        # refine by bisection
        a, b = hs[k], hs[k + 1]
        for _ in range(25):
            m = 0.5 * (a + b)
            dm = Lgap(14, m) - Lgap(12, m)
            if np.sign(dm) == np.sign(Lgap(14, a) - Lgap(12, a)):
                a = m
            else:
                b = m
        hc = 0.5 * (a + b)
    print(f"  J2={J2}: h_c (L=12/14 crossing of L*gap) = {hc:.5f}" + ("   (exact 1)" if J2 == 0 else ""))
    Ss, Lss = [], (8, 10, 12, 14, 16, 18)
    for L in Lss:
        E, V = low(tfim_nnn(L, J2=J2, h=hc), k=1)
        Ss.append(half_entropy(V[:, 0], L))
    Ss = np.array(Ss); x = np.log(np.array(Lss) / np.pi)
    c_pairs = [3 * (Ss[i + 1] - Ss[i]) / (x[i + 1] - x[i]) for i in range(len(Lss) - 1)]
    print(f"     S(L/2) PBC for L={Lss}: {np.round(Ss, 5)}")
    print(f"     c from successive pairs (S = (c/3) log(L/pi) + s0): {np.round(c_pairs, 4)}")
    # coupling-independence check slightly off h_c
    for dh in (+0.05, -0.05):
        S2 = [half_entropy(low(tfim_nnn(L, J2=J2, h=hc + dh), k=1)[1][:, 0], L) for L in (14, 18)]
        print(f"     h_c{dh:+.2f}: c_eff(14->18) = {3 * (S2[1] - S2[0]) / np.log(18 / 14):.4f}")


def xxz_sector(L, Delta, Sz):
    """H = sum (SxSx + SySy + Delta SzSz) PBC, fixed magnetization sector"""
    nup = L // 2 + Sz
    states = np.array([s for s in range(2 ** L) if bin(s).count("1") == nup], dtype=np.int64)
    index = {s: i for i, s in enumerate(states)}
    rows, cols, vals = [], [], []
    for a, s in enumerate(states):
        e = 0.0
        for i in range(L):
            j = (i + 1) % L
            bi, bj = (s >> i) & 1, (s >> j) & 1
            e += Delta * 0.25 * (1 if bi == bj else -1)
            if bi != bj:
                t = s ^ ((1 << i) | (1 << j))
                rows.append(index[t]); cols.append(a); vals.append(0.5)
        rows.append(a); cols.append(a); vals.append(e)
    return sps.csr_matrix((vals, (rows, cols)), shape=(len(states),) * 2)


print("\n=== D. hostile at c = 1: XXZ, x_1 = L (E_{Sz=1} - E_0) / (2 pi v) vs Bethe ansatz ===")
print("  Bethe: v = pi sqrt(1-D^2) / (2 arccos D);  x_1 = (pi - arccos D) / (2 pi)  [= 1/4 at D = 0]")
for D in (-0.6, -0.3, 0.0, 0.3, 0.6, 0.9):
    v = np.pi * np.sqrt(1 - D ** 2) / (2 * np.arccos(D))
    x_ba = (np.pi - np.arccos(D)) / (2 * np.pi)
    xs = []
    for L in (12, 16, 20):
        E0 = spla.eigsh(xxz_sector(L, D, 0), k=1, which="SA", return_eigenvectors=False)[0]
        E1 = spla.eigsh(xxz_sector(L, D, 1), k=1, which="SA", return_eigenvectors=False)[0]
        xs.append(L * (E1 - E0) / (2 * np.pi * v))
    print(f"  Delta={D:+.1f}: x_1(L=12,16,20) = {np.round(xs, 4)}   Bethe x_1 = {x_ba:.4f}")
print("  => at c = 1 the exponent is a continuous function of a microscopic coupling: no rigidity")
