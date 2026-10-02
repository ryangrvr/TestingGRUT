"""INDEPENDENT REPRODUCTION of S2-D-arrow key numbers (NR-D1..D8, priority order from the owner). No S2-D-arrow code.
Routines: Krylov time evolution (scipy expm_multiply) for pure states; n-general partial traces written here;
eigvalsh entropies; QR-built product frame; fresh Clifford (H/S/CZ). Each phenomenon at the ORIGINAL parameters
(number comparison) and at a VARIANT parameter set (robustness)."""
import itertools as it
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply
from scipy.linalg import expm

def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)
I2 = np.eye(2); X = np.array([[0, 1], [1, 0]], float); Z = np.diag([1.0, -1.0]); Yc = np.array([[0, -1j], [1j, 0]])
def sop(n, d):
    out = sp.identity(1, format="csr", dtype=complex)
    for q in range(n): out = sp.kron(out, sp.csr_matrix(d.get(q, I2)), format="csr")
    return out
def mfi(n, hx=0.9045, hz=0.809):
    return sum(sop(n, {i: Z, i + 1: Z}) for i in range(n - 1)) + hx * sum(sop(n, {i: X}) for i in range(n)) + hz * sum(sop(n, {i: Z}) for i in range(n))
def prod(n, a):
    v = np.ones(1, complex)
    for _ in range(n): v = np.kron(v, [np.cos(a / 2), np.sin(a / 2)])
    return v
def rdm_pure(psi, n, keep):
    t = psi.reshape([2] * n); rest = [q for q in range(n) if q not in keep]
    A = np.transpose(t, list(keep) + rest).reshape(2 ** len(keep), -1); return A @ A.conj().T
def rdm_mixed(rho, n, keep):
    t = rho.reshape([2] * (2 * n)); cur = n
    for q in sorted(set(range(n)) - set(keep), reverse=True):
        t = np.trace(t, axis1=q, axis2=q + cur); cur -= 1
    return t.reshape(2 ** len(keep), 2 ** len(keep))
def S(r):
    e = np.linalg.eigvalsh((r + r.conj().T) / 2); e = e[e > 1e-13]; return float(-(e * np.log2(e)).sum())
def traj(H, psi, ts):
    return expm_multiply(-1j * H, psi, start=ts[0], stop=ts[-1], num=len(ts), endpoint=True)
def td(a, b): return 0.5 * np.abs(np.linalg.eigvalsh(a - b)).sum()
RES = {}

# ---------------- NR-D1 exact time-reversal mirror / anti-arrow
hdr("NR-D1 exact time-reversal mirror / anti-arrow (real H => complex conjugation is a Sigma-local symmetry)")
for label, n, a, sys_ in [("original params n=10 tilt 1.4 S={0,1}", 10, 1.4, [0, 1]), ("variant n=8 tilt 1.1 S={3}", 8, 1.1, [3])]:
    H = mfi(n); ts = np.linspace(0, 8, 81); P = traj(H, prod(n, a), ts)
    Sv = np.array([S(rdm_pure(p, n, sys_)) for p in P]); k = int(np.argmax(Sv)); tstar = ts[k]
    rev = np.conj(P[k]); R = traj(H, rev, np.linspace(0, tstar, k + 1))
    Sr = np.array([S(rdm_pure(p, n, sys_)) for p in R])
    mirror = np.abs(Sr - Sv[:k + 1][::-1]).max()
    print(f"  {label}: t* = {tstar:.1f}, S_S(t*) = {Sv[k]:.3f}; reversed state forward: S {Sr[0]:.3f} -> {Sr[-1]:.3f}; mirror error {mirror:.1e}")
    RES.setdefault("D1", []).append(Sr[-1] < 1e-6 and mirror < 1e-8 and Sr[0] > 0.5)

# ---------------- NR-D3/D6 same marginals, different correlations
hdr("NR-D3/D6 same marginals, different correlations (mixed, S = qubit 0, bath initially maximally mixed)")
for label, n, tau in [("original params n=8 tau=3", 8, 3.0), ("variant n=7 tau=2", 7, 2.0)]:
    H = mfi(n).toarray(); U = lambda t: expm(-1j * H * t)
    rA = np.kron(np.diag([1.0, 0]), np.eye(2 ** (n - 1)) / 2 ** (n - 1)).astype(complex)
    rt = U(tau) @ rA @ U(tau).conj().T
    C = np.conj(rt); sS, sE = rdm_mixed(C, n, [0]), rdm_mixed(C, n, list(range(1, n))); D = np.kron(sS, sE)
    same = np.allclose(rdm_mixed(D, n, [0]), sS) and np.allclose(rdm_mixed(D, n, list(range(1, n))), sE)
    def SS(r, t): u = U(t); return S(rdm_mixed(u @ r @ u.conj().T, n, [0]))
    c_traj = [SS(C, t) for t in (0, tau / 3, 2 * tau / 3, tau)]; d_traj = [SS(D, t) for t in (0, tau / 3, 2 * tau / 3, tau)]
    print(f"  {label}: marginals identical {same}; bath entropy {S(sE):.3f} of {n-1}; correlated C: S_S " + " ".join(f"{x:.3f}" for x in c_traj)
          + " | product-of-marginals D: S_S " + " ".join(f"{x:.3f}" for x in d_traj))
    RES.setdefault("D3", []).append(same and c_traj[-1] < 1e-6 and min(d_traj) > 0.5)

# ---------------- NR-D4 / D5 collision models
hdr("NR-D4/D5 partial-swap collisions (theta = 0.6): fresh maximally mixed bath vs correlated / reused bath")
th = 0.6; SW = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], complex); Uc = expm(-1j * th * SW)
def channel_fresh(rho_s, sigma):
    """exact reduced map for a FRESH uncorrelated ancilla sigma (Stinespring, written independently)"""
    J = Uc @ np.kron(rho_s, sigma) @ Uc.conj().T
    return np.einsum('aibi->ab', J.reshape(2, 2, 2, 2))
def joint_run(anc, order):
    M = int(np.log2(anc.shape[0])); n = M + 1; out = {}
    for nm, s0 in (("a", np.diag([1.0, 0])), ("b", np.full((2, 2), 0.5))):
        rho = np.kron(s0, anc).astype(complex); seq = []
        for k in order:
            t = rho.reshape([2] * (2 * n)); perm = list(range(n)); perm[1], perm[k] = perm[k], perm[1]
            t = np.transpose(t, perm + [p + n for p in perm]).reshape(2 ** n, 2 ** n)
            G = np.kron(Uc, np.eye(2 ** (n - 2))); t = G @ t @ G.conj().T
            t = np.transpose(t.reshape([2] * (2 * n)), perm + [p + n for p in perm]).reshape(2 ** n, 2 ** n)
            rho = t; seq.append(rdm_mixed(rho, n, [0]))
        out[nm] = (rho, seq)
    tds = [td(a, b) for a, b in zip(out["a"][1], out["b"][1])]
    return out, tds
# fresh maximally mixed, via the exact channel (independent of the joint simulation)
ra, rb = np.diag([1.0, 0]).astype(complex), np.full((2, 2), 0.5, complex); tds, Ss = [], []
for _ in range(7):
    ra, rb = channel_fresh(ra, np.eye(2) / 2), channel_fresh(rb, np.eye(2) / 2); tds.append(td(ra, rb)); Ss.append(S(ra))
mono = all(tds[i + 1] <= tds[i] + 1e-12 for i in range(6))
# record in ancilla 1 (joint simulation, 7 fresh max-mixed ancillas)
M = 7; anc = np.eye(2 ** M) / 2 ** M
out, tds_j = joint_run(anc, list(range(1, M + 1)))
rec = td(rdm_mixed(out["a"][0], M + 1, [1]), rdm_mixed(out["b"][0], M + 1, [1]))
print(f"  fresh MAX-MIXED bath (channel): TD " + " ".join(f"{x:.3f}" for x in tds) + f" | monotone {mono} | S_sys {Ss[0]:.3f} -> {Ss[-1]:.3f};"
      f" joint-sim TD agrees: {np.allclose(tds, tds_j)}; record in ancilla 1: {rec:.3f}")
RES["D4"] = [mono and Ss[-1] > 0.99 and rec > 0.1 and np.allclose(tds, tds_j)]
cc = np.zeros((2 ** M, 2 ** M)); cc[0, 0] = cc[-1, -1] = 0.5
_, tds_c = joint_run(cc, list(range(1, M + 1))); mono_c = all(tds_c[i + 1] <= tds_c[i] + 1e-12 for i in range(len(tds_c) - 1))
print(f"  classically correlated bath (same max-mixed marginals): TD " + " ".join(f"{x:.3f}" for x in tds_c) + f" | monotone {mono_c}")
RES["D4"].append(not mono_c)
for Mr in (1, 2):
    outr, tdr = joint_run(np.eye(2 ** Mr) / 2 ** Mr, [1 + (i % Mr) for i in range(21)])
    Sa = [S(r) for r in outr["a"][1]]; monor = all(tdr[i + 1] <= tdr[i] + 1e-12 for i in range(20))
    print(f"  D5 reused bath M={Mr}: monotone TD {monor}; S_sys first {Sa[0]:.3f}, min over 21 collisions {min(Sa):.3f} (at {int(np.argmin(Sa))+1}), final {Sa[-1]:.3f}")
    RES.setdefault("D5", []).append(not monor)
RES["D5"].append(min(S(r) for r in joint_run(np.eye(2) / 2, [1] * 21)[0]["a"][1]) < 0.05)

# ---------------- NR-D7 record formation vs un-formation
hdr("NR-D7 record arrow vs its time reverse (S + 6 env, H = sum g_k Z_S Y_k; Theta = X_S K)")
for label, seed in (("seed A", 1), ("seed B", 2)):
    rng = np.random.default_rng(seed); n = 7; g = rng.uniform(0.8, 1.2, 6)
    H = sum(g[k] * sop(n, {0: Z, k + 1: Yc}) for k in range(6))
    p0 = np.zeros(2 ** n, complex); p0[0] = 1; p0 = (sop(n, {0: np.array([[1, 1], [1, -1]]) / np.sqrt(2)}) @ p0)
    tf = np.pi / 4 / g.mean(); ts = np.linspace(0, tf, 6)
    def mi(psi, k): return S(rdm_pure(psi, n, [0])) + S(rdm_pure(psi, n, [k])) - S(rdm_pure(psi, n, [0, k]))
    def red(psi):
        HS = S(rdm_pure(psi, n, [0])); return sum(1 for k in range(1, n) if HS > 0.1 and mi(psi, k) >= 0.9 * HS)
    F = traj(H, p0, ts); rev = sop(n, {0: X}) @ np.conj(F[-1]); B = traj(H, rev, ts)
    fwd = [red(p) for p in F]; bwd = [red(p) for p in B]
    print(f"  {label}: forward redundancy {fwd}; from Theta psi(t_f): {bwd}; <I(S:E_k)> fwd {np.mean([mi(F[-1], k) for k in range(1, n)]):.3f}, after reversal {np.mean([mi(B[-1], k) for k in range(1, n)]):.3f}")
    RES.setdefault("D7", []).append(fwd[0] == 0 and fwd[-1] >= 5 and bwd[0] >= 5 and bwd[-1] == 0)

# ---------------- NR-D8 Sigma dependence
hdr("NR-D8 same global state, inequivalent TPSs (n = 10, S = first two slots)")
n = 10; H = mfi(n); rng = np.random.default_rng(8)
st = traj(H, prod(n, 1.4), np.array([0.0, 4.0]))[-1]
A = np.c_[st, rng.normal(size=(2 ** n, 2 ** n - 1)) + 1j * rng.normal(size=(2 ** n, 2 ** n - 1))]; Q, Rr = np.linalg.qr(A); Q[:, 0] *= Rr[0, 0] / abs(Rr[0, 0])
Hd = np.array([[1, 1], [1, -1]]) / np.sqrt(2); C = np.eye(2 ** n, dtype=complex)
for _ in range(60):
    if rng.random() < 0.5:
        a, b = rng.choice(n, 2, replace=False); C = (sop(n, {a: np.diag([1, 0]), b: I2}) + sop(n, {a: np.diag([0, 1]), b: Z})) @ C
    else: C = sop(n, {int(rng.integers(n)): Hd}) @ C
ts = np.array([0, 0.5, 1, 2, 4, 8.0]); T = traj(H, st, ts)
for nm, U in (("computational", None), ("fresh Clifford frame", C), ("QR product frame (state -> |0..0>)", Q)):
    vals = [S(rdm_pure(p if U is None else U.conj().T @ p, n, [0, 1])) for p in T]
    print(f"  {nm:<36}: S_S(t) " + " ".join(f"{v:.3f}" for v in vals))
    RES.setdefault("D8", []).append(vals)
d8 = RES["D8"]; RES["D8"] = [d8[2][0] < 1e-6 and d8[2][1] > 0.3 and min(d8[1]) > 1.5]

# ---------------- NR-D2 typicality
hdr("NR-D2 typicality: Haar on an energy shell |E| < 1, n = 10, S = {0,1}")
H = mfi(n).toarray(); ev, V = np.linalg.eigh(H); sel = np.abs(ev) < 1.0; rng = np.random.default_rng(22)
for dt in (1.0, 0.3):
    d = []; S0 = []
    for _ in range(1500):
        v = rng.normal(size=sel.sum()) + 1j * rng.normal(size=sel.sum()); psi = V[:, sel] @ v / np.linalg.norm(v)
        a = S(rdm_pure(psi, n, [0, 1])); psit = V @ (np.exp(-1j * ev * dt) * (V.conj().T @ psi)); b = S(rdm_pure(psit, n, [0, 1]))
        d.append(b - a); S0.append(a)
    d = np.array(d); p = np.mean(d > 0)
    print(f"  dt={dt}: <S_S(0)> = {np.mean(S0):.3f} (max 2); <dS> = {d.mean():+.5f} +- {d.std()/np.sqrt(len(d)):.5f}; P(increase) = {p:.3f} +- {np.sqrt(.25/len(d)):.3f}")
    RES.setdefault("D2", []).append(abs(p - 0.5) < 3 * np.sqrt(0.25 / len(d)) and np.mean(S0) > 1.9)

hdr("SUMMARY")
claims = {"D1": "exact time-reversal mirror; reversed state runs an anti-arrow to S = 0",
          "D3": "same marginals: correlated state recoheres to purity, product-of-marginals does not (max-mixed-ish bath)",
          "D4": "fresh maximally mixed bath: Markovian contraction, entropy rise, records; correlated bath: non-monotone",
          "D5": "bath reuse destroys monotone (semigroup) behaviour; single reused ancilla re-purifies",
          "D7": "records form forward (0 -> >=5) and un-form from the reversed state (>=5 -> 0)",
          "D8": "same global state: product/rising in one TPS, flat high-entropy in another",
          "D2": "typicality: P(increase) ~ 1/2, states already near maximal"}
for k in ("D1", "D3", "D4", "D5", "D7", "D8", "D2"):
    v = RES[k]; print(f"NR-{k}: {'REPRODUCED' if all(v) else 'PARTIAL' if any(v) else 'NOT REPRODUCED'} ({sum(map(bool, v))}/{len(v)}) — {claims[k]}")
