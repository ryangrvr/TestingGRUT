"""SCOUT-2 S2-D-arrow numerics (D1-D10).  Theorem part: results/S2-D-arrow_D0_THEOREM.md.
Main model: real (time-reversal-symmetric) mixed-field Ising chain, H = sum Z_i Z_{i+1} + 0.9045 sum X_i + 0.809 sum Z_i
(chaotic, real symmetric -> complex conjugation K in the product basis is a Sigma-local antiunitary symmetry)."""
import itertools as it
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

rng = np.random.default_rng(314)
I2 = np.eye(2); X = np.array([[0, 1], [1, 0]], float); Z = np.diag([1., -1]); Yc = np.array([[0, -1j], [1j, 0]])
def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)
def op(n, d):
    out = np.eye(1)
    for q in range(n): out = np.kron(out, d.get(q, I2))
    return out
def mfi(n, hx=0.9045, hz=0.809):
    return sum(op(n, {i: Z, i + 1: Z}) for i in range(n - 1)) + hx * sum(op(n, {i: X}) for i in range(n)) + hz * sum(op(n, {i: Z}) for i in range(n))
def ent(p):
    p = p[p > 1e-14]; return float(-(p * np.log2(p)).sum())
def S_red(psi, n, keep):
    t = psi.reshape([2] * n); rest = [q for q in range(n) if q not in keep]
    m = np.transpose(t, list(keep) + rest).reshape(2 ** len(keep), -1)
    return ent(np.linalg.svd(m, compute_uv=False) ** 2)
def rho_red(rho, n, keep):
    t = rho.reshape([2] * (2 * n)); rest = [q for q in range(n) if q not in keep]
    m = np.transpose(t, list(keep) + rest + [x + n for x in keep] + [x + n for x in rest]).reshape(2 ** len(keep), 2 ** len(rest), 2 ** len(keep), 2 ** len(rest))
    return np.einsum('ajbj->ab', m)
def S_rho(r): return ent(np.linalg.eigvalsh((r + r.conj().T) / 2))
def td(a, b): return 0.5 * np.abs(np.linalg.eigvalsh(a - b)).sum()

n = 10; N = 2 ** n; SYS = [0, 1]
H = mfi(n); E, V = np.linalg.eigh(H)
def ev(psi, t): return V @ (np.exp(-1j * E * t) * (V.T @ psi))
def prod_state(angles):
    out = np.ones(1)
    for a in angles: out = np.kron(out, np.array([np.cos(a / 2), np.sin(a / 2)]))
    return out.astype(complex)
ts = [0, 0.5, 1, 2, 4, 8]

# ------------------------------------------------------------ D1
hdr("D1  subsystem entropy S(rho_S(t)), S = qubits {0,1} (max 2 bits), n = 10, every-state hostile")
p0 = prod_state(np.full(n, 1.4))
haar = rng.normal(size=N) + 1j * rng.normal(size=N); haar /= np.linalg.norm(haar)
traj = [S_red(ev(p0, t), n, SYS) for t in np.linspace(0, 8, 81)]
tstar = np.linspace(0, 8, 81)[int(np.argmax(traj))]
rev = np.conj(ev(p0, tstar))
states = {"product (real, tilt 1.4)": p0, "Haar": haar, f"time-reversed K psi_prod(t*={tstar:.1f})": rev,
          "energy eigenstate (mid)": V[:, N // 2].astype(complex)}
print("  state" + " " * 36 + "".join(f"t={t:<6}" for t in ts))
for nm, st in states.items():
    print(f"  {nm:<41}" + "".join(f"{S_red(ev(st, t), n, SYS):<8.3f}" for t in ts))
r = [S_red(ev(rev, s), n, SYS) for s in np.linspace(0, tstar, 9)]; f = [S_red(ev(p0, tstar - s), n, SYS) for s in np.linspace(0, tstar, 9)]
print(f"  mirror check (Corollary 2): max |S(rev(s)) - S(prod(t*-s))| = {max(abs(a - b) for a, b in zip(r, f)):.1e};"
      f" reversed state: S {r[0]:.3f} -> {r[-1]:.3f} (forward DEcrease)")
# recurrence: 3 qubits, S = qubit 0
n3 = 3; H3 = mfi(n3); E3, V3 = np.linalg.eigh(H3); c = V3.T @ prod_state(np.full(n3, 0.3)); w = np.abs(c) ** 2
best = (0, 0)
for k0 in range(0, 2_000_000, 200_000):
    tt = np.arange(k0, k0 + 200_000) * 0.05 + 50
    F = np.abs(np.exp(-1j * np.outer(tt, E3)) @ w)
    i = int(np.argmax(F))
    if F[i] > best[0]: best = (F[i], tt[i])
psi3 = V3 @ (np.exp(-1j * E3 * best[1]) * c)
print(f"  recurrence (3 qubits, Lemma 1): max fidelity over t in [50, 1e5] = {best[0]:.5f} at t = {best[1]:.2f};"
      f" S(qubit0) there = {S_red(psi3, n3, [0]):.4f} (started at 0; max reached early "
      f"{max(S_red(V3 @ (np.exp(-1j*E3*t)*c), n3, [0]) for t in np.linspace(0,10,101)):.3f})")

# ------------------------------------------------------------ D2
hdr("D2  typicality is not an arrow (Haar on the full space; and an energy shell)")
for label, gen in [("Haar on C^1024", lambda: (lambda v: v / np.linalg.norm(v))(rng.normal(size=N) + 1j * rng.normal(size=N))),
                   ("energy shell |E - 0| < 1 (Haar on the shell)", lambda: (lambda sel: (lambda v: V[:, sel] @ v / np.linalg.norm(v))(rng.normal(size=sel.sum()) + 1j * rng.normal(size=sel.sum())))(np.abs(E) < 1.0))]:
    S0, dS = [], []
    for _ in range(150):
        st = gen(); a = S_red(st, n, SYS); b = S_red(ev(st, 1.0), n, SYS); S0.append(a); dS.append(b - a)
    dS = np.array(dS)
    print(f"  {label:<44}: <S_S(0)> = {np.mean(S0):.3f} (max 2), <dS over dt=1> = {dS.mean():+.4f}, "
          f"increase {np.mean(dS > 0):.2f} / decrease {np.mean(dS < 0):.2f}")
print("  => typical states are already near-maximal; no entropy budget, no preferred direction: EQUILIBRIUM TYPICALITY != ARROW.")

# ------------------------------------------------------------ D3 / D6
hdr("D3/D6  same marginals, different S-E correlations (mixed states, n = 8, S = qubit 0, bath MAXIMALLY MIXED at start)")
n8 = 8; H8 = mfi(n8); E8, V8 = np.linalg.eigh(H8)
def U8(t): return (V8 * np.exp(-1j * E8 * t)) @ V8.T
rhoA = np.kron(np.diag([1., 0]), np.eye(2 ** 7) / 2 ** 7).astype(complex)
tau = 3.0
rt = U8(tau) @ rhoA @ U8(tau).conj().T
sS, sE = rho_red(rt, n8, [0]), rho_red(rt, n8, list(range(1, n8)))
C = np.conj(rt)                                    # correlated, time-reversed
D = np.kron(np.conj(sS), np.conj(sE))              # same marginals, no correlation
def report(nm, rho):
    out = []
    for t in [0, 1, 2, 3, 5]:
        r = U8(t) @ rho @ U8(t).conj().T
        rs = rho_red(r, n8, [0]); re = rho_red(r, n8, list(range(1, n8)))
        out.append((S_rho(rs), S_rho(rs) + S_rho(re) - S_rho(r), td(rs, np.eye(2) / 2)))
    print(f"  {nm:<44}" + "  ".join(f"{a:.3f}/{b:.3f}/{c_:.3f}" for a, b, c_ in out))
print("  columns t = 0,1,2,3,5 ; each entry S_S / I(S:E) / TD(rho_S, I/2)")
report("A product |0><0| (x) I/128", rhoA)
report("C correlated: K rho_A(tau) K", C)
report("D product of C's marginals", D)
for lam in [0.3, 0.7]:
    report(f"mixture lam*C + (1-lam)*D, lam={lam}", lam * C + (1 - lam) * D)
print(f"  bath entropy in C and D at t=0: {S_rho(rho_red(C, n8, list(range(1, n8)))):.3f} bits (max 7); marginals identical: "
      f"{np.allclose(rho_red(C, n8, [0]), rho_red(D, n8, [0])) and np.allclose(rho_red(C, n8, list(range(1,n8))), rho_red(D, n8, list(range(1,n8))))}")

# ------------------------------------------------------------ D4 / D5 collision models
hdr("D4/D5  collision models: partial swap theta = 0.6; ancilla preparation and freshness varied")
th = 0.6
SW = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], complex); Uc = expm(-1j * th * SW)
def collide(rho, nq, k):
    U = np.eye(1)
    # build U on (0,k) via permutation trick
    full = np.kron(Uc, np.eye(2 ** (nq - 2))).reshape([2] * (2 * nq))
    perm = list(range(nq)); perm[1], perm[k] = perm[k], perm[1]
    full = np.transpose(full, perm + [p + nq for p in perm]).reshape(2 ** nq, 2 ** nq)
    return full @ rho @ full.conj().T
def run(anc, order, label):
    M = int(np.log2(anc.shape[0])); nq = M + 1
    res = {}
    for nm, s0 in [("a", np.diag([1., 0])), ("b", np.full((2, 2), 0.5))]:
        rho = np.kron(s0, anc).astype(complex); tds = []; ent_ = []; recs = []
        for k in order:
            rho = collide(rho, nq, k)
            ent_.append(S_rho(rho_red(rho, nq, [0]))); recs.append(rho_red(rho, nq, [k]))
        res[nm] = (rho, ent_, recs)
    tdS = [td(rho_red(collide_seq_a, 1 + 0, [0]), 0) if False else None for collide_seq_a in []]
    # system TD per step
    ra = np.kron(np.diag([1., 0]), anc).astype(complex); rb = np.kron(np.full((2, 2), 0.5), anc).astype(complex); seq = []
    for k in order:
        ra, rb = collide(ra, nq, k), collide(rb, nq, k); seq.append(td(rho_red(ra, nq, [0]), rho_red(rb, nq, [0])))
    mono = all(seq[i + 1] <= seq[i] + 1e-12 for i in range(len(seq) - 1))
    rec_first = td(rho_red(ra, nq, [order[0]]), rho_red(rb, nq, [order[0]]))
    sysfin = rho_red(ra, nq, [0])
    print(f"  {label:<40} TD_sys: " + " ".join(f"{x:.3f}" for x in seq[::max(1, len(seq) // 8)]) +
          f" | monotone(Markov-like) {mono} | S_sys(a): {res['a'][1][0]:.3f}->{res['a'][1][-1]:.3f} | "
          f"rho_sys(a)[0,0] final {sysfin[0,0].real:.3f} | record in 1st ancilla at end (TD a vs b) {rec_first:.3f}")
M = 7
def prod_anc(s, M): 
    out = np.eye(1)
    for _ in range(M): out = np.kron(out, s)
    return out
fresh = list(range(1, M + 1))
run(prod_anc(np.diag([1., 0]), M), fresh, "fresh pure |0>")
run(prod_anc(np.diag([.8, .2]), M), fresh, "fresh thermal p0 = 0.8")
run(prod_anc(np.eye(2) / 2, M), fresh, "fresh MAXIMALLY MIXED")
cc = np.zeros((2 ** M, 2 ** M)); cc[0, 0] = cc[-1, -1] = 0.5
run(cc, fresh, "classically correlated (all-equal bits)")
g = np.zeros(2 ** M); g[0] = g[-1] = 1 / np.sqrt(2)
run(np.outer(g, g), fresh, "GHZ-entangled ancillas")
print("  -- freshness (D5): maximally mixed bath, 21 collisions --")
for Mr in [1, 2, 4]:
    run(prod_anc(np.eye(2) / 2, Mr), [1 + (i % Mr) for i in range(21)], f"finite bath M={Mr}, reused cyclically")
run(prod_anc(np.eye(2) / 2, 7), [1 + (i % 7) for i in range(21)], "finite bath M=7, reused cyclically")

# ------------------------------------------------------------ D7 records
hdr("D7  record arrow vs its time reverse (S + 6 env, H = sum g_k Z_S Y_k; Theta = X_S K is a symmetry)")
n7 = 7; gk = rng.uniform(0.8, 1.2, 6)
H7 = sum(gk[k] * op(n7, {0: Z, k + 1: Yc}) for k in range(6))
def U7(t): return expm(-1j * H7 * t)
def MI(psi, A, B):
    return S_red(psi, n7, A) + S_red(psi, n7, B) - S_red(psi, n7, sorted(A + B))
def redund(psi):
    HS = S_red(psi, n7, [0]); return sum(1 for k in range(1, n7) if HS > 0.1 and MI(psi, [0], [k]) >= 0.9 * HS)
p7 = prod_state([np.pi / 2] + [0] * 6); ts7 = np.linspace(0, np.pi / 4 / gk.mean(), 6)
fwd = [U7(t) @ p7 for t in ts7]
X0 = op(n7, {0: X}); revr = X0 @ np.conj(fwd[-1])
bwd = [U7(t) @ revr for t in ts7]
print("  t      " + "  ".join(f"{t:.3f}" for t in ts7))
print("  forward from |+>|0..>: <I(S:E_k)> " + "  ".join(f"{np.mean([MI(s,[0],[k]) for k in range(1,n7)]):.3f}" for s in fwd) + "   redundancy " + " ".join(str(redund(s)) for s in fwd))
print("  from Theta psi(t_f):    <I(S:E_k)> " + "  ".join(f"{np.mean([MI(s,[0],[k]) for k in range(1,n7)]):.3f}" for s in bwd) + "   redundancy " + " ".join(str(redund(s)) for s in bwd))
hs = (lambda v: v / np.linalg.norm(v))(rng.normal(size=2 ** n7) + 1j * rng.normal(size=2 ** n7))
print(f"  Haar global state: <I(S:E_k)> = {np.mean([MI(hs,[0],[k]) for k in range(1,n7)]):.3f}, redundancy {redund(hs)} (no records)")

# ------------------------------------------------------------ D8 Sigma dependence
hdr("D8  same global state and H, inequivalent TPSs (n = 10, S = first two slots)")
def householder_to_zero(psi):
    e = np.zeros_like(psi); e[0] = 1
    ph = psi[0] / abs(psi[0]) if abs(psi[0]) > 1e-12 else 1
    v = psi - ph * e; v /= np.linalg.norm(v)
    return np.eye(len(psi)) - 2 * np.outer(v, v.conj())   # maps psi -> ph*e (up to phase)
def cliff(nq, depth):
    Hd = np.array([[1, 1], [1, -1]]) / np.sqrt(2); U = np.eye(2 ** nq, dtype=complex)
    for _ in range(depth):
        if rng.random() < 0.5:
            a, b = rng.choice(nq, 2, replace=False)
            U = (op(nq, {a: np.diag([1., 0])}) + op(nq, {a: np.diag([0., 1]), b: X})) @ U
        else:
            U = op(nq, {int(rng.integers(nq)): Hd}) @ U
    return U
st = ev(p0, 4.0)        # entangled in the computational TPS
Wh = householder_to_zero(st)
frames = {"computational": np.eye(N), "Clifford frame": cliff(n, 40), "Householder frame (state -> product)": Wh}
for nm, W in frames.items():
    traj = [S_red(W @ ev(st, t), n, SYS) for t in ts]
    print(f"  {nm:<38} S_S(t) at t = {ts}: " + " ".join(f"{x:.3f}" for x in traj))
c_loc = lambda W: (lambda Hl: np.abs(Hl - np.diag(np.diag(Hl))).sum() / np.abs(Hl).sum())(W @ H @ W.conj().T)
print("  in the Householder frame the SAME state is product (S_S(0) = 0, a 'low-entropy environment'), but H written in")
print("  that frame is non-local; low correlation is a Sigma-relative notion -> arrow price is jointly H + Sigma.")

# ------------------------------------------------------------ D9 access
hdr("D9  access: full state vs local marginal vs coarse macro-variable (Z-magnetization histogram)")
def coarse(psi):
    pz = np.abs(psi) ** 2; w = np.array([bin(i).count("1") for i in range(N)])
    return ent(np.bincount(w, weights=pz, minlength=n + 1))
q0 = prod_state(np.full(n, 1.4)); q1 = prod_state(np.full(n, 1.41))
for nm, a, b in [("product start", q0, q1), ("time-reversed start", np.conj(ev(q0, 4.0)), np.conj(ev(q1, 4.0)))]:
    row = []
    for t in [0, 1, 2, 4]:
        A_, B_ = ev(a, t), ev(b, t)
        row.append(f"t={t}: global TD {np.sqrt(1 - abs(np.vdot(A_, B_))**2):.4f}, S_S {S_red(A_, n, SYS):.3f}, coarse S {coarse(A_):.3f}")
    print(f"  {nm}:\n     " + "\n     ".join(row))

# ------------------------------------------------------------ D10 Janus and the (H, Sigma)-fixed special state
hdr("D10  Janus: entropy on both sides of a special state; a special state fixed by (H, Sigma)")
sym = max(abs(S_red(ev(p0, t), n, SYS) - S_red(ev(p0, -t), n, SYS)) for t in np.linspace(0, 6, 13))
print(f"  real product state: max |S_S(t) - S_S(-t)| over t in [0,6] = {sym:.1e}  -> two-sided (Janus) arrow; orientation NOT selected")
def energy(angles): 
    s = prod_state(angles); return float(np.real(np.vdot(s, H @ s)))
sols = []
for _ in range(12):
    r_ = minimize(energy, rng.uniform(0, 2 * np.pi, n), method="BFGS")
    sols.append((round(r_.fun, 6), np.mod(r_.x, 2 * np.pi)))
emin = min(s_[0] for s_ in sols); mins = [s_ for s_ in sols if abs(s_[0] - emin) < 1e-5]
mf = prod_state(mins[0][1])
distinct = []
for m_ in mins:
    v_ = prod_state(m_[1])
    if all(abs(np.vdot(v_, d_)) < 0.999 for d_ in distinct): distinct.append(v_)
P_ref = np.zeros((N, N))
for i_ in range(N): P_ref[int(format(i_, f"0{n}b")[::-1], 2), i_] = 1
spread = len(distinct)
refl_related = len(distinct) == 2 and abs(np.vdot(P_ref @ distinct[0], distinct[1])) > 0.999
print(f"  mean-field (minimal-energy product) state: E = {emin:.4f} (ground {E[0]:.4f}); found in {len(mins)}/12 restarts, "
      f"distinct minimizers {spread} (reflection-related pair: {refl_related if spread == 2 else 'n/a'});  energy variance {float(np.real(np.vdot(mf, H@H@mf))) - emin**2:.3f}")
for j_, d_ in enumerate(distinct):
    print(f"  minimizer {j_}: S_S(t), t = -4..4: " + " ".join(f"{S_red(ev(d_, t), n, SYS):.3f}" for t in [-4, -2, -1, 0, 1, 2, 4])
          + f";  half-chain S(t=4) = {S_red(ev(d_, 4), n, list(range(n//2))):.3f}")
print("  stationarity of an H-only state (Proposition 2): ground state S_S(t) = " + " ".join(f"{S_red(ev(V[:,0].astype(complex), t), n, SYS):.3f}" for t in [0, 2, 4]))
# alternative (H, Sigma)-covariant selector functionals
alt = {}
for nm, fn in [("max-energy product state", lambda a: -energy(a)),
               ("min energy-variance product state", lambda a: (lambda s_: float(np.real(np.vdot(s_, H @ H @ s_)) - np.real(np.vdot(s_, H @ s_)) ** 2))(prod_state(a)))]:
    best_ = min((minimize(fn, rng.uniform(0, 2 * np.pi, n), method="BFGS") for _ in range(8)), key=lambda r_: r_.fun)
    v_ = prod_state(best_.x); alt[nm] = v_
    print(f"  alternative selector '{nm}': E = {np.real(np.vdot(v_, H @ v_)):.3f}; S_S(t), t=-4..4: "
          + " ".join(f"{S_red(ev(v_, t), n, SYS):.3f}" for t in [-4, -2, -1, 0, 1, 2, 4]) + f"; overlap with mean-field {abs(np.vdot(v_, mf)):.3f}")
print("  => several (H, Sigma)-covariant extremal principles each fix a different special Janus state: the boundary condition")
print("     is fixed by (H, Sigma) ONLY after a selector functional is chosen (preference-priced, as in S2-Sigma).")
