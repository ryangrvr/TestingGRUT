"""INDEPENDENT REPRODUCTION of S2-H2 priorities (NR-H1..H5). No S2-H2 code.
Routines: fixed-step RK4 (own) instead of solve_ivp; random-sequential Glauber (own) instead of checkerboard; exact
rational arithmetic for invertibility; fresh models / seeds."""
from fractions import Fraction as Fr
import numpy as np
from scipy.linalg import expm

def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)
rng = np.random.default_rng(2718); RES = {}

# ---------------- NR-H1 structural positive controls
hdr("NR-H1 structural positive controls: symmetric contraction; doubly stochastic primitive chain")
c, ph = 0.7, 1.3; R = np.array([[np.cos(ph), -np.sin(ph)], [np.sin(ph), np.cos(ph)]])
starts = [np.array(v, float) for v in ([5e5, -2e5], [0, 0], [-3, 9], [1e-12, 0])]
fin = []
for x in starts:
    for _ in range(300): x = c * R @ x
    fin.append(np.linalg.norm(x))
print(f"  F(x) = 0.7 R x (equivariant about 0, no offset): max |x_300| over starts incl. 5e5 = {max(fin):.1e}  -> attractor 0 fixed by symmetry")
perms = [np.eye(6)[rng.permutation(6)] for _ in range(4)]; P = 0.4 * np.eye(6) + 0.15 * sum(perms)
prim = (np.linalg.matrix_power(P, 12) > 0).all(); dev = 0
for p0 in [np.eye(6)[i] for i in range(6)] + [rng.dirichlet(np.ones(6)) for _ in range(5)]:
    dev = max(dev, np.abs(p0 @ np.linalg.matrix_power(P, 80) - 1 / 6).max())
Pp = np.abs(P + 0.03 * rng.normal(size=P.shape)); Pp /= Pp.sum(1, keepdims=True)
w, V = np.linalg.eig(Pp.T); pi = np.real(V[:, np.argmin(abs(w - 1))]); pi /= pi.sum()
print(f"  doubly stochastic primitive 6-state chain (row/col sums {P.sum(1).max():.3f}/{P.sum(0).max():.3f}; primitive {prim}): max deviation from uniform after 80 steps over 11 starts = {dev:.1e}")
print(f"  perturbed (row-renormalized noise 0.03): stationary pi = {np.round(pi, 4)} (max |pi - 1/6| = {np.abs(pi - 1/6).max():.4f}) -> uniformity is tuned")
RES["H1"] = [max(fin) < 1e-30, prim and dev < 1e-12, np.abs(pi - 1 / 6).max() > 1e-3]

# ---------------- NR-H2 explicit-target relocation
hdr("NR-H2 explicit-target relocation: the attractor value tracks a state-valued parameter of the law")
for b in ([1.0, 0.0], [-2.0, 3.0]):
    b = np.array(b); xs = np.linalg.solve(np.eye(2) - c * R, b); x = np.array([40.0, -40.0])
    for _ in range(300): x = c * R @ x + b
    print(f"  F(x) = 0.7 R x + b, b = {b}: x_300 = {np.round(x, 6)}; (I - cR)^-1 b = {np.round(xs, 6)}")
def rk4(f, x, T, h=1e-3):
    for _ in range(int(T / h)):
        k1 = f(x); k2 = f(x + h / 2 * k1); k3 = f(x + h / 2 * k2); k4 = f(x + h * k3); x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return x
mins = []
for a in (0.3, -1.1):
    ends = [rk4(lambda x: -(x ** 3 + x - a), x0, 40.0) for x0 in (-8.0, 0.0, 8.0)]
    root = np.real([r for r in np.roots([1, 0, 1, -a]) if abs(r.imag) < 1e-9][0])
    print(f"  gradient flow V = x^4/4 + x^2/2 - a x, a = {a}: endpoints from -8, 0, 8: {np.round(ends, 6)}; unique minimum (root of x^3+x-a) = {root:.6f}")
    mins.append(max(abs(e - root) for e in ends) < 1e-6)
RES["H2"] = mins

# ---------------- NR-H3 all-to-all vs local coupling
hdr("NR-H3 alignment: all-to-all vs nearest-neighbour ring (identical oscillators, own RK4)")
def kuramoto(th, A, T=200.0, h=0.05):
    deg = A.sum(1)
    f = lambda t_: (A * np.sin(t_[None, :] - t_[:, None])).sum(1) / deg
    return rk4(f, th, T, h)
Nk = 40; comp = np.ones((Nk, Nk)) - np.eye(Nk); ring = np.zeros((Nk, Nk))
for i in range(Nk): ring[i, (i + 1) % Nk] = ring[i, (i - 1) % Nk] = 1
order = lambda th: abs(np.exp(1j * th).mean())
def winding(th): d = np.angle(np.exp(1j * (np.roll(th, -1) - th))); return int(round(d.sum() / (2 * np.pi)))
rc = [order(kuramoto(rng.uniform(0, 2 * np.pi, Nk), comp)) for _ in range(20)]
wr = [winding(kuramoto(rng.uniform(0, 2 * np.pi, Nk), ring, T=600.0)) for _ in range(30)]
tw = sum(1 for w_ in wr if w_ != 0)
print(f"  complete graph: 20 random starts, final r in [{min(rc):.6f}, {max(rc):.6f}]")
print(f"  ring: 30 random starts, twisted (winding != 0) finals: {tw}/30; winding histogram {dict(zip(*np.unique(wr, return_counts=True)))}")
RES["H3"] = [min(rc) > 0.9999, tw >= 5]

# ---------------- NR-H4 frozen stripes
hdr("NR-H4 zero-temperature Ising frame field: 2D lattice (random-sequential updates) vs complete graph")
def zeroT_sequential(s, sweeps, rng):
    L = s.shape[0]
    for _ in range(sweeps):
        idx = rng.integers(0, L, size=(L * L, 2))
        for i, j in idx:
            h = s[(i + 1) % L, j] + s[(i - 1) % L, j] + s[i, (j + 1) % L] + s[i, (j - 1) % L]
            if h != 0: s[i, j] = np.sign(h)
            elif rng.random() < 0.5: s[i, j] = -s[i, j]
    return s
L = 24; runs = 60; froz = 0; straight = 0; uniform = 0
for _ in range(runs):
    s = np.where(rng.random((L, L)) < 0.5, 1, -1)
    s = zeroT_sequential(s, 1500, rng)
    if abs(s.mean()) == 1: uniform += 1; continue
    s2 = zeroT_sequential(s.copy(), 500, rng)
    rows = (np.abs(s2.sum(1)) == L).all(); cols = (np.abs(s2.sum(0)) == L).all()
    if (s2 == s).all() or rows or cols: froz += 1
    if rows or cols: straight += 1
print(f"  2D L={L}, {runs} random starts, 1500 sweeps: uniform {uniform/runs:.2f}; non-uniform {1-uniform/runs:.2f}, of which straight stripes {straight}/{runs-uniform}")
Nm = 201; cfg = np.where(rng.random((100, Nm)) < 0.5, 1, -1)
for _ in range(40):
    for i in rng.permutation(Nm):
        h = cfg.sum(1) - cfg[:, i]; cfg[:, i] = np.where(h != 0, np.sign(h), cfg[:, i])
print(f"  complete graph N={Nm}: 100 starts, aligned finals {np.mean(np.abs(cfg.mean(1)) == 1):.2f}; '+' orientation {np.mean(cfg.mean(1) == 1):.2f} (gauge)")
RES["H4"] = [0.15 < 1 - uniform / runs < 0.7 and straight == runs - uniform, np.mean(np.abs(cfg.mean(1)) == 1) == 1.0]

# ---------------- NR-H5 invertibility / information preservation
hdr("NR-H5 exact invertibility / information preservation where claimed")
cr, sr, cc = Fr(5, 13), Fr(12, 13), Fr(7, 10); bb = (Fr(3, 2), Fr(-1, 3))
F_ = lambda x: (cc * (cr * x[0] - sr * x[1]) + bb[0], cc * (sr * x[0] + cr * x[1]) + bb[1])
Fi = lambda y: (lambda u, v: (cr * u + sr * v, -sr * u + cr * v))((y[0] - bb[0]) / cc, (y[1] - bb[1]) / cc)
x0 = (Fr(2), Fr(-7, 5)); y = x0
for _ in range(150): y = F_(y)
z = y
for _ in range(150): z = Fi(z)
print(f"  affine contraction, exact rationals, 150 steps forward then back: recovered x0 exactly: {z == x0}")
f = lambda x: -(x ** 3 + x - 0.7); xT = rk4(f, 2.5, 2.0, 1e-4); xb = rk4(lambda x: -f(x), xT, 2.0, 1e-4)
print(f"  gradient flow T = 2 forward/back (RK4 h = 1e-4): x0 = 2.5 recovered as {xb:.8f}")
th = 0.6; SW = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], complex); Uc = expm(-1j * th * SW)
na = 6; a = np.zeros(2 ** (na + 1), complex); a[0] = 1; b = np.zeros_like(a); b[2 ** na] = 1 / np.sqrt(2); b[0] = 1 / np.sqrt(2)
def coll(psi, k):
    t = psi.reshape([2] * (na + 1)); t = np.moveaxis(t, [0, k], [0, 1]).reshape(4, -1); t = (Uc @ t).reshape([2, 2] + [2] * (na - 1))
    return np.moveaxis(t, [0, 1], [0, k]).reshape(-1)
g = []
for k in range(1, na + 1): a, b = coll(a, k), coll(b, k); g.append(np.sqrt(max(0, 1 - abs(np.vdot(a, b)) ** 2)))
print(f"  unitary dilation (6 fresh |0> ancillas): global distinguishability per collision {np.round(g, 9)} (constant)")
P2 = np.array([[0.6, 0.3, 0.1], [0.2, 0.5, 0.3], [0.25, 0.25, 0.5]]); p0 = np.array([0.7, 0.2, 0.1])
rec = np.linalg.solve(np.linalg.matrix_power(P2, 3).T, p0 @ np.linalg.matrix_power(P2, 3))
print(f"  Markov distribution map P^3 (det P = {np.linalg.det(P2):.3f} != 0): exact inversion recovers p0 to {np.abs(rec - p0).max():.1e}")
RES["H5"] = [z == x0, abs(xb - 2.5) < 1e-6, np.ptp(g) < 1e-12, np.abs(rec - p0).max() < 1e-10]

hdr("SUMMARY")
claims = {"H1": "structural attractors: symmetric contraction -> 0; doubly stochastic -> uniform from every start; uniformity tuned",
          "H2": "explicit-target laws relocate the state value into D (attractor tracks b, a)",
          "H3": "all-to-all alignment synchronizes every random start; local ring keeps twisted basins",
          "H4": "2D zero-T Ising freezes a finite fraction in straight stripes; complete graph aligns every start (mod Z2)",
          "H5": "bijective / finite-time / unitary / invertible-Markov laws preserve fine-grained information"}
for k, v in RES.items():
    print(f"NR-{k}: {'REPRODUCED' if all(v) else 'PARTIAL' if any(v) else 'NOT REPRODUCED'} ({sum(map(bool, v))}/{len(v)}) — {claims[k]}")
