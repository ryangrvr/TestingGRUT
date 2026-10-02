"""SCOUT-2 S2-H2: state / basin selection by dynamics.  Pre-registered (probes/PROBE_CHARTERS.md, S2-H2).
Three notions kept separate:  A attractor uniqueness | B global reachability (EVERY admissible state) |
C information erasure (from the EXACT final microstate, including all degrees of freedom in the model)."""
from fractions import Fraction as Fr
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

rng = np.random.default_rng(42)
def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)

# ---------------------------------------------------------------- H2-1 contraction
hdr("H2-1  deterministic contraction on R^2:  F(x) = c R(phi) x + b,  c = 0.8")
c, phi = 0.8, 0.9
R = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
for bname, b in [("b = 0 (law is rotation-equivariant about 0)", np.zeros(2)), ("b = (1.3, -0.4) (explicit offset)", np.array([1.3, -0.4]))]:
    xstar = np.linalg.solve(np.eye(2) - c * R, b)
    starts = [np.array([1e6, -3e5]), np.array([0.0, 0.0]), np.array([-7.0, 2.0]), xstar + 1e-12]
    finals = []
    for x in starts:
        for _ in range(400):
            x = c * R @ x + b
        finals.append(x)
    print(f"  {bname}: fixed point x* = {np.round(xstar, 6)};  max |x_400 - x*| over starts incl. 1e6: "
          f"{max(np.linalg.norm(f - xstar) for f in finals):.1e}")
print("  A: unique fixed point (Banach).  B: EVERY start in R^2 converges; no measure, no exclusions.")
# C: exact invertibility (rational arithmetic, rational rotation via Pythagorean triple)
cr, sr = Fr(3, 5), Fr(4, 5); cc = Fr(4, 5); bb = (Fr(13, 10), Fr(-2, 5))
def Fx(x): return (cc * (cr * x[0] - sr * x[1]) + bb[0], cc * (sr * x[0] + cr * x[1]) + bb[1])
def Finv(y):
    u, v = (y[0] - bb[0]) / cc, (y[1] - bb[1]) / cc
    return (cr * u + sr * v, -sr * u + cr * v)
x0a, x0b = (Fr(1), Fr(0)), (Fr(1), Fr(1, 10 ** 6))
ya, yb = x0a, x0b
for _ in range(200):
    ya, yb = Fx(ya), Fx(yb)
ra = ya
for _ in range(200):
    ra = Finv(ra)
dist = float(abs(ya[0] - yb[0]) + abs(ya[1] - yb[1]))
print(f"  C (exact rational arithmetic, n = 200): |x_a - x_b|_1 at n=200 = {dist:.2e} (initial 1e-6); exact inversion"
      f" recovers x0_a: {ra == x0a}")
xf = np.array([1.0, 0.0]); bf = np.array([1.3, -0.4]); Rr = np.array([[0.6, -0.8], [0.8, 0.6]])
for n in [20, 60, 120, 200]:
    y = xf.copy()
    for _ in range(n): y = 0.8 * Rr @ y + bf
    for _ in range(n): y = np.linalg.solve(0.8 * Rr, y - bf)
    print(f"     float64 inversion after n = {n:3d}: recovery error {np.linalg.norm(y - xf):.1e}")
print("  => the affine contraction is a BIJECTION: fine-grained information is NOT erased; it is pushed below any")
print("     finite resolution (A-priced forgetting).  Genuine erasure needs a non-injective law:")
xa, xb = 0.37, -0.37
fa, fb = 0.8 * abs(xa) + 0.5, 0.8 * abs(xb) + 0.5
print(f"     F(x) = 0.8|x| + 0.5 (Lipschitz 0.8, non-injective): F(0.37) = {fa}, F(-0.37) = {fb}  -> merged EXACTLY in one step")

# ---------------------------------------------------------------- H2-1' dilation
hdr("H2-1'  the contraction realized unitarily: qubit system + fresh ancillas |0>, partial-swap collisions")
nA = 8; theta = 0.6
SW = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], complex)
Ucol = expm(-1j * theta * SW)
def collide(psi, k):
    # system = qubit 0, ancilla k (1..nA); total qubits nA+1
    n = nA + 1
    t = psi.reshape([2] * n)
    t = np.moveaxis(t, [0, k], [0, 1]).reshape(4, -1)
    t = (Ucol @ t).reshape([2, 2] + [2] * (n - 2))
    return np.moveaxis(t, [0, 1], [0, k]).reshape(-1)
def sys_rho(psi):
    m = psi.reshape(2, -1); return m @ m.conj().T
def td(r1, r2): return 0.5 * np.abs(np.linalg.eigvalsh(r1 - r2)).sum()
anc = np.zeros(2 ** nA); anc[0] = 1
sA = np.kron(np.array([0, 1], complex), anc); sB = np.kron(np.array([1, 1j]) / np.sqrt(2), anc)
print("   k   TD(system)   TD(global = sqrt(1-|<a|b>|^2))   <Z_sys> for both")
for k in range(nA + 1):
    if k:
        sA, sB = collide(sA, k), collide(sB, k)
    g = np.sqrt(max(0.0, 1 - abs(np.vdot(sA, sB)) ** 2))
    zA = np.real(np.trace(sys_rho(sA) @ np.diag([1, -1]))); zB = np.real(np.trace(sys_rho(sB) @ np.diag([1, -1])))
    print(f"  {k:2d}   {td(sys_rho(sA), sys_rho(sB)):.4f}       {g:.6f}                        {zA:+.3f} {zB:+.3f}")
print("  => system converges toward the ANCILLA state |0>; global distinguishability is conserved exactly;")
print("     the attractor value is the bath's PREPARATION (H relocated into the environment's state).")

# ---------------------------------------------------------------- H2-2 gradient flows
hdr("H2-2  gradient flows  dx/dt = -V'(x)")
def flow(Vp, x0, T):
    s = solve_ivp(lambda t, x: -Vp(x), (0, T), [x0], rtol=1e-12, atol=1e-14)
    return s.y[0, -1]
cases = {
    "a  unique min, V = x^4/4 + x^2/2 - 0.7x": lambda x: x ** 3 + x - 0.7,
    "a0 unique min, symmetric (a = 0)": lambda x: x ** 3 + x,
    "b  double well + tilt eps = 0.05": lambda x: x ** 3 - x + 0.05,
    "c  symmetric double well": lambda x: x ** 3 - x,
}
starts = [-50.0, -1.0, -0.3, -1e-9, 0.0, 1e-9, 0.3, 1.0, 50.0]
for name, Vp in cases.items():
    fin = [flow(Vp, x0, 60.0) for x0 in starts]
    print(f"  {name}:")
    print("     x0 :   " + " ".join(f"{x:>8.1e}" for x in starts))
    print("     x(60): " + " ".join(f"{x:>8.4f}" for x in fin))
print("  a / a0: every start reaches the unique minimum (strictly convex V): A + B hold for EVERY x0.")
print("  b: two stable minima persist under the tilt; the unstable stationary point never leaves -> BASIN DATA SURVIVES.")
print("  c: outcomes +-1 are gauge-related by x -> -x, BUT x0 = 0 stays at 0 forever (an admissible stationary state).")
# C for gradient flow: finite-time flow is a diffeomorphism
Vp = cases["a  unique min, V = x^4/4 + x^2/2 - 0.7x"]
for T in [2.0, 8.0, 20.0]:
    xT = flow(Vp, 3.0, T)
    s = solve_ivp(lambda t, x: Vp(x), (0, T), [xT], rtol=1e-13, atol=1e-15)
    print(f"  C: forward then backward integration, T = {T:4.1f}: recovered x0 = {s.y[0,-1]:.6f} (true 3.0)")
print("  => finite-time gradient flow is invertible (ODE uniqueness); convergence is asymptotic; information is lost")
print("     only below finite resolution / as T -> infinity (A-priced), not from the exact microstate at finite T.")

# ---------------------------------------------------------------- H2-3 Markov / channel
hdr("H2-3  primitive Markov chains and a primitive qubit channel")
def primitive_ds(n):
    # convex combination of permutation matrices (Birkhoff) -> doubly stochastic
    P = 0.35 * np.eye(n)
    for _ in range(4):
        P += 0.1625 * np.eye(n)[rng.permutation(n)]
    return P
Pd = primitive_ds(5)
Pg = rng.uniform(0.05, 1, (5, 5)); Pg /= Pg.sum(1, keepdims=True)
for name, P in [("doubly stochastic (Birkhoff mixture)", Pd), ("generic row-stochastic", Pg)]:
    prim = (np.linalg.matrix_power(P, 10) > 0).all()
    w, V = np.linalg.eig(P.T); pi = np.real(V[:, np.argmin(abs(w - 1))]); pi /= pi.sum()
    print(f"  {name}: primitive {prim};  stationary pi = {np.round(pi, 4)};  |lambda_2| = {sorted(abs(w))[-2]:.3f};  det P = {np.linalg.det(P):.2e}")
    for p0 in [np.eye(5)[0], np.eye(5)[4], np.ones(5) / 5]:
        pt = p0 @ np.linalg.matrix_power(P, 60)
        print(f"     from {np.round(p0, 2)}:  p_60 = {np.round(pt, 4)}")
    for t in [3, 10, 30]:
        Pt = np.linalg.matrix_power(P, t); pt = np.eye(5)[0] @ Pt
        rec = np.linalg.lstsq(Pt.T, pt, rcond=None)[0]
        print(f"     C: invert P^{t:<2d} on the float64 distribution: recovery error {np.abs(rec - np.eye(5)[0]).max():.1e}  (cond {np.linalg.cond(Pt):.1e})")
print("  => unique stationary distribution from EVERY initial distribution (A + B).  At distribution level P^t is")
print("     injective here (det != 0): the exact distribution still encodes p_0; forgetting is resolution-priced.")
print("     Stochasticity (probability) is supplied in the transition law: PROBABILITY D-PRICED.")
Pp = Pd + 0.02 * rng.normal(size=(5, 5)); Pp = np.abs(Pp); Pp /= Pp.sum(1, keepdims=True)
w, V = np.linalg.eig(Pp.T); pi = np.real(V[:, np.argmin(abs(w - 1))]); pi /= pi.sum()
print(f"  H2-6 robustness: perturb the doubly stochastic chain (row-renormalized noise 0.02): pi = {np.round(pi, 4)}")
print("     primitivity (uniqueness) survives (open condition); UNIFORMITY of pi does not (tuned constraint).")
# qubit channel
g_ad, p_dp = 0.3, 0.2
K = [np.array([[1, 0], [0, np.sqrt(1 - g_ad)]]), np.array([[0, np.sqrt(g_ad)], [0, 0]])]
def chan(r):
    r = sum(k @ r @ k.conj().T for k in K)
    return (1 - p_dp) * r + p_dp * np.diag(np.diag(r))
ra = np.array([[0, 0], [0, 1]], complex); rb = 0.5 * np.array([[1, 1j], [-1j, 1]])
for t in range(0, 41, 10):
    if t:
        for _ in range(10): ra, rb = chan(ra), chan(rb)
    print(f"  amplitude-damping+dephasing qubit channel, step {t:2d}: TD = {td(ra, rb):.4f};  rho_a[0,0] = {ra[0,0].real:.4f}")
print("  => unique fixed point |0><0| = the zero-temperature bath state (Stinespring: same as H2-1').")

# ---------------------------------------------------------------- H2-4 unitary hostile
hdr("H2-4  unitary hostile: random Hamiltonian, 6 qubits")
n = 6; N = 2 ** n
A = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N)); H = (A + A.conj().T) / 2
a = rng.normal(size=N) + 1j * rng.normal(size=N); a /= np.linalg.norm(a)
b = rng.normal(size=N) + 1j * rng.normal(size=N); b /= np.linalg.norm(b)
def red(psi, q=0):
    t = np.moveaxis(psi.reshape([2] * n), q, 0).reshape(2, -1); return t @ t.conj().T
ev = np.linalg.eigvals(expm(-1j * H * 0.3))
print(f"  eigenvalues of U(0.3): |lambda| in [{abs(ev).min():.12f}, {abs(ev).max():.12f}]  (no attracting fixed point)")
for t in [0, 1, 5, 20, 100]:
    U = expm(-1j * H * t); at, bt = U @ a, U @ b
    print(f"   t={t:4d}: global TD = {np.sqrt(1 - abs(np.vdot(at, bt))**2):.12f};  qubit-0 TD = {td(red(at), red(bt)):.4f}")
print("  => NO FINE-GRAINED H SELECTION: exact distinguishability conserved; apparent equilibration only in reduced /")
print("     coarse observables (cf. S2-3b).  Scope: unitary / Hamiltonian / isometric dynamics.")

# ---------------------------------------------------------------- H2-7 alignment / frame field
hdr("H2-7  shared reference frames: dissipative alignment (U(1)) and Ising frame field (Z2: J <-> -J)")
def kuramoto(theta0, adj, K=1.0, omega=None, T=200.0):
    N_ = len(theta0); deg = adj.sum(1)
    om = np.zeros(N_) if omega is None else omega
    f = lambda t, th: om + K * (adj * np.sin(th[None, :] - th[:, None])).sum(1) / deg
    s = solve_ivp(f, (0, T), theta0, rtol=1e-8, atol=1e-10)
    return s.y[:, -1]
def order(th): return abs(np.exp(1j * th).mean())
def winding(th):
    d = np.angle(np.exp(1j * (np.roll(th, -1) - th))); return int(round(d.sum() / (2 * np.pi)))
Nk = 40
comp = np.ones((Nk, Nk)) - np.eye(Nk)
ring = np.zeros((Nk, Nk))
for i in range(Nk): ring[i, (i + 1) % Nk] = ring[i, (i - 1) % Nk] = 1
rs = [order(kuramoto(rng.uniform(0, 2 * np.pi, Nk), comp)) for _ in range(20)]
print(f"  complete graph, identical frames, 20 random starts: final order r in [{min(rs):.6f}, {max(rs):.6f}]")
anti = np.r_[np.zeros(Nk // 2 + 3), np.full(Nk // 2 - 3, np.pi)]
print(f"  complete graph, admissible stationary start (two antipodal clusters 23/17): r = {order(anti):.3f} -> "
      f"{order(kuramoto(anti, comp)):.3f};  with 1e-6 kick -> {order(kuramoto(anti + 1e-6 * rng.normal(size=Nk), comp)):.3f}")
wins = []
for _ in range(30):
    th = kuramoto(rng.uniform(0, 2 * np.pi, Nk), ring, T=600.0); wins.append((winding(th), order(th)))
wc = {}
for w_, r_ in wins: wc[w_] = wc.get(w_, 0) + 1
print(f"  ring (nearest-neighbour), 30 random starts: final winding numbers {dict(sorted(wc.items()))};  "
      f"r of non-synced finals: {sorted({round(r_, 3) for w_, r_ in wins if w_ != 0})[:5]}")
print("     twisted states (winding q != 0) are stable attractors -> BASIN DATA SURVIVES on the ring.")
for sig in [0.0, 0.2, 0.5, 1.0]:
    om = rng.normal(0, sig, Nk)
    print(f"  H2-6 complete graph with frequency disorder sigma = {sig}: r = {order(kuramoto(rng.uniform(0, 2*np.pi, Nk), comp, omega=om, T=300)):.3f}")
print("     disorder makes order a REGIME of D (K > K_c ~ 2 sqrt(2/pi) sigma), not structural.")

# Ising frame field
def glauber(s, beta, sweeps, rngl):
    L = s.shape[-1]
    ii, jj = np.indices((L, L)); masks = [((ii + jj) % 2 == p) for p in (0, 1)]
    for _ in range(sweeps):
        for m in masks:
            h = np.roll(s, 1, -1) + np.roll(s, -1, -1) + np.roll(s, 1, -2) + np.roll(s, -1, -2)
            if np.isinf(beta):
                new = np.sign(h); tie = h == 0
                new[tie] = np.where(rngl.random(new[tie].shape) < 0.5, 1, -1)
            else:
                p = 1 / (1 + np.exp(-2 * beta * h))
                new = np.where(rngl.random(s.shape) < p, 1, -1)
            s = np.where(m, new, s)
    return s
L = 32; S = 120
s0 = np.where(rng.random((S, L, L)) < 0.5, 1, -1)
sT0 = glauber(s0.copy(), np.inf, 3000, rng)
mag = np.abs(sT0.mean((1, 2)))
print(f"  2D Ising frame field, ZERO-temperature Glauber (pure energy descent), L = {L}, {S} random starts, 3000 sweeps:")
print(f"     reached a uniform aligned state (|m| = 1): {np.mean(mag == 1):.2f};  frozen striped states: {np.mean(mag < 1):.2f}")
for beta in [1 / 1.8, 1 / 3.2]:
    sb = glauber(np.where(rng.random((60, L, L)) < 0.5, 1, -1), beta, 3000, rng)
    corr = np.mean(sb[:, :, 0] * sb[:, :, L // 2])
    print(f"  finite-T Glauber T = {1/beta:.1f}: <|m|> = {np.abs(sb.mean((1,2))).mean():.3f};  frame correlation 2<f_A f_B> at distance {L//2}: {2*corr:.3f}")
Nmf = 201
cfg = np.where(rng.random((200, Nmf)) < 0.5, 1, -1)
for _ in range(50):
    for i in rng.permutation(Nmf):
        h = cfg.sum(1) - cfg[:, i]
        cfg[:, i] = np.where(h != 0, np.sign(h), cfg[:, i])
print(f"  complete-graph Ising, zero-T, N = {Nmf} (odd), 200 starts: uniform aligned final states {np.mean(np.abs(cfg.mean(1)) == 1):.2f};"
      f" orientation + : {np.mean(cfg.mean(1) == 1):.2f}")

hdr("H2-7 checks: exact stationarity of the antipodal start; frozenness of zero-T striped states")
f_anti = (comp * np.sin(anti[None, :] - anti[:, None])).sum(1) / comp.sum(1)
J = comp * np.cos(anti[None, :] - anti[:, None]); J = J - np.diag(J.sum(1)); J /= (Nk - 1)
print(f"  antipodal 23/17 config: max |d theta/dt| = {np.abs(f_anti).max():.1e} (float sin(pi) = {np.sin(np.pi):.1e}; exactly 0 in exact arithmetic);"
      f" largest Jacobian eigenvalue = {np.linalg.eigvalsh((J + J.T) / 2).max():+.3f} (unstable)")
print("  => the earlier 'r -> 1 without kick' is float roundoff acting as the kick. In exact dynamics this admissible")
print("     state never synchronizes. The set of such non-synchronized equilibria is closed, nowhere dense and Lebesgue-null.")
stuck = sT0[mag < 1]
later = glauber(stuck.copy(), np.inf, 3000, rng)
mag2 = np.abs(later.mean((1, 2)))
def straight(c):
    rows = (np.abs(c.sum(1)) == c.shape[1]).all(); cols = (np.abs(c.sum(0)) == c.shape[0]).all()
    return rows or cols
print(f"  non-uniform zero-T states after 3000 more sweeps: {len(stuck)} -> still non-uniform {np.sum(mag2 < 1)};"
      f" straight horizontal/vertical stripes: {sum(straight(c) for c in later[mag2 < 1])}")
