"""SCOUT-2 S2-3b: horocycle flow on X2 = SL(2,R)/SL(2,Z) (unimodular lattices in R^2).

Pre-registered firewall (probes/PROBE_CHARTERS.md, S2-3b): three notions measured SEPARATELY
  A. unique invariant measure
  B. time-average universality (single trajectories)
  C. preparation forgetting (ensembles, relative to an observable class)
and fine-grained information conservation under pullback.

Lattice L = span(b1, b2), det = 1.  Horocycle flow: h_t (a, b) = (a + t b, b).
Invariant coordinates after Gauss reduction: shape z = x + i y (y = 1/|v1|^2, |x| <= 1/2),
phi = arg(v1) mod pi.  Haar = (3/pi) dx dy / y^2 on the modular fundamental domain x uniform phi.
Analytic check: Haar P(|v1|^2 < s) = 3 s / pi for s <= 1.
"""
import numpy as np

rng = np.random.default_rng(20261001)
S = 0.5
HAAR_F = 3 * S / np.pi  # 0.477465


def reduce(b1, b2, iters=200):
    b1 = b1.copy(); b2 = b2.copy()
    for _ in range(iters):
        n1 = (b1 ** 2).sum(-1); n2 = (b2 ** 2).sum(-1)
        sw = n2 < n1
        if sw.any():
            t = b1[sw].copy(); b1[sw] = b2[sw]; b2[sw] = t
            n1 = (b1 ** 2).sum(-1)
        m = np.round((b1 * b2).sum(-1) / n1)
        if not (m != 0).any():
            n2 = (b2 ** 2).sum(-1)
            if not (n2 < n1).any():
                break
        b2 = b2 - m[:, None] * b1
    return b1, b2


def coords(b1, b2):
    v1, v2 = reduce(b1, b2)
    det = v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0]
    v2 = v2 * np.sign(det)[:, None]
    n1 = (v1 ** 2).sum(-1)
    y = 1.0 / n1
    x = (v1 * v2).sum(-1) / n1
    phi = np.mod(np.arctan2(v1[:, 1], v1[:, 0]), np.pi)
    return x, y, phi


def build(x, y, phi):
    c, s = np.cos(phi), np.sin(phi)
    r = 1 / np.sqrt(y)
    b1 = np.stack([r * c, r * s], -1)
    b2 = np.stack([r * (c * x - s * y), r * (s * x + c * y)], -1)
    return b1, b2


def haar(n):
    out = []
    while sum(len(o) for o in out) < n:
        x = rng.uniform(-.5, .5, 2 * n)
        y0 = np.sqrt(1 - x ** 2)
        keep = rng.uniform(size=2 * n) < (1 / y0) / (2 / np.sqrt(3))
        x, y0 = x[keep], y0[keep]
        y = y0 / rng.uniform(size=len(x))
        out.append(np.stack([x, y], -1))
    xy = np.concatenate(out)[:n]
    return build(xy[:, 0], xy[:, 1], rng.uniform(0, np.pi, n))


def flow(b1, b2, t):
    t = np.asarray(t, float)
    h = lambda v: np.stack([v[..., 0] + t * v[..., 1], v[..., 1]], -1)
    return h(b1), h(b2)


def f_obs(b1, b2):
    x, y, phi = coords(b1, b2)
    return (1 / y < S).astype(float), np.cos(2 * phi)


print("=" * 78)
print("S2-3b  horocycle flow on X2 = SL(2,R)/SL(2,Z)")
print("=" * 78)

# ---------- sanity: Haar sampler vs analytic ----------
N = 400_000
b1, b2 = haar(N)
F, G = f_obs(b1, b2)
print(f"[sanity] Haar P(|v1|^2<{S}) = {F.mean():.4f}  analytic 3s/pi = {HAAR_F:.4f};  <cos 2phi> = {G.mean():+.4f}")
F1, _ = f_obs(*flow(b1, b2, 7.3))
print(f"[sanity] h_t-invariance of Haar (t=7.3): P = {F1.mean():.4f}")

# ---------- A: unique measure? ----------
print("\n--- A. UNIQUE INVARIANT MEASURE? ---")
for a in [1.0, 0.5, 0.2]:
    p1 = np.array([[a, 0.0]]); p2 = np.array([[0.0, 1 / a]])
    ts = np.linspace(0, a * a, 9)[:-1]
    vals = [f_obs(*flow(p1, p2, t))[0][0] for t in ts]
    x, y, phi = coords(*flow(p1, p2, a * a))
    _, y0, _ = coords(p1, p2)
    print(f"  lattice diag({a},{1/a})Z^2: h_(a^2) returns (y {y0[0]:.4f} -> {y[0]:.4f}); period a^2 = {a*a:.3f}; "
          f"orbit average of 1[|v1|^2<0.5] = {np.mean(vals):.3f}  (Haar {HAAR_F:.3f})")
print("  => a one-parameter family of periodic orbits, each carrying its own invariant probability measure;")
print("     their averages (0, 1, ...) differ from Haar. A FAILS on non-compact X2 (Dani: Haar + periodic-orbit measures).")

# ---------- B: time-average universality ----------
print("\n--- B. TIME-AVERAGE UNIVERSALITY (Birkhoff averages of 1[|v1|^2<0.5], target 0.4775) ---")
def rot(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, s]]), np.array([[-s, c]])
hb1, hb2 = haar(1)
starts = {
    "Haar-random lattice": (hb1, hb2),
    "Z^2 rotated by sqrt(2)-1": rot(np.sqrt(2) - 1),
    "Z^2 rotated by 1e-2": rot(1e-2),
    "Z^2 rotated by 1e-3": rot(1e-3),
    "Z^2 (periodic, horizontal vector)": rot(0.0),
}
Ts = [1e2, 1e3, 1e4, 1e5]
dt = 0.02
print("  " + " " * 36 + "".join(f"T={T:<9.0e}" for T in Ts))
for name, (p1, p2) in starts.items():
    acc = 0.0; done = 0.0; row = []
    tgrid_all = np.arange(0, Ts[-1], dt) + dt / 2
    sums = np.zeros(len(Ts))
    for k0 in range(0, len(tgrid_all), 250_000):
        tg = tgrid_all[k0:k0 + 250_000]
        q1 = np.repeat(p1, len(tg), 0); q2 = np.repeat(p2, len(tg), 0)
        val = f_obs(*flow(q1, q2, tg[:, None][:, 0]))[0] if False else None
        a1 = np.stack([q1[:, 0] + tg * q1[:, 1], q1[:, 1]], -1)
        a2 = np.stack([q2[:, 0] + tg * q2[:, 1], q2[:, 1]], -1)
        val = f_obs(a1, a2)[0]
        for i, T in enumerate(Ts):
            sums[i] += val[tg < T].sum()
    row = [sums[i] / (T / dt) for i, T in enumerate(Ts)]
    print(f"  {name:<36}" + "".join(f"{r:<11.4f}" for r in row))
print("  => B holds for every non-periodic start tested (Dani-Smillie: every non-periodic orbit equidistributes),")
print("     with point-dependent, non-uniform convergence time (~1/angle near the periodic set); fails on periodic orbits.")

# ---------- C-i: decay of correlations (Haar) ----------
print("\n--- C-i. DECAY OF CORRELATIONS under Haar (N = 400k; MC noise ~ 1e-3) ---")
Fc = F - F.mean(); Gc = G - G.mean()
print("   t      C_ff(t)    C_gg(t)    C_fg(t)")
for t in [0, 0.5, 1, 2, 5, 10, 20, 50, 100, 200]:
    Ft, Gt = f_obs(*flow(b1, b2, t))
    print(f"  {t:5}  {np.mean(Fc*(Ft-Ft.mean())):+.4f}    {np.mean(Gc*(Gt-Gt.mean())):+.4f}    {np.mean(Fc*(Gt-Gt.mean())):+.4f}")

# ---------- C-ii: weak convergence of preparations ----------
print("\n--- C-ii. WEAK CONVERGENCE of accessible observable <1[|v1|^2<0.5]> under pushforward ---")
M = 200_000
def blob(xr, yr, pr, n=M):
    return build(rng.uniform(*xr, n), rng.uniform(*yr, n), rng.uniform(*pr, n))
B1 = ((0.0, 0.1), (1.4, 1.6), (0.3, 0.4))
B2 = ((-0.3, -0.2), (2.5, 3.0), (1.2, 1.3))
P1 = blob(*B1); P2 = blob(*B2)
dirac = rot(np.sqrt(2) - 1)
per = (np.array([[0.5, 0.0]]), np.array([[0.0, 2.0]]))
tlist = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500]
print("   t     blob1    blob2    Dirac(irrational)  periodic diag(.5,2)")
for t in tlist:
    a = f_obs(*flow(*P1, t))[0].mean(); b = f_obs(*flow(*P2, t))[0].mean()
    d = f_obs(*flow(*dirac, t))[0][0]; p = f_obs(*flow(*per, t))[0][0]
    print(f"  {t:4}   {a:.4f}   {b:.4f}   {d:.0f}                  {p:.0f}")
print(f"  Haar value {HAAR_F:.4f}")

# coarse partition TV distance between the two pushed ensembles
def cells(b1, b2, nb=10):
    x, y, phi = coords(b1, b2)
    u = np.clip(1 - 1 / y, 0, 1 - 1e-12)          # cusp coordinate in [0,1) (Haar-nonuniform, fine)
    ix = np.clip(((x + .5) * nb).astype(int), 0, nb - 1)
    iu = (u * nb).astype(int); ip = np.clip((phi / np.pi * nb).astype(int), 0, nb - 1)
    return (ix * nb + iu) * nb + ip
print("\n  coarse partition (10x10x10 cells in x, 1-1/y, phi): TV distance between pushed blob1, blob2")
print("   t     TV_coarse   [sampling floor for identical ensembles ~ %.3f]" % 0.0)
H0 = cells(*haar(M)); H0b = cells(*haar(M))
floor = 0.5 * np.abs(np.bincount(H0, minlength=1000) - np.bincount(H0b, minlength=1000)).sum() / M
for t in tlist:
    c1 = np.bincount(cells(*flow(*P1, t)), minlength=1000); c2 = np.bincount(cells(*flow(*P2, t)), minlength=1000)
    print(f"  {t:4}   {0.5*np.abs(c1-c2).sum()/M:.4f}")
print(f"  sampling floor (two independent Haar ensembles): {floor:.4f}")

# ---------- C-iii: fine-grained conservation under pullback ----------
print("\n--- C-iii. FINE-GRAINED INFORMATION under pullback: observable 1_{B1} o h_{-t} ---")
def inbox(b1, b2, box):
    x, y, phi = coords(b1, b2)
    (x0, x1), (y0, y1), (p0, p1) = box
    return (x >= x0) & (x <= x1) & (y >= y0) & (y <= y1) & (phi >= p0) & (phi <= p1)
print("   t      P(B1 | prep1)  P(B1 | prep2)  contrast   max round-trip |dy|")
for t in tlist:
    q1 = flow(*P1, t); q2 = flow(*P2, t)
    # the agent holds only the time-t lattice: reduce it, then apply the pullback h_{-t}
    r1 = reduce(*q1); r2 = reduce(*q2)
    back1 = flow(*r1, -t); back2 = flow(*r2, -t)
    m1 = inbox(*back1, B1).mean(); m2 = inbox(*back2, B1).mean()
    _, yb, _ = coords(*back1); _, y0, _ = coords(*P1)
    print(f"  {t:4}    {m1:.4f}         {m2:.4f}        {m1-m2:.4f}     {np.abs(yb-y0).max():.1e}")

print("\n  complexity of the pulled-back observable: # of 20^3 coarse cells occupied by h_t(B1) (N=200k)")
for t in [0, 1, 2, 5, 10, 20, 50, 100, 200, 500]:
    print(f"  t={t:4}: {len(np.unique(cells(*flow(*P1, t), nb=20))):5d} cells")
print("  contrast: doubling map T^n of an interval of width 2^-12, cells of width 2^-12: 2^n (exponential)")
