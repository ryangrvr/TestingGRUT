"""A3: classes (I-z, I-q) and v_b for the registered families, fermions (SF-1) and free bosons.
P  : fixed q, L->inf along the family sequence (exact finite-L lower edge), then q->0.
P' : omega_1(L) = lower edge at q=2pi/L, exponent from L->inf.
Closed forms derived in V0_3_P02_REPRODUCTION.md are checked against brute-force finite-L values."""
import numpy as np
from v0_3_lib import *

eps = eps_from_hops({1: -1.0, -1: -1.0})


def closed_lower_edge(nu, q):
    """thermodynamic SF-1 lower edge, derived: 4 sin(q/2) |sin(pi*min(nu,1-nu) - q/2)|"""
    kb = np.pi * min(nu, 1 - nu)
    return 4 * np.sin(q / 2) * np.abs(np.sin(kb - q / 2))


def omega1_closed(L, N):
    return 4 * np.sin(np.pi / L) * np.sin(np.pi * N / L)


def nearest_odd_sqrt(L):
    s = np.sqrt(L)
    lo = int(np.floor(s))
    cands = [n for n in (lo - 1, lo, lo + 1, lo + 2) if n % 2 == 1 and n > 0]
    d = [abs(n - s) for n in cands]
    return cands[int(np.argmin(d))]   # ties -> lower (only at even perfect squares)


families = {
    'D':    (lambda L: L // 2,     lambda j: 4 * j + 2,  0.5),
    'D-1/4': (lambda L: L // 4,    lambda j: 8 * j + 4,  0.25),
    'D-3/4': (lambda L: 3 * L // 4, lambda j: 8 * j + 4, 0.75),
    'E':    (lambda L: 1,          lambda j: 2 * j,      0.0),
    'E-3':  (lambda L: 3,          lambda j: 2 * j,      0.0),
    'Ebar': (lambda L: L - 1,      lambda j: 2 * j,      1.0),
    'C-6':  (nearest_odd_sqrt,     lambda j: 2 * j,      0.0),
}

qtest = [np.pi / 8, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]
print("=" * 100)
print("FERMIONS (SF-1).  Exact finite-L lower edge at q_L (nearest grid, ties->lower m) vs derived closed form")
summary = []
for name, (Nrule, Lseq, nu) in families.items():
    print(f"\n--- family {name}: nu_inf={nu}")
    # convergence of P at fixed q
    Ls = [Lseq(j) for j in (25, 250, 2500, 25000, 250000)]
    for q in qtest:
        vals = []
        for L in Ls:
            N = Nrule(L)
            assert N % 2 == 1, (name, L, N)
            occ, uniq, gap = fermion_occupation(eps, L, N)
            assert uniq
            m = nearest_m(L, q)
            vals.append(ff_lower_edge(eps, L, occ, m))
        print(f"  q={q:.6f}: w-_L(q_L) for L={Ls}: " + ", ".join(f"{v:.8f}" for v in vals)
              + f"   closed form w-(q)={closed_lower_edge(nu, q):.8f}")
    # P exponent: from the limit function
    f = lambda q: closed_lower_edge(nu, q)
    le = local_exponent(f, [1e-2, 1e-4, 1e-6, 1e-8])
    for (q, w, s, r) in le:
        print(f"  P: q={q:.0e}  w-={w:.6e}  dlog w/dlog q={s:.6f}  log w/log q={r:.6f}")
    zP = round(le[-1][2], 3)
    # soft set from closed form zeros on (0, pi]
    qs = np.linspace(1e-9, np.pi, 200001)
    w = f(qs)
    zs = [0.0] + [qs[t] for t in range(1, len(qs) - 1) if w[t] <= w[t - 1] and w[t] <= w[t + 1] and w[t] < 1e-4]
    if w[-1] < 1e-4 and w[-1] <= w[-2]:
        zs.append(qs[-1])
    zs = sorted(set(round(z, 4) for z in zs))
    # P'
    Lp = [Lseq(j) for j in (250, 2500, 25000, 250000, 1000000)]
    w1 = []
    for L in Lp:
        N = Nrule(L)
        occ, uniq, gap = fermion_occupation(eps, L, N)
        wb = ff_lower_edge(eps, L, occ, 1)
        assert abs(wb - omega1_closed(L, N)) < 1e-12 * max(1, abs(wb)) + 1e-15, (wb, omega1_closed(L, N))
        w1.append(wb)
    sl = [-np.log(w1[t + 1] / w1[t]) / np.log(Lp[t + 1] / Lp[t]) for t in range(len(Lp) - 1)]
    rat = [-np.log(w) / np.log(L) for w, L in zip(w1, Lp)]
    print(f"  P': L={Lp}  w1(L) brute force == 4 sin(pi/L) sin(pi N/L) to 1e-12")
    print(f"      local -dlog w1/dlog L: " + ", ".join(f"{s:.5f}" for s in sl)
          + "   -log w1/log L: " + ", ".join(f"{r:.4f}" for r in rat))
    vb = abs(2 * np.sin(np.pi * nu))
    zPp = round(sl[-1], 3)
    summary.append((name, nu, vb, zP, zPp, zs, len(zs)))

print("\n" + "=" * 100)
print("C-6 under P: limit function is 4 sin^2(q/2) (k_F(L)->0), check at large L:")
for q in (np.pi / 8, np.pi / 2, np.pi):
    row = []
    for L in (10**3, 10**4, 10**5, 10**6, 4 * 10**6):
        N = nearest_odd_sqrt(L)
        occ, u, g = fermion_occupation(eps, L, N)
        row.append(ff_lower_edge(eps, L, occ, nearest_m(L, q)))
    print(f"  q={q:.5f}: " + ", ".join(f"{v:.7f}" for v in row) + f"  -> 4sin^2(q/2)={4*np.sin(q/2)**2:.7f}")
print("C-6 two-scale analysis. Exact finite-L lower edge (closed form, verified below):")
print("  w-_L(m) = 4 sin(pi m/L) sin(pi(N+1-m)/L) for m<=N,  4 sin(pi m/L) sin(pi(m-N+1)/L) for m>N")
print("  small-q scaling form: w- ~= q |q - 2k_F(L)|,  k_F(L) ~ pi L^{-1/2}: two scales q and k_F(L).")
print("  Path q_L = 2pi m/L with m = max(1, round(c L^{1-a})); successive slopes dlog w/dlog q_L over L=10^4..10^12:")
def wL(L, N, m):
    return 4*np.sin(np.pi*m/L)*(np.sin(np.pi*(N+1-m)/L) if m <= N else np.sin(np.pi*(m-N+1)/L))
for a, c in ((0.25, 0.05), (0.5, 0.05), (0.5, 3.0), (0.5, 1.0), (0.6, 0.05), (0.75, 1.0), (0.9, 1.0), (1.0, 1.0)):
    pts = []
    for e in range(4, 13, 2):
        L = 2*int(10**e/2)
        N = nearest_odd_sqrt(L)
        m = max(1, int(round(c*L**(1-a))))
        pts.append((2*np.pi*m/L, wL(L, N, m), m, N))
    sl = [np.log(pts[t+1][1]/pts[t][1])/np.log(pts[t+1][0]/pts[t][0]) for t in range(len(pts)-1)]
    pred = 2.0 if a <= 0.5 else 1 + 1/(2*a)
    print(f"  a={a:4.2f} c={c:4.2f}: slopes " + ", ".join(f"{s:.4f}" for s in sl)
          + f"   (m/N at largest L = {pts[-1][2]/pts[-1][3]:.3g})  predicted {pred:.4f}")
print("  a=1/2 with c ~ 1 (m ~ N, i.e. q ~ 2k_F(L)) rides the moving soft line; exponent then not 2 (see c=1.0 row).")
# verify the finite-L closed form used above against brute force
bad = 0
for L in (1000, 4096, 10000):
    N = nearest_odd_sqrt(L)
    occ, u, g = fermion_occupation(eps, L, N)
    for m in range(1, L // 2 + 1):
        bf = ff_lower_edge(eps, L, occ, m)
        cf = 4 * np.sin(np.pi * m / L) * (np.sin(np.pi * (N + 1 - m) / L) if m <= N else np.sin(np.pi * (m - N + 1) / L))
        if abs(bf - cf) > 1e-12:
            bad += 1
print(f"  finite-L closed form vs brute force (C-6, L=1000,4096,10000, all m<=L/2): mismatches={bad}")
print("General N ~ L^alpha: w1 = 4 sin(pi/L) sin(pi N/L) ~ 4 pi^2 L^{alpha-2} -> z_P' = 2-alpha (alpha<1); "
      "C-6 (alpha=1/2) z_P'=3/2 while z_P=2.")

print("\n" + "=" * 100)
print("FREE BOSONS: GS = condensate in k=0 for every N (unique min of -2cos k); exact finite-L lower edge")
for name, (Nrule, Lseq, nu) in families.items():
    vals = []
    for L in (Lseq(25), Lseq(2500), Lseq(250000)):
        N = Nrule(L)
        nocc = np.zeros(L, int); nocc[0] = N
        m = nearest_m(L, np.pi / 2)
        sup = fb_support(eps, L, nocc, m)
        vals.append(min(e for e, W in sup))
    w1 = []
    Lp = [Lseq(2500), Lseq(250000)]
    for L in Lp:
        N = Nrule(L)
        nocc = np.zeros(L, int); nocc[0] = N
        w1.append(min(e for e, W in fb_support(eps, L, nocc, 1)))
    zpp = -np.log(w1[1] / w1[0]) / np.log(Lp[1] / Lp[0])
    print(f"  {name:6s}: w-_L(pi/2) = " + ", ".join(f"{v:.8f}" for v in vals) + f" (->2)   z_P' = {zpp:.5f}   "
          f"w-(q)=4sin^2(q/2) => z_P=2, Q_soft={{0}}, class (2,1)")

print("\n" + "=" * 100)
print("SUMMARY (fermions, SF-1):  family | nu | v_b=|eps'(pi nu)| | z_P | z_P' | Q_soft reps | I-q")
for s in summary:
    print(f"  {s[0]:6s} | {s[1]:.2f} | {s[2]:.6f} | {s[3]} | {s[4]} | {s[5]} | {s[6]}")
print("NB: for C-6 the P column is from the limit function 4sin^2(q/2) (z_P=2, I-q=1); its v_b row uses nu_inf=0.")
