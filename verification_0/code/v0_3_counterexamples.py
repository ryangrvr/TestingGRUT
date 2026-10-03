"""Counterexamples / hostile checks for the generality claims of the occupation-edge theorem (TARGET 1 & 3, audit 1-8).
Every model is a free, number-conserving 1D ring H = sum_j sum_r h_r a_j^dag a_{j+r}; eps(k) = sum_r h_r e^{ikr};
rho_q = sum_j e^{-iqj} n_j = sum_k a_{k-q}^dag a_k (q>0: k -> k-q).  Small-L exact many-body ED cross-checks the
free-particle support formula; large-L / continuum evaluations give the P and P' exponents and the soft sets."""
import numpy as np
from math import factorial
from v0_3_lib import *

np.set_printoptions(linewidth=160)


def ed_vs_formula_fermion(hops, L, N, xs=None, nsites=None, bonds=None, eps=None, occ=None, label=''):
    nsites = nsites or L
    bonds = bonds or ring_bonds(L, hops)
    xs = list(range(L)) if xs is None else xs
    B, idx, H = hamiltonian(nsites, N, 'fermion', bonds)
    ev, V = np.linalg.eigh(H)
    deg = int(np.sum(ev - ev[0] < 1e-9))
    ok = True
    for m in range(L):
        s_ed = ed_support(ev, V, V[:, 0], ev[0], rho_diag(B, xs, 2 * np.pi * m / L), N)
        s_ff = ff_support(eps, L, occ, m) if eps is not None else None
        if s_ff is not None and not same_support(s_ed, s_ff):
            ok = False
    return deg, ev[0], ok


def thermo_report(eps, O, name, qs=(1e-2, 1e-3, 1e-4, 1e-5)):
    f = lambda q: thermo_lower_edge(eps, O, q)
    le = local_exponent(f, qs)
    print(f"  [{name}] P: " + "; ".join(f"q={q:.0e}: w-={w:.3e}, dlogw/dlogq={s:.4f}" for q, w, s, r in le))
    soft, _, _ = soft_set(eps, O, nq=1501)
    print(f"  [{name}] Q_soft reps in [0,pi] = {[round(float(s), 6) for s in soft]}  ->  I-q = {len(soft)}")
    return le[-1][2], soft


def Pprime(eps, Lseq, Nrule, name):
    w1, Ls = [], []
    for L in Lseq:
        N = Nrule(L)
        occ, uniq, gap = fermion_occupation(eps, L, N)
        assert uniq, (name, L, N)
        w1.append(ff_lower_edge(eps, L, occ, 1)); Ls.append(L)
    sl = [-np.log(w1[t + 1] / w1[t]) / np.log(Ls[t + 1] / Ls[t]) for t in range(len(Ls) - 1)]
    print(f"  [{name}] P': L={Ls}  successive -dlog w1/dlog L = " + ", ".join(f"{s:.4f}" for s in sl))
    return sl[-1]


print("=" * 110)
print("CE-1  flat (quartic) minimum: eps(k) = -4cos k + cos 2k  =>  eps(k)-eps(0) = 8 sin^4(k/2)")
hq = {1: -2.0, -1: -2.0, 2: 0.5, -2: 0.5}
eq = eps_from_hops(hq)
kk = np.linspace(-np.pi, np.pi, 7)
print("  identity check max|eps(k)-eps(0)-8sin^4(k/2)| =", np.max(np.abs(eq(kk) - eq(0) - 8 * np.sin(kk / 2) ** 4)))
for L, N in ((6, 3), (10, 3), (10, 1)):
    B, idx, H = hamiltonian(L, N, 'boson', ring_bonds(L, hq))
    ev, V = np.linalg.eigh(H)
    deg = int(np.sum(ev - ev[0] < 1e-9))
    lines = []
    ok = True
    for m in range(1, L):
        s = ed_support(ev, V, V[:, 0], ev[0], rho_diag(B, range(L), 2 * np.pi * m / L), N)
        ok &= len(s) == 1 and abs(s[0][0] - 8 * np.sin(np.pi * m / L) ** 4) < 1e-9 and abs(s[0][1] - 1) < 1e-9
    print(f"  bosons L={L} N={N}: GS deg={deg}  E0={ev[0]:+.10f} (N eps(0)={N*eq(0):+.1f})  support = single line "
          f"8 sin^4(q/2) with weight 1 for all m != 0: {ok}")
f = lambda q: 8 * np.sin(q / 2) ** 4
print("  bosons P: w-(q) = 8 sin^4(q/2): " + "; ".join(f"q={q:.0e}: dlogw/dlogq={s:.4f}" for q, w, s, r in local_exponent(f, (1e-2, 1e-4))))
print("  bosons P': w1 = 8 sin^4(pi/L): exponent 4.  => free-boson class (4,1), NOT (2,1).  Fermion E (N=1) identical.")
zp = Pprime(eq, [2 * j for j in (25, 100, 400)], lambda L: 1, "fermion E, quartic band")

print("=" * 110)
print("CE-2  degenerate minima: eps(k) = 2cos k + cos 2k, minima at k0 = +-2pi/3 (eps=-3/2); 3 | L")
hd = {1: 1.0, -1: 1.0, 2: 0.5, -2: 0.5}
ed = eps_from_hops(hd)
for L, N in ((6, 2), (6, 3), (12, 3)):
    B, idx, H = hamiltonian(L, N, 'boson', ring_bonds(L, hd))
    ev, V = np.linalg.eigh(H)
    deg = int(np.sum(ev - ev[0] < 1e-9))
    print(f"  bosons L={L} N={N}: E0={ev[0]:+.10f} (N*(-1.5)={-1.5*N})  GS degeneracy = {deg}  (expected N+1 = {N+1})")
    # reference states inside the GS manifold: |n+, n-> Fock states and a superposition
    mp, mm = L // 3, (2 * L) // 3          # k0 = 2pi/3 (index L/3), -k0 = 4pi/3 (index 2L/3)

    def fock(npl, nmi):
        vec, Bprev = np.array([1.0 + 0j]), [tuple([0] * L)]
        cnt = 0
        for (mk, nk) in ((mp, npl), (mm, nmi)):
            for _ in range(nk):
                cnt += 1
                Bn = basis(L, cnt, 'boson'); idn = {b: t for t, b in enumerate(Bn)}
                out = np.zeros(len(Bn), complex)
                for t, b in enumerate(Bprev):
                    if vec[t] == 0: continue
                    for j in range(L):
                        n = list(b); a = np.sqrt(n[j] + 1); n[j] += 1
                        out[idn[tuple(n)]] += vec[t] * a * np.exp(1j * 2 * np.pi * mk / L * j) / np.sqrt(L)
                vec, Bprev = out, Bn
        return vec / np.linalg.norm(vec)
    refs = {f"all in +k0 (n+={N})": fock(N, 0), f"all in -k0 (n-={N})": fock(0, N),
            f"n+={N//2+N%2}, n-={N//2}": fock(N // 2 + N % 2, N // 2)}
    if N == 2:
        refs["(|2,0> - |0,2>)/sqrt2"] = (fock(2, 0) - fock(0, 2)) / np.sqrt(2)
    for nm, psi in refs.items():
        E = np.real(np.vdot(psi, H @ psi))
        row = []
        for m in range(1, L // 2 + 1):
            s = ed_support(ev, V, psi, E, rho_diag(B, range(L), 2 * np.pi * m / L), N)
            w0 = sum(W for e, W in s if abs(e) < 1e-9)
            row.append(f"m={m}:min w={min(e for e, W in s):+.4f},W(w=0)={w0:.3f}")
        print(f"    ref {nm:24s} E-E0={E-ev[0]:.1e}: " + " | ".join(row))
print("  => q* = 2pi/3 (= 2k0 mod quotient) is soft or not depending on WHICH ground state is declared;")
print("     the 'unique condensate' hypothesis of TARGET 1 fails and (I-z,I-q) = (2,1) or (2,2) is state-dependent.")
occ, uniq, gap = fermion_occupation(ed, 12, 1)
print(f"  fermions, same band, N=1, L=12: unique sector GS? {uniq} (gap {gap:.1e}) -> family E ill-defined here.")

print("=" * 110)
print("CE-3  symmetric band with a stationary inflection at the D-type edge: eps(k) = -cos^3(k)/3")
hc = {1: -1 / 8, -1: -1 / 8, 3: -1 / 24, -3: -1 / 24}
ec = eps_from_hops(hc)
print("  identity check max|eps+cos^3/3| =", np.max(np.abs(ec(kk) + np.cos(kk) ** 3 / 3)),
      " eps'(pi/2) = sin cos^2 = 0, eps''(pi/2)=0, eps'''(pi/2) = 2")
for L in (10, 14):
    N = L // 2
    occ, uniq, gap = fermion_occupation(ec, L, N)
    deg, E0, ok = ed_vs_formula_fermion(hc, L, N, eps=ec, occ=occ)
    print(f"  fermion ED L={L} N={N}: GS deg={deg}, unique by levels={uniq}, ED support == form-factor support: {ok}")
O = [(-np.pi / 2, np.pi / 2)]
z3, s3 = thermo_report(ec, O, "fermion D, cubic-inflection band (v_b=0)", qs=(1e-1, 1e-2, 1e-3, 3e-4))
print("  predicted: w-(q) ~ (eps'''/24) q^3 = q^3/12 ->", f"{1e-3**3/12:.3e} at q=1e-3")
zp3 = Pprime(ec, [4 * j + 2 for j in (25, 100, 400)], lambda L: L // 2, "fermion D, cubic band")
print("  => class (3,2) with v_b = 0, while cosine-band E/E-3/Ebar have v_b = 0 and class (2,1): class is NOT a function of v_b.")

print("=" * 110)
print("CE-4  asymmetric band (broken reflection): eps(k) = sin(k)/2 - sin(2k)/4 ; mirror eps_m(k) = eps(-k)")
ha = {1: -0.25j, -1: 0.25j, 2: 0.125j, -2: -0.125j}
ea = eps_from_hops(ha)
print("  identity check:", np.max(np.abs(ea(kk) - (np.sin(kk) / 2 - np.sin(2 * kk) / 4))),
      " eps(0)=eps(pi)=0; eps'(0)=eps''(0)=0, eps'''(0)=3/2; eps'(pi) = -1")
hm = {r: np.conj(h) for r, h in ha.items()}
em = eps_from_hops(hm)
for L in (8, 12):
    N = L // 2 - 1
    for nm, hh, ee in (("eps", ha, ea), ("mirror", hm, em)):
        occ, uniq, gap = fermion_occupation(ee, L, N)
        B, idx, H = hamiltonian(L, N, 'fermion', ring_bonds(L, hh))
        evs, V = np.linalg.eigh(H)
        ok_minus, ok_plus = True, True
        for m in range(L):
            s_ed = ed_support(evs, V, V[:, 0], evs[0], rho_diag(B, range(L), 2 * np.pi * m / L), N)
            ok_minus &= same_support(s_ed, ff_support(ee, L, occ, m))
            ok_plus &= same_support(s_ed, ff_support(ee, L, occ, (-m) % L))
        print(f"  ED L={L} N={N} [{nm}]: GS deg={int(np.sum(evs-evs[0]<1e-9))}; ED support matches k->k-q formula: "
              f"{ok_minus}; matches k->k+q formula: {ok_plus}")
za, sa = thermo_report(ea, [(-np.pi, 0.0)], "eps (occupied (-pi,0))")
zm, sm = thermo_report(em, [(0.0, np.pi)], "mirror (occupied (0,pi))")
Pprime(ea, [4 * j for j in (25, 100, 400)], lambda L: L // 2 - 1, "eps")
Pprime(em, [4 * j for j in (25, 100, 400)], lambda L: L // 2 - 1, "mirror")
print("  => same nu=1/2, same set of edges {0,pi}, same |v| data; z_P = 1 or 3 depending on chirality; occupied set is")
print("     NOT [-pi nu, pi nu] (k_b = pi nu fails); P restricted to q in (0,pi] sees only one orientation of edges.")

print("=" * 110)
print("CE-5  multiple disconnected pockets (exclusion): eps(k) = -2cos 3k - 0.5 cos k, mu = -1.0")
h3 = {3: -1.0, -3: -1.0, 1: -0.25, -1: -0.25}
e3 = eps_from_hops(h3)
mu = -1.0
O3 = occupied_intervals(e3, mu)
nu3 = sum(b - a for a, b in O3) / (2 * np.pi)
fp = sorted([np.mod(x + np.pi, 2 * np.pi) - np.pi for iv in O3 for x in iv])
vel = [(e3(x + 1e-6) - e3(x - 1e-6)) / 2e-6 for x in fp]
print(f"  occupied intervals: {[(round(a,5), round(b,5)) for a,b in O3]}  nu={nu3:.6f}")
print(f"  Fermi points: {[round(x,5) for x in fp]}  velocities eps'(k_F): {[round(v,5) for v in vel]}")
diffs = set()
for x in fp:
    for y in fp:
        d = np.mod(x - y, 2 * np.pi); d = min(d, 2 * np.pi - d)
        if all(abs(d - z) > 1e-6 for z in diffs):
            diffs.add(round(d, 7))
print(f"  predicted Q_soft = pairwise Fermi-point differences mod quotient: {sorted(diffs)} -> {len(diffs)}")
z5, s5 = thermo_report(e3, O3, "three pockets")
# finite-L sector with this filling: unique odd-N closed shell?
for L in (12, 600, 6000):
    k = kgrid(L); N = int(np.sum(e3(k) < mu))
    occ, uniq, gap = fermion_occupation(e3, L, N)
    msg = f"  L={L}: N=#(eps<mu)={N} odd={N%2==1} unique={uniq}"
    if L == 12 and uniq:
        deg, E0, ok = ed_vs_formula_fermion(h3, L, N, eps=e3, occ=occ)
        msg += f"  ED GS deg={deg}, ED support == form-factor: {ok}"
    print(msg)
print("  => z=1 but I-q = 10 (six Fermi points, not 2); there is no single k_b or v_b.")

print("=" * 110)
print("CE-6  two bands: two-leg ladder, rung t_perp=1: eps_-(k) = -2cos k - 1, eps_+(k) = -2cos k + 1; rho_q summed over legs")
def ladder_bonds(L, tp=1.0):
    b = []
    for leg in (0, 1):
        for j in range(L):
            i, l = leg * L + j, leg * L + (j + 1) % L
            b += [(i, l, -1.0), (l, i, -1.0)]
    for j in range(L):
        b += [(j, L + j, -tp), (L + j, j, -tp)]
    return b
def ladder_support(L, N, m):
    k = kgrid(L)
    lev = np.concatenate([-2 * np.cos(k) - 1, -2 * np.cos(k) + 1])
    order = np.argsort(lev, kind='stable'); occ = np.zeros(2 * L, bool); occ[order[:N]] = True
    ws = []
    for band in (0, 1):
        for a in range(L):
            b = (a - m) % L
            if occ[band * L + a] and not occ[band * L + b]:
                ws.append(lev[band * L + b] - lev[band * L + a])
    if m % L == 0:
        return [(0.0, float(N))]
    return group(np.array(ws), np.ones(len(ws)) / N)
for L, N in ((6, 4), (6, 8)):
    B, idx, H = hamiltonian(2 * L, N, 'fermion', ladder_bonds(L))
    ev, V = np.linalg.eigh(H)
    xs = list(range(L)) * 2
    ok = all(same_support(ed_support(ev, V, V[:, 0], ev[0], rho_diag(B, xs, 2 * np.pi * m / L), N),
                          ladder_support(L, N, m)) for m in range(L))
    print(f"  ED ladder L={L} N={N}: GS deg={int(np.sum(ev-ev[0]<1e-9))}, ED support == two-band intraband formula: {ok}")
mu = 0.5
km, kp = np.arccos(-(mu + 1) / 2), np.arccos((1 - mu) / 2)
print(f"  mu={mu}: k_-={km:.6f} (v={2*np.sin(km):.5f}), k_+={kp:.6f} (v={2*np.sin(kp):.5f}); predicted Q_soft = "
      f"{{0, 2k_-, 2k_+}} mod quotient = {sorted({0.0, round(min(2*km % (2*np.pi), 2*np.pi-2*km % (2*np.pi)),6), round(2*kp,6)})}")
def ladder_thermo(q):
    w1 = thermo_lower_edge(lambda k: -2 * np.cos(k) - 1, [(-km, km)], q)
    w2 = thermo_lower_edge(lambda k: -2 * np.cos(k) + 1, [(-kp, kp)], q)
    return min(w1, w2)
le = local_exponent(ladder_thermo, (1e-3, 1e-5))
print("  P: " + "; ".join(f"q={q:.0e}: dlogw/dlogq={s:.4f}" for q, w, s, r in le))
qs = np.linspace(1e-6, np.pi, 3001); w = np.array([ladder_thermo(q) for q in qs])
mins = [round(qs[t], 3) for t in range(1, len(qs) - 1) if w[t] <= w[t - 1] and w[t] <= w[t + 1] and w[t] < 5e-3]
print(f"  numerical zeros of w-(q) on (0,pi]: {[0.0] + mins} -> I-q = {1 + len(mins)}; z=1 with two different velocities")

print("=" * 110)
print("CE-7  d>1 scope of TARGET 2: 4x4 periodic square lattice, N=5 (closed shell for fermions)")
Lx = 4
def sq_bonds():
    b = []
    for x in range(Lx):
        for y in range(Lx):
            s = x + Lx * y
            for (dx, dy) in ((1, 0), (0, 1)):
                t = (x + dx) % Lx + Lx * ((y + dy) % Lx)
                b += [(s, t, -1.0), (t, s, -1.0)]
    return b
xs_x = [s % Lx for s in range(Lx * Lx)]
for kind in ('hcb', 'fermion'):
    B, idx, H = hamiltonian(Lx * Lx, 5, kind, sq_bonds())
    ev, V = np.linalg.eigh(H.real)
    deg = int(np.sum(ev - ev[0] < 1e-9))
    s = ed_support(ev, V, V[:, 0], ev[0], rho_diag(B, xs_x, np.pi / 2), 5) if deg == 1 else None
    print(f"  {kind:8s}: E0={ev[0]:+.10f} GS deg={deg}; support at q=(pi/2,0): "
          + (", ".join(f"({e:.4f},{W:.4f})" for e, W in s[:6]) + (" ..." if s and len(s) > 6 else "") if s else "n/a"))
print("  free-fermion E0 = -4 + 4*(-2) = -12.  HCB != free fermions in 2D (JW strings do not cancel).")
