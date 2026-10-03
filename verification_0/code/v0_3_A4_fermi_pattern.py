"""A4 (diagnostic): free-boson Fock state with ONE boson in each of the N lowest levels ('Fermi pattern').
(a) small-L exact many-body check (real-space ED) that it is an eigenstate, not the GS, and that its grouped rho_q
    support equals {eps(k-q)-eps(k): k in O} (k-q may be occupied: bosonic 'anything admissible'), incl. w<0;
(b) the quoted numbers at L=1026: E-E_gs and min w at q=pi/2, pi/8, pi/32; which sector(s) reproduce them."""
import numpy as np
from v0_3_lib import *

hops = {1: -1.0, -1: -1.0}
eps = eps_from_hops(hops)

# ---------------- (a) small-L ED
def create_k(vec, Bfrom, Bto_idx, L, k):
    """b_k^dag = L^{-1/2} sum_j e^{ikj} b_j^dag acting on vector over basis Bfrom -> basis Bto"""
    out = np.zeros(len(Bto_idx), complex)
    for t, b in enumerate(Bfrom):
        if vec[t] == 0:
            continue
        for j in range(L):
            n = list(b); amp = np.sqrt(n[j] + 1); n[j] += 1
            out[Bto_idx[tuple(n)]] += vec[t] * amp * np.exp(1j * k * j) / np.sqrt(L)
    return out

print("(a) small-L exact checks")
for L, N in ((10, 5), (10, 3), (6, 3)):
    occ, uniq, gap = fermion_occupation(eps, L, N)
    ms = np.nonzero(occ)[0]
    vec, Bprev = np.array([1.0 + 0j]), [tuple([0] * L)]
    for n_, m in enumerate(ms, start=1):
        Bn = basis(L, n_, 'boson'); idn = {b: t for t, b in enumerate(Bn)}
        vec = create_k(vec, Bprev, idn, L, 2 * np.pi * m / L); Bprev = Bn
    vec /= np.linalg.norm(vec)
    B, idx, H = hamiltonian(L, N, 'boson', ring_bonds(L, hops))
    assert B == Bprev
    E = np.real(np.vdot(vec, H @ vec)); res = np.linalg.norm(H @ vec - E * vec)
    ev, V = np.linalg.eigh(H)
    nocc = occ.astype(int)
    allok = True; minw = []
    for m in range(L):
        s_ed = ed_support(ev, V, vec, E, rho_diag(B, range(L), 2 * np.pi * m / L), N)
        s_fb = fb_support(eps, L, nocc, m)
        allok &= same_support(s_ed, s_fb)
        minw.append(min(e for e, W in s_ed))
    print(f"  L={L} N={N}: O(m)={list(ms)}  E={E:+.10f}  ||H psi - E psi||={res:.1e}  E_gs={ev[0]:+.10f} "
          f"(N eps(0)={-2*N})  E-E_gs={E-ev[0]:.10f}  ED support == boson-formula support (all m, with weights): {allok}")
    print(f"     min w per m: " + ", ".join(f"{w:+.6f}" for w in minw))

# ---------------- (b) L = 1026
L = 1026
k = kgrid(L); e = eps(k)
def stats(N):
    occ, uniq, gap = fermion_occupation(eps, L, N)
    dE = e[occ].sum() - N * (-2.0)
    nocc = occ.astype(int)
    mins = []
    for q in (np.pi / 2, np.pi / 8, np.pi / 32):
        m = nearest_m(L, q)
        mins.append((m, min(x for x, W in fb_support(eps, L, nocc, m))))
    return uniq, dE, mins

print("\n(b) L=1026")
for N in (513, 257, 769, 1, 3, 1025, 31, 33):
    uniq, dE, mins = stats(N)
    print(f"  N={N:4d} (pattern unique={uniq}): E-E_gs={dE:9.4f}  min w: " +
          ", ".join(f"q=2pi*{m}/L -> {w:+.6f}" for m, w in mins))
print("  closed form for N=L/2: E-E_gs = 2N - 2/sin(pi/L) =", f"{2*513 - 2/np.sin(np.pi/L):.6f}")
print("  closed-form min w (N=L/2, finite L): -4 sin(pi m/L) cos(pi(m+1)/L);  thermodynamic: -2 sin q")
for q in (np.pi / 2, np.pi / 8, np.pi / 32):
    m = nearest_m(L, q)
    print(f"    q={q:.6f}: m={m}  finite={-4*np.sin(np.pi*m/L)*np.cos(np.pi*(m+1)/L):+.6f}  -2 sin(q)={-2*np.sin(q):+.6f}")
print("  sectors N (1..1025) with |E-E_gs - 372.8| < 0.05:")
hits = []
for N in range(1, L):
    occ, uniq, gap = fermion_occupation(eps, L, N)
    dE = e[np.argsort(e, kind='stable')[:N]].sum() + 2 * N
    if abs(dE - 372.8) < 0.05:
        hits.append((N, uniq, dE))
print("   ", hits)
