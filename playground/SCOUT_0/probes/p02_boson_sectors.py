# P-02: statistics-swap of SF-1 (requires numpy, sympy). See P02_SF1_CONTROL_MAP.md, P02_RESULT.md.
# Main model: free bosons, H_B = -sum (b_j^dag b_{j+1} + h.c.), periodic ring, conserved N, sector ground state,
# readout rho_q = sum_j e^{-iqj} n_j, support of S(q, w) (weight > 1e-12, grouped over degenerate eigenspaces).
import itertools, math
import numpy as np, sympy as sp

TOL_W, TOL_E = 1e-12, 1e-9

def support(H, basis_nj, L, gs_tol=1e-9):
    """Exact ground state of the sector block H, and the grouped support of rho_q|0> for every q = 2 pi m/L."""
    E, V = np.linalg.eigh(H)
    gap = E[1] - E[0] if len(E) > 1 else np.inf
    g = V[:, 0]
    groups, start = [], 0
    for i in range(1, len(E)+1):
        if i == len(E) or E[i] - E[start] > 1e-8:
            groups.append((E[start], V[:, start:i])); start = i
    out = {}
    for m in range(1, L):
        q = 2*np.pi*m/L
        rho = np.array([sum(np.exp(-1j*q*j)*n[j] for j in range(L)) for n in basis_nj])   # diagonal in site basis
        v = rho*g
        pts = []
        for Eg, Vg in groups:
            w = np.sum(np.abs(Vg.conj().T @ v)**2)
            if w > TOL_W: pts.append(round(Eg - E[0], 9))
        out[m] = sorted(pts)
    return E[0], gap, g, out

def boson_block(L, N):
    basis = [c for c in itertools.product(range(N+1), repeat=L) if sum(c) == N]
    idx = {c: i for i, c in enumerate(basis)}
    H = np.zeros((len(basis), len(basis)))
    for c, i in idx.items():
        for a in range(L):
            for (s, t) in ((a, (a+1) % L), ((a+1) % L, a)):        # b_s^dag b_t
                if c[t] == 0: continue
                d = list(c); amp = math.sqrt(d[t]); d[t] -= 1; amp *= math.sqrt(d[s]+1); d[s] += 1
                H[idx[tuple(d)], i] -= amp
    return H, basis

def hcb_block(L, N):   # hard-core bosons: n_j in {0,1}, commuting operators, NO Jordan-Wigner signs
    basis = [c for c in itertools.product((0, 1), repeat=L) if sum(c) == N]
    idx = {c: i for i, c in enumerate(basis)}
    H = np.zeros((len(basis), len(basis)))
    for c, i in idx.items():
        for a in range(L):
            for (s, t) in ((a, (a+1) % L), ((a+1) % L, a)):
                if c[t] == 1 and c[s] == 0:
                    d = list(c); d[t] = 0; d[s] = 1; H[idx[tuple(d)], i] -= 1
    return H, basis

eps = lambda k: -2*math.cos(k)
def ff_support(L, N):  # free-fermion I-FF particle-hole set (SF-1 identity), odd N symmetric fill
    M = (N-1)//2; occ = {(j % L) for j in range(-M, M+1)}
    out = {}
    for m in range(1, L):
        out[m] = sorted({round(eps(2*math.pi*((k+m) % L)/L) - eps(2*math.pi*k/L), 9) for k in occ if (k+m) % L not in occ})
    return out
def bose_closed(L):   # free-boson closed form: one line, w = eps(q) - eps(0) = 4 sin^2(q/2)
    return {m: [round(4*math.sin(math.pi*m/L)**2, 9)] for m in range(1, L)}

print("=== V-BFOCK / A-BOSE: exact free-boson sector blocks ===")
for L in (6, 10):
    for N in (1, 3, 5):
        H, B = boson_block(L, N)
        E0, gap, g, sup = support(H, B, L)
        # A-BOSE: E0 = N*eps(0) with a positive gap <=> unique ground state = (b_{k=0}^dag)^N|0>/sqrt(N!)
        ok = sup == bose_closed(L)
        print(f"  L={L:2d} N={N}: dim={len(B):5d}  E0={E0:+.10f} (N*eps(0)={-2*N:+d})  gap={gap:.4f}  support == {{4 sin^2(q/2)}} for all q: {ok}")

print("\n=== HCB diagnostic (exclusion WITHOUT antisymmetry; NEW ASSUMPTION) vs free-fermion I-FF ===")
for L, Ns in ((6, (1, 3, 5)), (10, (1, 3, 5, 9)), (12, (1, 3, 9, 11, 6))):
    for N in Ns:
        H, B = hcb_block(L, N)
        E0, gap, g, sup = support(H, B, L)
        if N % 2:
            print(f"  L={L:2d} N={N:2d} (odd):  gap={gap:.4f}  support == free-fermion SF-1 set: {sup == ff_support(L, N)}")
        else:
            print(f"  L={L:2d} N={N:2d} (even): gap={gap:.2e}  (JW string gives an antiperiodic fermion ring: degenerate/twisted; reported, not used)")

print("\n=== Closed-form classes, prescriptions P and P' (all SF-1 families) ===")
q, Ls = sp.symbols('q L', positive=True)
w_b = 4*sp.sin(q/2)**2
zP = sp.limit(sp.log(w_b)/sp.log(q), q, 0, '+')
w1 = w_b.subs(q, 2*sp.pi/Ls); zPp = sp.limit(-sp.log(w1)/sp.log(Ls), Ls, sp.oo)
print(f"  free bosons, EVERY family (D, E, D-1/4, D-3/4, E-3, E-bar, C-6): w^-(q) = {w_b}, z_P = {zP}, z_P' = {zPp}, soft set {{0}} -> class (2, 1)")
print("  SF-1 fermions (record): D, D-1/4, D-3/4 -> (1, 2); E, E-3, E-bar -> (2, 1); C-6 multiscale")

print("\n=== Sector invariant -> law class: occupation edge k_b against the dispersion ===")
nu = sp.symbols('nu', nonnegative=True); k = sp.symbols('k')
for name, kb in (("fermions (exclusion)", sp.pi*nu), ("free bosons", sp.Integer(0))):
    v = sp.diff(-2*sp.cos(k), k).subs(k, kb)
    print(f"  {name:22}: k_b^inf(nu) = {kb};  edge velocity eps'(k_b) = {sp.simplify(v)}")
print("  rule: z = 1 iff eps'(k_b^inf) != 0 (occupation edge in band interior); z = 2 iff k_b^inf at a band extremum")
for nuv, lab in ((sp.Rational(1, 2), "D"), (sp.Rational(1, 4), "D-1/4"), (sp.Rational(3, 4), "D-3/4"), (0, "E/E-3"), (1, "E-bar")):
    v = sp.simplify(2*sp.sin(sp.pi*nuv))
    print(f"    fermion {lab:6} nu={str(nuv):4}: eps'(k_F) = {v} -> z = {1 if v != 0 else 2}   (record: {'1' if lab.startswith('D') else '2'})")

print("\n=== FPB diagnostic: bosons placed in the Fermi pattern (exact eigenstate, NOT the sector ground state) ===")
L, N = 1026, 513; M = (N-1)//2
occ = [2*math.pi*j/L for j in range(-M, M+1)]
E_fpb = sum(eps(kk) for kk in occ); E_gs = N*eps(0)
for qq in (math.pi/2, math.pi/8, math.pi/32):
    w = [eps(kk+qq) - eps(kk) for kk in occ]      # no Pauli blocking: transitions into occupied levels allowed
    print(f"  q={qq:.4f}: min support w = {min(w):+.4f} (negative: de-excitations present); E_FPB - E_gs = {E_fpb-E_gs:.1f} > 0")
print("  -> the Fermi pattern needs exclusion twice: to BE the ground state, and to BLOCK in-sea transitions")
