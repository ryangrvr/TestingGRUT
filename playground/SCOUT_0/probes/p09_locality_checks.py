# P-09: locality-of-representation, one reconnaissance pass. Pure Python (no numpy).
# Substrate = the recorded anchor bath block K_b (build_K(24,0), (1:,1:) block, N = 23; pin 0.3, unit springs).
# Abstract small-mode models (2-3 modes) are labelled as mathematics. See P09_RESULT.md.
import math

N, PIN = 23, 0.3

def build_K(n):
    K = [[0.0]*n for _ in range(n)]
    for i in range(n):
        K[i][i] = PIN
    for i in range(n-1):
        K[i][i+1] -= 1; K[i+1][i] -= 1; K[i][i] += 1; K[i+1][i+1] += 1
    return K

def matmul(A, B):
    return [[sum(A[i][t]*B[t][j] for t in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def T(A): return [list(r) for r in zip(*A)]
def maxdiff(A, B): return max(abs(A[i][j]-B[i][j]) for i in range(len(A)) for j in range(len(A[0])))
def eye(n): return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

K24 = build_K(24)
Kb = [row[1:] for row in K24[1:]]

def graph_support(M, tol=1e-12):
    return {(i, j) for i in range(len(M)) for j in range(len(M)) if i != j and abs(M[i][j]) > tol}
G = graph_support(Kb)
print("=== C0 substrate: K_b support graph = nearest-neighbour path ===")
print("  edges:", len(G)//2, " all |i-j|=1:", all(abs(i-j) == 1 for i, j in G))

# closed-form eigenpairs (L0-1e pre-freeze review): lam_k = 2.3 - 2cos((2k-1)pi/47), v_k(i) ~ sin((2k-1) i pi/47)
lam, V = [], []
for k in range(1, N+1):
    th = (2*k-1)*math.pi/47
    v = [math.sin(th*i) for i in range(1, N+1)]
    nv = math.sqrt(sum(x*x for x in v)); v = [x/nv for x in v]
    lam.append(2.3 - 2*math.cos(th)); V.append(v)
res = max(abs(sum(Kb[i][j]*V[k][j] for j in range(N)) - lam[k]*V[k][i]) for k in range(N) for i in range(N))
print("  closed-form eigen-residual:", f"{res:.1e}", " min lam:", f"{min(lam):.4f}", "(K_b PD, accretive)")

def fn(f):  # f(K_b) by spectral calculus
    return [[sum(f(lam[k])*V[k][i]*V[k][j] for k in range(N)) for j in range(N)] for i in range(N)]

print("\n=== C1 cotangent lift H = p.f(x), f(x) = -K_b x - 4 beta x^3: generator is graph-local ===")
beta, x = 0.7, [math.sin(1.3*i+0.2) for i in range(N)]
# mixed Hessian d2H/dp_i dx_j = -(K_b)_ij - 12 beta x_i^2 delta_ij ; d2H/dp dp = 0 ; d2H/dx dx = -12 beta * diag(p_i x_i) * 2
Hpx = [[-Kb[i][j] - (12*beta*x[i]**2 if i == j else 0.0) for j in range(N)] for i in range(N)]
print("  off-diagonal support of d2H/dp dx equals K_b graph:", graph_support(Hpx) == G)
print("  site Poisson algebras {x_i,p_i} commute across sites (canonical bracket) -> commuting site net: by construction")

print("\n=== C2 Sz.-Nagy / one-particle dilation of e^{-K_b t}: coupling locality is NOT a unitary invariant ===")
# Any unitary (quantum-noise / Langevin) dilation with no system Hamiltonian needs channel map C with C^T C = 2 K_b.
Cmin = fn(lambda l: math.sqrt(2*l))          # minimal / canonical channel basis: C = (2K_b)^{1/2}, N channels
print("  C_min^T C_min == 2K_b:", f"{maxdiff(matmul(T(Cmin), Cmin), [[2*v for v in r] for r in Kb]):.1e}")
decay = [abs(Cmin[0][j]) for j in range(N)]
print("  minimal coupling is DENSE: |C_min[1,j]| for j=1,2,3,6,12,23:",
      ", ".join(f"{decay[j]:.2e}" for j in (0, 1, 2, 5, 11, 22)))
nz = sum(1 for i in range(N) for j in range(N) if abs(Cmin[i][j]) > 1e-12)
print("  nonzero entries:", nz, "of", N*N, " (K_b has", N + len(G), ")")
rates = [math.log(abs(Cmin[0][j])/abs(Cmin[0][j+1])) for j in (4, 8, 12)]
print("  exponential decay rate per site (quasi-local, not strictly local):", ", ".join(f"{r:.3f}" for r in rates))
# local (non-minimal) channel basis: K_b = B B^T with B = [sqrt(pin) e_i | e_1 (cut spring) | e_i - e_{i+1}]
cols = [[math.sqrt(PIN) if r == i else 0.0 for r in range(N)] for i in range(N)]
cols.append([1.0 if r == 0 else 0.0 for r in range(N)])
for i in range(N-1):
    cols.append([1.0 if r == i else (-1.0 if r == i+1 else 0.0) for r in range(N)])
B = T(cols)                                    # N x 46
Cloc = [[math.sqrt(2)*v for v in c] for c in cols]   # 46 x N, each channel touches <= 2 adjacent sites
print("  local basis: B B^T == K_b:", f"{maxdiff(matmul(B, T(B)), Kb):.1e}",
      " channels:", len(Cloc), " max sites per channel:", max(sum(1 for v in c if v != 0) for c in Cloc))
Wiso = matmul(Cloc, fn(lambda l: 1/math.sqrt(2*l)))      # polar factor: Cloc = W Cmin
print("  W isometry (W^T W = I_23):", f"{maxdiff(matmul(T(Wiso), Wiso), eye(N)):.1e}",
      "  W C_min == C_loc:", f"{maxdiff(matmul(Wiso, Cmin), Cloc):.1e}")
print("  => same dilation (minimal part) up to a channel-space isometry; strict locality depends on the bath factorization")
# finite-time maps are dense for every lift alike (shared descended data): -K_b is Metzler and irreducible
# (connected path), so e^{-K_b t} > 0 entrywise for every t > 0; shown at t = 20 where it clears roundoff.
Et = fn(lambda l: math.exp(-l*20.0))
print("  shared by all lifts: e^{-K_b t} is entrywise positive (irreducible Metzler); t=20: [.]_{1,23} =",
      f"{Et[0][22]:.2e}", " min entry:", f"{min(min(r) for r in Et):.2e}")

print("\n=== C3 Lambda-B: locality rides on the doubling complex structure J (abstract, 2 sites) ===")
# real symplectic form induced by J on R^{2N} (x_1..x_N, pi_1..pi_N): sigma(u,v) = <u, J v>; [Phi(u),Phi(v)] = i sigma(u,v)
def Jstd(n):
    J = [[0.0]*(2*n) for _ in range(2*n)]
    for i in range(n):
        J[i][n+i] = -1.0; J[n+i][i] = 1.0
    return J
def sigma(J, u, v): return sum(u[i]*sum(J[i][j]*v[j] for j in range(len(v))) for i in range(len(u)))
n = 2; e = lambda k: [1.0 if r == k else 0.0 for r in range(2*n)]
J0 = Jstd(n)
print("  recorded site-local doubling (pi_i conjugate to x_i): sigma(x1,x2) =", sigma(J0, e(0), e(1)),
      " sigma(x1,pi2) =", sigma(J0, e(0), e(3)), "-> disjoint-site Weyl algebras commute")
c, s = math.cos(0.6), math.sin(0.6)                    # orthogonal O rotating x_2 into pi_1
O = eye(2*n); O[1][1], O[1][2], O[2][1], O[2][2] = c, -s, s, c
J1 = matmul(matmul(O, J0), T(O))
print("  another metric-compatible J (J^2 = -1:", f"{maxdiff(matmul(J1, J1), [[-v for v in r] for r in eye(2*n)]):.0e})"
      ": sigma(x1,x2) =", f"{sigma(J1, e(0), e(1)):.4f}", "-> site-1 and site-2 coordinate images no longer commute")
print("  => Lambda-B tensor-product locality holds iff the (priced, V-7) doubling J is chosen site-local")

print("\n=== C4 Lambda-F: graded vs tensor-product locality (abstract, 3 sites, Jordan-Wigner) ===")
def kron(A, Bm):
    n1, n2 = len(A), len(Bm)
    return [[A[i//n2][j//n2]*Bm[i % n2][j % n2] for j in range(n1*n2)] for i in range(n1*n2)]
def kr(*ms):
    out = [[1.0]]
    for m in ms: out = kron(out, m)
    return out
I2, Z, sm = eye(2), [[1.0, 0.0], [0.0, -1.0]], [[0.0, 1.0], [0.0, 0.0]]
a = [kr(*([Z]*i + [sm] + [I2]*(2-i))) for i in range(3)]
ad = [T(m) for m in a]
add = lambda A, Bm, s=1: [[A[i][j]+s*Bm[i][j] for j in range(len(A))] for i in range(len(A))]
nrm = lambda A: max(abs(v) for r in A for v in r)
anti = nrm(add(matmul(a[0], a[2]), matmul(a[2], a[0])))
comm = nrm(add(matmul(a[0], a[2]), matmul(a[2], a[0]), -1))
hop = add(matmul(ad[0], a[1]), matmul(ad[1], a[0]))
ev = nrm(add(matmul(hop, a[2]), matmul(a[2], hop), -1))
print(f"  odd generators, disjoint sites: {{a1,a3}} = {anti:.0f}, [a1,a3] norm = {comm:.0f}  (tensor-product locality FAILS)")
print(f"  even hop(1,2) vs odd a3: commutator norm = {ev:.0f}  (graded locality HOLDS; even net commutes, V-6/LS-8)")
print("  => Lambda-F passes iff locality is read as graded (record: D-6 CRITERION, ruling 02 s3)")
