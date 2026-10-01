# P-08 CORRECTION 01: algebraic checks against the RECORDED source-observable maps of the four lifts.
# Pure Python (complex floats + math); no repo data imported. Small explicit models only (1-3 modes),
# labelled as mathematics. See P08_CORRECTION_01.md.
import math, cmath

def matmul(A, B):
    n, m, k = len(A), len(B[0]), len(B)
    return [[sum(A[i][t]*B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]
def matadd(A, B, s=1):
    return [[A[i][j] + s*B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def dagger(A):
    return [[A[j][i].conjugate() for j in range(len(A))] for i in range(len(A[0]))]
def eye(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
def kron(A, B):
    n, m = len(A), len(B)
    return [[A[i//m][j//m]*B[i%m][j%m] for j in range(n*m)] for i in range(n*m)]
def maxabs(A, skip_last=0):
    n = len(A) - skip_last
    return max(abs(A[i][j]) for i in range(n) for j in range(n))
def annihilator(nmax):
    a = [[0.0]*(nmax+1) for _ in range(nmax+1)]
    for n in range(1, nmax+1):
        a[n-1][n] = math.sqrt(n)
    return a

print("=== C1 FERMIONIC: finite-dimensional observable algebra => no injective linear map from R[x] ===")
N = 23
def binom(n, k): return math.comb(n, k)
def smallest_d(dim):
    d = 0
    while binom(N + d, d) <= dim:
        d += 1
    return d
# Three recorded fermionic algebras, all finite-dimensional (distinguish Hilbert-space vs algebra dimension):
#   real Clifford Cl(R^23): real dim 2^23;  exterior algebra Lambda(R^23): 2^23;
#   fermionic Fock space Lambda(C^23): Hilbert dim 2^23  =>  CAR algebra = B(Lambda C^23) = M_{2^23}(C): complex dim 2^46 (real 2^47);
#   even (parity-superselected) subalgebras: 2^22 resp. 2^45 — still finite.
for name, dim in [("Cl(R^23) / Lambda(R^23) [real dim]", 2**N), ("CAR(C^23) = M_{2^23}(C) [complex dim]", 2**(2*N)), ("CAR(C^23) as a real space", 2**(2*N+1))]:
    d = smallest_d(dim)
    print(f"{name} = {dim}: smallest d with dim R[x_1..x_23]_(deg<=d) = C(23+d,d) > that is d = {d} ({binom(N+d,d)})")
print("=> in every recorded fermionic algebra, polynomials of bounded degree already cannot inject; R[x] (infinite-dim) a fortiori; independent of the coordinate map")
# odd-image lemma, concrete 2-mode CAR (Jordan-Wigner), 4x4
sz = [[1,0],[0,-1]]; sm = [[0,1],[0,0]]; I2 = eye(2)
a1 = kron(sm, I2); a2 = kron(sz, sm)
anti = matadd(matmul(a1,a2), matmul(a2,a1))
print(f"CAR 2-mode: ||{{a1,a2}}|| = {maxabs(anti):.1e}, ||a1 a2|| = {maxabs(matmul(a1,a2)):.3f}  (nonzero) => a1a2 = -a2a1 != a2a1: no commutative image")
g1 = matadd(a1, dagger(a1)); g2 = matadd(a2, dagger(a2))      # Majorana-type, g_i^2 = 1
print(f"Majorana: ||g1^2 - I|| = {maxabs(matadd(matmul(g1,g1), eye(4), -1)):.1e} => pi(x_1^2 - 1) = 0, non-injective; ||g1 g2|| = {maxabs(matmul(g1,g2)):.3f} (invertible) and {{g1,g2}} = {maxabs(matadd(matmul(g1,g2),matmul(g2,g1))):.1e}")

print()
print("=== C2 BOSONIC, doubling complex structure: site-coordinate images COMMUTE ===")
nmax = 6
a = annihilator(nmax); I = eye(nmax+1)
A1 = kron(a, I); A2 = kron(I, a)
Phi1 = [[(A1[i][j] + dagger(A1)[i][j])/math.sqrt(2) for j in range(len(A1))] for i in range(len(A1))]
Phi2 = [[(A2[i][j] + dagger(A2)[i][j])/math.sqrt(2) for j in range(len(A2))] for i in range(len(A2))]
comm = matadd(matmul(Phi1, Phi2), matmul(Phi2, Phi1), -1)
print(f"doubling (e_1, e_2 real, orthogonal; separate modes): ||[Phi(e1),Phi(e2)]|| = {maxabs(comm):.1e}  => commuting => functional calculus f(Phi(e1),Phi(e2)) is a unital hom")
# non-doubling J on R^2 (N even): e_2 = J e_1 = i e_1 -> Phi(e_2) = (a - a^dag)/(i sqrt2); commutator = i
nmax1 = 30
a = annihilator(nmax1)
P1 = [[(a[i][j] + dagger(a)[i][j])/math.sqrt(2) for j in range(nmax1+1)] for i in range(nmax1+1)]
P2 = [[(a[i][j] - dagger(a)[i][j])/(1j*math.sqrt(2)) for j in range(nmax1+1)] for i in range(nmax1+1)]
c12 = matadd(matmul(P1, P2), matmul(P2, P1), -1)
print(f"non-doubling J on R^2 (e2 = J e1): [Phi(e1),Phi(e2)] = {c12[0][0]:.3f} (expected i; truncation only at top level) => images do NOT commute")

print()
print("=== C3 BOSONIC dichotomy: coherent state reproduces <Phi> but not <Phi^2>; Wick does, but is not multiplicative ===")
nmax2 = 60
a = annihilator(nmax2); ad = dagger(a)
x = 0.9; alpha = x/math.sqrt(2)                    # <alpha|Phi|alpha> = sqrt2 Re alpha = x
coh = [math.exp(-abs(alpha)**2/2) * alpha**n / math.sqrt(math.factorial(n)) for n in range(nmax2+1)]
Phi = [[(a[i][j] + ad[i][j])/math.sqrt(2) for j in range(nmax2+1)] for i in range(nmax2+1)]
Phi2 = matmul(Phi, Phi)
def expv(M): return sum(coh[i].conjugate()*M[i][j]*coh[j] for i in range(len(M)) for j in range(len(M))).real
e1, e2 = expv(Phi), expv(Phi2)
print(f"x = {x}: <Phi> = {e1:.6f}, <Phi^2> = {e2:.6f}, x^2 = {x*x:.6f}, variance = {e2 - e1*e1:.6f} (= 1/2 vacuum variance)")
print(f"Wick: <:Phi^2:> = <Phi^2> - 1/2 = {e2-0.5:.6f} = x^2 (reproduces) but :Phi:^2 = Phi^2 != :Phi^2: = Phi^2 - 1/2 (not multiplicative)")
# The RECORDED quadratic image (evaluation section 3, R-b column) is dGamma(P), not (1/2) sum P_ij Phi_i Phi_j.
# One mode, P = 1: dGamma = a^dag a ; (1/2)Phi^2 = (1/4)(a^2 + a^dag^2 + 2 a^dag a + 1). They differ as operators,
# but dGamma(P) is expectation-reproducing in coherent states (<x|a^dag a|x> = x^2/2 = V(x)).
NN = matmul(ad, a)
halfPhi2 = [[0.5*Phi2[i][j] for j in range(nmax2+1)] for i in range(nmax2+1)]
print(f"recorded quadratic image dGamma(1) = a^dag a: <x|dGamma|x> = {expv(NN):.6f} = x^2/2 = {x*x/2:.6f} (reproduces V); <x|(1/2)Phi^2|x> = {expv(halfPhi2):.6f} = x^2/2 + 1/4")
print(f"||dGamma(1) - (1/2)Phi^2|| (interior block) = {maxabs(matadd(NN, halfPhi2, -1), skip_last=2):.3f} != 0 => the recorded identification is NOT the restriction of a homomorphism (it is the expectation-reproducing choice)")
# LS-5: Mehler/OU on the commutative image is not multiplicative
T = 0.6
print(f"Mehler P_T on the Lagrangian image (1 mode, T={T}): P_T(x^2) - (P_T x)^2 = 1 - T^2 = {1-T*T:.3f} != 0 => dynamics not intertwined with Koopman")

print()
print("=== C2b BOSONIC, NON-* image pi_a(x_i) = sqrt2 a_i: satisfies ALL THREE sub-readings (price: not self-adjoint) ===")
# a|z> = z|z>, coherent z = x/sqrt2  =>  <x|sqrt2 a|x> = x and <x|(sqrt2 a)^n|x> = x^n EXACTLY (joint eigenvector).
sq2a = [[math.sqrt(2)*a[i][j] for j in range(nmax2+1)] for i in range(nmax2+1)]
p1 = expv(sq2a); p2 = expv(matmul(sq2a, sq2a)); p3 = expv(matmul(matmul(sq2a, sq2a), sq2a))
print(f"x={x}: <pi_a(x)> = {p1:.6f} = x ; <pi_a(x^2)> = {p2:.6f} = x^2 = {x*x:.6f} ; <pi_a(x^3)> = {p3:.6f} = x^3 = {x**3:.6f}")
print("  => coherent states are JOINT EIGENVECTORS of the a_i: expectation-reproducing for EVERY polynomial (no vacuum variance)")
print("  => and quasi-free Heisenberg acts by substitution a -> T a, so f(a) -> f(Ta): multiplicative AND Koopman-intertwined")
print(f"  PRICE: pi_a(x_i)* = sqrt2 a_i^dag != pi_a(x_i)  (not a *-representation) -- the bosonic analogue of the recorded Lambda-F odd-image price LS-8")
# *-case, correct proof: coherent-state non-orthogonality rules out EVERY symmetric generator image at once
xp = 0.3; zp = xp/math.sqrt(2)
cohp = [math.exp(-abs(zp)**2/2) * zp**n / math.sqrt(math.factorial(n)) for n in range(nmax2+1)]
ov = sum(coh[i].conjugate()*cohp[i] for i in range(nmax2+1))
print(f"*-case: <x|x'> = {ov.real:.6f} = exp(-|x-x'|^2/4) = {math.exp(-(x-xp)**2/4):.6f} != 0 => |x> cannot be an eigenvector of a symmetric pi(x_i) for every x")
print(f"normalization note: in the Phi=(a+a^dag)/sqrt2 convention Var_vac = 1/2 and the quasi-free Mehler defect on the Phi-image is (1-T^2)/2 = {(1-T*T)/2:.3f}; LS-5's '1-T^2' = {1-T*T:.3f} is the unit-variance L^2(gamma) convention")

print()
print("=== C4 SZ.-NAGY ===")
print("--- C4a: under the RECORDED Poisson/function-algebra reading (L0_LIFT_SELECTION_01 section 1: Lambda-H row algebra = 'Poisson'),")
print("---      pi(f) = f o P_H (pullback along the projection K -> H, dual to the inclusion H in K) IS faithful ---")
nK, nH = 5, 3
def P_H(v): return v[:nH]                                  # orthogonal projection K -> H
vx = [0.41, -0.77, 0.52, 1.9, -0.6]                        # a point of K with NONZERO bath part
vx0 = [0.41, -0.77, 0.52, 0.0, 0.0]                        # the recorded zero-bath point x (+) 0
polys = {'1': lambda u: 1.0, 'x1 (odd)': lambda u: u[0], 'x1 x2': lambda u: u[0]*u[1], 'x1^2 x3 - 2x2': lambda u: u[0]**2*u[2] - 2*u[1]}
for nm, fpoly in polys.items():
    lhs = fpoly(P_H(vx0)); rhs = fpoly([0.41, -0.77, 0.52])
    print(f"  delta_(x+0)(pi({nm})) = {lhs:+.6f}  vs  f(x) = {rhs:+.6f}   match={abs(lhs-rhs)<1e-12}")
print("  multiplicative: (fg) o P_H = (f o P_H)(g o P_H) identically; injective: P_H surjective => f o P_H = 0 forces f = 0")
print("  => survives 2-map AND 2-state under this reading, INCLUDING odd polynomials (cf. cotangent's f o pi, same mechanism)")
print("--- C4b: under a QUANTUM operator / vector-state reading it fails instead ---")
dimBH = nH*nH if False else 23*23
print(f"  recorded operator class B(H)(+)0 (compressions): dim B(H) = 23^2 = {dimBH} < dim R[x_1..x_23]_(deg<=3) = C(26,3) = {binom(26,3)} => no injective linear map")
# recorded R-b image V_P = (1/2) x^T P x -> P (+) 0 is not multiplicative: commuting source quadratics -> non-commuting images
E11 = [[1,0,0],[0,0,0],[0,0,0]]; E12s = [[0,1,0],[1,0,0],[0,0,0]]; E22 = [[0,0,0],[0,1,0],[0,0,0]]
c = matadd(matmul(E11,E12s), matmul(E12s,E11), -1)
pr = matmul([[2*E11[i][j] for j in range(3)] for i in range(3)], [[2*E22[i][j] for j in range(3)] for i in range(3)])
print(f"  recorded quadratic image P(+)0: ||[P_(x1^2), P_(x1x2)]|| = {max(abs(c[i][j]) for i in range(3) for j in range(3)):.3f} != 0 though x1^2, x1x2 commute;")
print(f"    and image(x1^2)*image(x2^2) = {max(abs(pr[i][j]) for i in range(3) for j in range(3)):.3f} = 0 though x1^2 x2^2 != 0 => not multiplicative")
print("--- C4c: vector-state expectations are quadratic forms => cannot reproduce odd polynomials ---")
import random
random.seed(7)
n = 5
Aop = [[complex(random.uniform(-1,1), random.uniform(-1,1)) for _ in range(n)] for _ in range(n)]
v = [random.uniform(-1,1) for _ in range(3)] + [0.0, 0.0]      # x (+) 0 in K = H (+) H^perp
def qf(A, v): return sum(v[i]*A[i][j]*v[j] for i in range(n) for j in range(n))
print(f"<x(+)0, A (x(+)0)> = {qf(Aop, v):.6f}; with x -> -x: {qf(Aop, [-t for t in v]):.6f}  (even) ; target f(x) = x_1 = {v[0]:.6f} -> -x_1 = {-v[0]:.6f} (odd) => no operator A reproduces x_1")
print("Also: the recorded coordinate 'image' <e_i, .> is a linear functional on K, not an operator; the inclusion H in K maps vectors, not the algebra R[x].")

print()
print("=== C5 COTANGENT: pullback along pi: T*R^N -> R^N is an injective unital homomorphism; H = p^T f projects to the source flow ===")
# 1-d check: H = p f(x); Hamilton: xdot = dH/dp = f(x) (source flow), pdot = -p f'(x).
f = lambda x: -1.3*x - 4*0.2*x**3
x0, p0, h = 0.7, 2.0, 1e-3
# one explicit RK4 step on (x,p) vs source flow on x alone
def rk4(F, y, h):
    k1 = F(y); k2 = F([y[i]+h/2*k1[i] for i in range(len(y))]); k3 = F([y[i]+h/2*k2[i] for i in range(len(y))]); k4 = F([y[i]+h*k3[i] for i in range(len(y))])
    return [y[i] + h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]) for i in range(len(y))]
Fc = lambda y: [f(y[0]), -y[1]*(-1.3 - 2.4*y[0]**2)]
Fs = lambda y: [f(y[0])]
xc = rk4(Fc, [x0, p0], h)[0]; xs = rk4(Fs, [x0], h)[0]
print(f"x-projection of the cotangent flow vs source flow after one step: |diff| = {abs(xc-xs):.1e} (identical; independent of p0) => (f o pi) o Phi_t = (f o phi_t) o pi")
print("ALL P-08B RECORDED-MAP CHECKS COMPLETE")
