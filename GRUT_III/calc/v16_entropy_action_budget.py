#!/usr/bin/env python3
"""
V16 -- Entropy-action budget: why bounded-degree graph models cannot select geometry, and why causal-set-type
nonlocal actions can at least compete.  (Synthesis 13 Sec. 2; numbers for the codex fixed-bond / CP12 obstructions.)

CLAIM (heuristic, made quantitative).  A weight w(c) can make a "good" (geometric) class dominate a "bad" class only if the
weights can fight the bad class's excess log-count  Delta L = ln|bad| - ln|good|.
 * Boltzmann weights e^{-beta E} with extensive E (|E| <= C N): the codex theorem gives P(good) <= exp(-Delta L + O(N)).
 * Lorentzian weights e^{i lambda S}: a phase can cancel, but for an action that is a sum of ~n roughly independent local
   terms the class sum is generically |bad| x exp(-n g(lambda)) with g = O(1) -- a product of local characteristic
   functions.  g ~ lambda^2 Var(S)/(2n) only at weak coupling (the Gaussian regime); at strong coupling g is smaller, and it
   can vanish (for the degree action below, S is always even, so at lambda = pi there is no suppression at all).  So
   cancellation needs n >~ Delta L, generically.  This is a heuristic, not a theorem: exact cancellations at tuned couplings
   are not excluded, and the phase also acts on the good class, which is not evaluated here.
For N vertices and K = alpha N edges, Delta L >= (alpha - 1) N ln N + O(N) for any family with <= N! e^{O(N)} labelled members
(lattice templates), while a LOCAL action on a sparse graph has n = O(N) terms: Var(S) = O(N).  So neither Boltzmann nor
Lorentzian weights with a fixed-strength local action can select the family; the deficit grows like ln N.
Layered orders (ln count ~ (N^2/4) ln 2 for the Kleitman-Rothschild three-layer class): an action counting links/relations
(causal-set theory) has ~N^2 independent terms, so its phase can suppress the entropically dominant, balanced layered orders by
exp(-c N^2): for a fixed balanced, naturally labelled bipartition of two-level orders the link-action sum is exactly
(1 + e^{i theta})^{N^2/4} (checked by brute force below), absolutely suppressed only for |cos(theta/2)| < 1/2.  Literature (checked, 13): Loomis & Carlip 2018 (CQG 35, 024002) -- two-level (bilayer) orders, leading
order in N^2; Mathur, Singh & Surya 2021 (CQG 38, 045017) -- k-layer / KR-type orders with the link action; Carlip, Carlip &
Surya 2023-24 (CQG 40, 095004; CQG 41, 145005) -- KR / layered orders with the full BDG action.  CAVEATS: the sum over ALL
two-level orders is not exponentially small in N^2 (unbalanced orders dominate; workflow computation recorded in 13), results
depend on which subset is summed, and no paper compares these sums with the manifold-like contribution.  So the N^2 budget is
NECESSARY for selection and has been shown sufficient to suppress the dominant layered classes, not to select manifolds.

CHECKS.
 16.1 exact: ln C(M, alpha N) - ln N! = (alpha - 1) N ln N + O(N) (Stirling, numerically for N to 10^6).
 16.2 Monte Carlo over uniform fixed-K graphs (alpha = 2): the valence action S = sum_v (d_v - 4)^2 has Var(S)/N -> const
      (independent-Poisson reference 36 - 4 = 32), so at weak coupling ln |<e^{i lambda S}>| = -lambda^2 Var(S)/2 = -O(N)
      (checked at lambda = 0.05); at strong coupling the per-vertex rate g(lambda) is non-Gaussian and smaller (independent-
      Poisson comparison), and at lambda = pi it is zero because S is always even.  In every case the suppression is
      exp(-O(N)) or none -- never exp(-N ln N).  The ratio (alpha - 1) N ln N / Var(S) grows like ln N.
 16.3 balanced bilayer (two-level) orders with a fixed, naturally labelled bipartition: the link-action sum is
      (1 + e^{i theta})^{N^2/4} (brute force at small N), an exponent in N^2 matching the (N^2/4) ln 2 entropy; absolute
      suppression (|class sum| < 1) requires |cos(theta/2)| < 1/2, i.e. theta mod 2 pi in (2 pi/3, 4 pi/3).
"""
import math, random
import numpy as np
ok_all = True
def check(label, cond):
    global ok_all
    print(f"  {label}:", "PASS" if cond else "FAIL"); ok_all &= bool(cond)

print("== 16.1 excess log-count of sparse graphs over any N! e^{O(N)} family ==")
alpha = 2
rows = []
for N in [100, 1000, 10000, 100000, 1000000]:
    M = N*(N - 1)//2; K = alpha*N
    lnC = math.lgamma(M + 1) - math.lgamma(K + 1) - math.lgamma(M - K + 1)
    excess = lnC - math.lgamma(N + 1)
    rows.append((N, excess/(N*math.log(N))))
    print(f"    N = {N:>8}:  [ln C(M, 2N) - ln N!]/(N ln N) = {excess/(N*math.log(N)):.5f}   (-> alpha - 1 = 1 as N -> inf; O(N)/(N ln N) corrections)")
cres = [(r - 1)*math.log(N) for N, r in rows]
print("    (excess/(N ln N) - 1) x ln N =", ", ".join(f"{c:.4f}" for c in cres), " -> constant: the O(N) term")
check("excess = (alpha - 1) N ln N + c N with c constant (residual x ln N stable to 2% for N >= 10^4)", abs(cres[-1]/cres[-3] - 1) < 0.02 and abs(rows[-1][1] - 1) < 0.02)

print("== 16.2 a local action on sparse graphs has only O(N) fluctuation budget ==")
rng = np.random.default_rng(12345)
def sample_degrees(N, K):
    # uniform K-edge simple graph on N labelled vertices: sample K distinct pairs
    edges = set()
    while len(edges) < K:
        i, j = rng.integers(0, N, size=2)
        if i != j:
            edges.add((min(i, j), max(i, j)))
    d = np.zeros(N, dtype=np.int64)
    for i, j in edges:
        d[i] += 1; d[j] += 1
    return d
varN = []
for N in [100, 400, 1600, 6400]:
    S = np.array([np.sum((sample_degrees(N, alpha*N) - 4)**2) for _ in range(400)], dtype=float)
    v = S.var(ddof=1)
    varN.append((N, v/N))
    print(f"    N = {N:>5}: <S> = {S.mean():10.1f}   Var(S)/N = {v/N:7.3f}   entropy excess / Var(S) = {(alpha - 1)*N*math.log(N)/v:7.3f}")
vals = [x[1] for x in varN]
print("    (independent-Poisson reference: Var((d-4)^2) = 2*4^2 + 4 = 36 per vertex; the fixed-K constraint lowers it slightly)")
check("Var(S)/N converges to an N-independent value below the Poisson bound 36 (N >= 400 within 15%)", max(vals[1:])/min(vals[1:]) < 1.15 and max(vals) < 36.5)
ratios = [(alpha - 1)*N*math.log(N)/(r*N) for N, r in varN]
check("entropy excess / Var(S) grows with N (like ln N)", all(ratios[i + 1] > ratios[i] for i in range(len(ratios) - 1)))
# direct characteristic function at small N: ln|<e^{i lam S}>| tracks -lam^2 Var_N(S)/2 at the same N (O(N), not O(N ln N))
lam = 0.05
gauss_ratio = []
for N in [20, 40, 80]:
    Ss = np.array([np.sum((sample_degrees(N, alpha*N) - 4)**2) for _ in range(40000)], dtype=float)
    lnphi = math.log(abs(np.exp(1j*lam*Ss).mean())); pred = -lam**2*Ss.var(ddof=1)/2
    gauss_ratio.append(lnphi/pred)
    print(f"    N = {N:>3}: ln|<exp(i 0.05 S)>| = {lnphi:8.4f}   -lambda^2 Var_N(S)/2 = {pred:8.4f}   ratio {lnphi/pred:.3f}   (MC noise floor ln(1/200) = {math.log(1/200):.2f})")
check("weak coupling (lambda = 0.05): ln|<e^{i lambda S}>| = -lambda^2 Var(S)/2 to 10% at each N (Gaussian regime; Var = O(N))", all(abs(r - 1) < 0.10 for r in gauss_ratio))
# strong coupling: S = sum (d-4)^2 = sum d^2 - 8 sum d + 16 N, and sum d^2 = sum d (mod 2) = 2K (mod 2), so S is even
parity_ok = all(int(np.sum((sample_degrees(N, alpha*N) - 4)**2)) % 2 == 0 for N in [20, 40, 80, 160] for _ in range(200))
check("S is always even, so <e^{i pi S}> = 1 exactly: no suppression at lambda = pi (the Gaussian law fails at strong coupling)", parity_ok)
# per-vertex rate for independent Poisson(4) degrees: g(lambda) = -ln|E e^{i lambda (d-4)^2}| vs the Gaussian lambda^2 Var/2
def g_poisson(lam, mu=4.0, dmax=60):
    tot = 0j
    for d in range(dmax):
        tot += math.exp(-mu + d*math.log(mu) - math.lgamma(d + 1))*complex(math.cos(lam*(d - 4)**2), math.sin(lam*(d - 4)**2))
    return -math.log(abs(tot))
for lamb in [0.05, 1.0, math.pi/2, math.pi]:
    print(f"    independent Poisson(4): per-vertex suppression rate g({lamb:.3f}) = {g_poisson(lamb):8.4f}   Gaussian lambda^2*36/2 = {lamb**2*18:8.3f}")
check("strong coupling: the per-vertex rate is O(1) and below the Gaussian value (g(1) = 1.87 < 18)", g_poisson(1.0) < 18)
print("    note: independent Poisson degrees give g(pi) = 8 (each vertex contributes (-1)^d, E = e^-8), but the true fixed-K ensemble has")
print("    sum d = 2K even, so its phase at lambda = pi is exactly 1: a global constraint the per-vertex approximation misses. Strong-")
print("    coupling suppression rates therefore need the true ensemble; the O(N)-not-O(N ln N) conclusion is unaffected.")

print("== 16.3 balanced bilayer (two-level) orders: an N^2 action budget ==")
for N in [20, 100, 1000]:
    pairs = (N//2)**2
    for lamb in [math.pi/2, math.pi]:
        ln_supp = pairs*math.log(abs(math.cos(lamb/2))) if abs(math.cos(lamb/2)) > 0 else -math.inf
        print(f"    N = {N:>5}, lambda = {lamb:.4f}: ln|<e^(i lambda S)>| = {ln_supp:.1f}   vs entropy of the class ln 2^(N^2/4) = {pairs*math.log(2):.1f}")
# brute force: all 2^{(N/2)^2} relation sets across a fixed balanced bipartition (N = 6: 512 orders), sum of e^{i theta N0}
from itertools import product
ok_bf = True
for Nb in [4, 6]:
    pairs_b = (Nb//2)**2
    for theta in [0.7, 2.4, math.pi]:
        Zbf = sum(complex(math.cos(theta*sum(bits)), math.sin(theta*sum(bits))) for bits in product((0, 1), repeat=pairs_b))
        Zcl = (1 + complex(math.cos(theta), math.sin(theta)))**pairs_b
        ok_bf &= abs(Zbf - Zcl) < 1e-9*max(1, abs(Zcl))
check("brute force over all relation sets of a fixed balanced bipartition (N = 4, 6): sum e^{i theta N0} = (1 + e^{i theta})^{N^2/4}", ok_bf)
check("absolute suppression (|class sum| < 1) requires |cos(theta/2)| < 1/2: true at theta = pi (exact zero), false at theta = pi/2", abs(math.cos(math.pi/2)) < 0.5 and abs(math.cos(math.pi/4)) > 0.5)

print("\nCONCLUSION (heuristic, not a theorem): on sparse graphs with alpha > 1 the entropy excess of non-geometric configurations is")
print("(alpha - 1) N ln N, while a fixed-strength local action suppresses a class sum by at most exp(-O(N)) generically (Gaussian at weak")
print("coupling, weaker or none at strong coupling). Boltzmann energy cannot select a lattice family (theorem); a local action's phase")
print("is not expected to. An exponentially bounded configuration space, or a nonlocal action with an N^2 budget against N^2 entropy")
print("(causal sets: balanced layered orders), is what a geometry-selection mechanism needs.")
print("\nV16 OVERALL:", "PASS" if ok_all else "FAIL")
raise SystemExit(0 if ok_all else 1)
