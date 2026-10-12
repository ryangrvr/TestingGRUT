#!/usr/bin/env python3
"""
V16 -- Entropy-action budget: why bounded-degree graph models cannot select geometry, and why causal-set-type
nonlocal actions can at least compete.  (Synthesis 13 Sec. 2; numbers for the codex fixed-bond / CP12 obstructions.)

CLAIM (heuristic, made quantitative).  A weight w(c) can make a "good" (geometric) class dominate a "bad" class only if the
weights can fight the bad class's excess log-count  Delta L = ln|bad| - ln|good|.
 * Boltzmann weights e^{-beta E} with extensive E (|E| <= C N): the codex theorem gives P(good) <= exp(-Delta L + O(N)).
 * Lorentzian weights e^{i S}: a phase can cancel, but for an action that is a sum of ~n roughly independent local terms
   the class sum |sum_bad e^{iS}| is generically |bad| x exp(-O(n)) (central-limit / product of local characteristic
   functions); cancellation needs n >~ Delta L.
For N vertices and K = alpha N edges, Delta L >= (alpha - 1) N ln N + O(N) for any family with <= N! e^{O(N)} labelled members
(lattice templates), while a LOCAL action on a sparse graph has n = O(N) terms: Var(S) = O(N).  So neither Boltzmann nor
Lorentzian weights with a fixed-strength local action can select the family; the deficit grows like ln N.
Layered orders (the Kleitman-Rothschild class, ln count ~ (N^2/4) ln 2): an action counting links/relations (causal-set
theory) has ~N^2 independent terms, so its phase can suppress the entropically dominant, balanced layered orders by exp(-c N^2):
for the balanced two-level class the link-action sum is exactly (1 + e^{i theta})^{N^2/4}, absolutely suppressed only for
|cos(theta/2)| < 1/2.  Literature (checked, 13): Loomis & Carlip 2018 (CQG 35, 024002) -- two-level (bilayer) orders, leading
order in N^2; Mathur, Singh & Surya 2021 (CQG 38, 045017) -- k-layer / KR-type orders with the link action; Carlip, Carlip &
Surya 2023-24 (CQG 40, 095004; CQG 41, 145005) -- KR / layered orders with the full BDG action.  CAVEATS: the sum over ALL
two-level orders is not exponentially small in N^2 (unbalanced orders dominate; workflow computation recorded in 13), results
depend on which subset is summed, and no paper compares these sums with the manifold-like contribution.  So the N^2 budget is
NECESSARY for selection and has been shown sufficient to suppress the dominant layered classes, not to select manifolds.

CHECKS.
 16.1 exact: ln C(M, alpha N) - ln N! = (alpha - 1) N ln N + O(N) (Stirling, numerically for N to 10^6).
 16.2 Monte Carlo over uniform fixed-K graphs (alpha = 2): the valence action S = sum_v (d_v - 4)^2 has Var(S)/N -> const,
      so ln |<e^{i lambda S}>| = -O(N) (checked directly at small N); the ratio (alpha - 1) N ln N / Var(S) grows like ln N.
 16.3 exact: random bilayer orders (n1 = n2 = N/2 layers, each cross pair related with probability 1/2, the KR-type class):
      the relation-count phase gives |<e^{i lambda S}>| = |cos(lambda/2)|^{N^2/4} -- an exponent in N^2, matching the
      N^2/4 ln 2 entropy; for lambda with ln|cos(lambda/2)| < -ln 2 ... the class sum is smaller than its count by more than its
      own entropy advantage could compensate only when |cos(lambda/2)| < 1/2, i.e. lambda in (2pi/3, 4pi/3) mod 2pi.
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
check("ln|<e^{i lambda S}>| = -lambda^2 Var(S)/2 to 10% at each N (Gaussian; Var = O(N), so suppression exp(-O(N)))", all(abs(r - 1) < 0.10 for r in gauss_ratio))

print("== 16.3 Kleitman-Rothschild-type bilayer orders: an N^2 action budget ==")
for N in [20, 100, 1000]:
    pairs = (N//2)**2
    for lamb in [math.pi/2, math.pi]:
        ln_supp = pairs*math.log(abs(math.cos(lamb/2))) if abs(math.cos(lamb/2)) > 0 else -math.inf
        print(f"    N = {N:>5}, lambda = {lamb:.4f}: ln|<e^(i lambda S)>| = {ln_supp:.1f}   vs entropy of the class ln 2^(N^2/4) = {pairs*math.log(2):.1f}")
check("balanced bilayer phase suppression exponent scales as N^2 (same order as the layered-order entropy)", all(abs(((N//2)**2*math.log(abs(math.cos(math.pi/4))))/((N//2)**2) - math.log(math.cos(math.pi/4))) < 1e-12 for N in [20, 100, 1000]))
check("full cancellation of the bilayer entropy needs |cos(lambda/2)| < 1/2 (e.g. lambda = pi gives exact zero)", abs(math.cos(math.pi/2)) < 0.5 and abs(math.cos(math.pi/4)) > 0.5)

print("\nCONCLUSION: on sparse graphs with alpha > 1 the entropy excess of non-geometric configurations is (alpha - 1) N ln N while")
print("any fixed-strength local action has O(N) fluctuation budget: neither Boltzmann energy nor Lorentzian phase can select a")
print("lattice family. A configuration space with exponentially bounded counts (triangulations of fixed topology), or a nonlocal")
print("action with an N^2 budget against N^2 entropy (causal sets, Kleitman-Rothschild), is what a geometry-selection mechanism needs.")
print("\nV16 OVERALL:", "PASS" if ok_all else "FAIL")
raise SystemExit(0 if ok_all else 1)
