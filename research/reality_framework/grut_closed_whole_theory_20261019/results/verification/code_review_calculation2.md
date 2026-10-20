# Adversarial review: `marginal_compatibility.py` (GRUT closed-whole calc 2)

Target: `/home/user/TestingGRUT/research/reality_framework/grut_closed_whole_theory_20261019/code/marginal_compatibility.py`
(sha256 `1aad4b06...eb2cec`, matches `results/marginal_compatibility.json` provenance) and that JSON.
All review scripts and raw outputs: `/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/review2/`
(`check_a_ptrace_axis.py`, `check_b_c3_analytic.py`, `check_c_c4_thresholds.py`, `check_c2_clarabel_panic.py`,
`check_d_c5_lp.py`, `check_e_notes.py`, outputs `out_*.txt`, reproduction run under `repro/`).
Nothing under `/home/user/TestingGRUT` was modified.

## Verdict

**sound_with_minor_issues.** Every computed number in the JSON is correct and independently reproduced
(second solver, analytic closed forms, explicit optimal states). Defects are in prose: one false number inside the
C2 "analytic_proof" field (fidelity 1/2 should be 1/4; conclusion unaffected), one false pair of numbers and a
misleading "W-type" hint in a code comment, loose "Werner" naming for the isotropic state, and a cycle note that can
be read as stronger than the demonstrated fact. One robustness caveat about relying on bare solver status.

## Reproduction

```
cp marginal_compatibility.py <scratch>/repro/code/ && python3 marginal_compatibility.py
```
Output JSON identical to the repo JSON except `wall_seconds` (0.28 vs 0.29); script sha256 identical.

## (a) cvxpy `partial_trace` axis convention -- CONFIRMED

`check_a_ptrace_axis.py`: random 3-qubit density matrix vs an explicit index-loop partial trace (no reshape tricks),
plus a product state A=|0>, B=|+>, C=|+i> whose three reductions are mutually distinguishable.

| call | equals numpy keep | deviation | deviation from the other two reductions |
|---|---|---|---|
| `axis=2` | [0,1] = AB (C traced) | 0.0 | 1.65e-1, 8.6e-2 |
| `axis=0` | [1,2] = BC (A traced) | 0.0 | 1.65e-1, 1.63e-1 |
| `axis=1` | [0,2] = AC (B traced) | 0.0 | 8.6e-2, 1.63e-1 |

The script's own numpy `ptrace(rho, keep, n)` agrees with the loop implementation (0.0) and with the product-state
expectations (<=1.1e-16, including `keep=[1]` -> B). `tr_C/tr_A/tr_B` and the C3/C6 numpy reductions are mutually
consistent, with A = most significant qubit (index 4a+2b+c), as the C5 LP also assumes.

## (b) C3: max F_AB + F_BC = 3/2, symmetric max = 3/4 -- CONFIRMED (analytically and by random restarts)

`check_b_c3_analytic.py`. The objective is Tr(M rho) with M = Phi+_AB (x) I_C + I_A (x) Phi+_BC, so the maximum over
states is lambda_max(M), attained on a pure state.

- `eigh(M)`: eigenvalues {0,0,0,0, 1/2,1/2, 3/2,3/2}; lambda_max = 1.5 exactly, multiplicity 2.
- Reason: M = P + Q, two rank-2 projectors; the overlap matrix between their ranges is (1/2) I_2 (cosines of both
  principal angles = 0.5, computed), so lambda_max = 1 + 1/2 = 3/2 and every top eigenvector has
  <P> = <Q> = (1 + cos theta)/2 = 3/4. Hence the symmetric optimum is exactly 3/4 (it cannot exceed (3/2)/2).
- Explicit optimisers: (2|000> + |110> + |011>)/sqrt6 and (2|111> + |001> + |100>)/sqrt6, each giving
  F_AB = F_BC = 0.75 (computed).
- Random-restart L-BFGS over pure states (16 real parameters, 300 starts): max = 1.5000000000, all 300 starts within
  1e-6; at the best point (F_AB, F_BC) = (0.75000002, 0.74999998). Symmetric surrogate (soft-min, 300 starts):
  max min(F_AB,F_BC) = 0.75000000.
- Comparison states: W gives (1/6, 1/6), sum 1/3 (and at most 2/3 per pair with the best rotated target, sum 4/3 < 3/2);
  GHZ gives (1/2, 1/2); Phi+_AB (x) |0>_C gives (1, 1/4), sum 5/4.

JSON values 1.49999999899, 0.750000000106, [0.7499999995, 0.7499999995] are correct to solver tolerance.

## (c) C4: chain threshold 2/3, triangle threshold 1/3 -- CONFIRMED (SCS + CLARABEL + closed forms + explicit states)

`check_c_c4_thresholds.py`. Note: plain feasibility SDPs are a poor bisection oracle here -- SCS returns
`optimal_inaccurate` on both sides of the boundary (e.g. chain p=0.6657 and 0.6677 both `optimal_inaccurate`), and
CLARABEL **panics** (pyo3 `PanicException: Eigval error: Eigen(1)`, not a `SolverError`) deterministically at
p = 0.6666669846 and p = 0.6666666865 (`check_c2_clarabel_panic.py`, 3/3 repeats). So the bisection was run on the
always-feasible surrogate v(p) = min over states of the max-abs marginal violation (v(p)=0 iff compatible; v convex,
feasible set an interval).

| problem | solver | bisection on v(p) <= 1e-6 | direct max p | target |
|---|---|---|---|---|
| chain {AB,BC} | SCS (eps 1e-9) | 0.6666691 < p* <= 0.6666694 | 0.666666667 (err 7e-12) | 2/3 |
| chain | CLARABEL | same bracket | 0.666666674 (err 7e-9) | 2/3 |
| triangle {AB,BC,CA} | SCS | 0.3333373 < p* <= 0.3333375 | 0.333333333 (err 1e-11) | 1/3 |
| triangle | CLARABEL | same bracket | 0.333333336 (err 2e-9) | 1/3 |

(The ~3e-6 bracket offset is the v<=1e-6 tolerance: v grows ~0.375/unit p past the chain threshold, 0.25/unit for
the triangle, e.g. v(2/3+1e-3) = 3.751e-4, v(1/3+1e-3) = 2.500e-4.) CLARABEL's returned optimal states have
marginal deviations 6.5e-15 / 5.4e-15 from W(p*) and min eigenvalue -1.5e-9 / -8.6e-10. At p = 2/3 exactly the
chain admits a rank-2 extension (eigenvalues 0.5, 0.5); min violation 2.5e-11.

Analytic / literature cross-checks:
- **Chain = symmetric extension of W(p) (B shared; W is swap-symmetric).** Chen, Ji, Kribs, Luetkenhaus, Zeng,
  *Symmetric extension of two-qubit states*, PRA 90, 032318 (2014), arXiv:1310.3530 -- **fetched**; the abstract
  states the iff condition Tr(rho_B^2) >= Tr(rho_AB^2) - 4 sqrt(det rho_AB). For W(p): Tr rho_B^2 = 1/2,
  Tr rho_AB^2 = (1+3p^2)/4, det = (1+3p)(1-p)^3/256, so the condition is sqrt((1+3p)(1-p)^3) >= 3p^2 - 1, with
  equality exactly at p = 2/3 ((1+3p)(1-p)^3 = 1/9). Computed g(p): +2.2e-2 at 0.65, +2e-16 at 2/3, -4.5e-2 at 0.7.
- Johnson & Viola, *Compatible quantum correlations: on extension problems for Werner and isotropic states*,
  PRA 88, 032323 (2013), arXiv:1305.1342 -- **fetched (abstract only)**; it states iff conditions for 1-n sharability
  of isotropic states but the closed form is not in the abstract. My recollection F <= (1/d)(1 + (d-1)/k) gives
  F <= 3/4, i.e. p <= 2/3 for d=2, k=2, consistent, but this formula was **not** verified from the fetched page.
- **Triangle: independent analytic derivation (no reference checked).** The constraint set is S_3-invariant, so
  WLOG rho commutes with permutations: rho = w4 rho_4 (+) w2 (rho_2 (x) I/2) on the S=3/2 irrep and the S=1/2 pair.
  On S=1/2, S_a^2 = 1/4 identically, so each pair correlator is -1/3 regardless of rho_2. On S=3/2 the pair
  correlators are (4<S_a^2> - 3)/6; the Werner constraints <XX> = <ZZ> = -<YY> = p with 2u + v = 15/4 give
  p = (4u-5)/(13-4u), increasing in u = <S_x^2> = <S_z^2>, and u <= 7/4 because S_y^2 >= 1/4 on spin 3/2. Hence
  p_max = 1/3, attained at w4 = 2/3. The explicit state
  rho* = (1/3) P_{S=1/2}/4 + (2/3) (|3/2, m_y=+1/2><.| + |3/2, m_y=-1/2><.|)/2 was built numerically:
  trace 1, min eig -3e-17, and max|rho*_XY - W(1/3)| = 1.1e-16 for all three pairs. Ingredient checks: on the
  m_y=+-1/2 mixture <S_x^2> = <S_z^2> = 1.75, <S_y^2> = 0.25, pair correlators (2/3, -1/3, 2/3); on P_{1/2}/4
  (-1/3, -1/3, -1/3).
- Literature-adjacent sanity: the SU(2)-twirlable Werner-proper triangle (p Psi- + (1-p)I/4 on all pairs) gives the
  same number from Sum_{i<j} <sigma_i.sigma_j> = 2(S(S+1) - 9/4) >= -3 -> p <= 1/3; the two triangle problems are
  not locally equivalent (see defect 3) but coincide numerically.

## (d) C5: LP construction, classical chain gluing, even-cycle parity -- CONFIRMED

`check_d_c5_lp.py`. The compiled equality system from `lp.get_problem_data(cp.CLARABEL)` (printed in `out_d.txt`)
has the normalisation row, 12 marginal rows and 8 nonnegativity rows. The 12 marginal rows equal, as a set, the
hand-built rows for P(A=a,B=b), P(B=a,C=b), P(C=a,A=b) with 0.5 off-diagonal / 0 diagonal (`True (12 rows each)`),
so the three `idx(...)` permutations in lines 181-183 do encode the AB, BC and CA marginals. Independent solvers:
HiGHS via cvxpy `infeasible`; scipy `linprog(method="highs")` on the hand-built system: status 2 infeasible.
Logical check: no bit triple satisfies a!=b, b!=c, c!=a (empty list).

Chain gluing: over 20 random consistent (P_AB, P_BC) pairs, P(a,b,c) = P(a,b)P(c|b) and the general formula
P(a,b)P(b,c)/P(b) reproduce both marginals to 1.1e-16, and the HiGHS chain LP is feasible for all 20; with
disagreeing P(B) (0.6/0.4 vs 0.5/0.5) it is infeasible (status 2). The script's chain construction is correct
(its P_BC is built from the same P(b), so overlap agreement holds by construction, which is the hypothesis of the
gluing theorem, not a loophole).

Parity point: 4-cycle A-B-C-D-A with **four** anticorrelated edges is **feasible** (HiGHS status 0; support
{0101: 0.5, 1010: 0.5}); the same 4-cycle with three anticorrelated + one correlated edge is **infeasible**; the
3-cycle with two anticorrelated + one correlated edge is **feasible**. So the obstruction is odd frustration around
a cycle, which a cycle permits but does not force.

## (e) Wording -- defects and assessments

1. **C2 `analytic_proof` (line 110-112, in the JSON): "rho_BC = (I/2)(x)rho_C, whose fidelity with Phi+ is 1/2 != 1"
   is numerically wrong.** <Phi+|(I/2 (x) rho_C)|Phi+> = Tr(rho_C)/4 = **1/4** for every rho_C (computed 0.25 for
   rho_C = |0><0|, I/2, random; `check_e_notes.py` item 1). The conclusion (!= 1, hence no extension) stands. 1/2 is
   the maximum fidelity of a *product* state with Phi+, which is a different fact. Fix: "... whose fidelity with Phi+
   is Tr(rho_C)/4 = 1/4 != 1".
2. **Line 129 comment: "the product of a singlet with a free qubit gives (1, 1/2)" is wrong** -- Phi+_AB (x) rho_C
   gives (F_AB, F_BC) = (1, 1/4) for any rho_C (computed). Same comment: "analytic candidate: the W-type / symmetric
   three-qubit states" is misleading -- the W state gives (1/6, 1/6), sum 1/3, nowhere near the optimum; the
   optimisers are (2|000>+|110>+|011>)/sqrt6 and its bit-flip (top eigenspace of M, both with F = 3/4). Code comment
   only, not in the JSON. Fix the comment to name the actual optimiser and (1, 1/4).
3. **"Werner" naming (C4, `werner()`, JSON note).** W(p) = p Phi+ + (1-p) I/4 is the *isotropic* state (invariant
   under U (x) U*, verified 5.6e-17; not under U (x) U, deviation 0.243); the Werner state proper is
   p Psi- + (1-p) I/4 (U (x) U invariant). Pairwise they are locally equivalent ((I (x) Y) W(p) (I (x) Y) = Werner
   proper, deviation 0), and for the **chain** a Y on B maps one problem onto the other, so 2/3 holds for both. For the
   **triangle** no consistent local unitary does this (A-B, B-C, C-A cannot all be U (x) U*), so the two triangle
   problems are distinct; both happen to give 1/3 (section (c)). The "entangled iff p > 1/3" statement is true for
   both. Fix: call it "isotropic (Phi+-Werner) state" once and note the chain equivalence.
4. **C5 `cycle_note`: "overlap agreement on a CYCLE is not sufficient"** is correct in the existential sense (on any
   non-acyclic hypergraph some pairwise-consistent family fails: Vorob'ev 1962 / running-intersection property; for
   2-element edges acyclic = forest, so "acyclic (tree)" is right), and "fails already for a chain in the quantum case"
   is what C2 shows. But the demonstrated mechanism is odd frustration, not cyclicity: the even 4-cycle with four
   anticorrelated edges extends. Suggested addition: "(the obstruction is an odd number of anticorrelated edges around
   the cycle; the same 4-cycle with four anticorrelated edges extends, with 1/2 on 0101 and 1/2 on 1010)". The note
   does not say chordal; no chordal claim to assess.
5. **`chain_note` "specifically quantum (monogamy), not a bookkeeping failure"** -- accurate: classically any two
   marginals agreeing on the shared variable extend (verified by LP for 20 random instances), the quantum Bell chain
   does not (C2, infeasible with both CLARABEL and SCS). "Necessary vs sufficient" language in C1, C2, C6 is accurate.
6. **C1 note "partial trace is a conditional expectation"** -- loose: the trace-preserving conditional expectation is
   X -> Tr_C(X) (x) I/d_C; composition of restrictions is simply Tr_CD = Tr_C o Tr_D. Also the reported 0.0 for
   Tr_C(Tr_D rho) vs Tr_CD rho compares two bitwise-identical floating-point computations (same trace order); the
   other order gives 2.8e-17. Harmless (the statement is a tautology) but the 0.0 is not evidence.
7. **Robustness (not a correctness defect).** All "infeasible" conclusions rest on CLARABEL's certificate with no
   second solver and no exception handling; CLARABEL panics (uncatchable as `SolverError`) within ~3e-7 of the C4
   boundary. Here every status was confirmed: C2 and C6 `infeasible` with SCS (eps 1e-9) and CLARABEL, C5 LP
   infeasible with HiGHS, C4 via min-violation values. Suggest recording the minimum marginal-violation value
   (a number, not a status) next to each infeasibility claim, and/or a second solver's status in the JSON.

## Minimal fixes (text only; no numbers change)

- line 111: `whose fidelity with Phi+ is 1/2 != 1` -> `whose fidelity with Phi+ is Tr(rho_C)/4 = 1/4 != 1`.
- line 129: `# analytic candidate: the W-type / symmetric three-qubit states; also the product of a singlet with a free qubit gives (1, 1/2)`
  -> `# optimum = top eigenspace of Phi+_AB(x)I + I(x)Phi+_BC (lambda_max = 3/2, F_AB = F_BC = 3/4), e.g. (2|000>+|110>+|011>)/sqrt6; W gives (1/6,1/6); Phi+_AB (x) rho_C gives (1, 1/4)`.
- C4 note: `W(p) = p Phi+ + (1-p) I/4` -> `W(p) = p Phi+ + (1-p) I/4 (isotropic state; locally equivalent to the Werner state pairwise, and for the chain globally via Y on B)`.
- cycle_note: append the odd-frustration / even-4-cycle sentence from item 4.
- optional: add `"second_solver_status"` (SCS) or `"min_marginal_violation"` fields for C2, C4, C5, C6.

## Confirmed claims (with the deciding check)

- cvxpy `partial_trace` axis k traces out subsystem k; `tr_C` keeps AB, `tr_A` keeps BC, `tr_B` keeps AC (dev 0.0 vs loop implementation).
- C3: max F_AB + F_BC = 3/2 = lambda_max(M), multiplicity 2; symmetric max 3/4; 300/300 random pure-state restarts reach 1.5.
- C4: chain p* = 2/3 (SCS 0.666666667, CLARABEL 0.666666674, Chen et al. criterion root exactly 2/3, reference fetched); triangle p* = 1/3 (SCS 0.333333333, CLARABEL 0.333333336, explicit rho* with marginals W(1/3) to 1e-16, spin-decomposition bound 1/3).
- C5: LP rows match the hand-built AB/BC/CA anticorrelation system; infeasible under CLARABEL, HiGHS (cvxpy) and scipy HiGHS; chain gluing exact for 20 random instances; even anticorrelated 4-cycle feasible, odd-frustrated 4-cycle infeasible.
- C2, C6: infeasible with both CLARABEL and SCS; consistent product control feasible with both.
- JSON reproduces bit-for-bit (except wall time) from the script; sha256 matches.
