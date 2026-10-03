# Family A — Hamiltonian / quantum-mereology selectors (target Σ)

**Frozen Σ witness pairs** (`SELECTOR_CHALLENGE.md` §2):
- **W-Σ1:** one K, two **isospectral inequivalent local nets** N₁ / N₂ (N = 12; Bridge B1, BI-01: NO BRIDGE).
- **W-Σ2:** two frames on one Pareto front of (H, ψ) (S2-Σ / S2-ΣH: CRITERION-PRICED).
- **W-Σ3:** groupings 32×2 / 16×4 / 8×8 tied in compatible cases (SCOUT-2 §K).

---

## A2 — Locality from the spectrum (Cotler–Penington–Ranard) — *the strongest Σ candidate; treated first*

| C | content |
|---|---|
| C0 | Given only the spectrum of H, a k-local tensor-product structure (on n qudits of dimension d) in which H is k-local is **generically unique** up to unitary equivalence "when such a description exists" |
| C1 | Σ (subsystem / local net) |
| C2 | W-Σ1 (N₁ vs N₂, same spectrum) |
| C3 | spec(H); the locality class (k, n, d); genericity |
| C4 | the unique k-local TPS (when it exists and is generic) |
| C5 | the claimed elimination: the spectrum + "k-local" is less information than a named factorization. **Frozen record:** "H is k-local" is the canonical **RENAMING** of the IR-01 CPR criterion (BI-09). E-01 / review accounting #1 grade this route **CONDITIONAL DERIVATION + CRITERION-PRICED**. The locality class is the selecting criterion, and it is already a supplied canonical item. **FAIL S1** |
| C6 | W-Σ1 is a frozen isospectral pair of local nets (Bridge class: oscillator networks, so not literally CPR's finite-qudit setting; whether it is one of CPR's "special cases … multiple dual local descriptions" is **not checked** against the source text). **On this frozen witness**, spectrum + locality cannot choose N₁ over N₂. CPR's genericity cannot act on a witness that is already given as a twin |
| C7 | **Twin test:** World A = (K, N₁), World B = (K, N₂). Same spectrum, both local of the same class, both admissible. The principle accepts both → **NONSELECTING on the frozen witness** |
| C8 | vary k (locality class) or the qudit dimension d → the output changes: load-bearing (criterion-priced) |
| C9 | invariant under relabelling; fine |
| C10 | consistent with the frozen CONDITIONAL DERIVATION (E-S2). Adds nothing beyond it |
| C11 | genericity theorem, abstract-verified (search) |
| C12 | **KILLED — S1** (locality criterion = a frozen supplied item). C7 / S2 failure on W-Σ1 confirms it |

## A3 — Spectrum-only / unitary-invariant TPS constructions (Loizeau–Sels; Soulas–Franzmann–Di Biagio)

| C | content |
|---|---|
| C0 | construct the subsystem structure from unitary invariants of H (the spectrum) only |
| C1 | Σ |
| C2 | W-Σ1 |
| C3 | spec(H) only (plus, in Loizeau–Sels, an initial state and a TPS among the "minimal ingredients", per snippet) |
| C4 | a TPS defined up to unitary equivalence |
| C5 | Σ is the factorization *relative to the given operator / observables*. A rule whose only input is unitary-invariant data outputs an object defined only up to unitary equivalence, which cannot distinguish two factorizations of the same operator |
| C6 | cannot pick N₁ vs N₂ (same unitary invariants) |
| C7 | **Twin test, structural (source-independent):** any map f(spec H) returns the same output in World A (N₁) and World B (N₂) → **NONSELECTING** |
| C8 | adding a state or a TPS as input (Loizeau–Sels snippet) relocates the selection into that input |
| C9–C10 | n/a (already dead) |
| C11 | snippets only. The Soulas et al. vs Stoica dispute (2025–26) is unresolved, but **C7 kills any purely unitary-invariant selector of Σ regardless of how that dispute ends** |
| C12 | **KILLED — S2** (nonselecting on a frozen twin) |

## A1 — Quasiclassical factorization (Carroll–Singh) [repaired: SSR1-01]

| C | content |
|---|---|
| C0 | an explicit algorithm (primary / source-text verified, owner): (1) for a bipartite decomposition with **fixed dimensions d_A, d_B**, (2) construct a Candidate Pointer Observable, (3) construct **prescribed** initially peaked product states, (4) compute entanglement-growth and pointer-predictability measures, (5) define the **Schwinger Entropy**, and (6) minimize it over candidate factorizations |
| C1 | Σ (the full frozen factorization target) |
| C2 | W-Σ1 / W-Σ2 / W-Σ3 (including the grouping-dimension ties 32×2 / 16×4 / 8×8) |
| C3 | H; fixed (d_A, d_B); the CPO construction; the prescribed candidate-state family; the Schwinger-Entropy definition |
| C4 | the factorization minimizing Schwinger Entropy *within the fixed-dimension class* |
| C5 | **This is more than announcing a preference: a definite algorithm exists.** It is still not an input-free selector for the frozen Σ target: (1) d_A, d_B are fixed before optimization, so part of the factorization class is supplied (cf. frozen W-Σ3, where the grouping dimensions themselves are the open choice); (2) the CPO / state construction uses a prescribed family of candidate states; (3) the worked example is restricted to the quantum-measurement-limit regime; (4) the paper calls the Schwinger Entropy suggestive, with a definition / purview that may need refinement; (5) varying the factor dimensions is explicitly left to future work; (6) no general theorem proves a unique global minimizer for arbitrary H. **FAIL S1 at the frozen Σ target scope**: it conditionally selects within a supplied factorization class and does not eliminate all Σ information |
| C6 | across the frozen dimension-tied witness W-Σ3 it does not act (the dimensions are input) |
| C7 | supporting evidence only: Adil et al. 2026 (LOC) find that **for a fixed Hamiltonian several tensor factorizations can admit quasiclassical descriptions**. This is **not** a proof that Carroll–Singh's exact Schwinger-Entropy objective has degenerate minima |
| C8 | (d_A, d_B) and the candidate-state prescription are load-bearing by construction |
| C11 | Carroll & Singh: **PRIMARY / SOURCE-TEXT VERIFIED** (owner inspection) for the algorithm and its explicit scope |
| C12 | **KILLED — S1** at the frozen Σ target scope (conditional selection within a supplied factorization class) |

## A4 — Observable-induced TPS (Zanardi 2001; Zanardi–Lidar–Lloyd 2004)

| C | content |
|---|---|
| C0 | the accessible observable algebra 𝒜 induces the TPS |
| C1 | Σ / A_partition |
| C3 | 𝒜 |
| C5 | 𝒜 is A_interface / A_partition data. The source itself says the TPS is "relative and observable induced". **Relocation** |
| C8 | vary 𝒜 → the TPS changes (by the source's own statement) |
| C12 | **KILLED — S5** (accessible-observable algebra). A **strengthened nonselection** reading is recorded (TPS is relative, by theorem) |

## A5 — Minimal-scrambling operational mereology (Zanardi–Dallas–Andreadakis–Lloyd 2024) [repaired: SSR1-02]

| C | content |
|---|---|
| C0 | select the subsystem partition (algebra / commutant pair) minimizing a short-time scrambling rate under H, **given operational constraints** |
| C1 | Σ / A_partition |
| C2 | W-Σ2 |
| C3 | H; **a supplied family of operationally admissible algebras / partitions** (in the examples: a supplied family S of candidate subsets, or a supplied adjoint orbit); the scrambling functional |
| C5 | the authors formulate selection as dynamics **plus operational constraints**. The Hamiltonian selects only *relative to the supplied family of operational possibilities*, and that family is A_interface / A_partition information |
| C8 | change the admissible family → the selected partition changes |
| C11 | publisher abstract: PRIMARY VERIFIED (owner); full-text mirror inspected (owner) for the operational-family wording |
| C12 | **KILLED — S5** (the supplied operational family of algebras / partitions = A_interface / A_partition). Supersedes the earlier "KILLED — S1 (objective)" |

## A6 — Space from Hilbert space (Cao–Carroll–Michalakis)

| C | content |
|---|---|
| C0 | geometry / dimension from mutual information between factors of a state |
| C1 | dimension / geometry |
| C3 | "a decomposition of Hilbert space into a tensor product of factors" (assumed); a redundancy-constrained state; a fitting procedure |
| C5 | the target's prerequisite Σ is an input, and so is the state |
| C12 | **KILLED — S5** (TPS + state supplied). Comparable to the frozen GS1 one-way (geometry given net + access) |

## A7 — Stoica non-uniqueness (control)

| C | content |
|---|---|
| C0 | under "Hilbert-space fundamentalism" (only ψ and H), any physically relevant emergent structure is not unique |
| C12 | **NO-GO / NONSELECTION RESULT** (snippet grade; contested by Soulas et al., defended by Stoica). Consistent with W-Σ1 and with the frozen B1 |
