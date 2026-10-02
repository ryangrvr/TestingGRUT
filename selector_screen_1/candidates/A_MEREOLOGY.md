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

## A1 — Quasiclassical factorization (Carroll–Singh)

| C | content |
|---|---|
| C0 | among factorizations ℋ = ℋ_S ⊗ ℋ_E, choose the one maximizing quasiclassical (robust, predictable) dynamics for given H |
| C1 | Σ |
| C2 | W-Σ2 / W-Σ1 |
| C3 | H; factor dimensions; a quasiclassicality functional (purity / pointer-entropy-based); test states; a time window |
| C4 | the "most quasiclassical" factorization |
| C5 | the functional, the dimensions, the test states and the window carry the selecting preference. SCOUT S2-Σ already graded objective-based Σ selection **CRITERION-PRICED** with Pareto conflicts. **FAIL S1** |
| C6 | Adil et al. (PRD 113, 103535, 2026; snippet): **several factorizations of one H admit a quasiclassical description**. So the criterion is not selecting even with its inputs |
| C7 | the twin exists (Adil et al.; frozen Pareto ties) |
| C8 | changing the functional / window changes the output (Zanardi et al. 2024 use a different objective) |
| C12 | **KILLED — S1** (objective supplied). Non-uniqueness evidence noted |

## A4 — Observable-induced TPS (Zanardi 2001; Zanardi–Lidar–Lloyd 2004)

| C | content |
|---|---|
| C0 | the accessible observable algebra 𝒜 induces the TPS |
| C1 | Σ / A_partition |
| C3 | 𝒜 |
| C5 | 𝒜 is A_interface / A_partition data. The source itself says the TPS is "relative and observable induced". **Relocation** |
| C8 | vary 𝒜 → the TPS changes (by the source's own statement) |
| C12 | **KILLED — S5** (accessible-observable algebra). A **strengthened nonselection** reading is recorded (TPS is relative, by theorem) |

## A5 — Minimal-scrambling operational mereology (Zanardi et al. 2024)

| C | content |
|---|---|
| C0 | choose the (algebra, commutant) pair minimizing a short-time scrambling rate under H |
| C3 | H; a family of candidate algebras; a scrambling objective |
| C5 | both the objective and the candidate family are supplied |
| C12 | **KILLED — S1** (objective). Not distinct from A1 in input accounting |

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
