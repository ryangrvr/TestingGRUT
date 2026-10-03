# DA0 CORRECTION LEDGER (additive)

**Rule:** logs are kept as emitted. Where document wording differs, this ledger takes precedence.

## DA0 CHARTER REPAIR 01 (owner review of `bc9d3d2`; applied before any computation)

| ID | item | correction |
|---|---|---|
| **DC1-01** | cycle decomposition is not canonical | J, the cycle space and edge affinities are canonical. Any particular fundamental-cycle decomposition is **BASIS-PRICED / DIAGNOSTIC ONLY**, unless a basis-independent selector is proved. A "canonical cycle" must be invariantly identified from D; otherwise the terminal is **FLUX STRUCTURE ONLY — NO CANONICAL CYCLE**. EP = ½ Σ J log[(π_x G_xy)/(π_y G_yx)] ≥ 0 only for mutually supported edges; one-way edges need an explicit extended convention |
| **DC1-02** | real-part gap is not a sufficient non-normal cut certificate | Use r(λ) = −Re λ. Spectral sets must be conjugation-closed. Report E_N, ‖E_N‖, K_N(Γ) and conditioning. Perturbative stability needs the contour-resolvent product to vanish. Exploding conditioning ⇒ **SPECTRAL SEPARATION PRESENT — NONNORMAL RECOVERY UNSTABLE**, not a recovered macrostructure |
| **DC1-03** | orientation bit vs supplied non-reversible structure | P vs P* is one supplied bit. A and J are part of supplied D and carry more. The positive grade is **CONDITIONAL DERIVATION OF DIRECTED MACROSTRUCTURE FROM SUPPLIED NONREVERSIBLE D**. Test: covariance under P → P*; matched families with the same π and S but different A may differ (**A-PRICED**). R1 becomes a matched pair |
| **DC1-04** | partition recovery and directedness are separate | C1-P (partition without supplied inputs) and C1-D (directed macro-generator, conditional on C1-P) are separate questions. C1-D without C1-P earns no directed endogenous partition. For Q vs Q*: covariant partition, sign-reversed current and affinities, no supplied cycle basis. Final information accounting: derived / supplied / not earned |

## C1 note (process; no repair)

These changes were made **before the logged C1 run**:
- `macro_cycle` crashed on a round-off-level macro current in R1b. A declared computational zero (ZERO_J = 1e-13) and a
  guard requiring a cycle of length ≥ 3 were added.
- The Kato bound was first coded in the crude form ρR₀·δ/(1 − δ), which is O(√(η/g)). It was replaced by the
  residue-expanded form 2eaη + (e + ρa)δ²/(1 − δ), which is O(η/g) and is the one stated in THEOREM C1-T. Both are valid.
- The macro-current error constant was corrected from 2·q_max·m to **4**·q_max·m, as in PROP C1-D. All logged m are 0.
- An audit mapping blind macro-cycle labels to hidden labels, the R1 matched-pair check, and the DC1-02 defect-ring probe
  were added.
- The ad-hoc defect-ring probe that preceded this was run once outside the log. Its values are reproduced in the logged
  probe.

## DA0 C1 REPAIR 01 (owner review of `feb6bc3`)

The C1 terminal is preserved: **NONREVERSIBLE PARTITION DERIVATION**, graded **CONDITIONAL DERIVATION OF DIRECTED
MACROSTRUCTURE FROM SUPPLIED NONREVERSIBLE D — A-PRICED**. The original log is preserved.

| ID | item | correction |
|---|---|---|
| **C1R-01** | rank-K proof in THEOREM C1-T | The ‖E − E₀‖ < 1 argument is withdrawn; the hypotheses do not make the explicit bound < 1 at finite N. Replaced by a **homotopy** G_t = G₀ + tG₁ on the fixed circle |z| = ρ. For all t, ℓ_t ≤ ℓ < ρ and δ_t ≤ δ < 1, so the circle stays in the resolvent set and the rank is constant: rank E₁ = rank E₀ = K. The projector bound is kept for convergence only |
| **C1R-02** | macro-current recovery condition | N4 is replaced by **N4a** (j_min,N = min nonzero \|J_𝔅\|) and **N4b** (q_max,N·m_N / j_min,N → 0). Sufficient specialisation when j_min ≥ c_J η: q_max·η/g² → 0. N3 alone (η/g → 0) does not suffice when g → 0, e.g. η = g^{3/2}. **Audit:** R1a, R2 and R5 satisfy N4a / N4b asymptotically (g = O(1), η ∝ 1/M); **R4 fails both** (j_min/η → 0; q_max·η/g² grows ∝ M) |
| **C1R-03** | zero currents | "Support and sign pattern are recovered" is withdrawn. Under N4b, non-degenerate nonzero entries have their signs and relative magnitudes recovered. Hidden zeros are shown only to converge to zero, and an exact zero is never inferred from a vanishing sequence. Structural zeros by construction (R2's non-adjacent basins; R5's pair (1, 3)) are family-specific facts |
| **C1R-04** | macro EP | General macro-EP recovery is removed from PROP C1-D (log traffic ratios can be ill-conditioned). Kept: blind EP_Q, the numerical values, EP(G*) = EP(G), EP_Q ≤ EP |
| **C1R-05** | time-reversal certification | "G* satisfies N1 – N4 with the same constants" is withdrawn: the decoupled operator of G* need not be G₀*. Kept: J(G*) = −J(G) exactly, and J_Q(G*) = −J_Q(G) for any fixed partition. If G and G* are each certified for the same asymptotic partition, the recovered partitions agree and the currents reverse. The numerical covariance is evidence |
| **C1R-06** | K2 comparator scope | The theorem-like wording is replaced by a comparator observation: in the ideal first-order Markov description with predictively distinct macro-rows, causal states identify the macrostates. Verdict unchanged: **NOT EXCEEDED** |
| **C1R-07** | multi-cycle control | **R5 added**: K = 4, complete macro graph (cycle-space dim 3), supplied superposition of two independent circulations. No spanning tree, cycle basis or decomposition is used. Partition recovered blindly (misclassified 0 at M = 10 … 80); canonical J_Q reported as a matrix; J_Q is not a multiple of any single cycle (5 support edges on 4 vertices); time-reversal covariance holds to ≤ 1.4e-16. Terminal: **NONREVERSIBLE PARTITION DERIVATION + CANONICAL MACRO FLUX; NO CANONICAL CYCLE DECOMPOSITION** |
| **C1R-08** | non-normal instability | Remains an **untriggered possible terminal**; its absence is not claimed in general. N2 stays priced. No pathological family is required before C2 |

**Applied to:** `C1_NONREVERSIBLE_MACROSTRUCTURE.md`, `DA0_ZOOM_OUT_01.md`, `DA0_STATUS.md`. New numerics:
`c1/c1_repair.py`, `c1/c1_repair.log`.

## Owner ruling on `1c3e9ad` (no repair)

- **C1 ACCEPTED AFTER REPAIR 01. C2 APPROVED TO CHARTER.**
- **Wording discipline (standing).** At the logged finite sizes, the N4b sufficient quantity q_max·η/g² is still > 1 for
  R1a and R5 and ≈ 1 for R2 at the largest M. Their certification is **asymptotic, from the analytic 1/M scaling**, not
  from having entered a numerically small-error regime.

## DA0 C2 CHARTER REPAIR 01 (owner review of `92393d1`; before any C2 computation)

| ID | correction |
|---|---|
| **C2R-01** | Bounded description alone does not separate emergent from encoded architecture. **C2-F3′** requires both a DESCRIPTION TEST and an ARCHITECTURE-EXPLICITNESS TEST (no tree / level / basin / block / layer-count / level-indexed coupling / isomorphic latent variable). Otherwise the family is ARCHITECTURE EXPLICIT IN FAMILY, and recovery is graded ARCHITECTURE SUPPLIED BY FAMILY |
| **C2R-02** | GREM / CREM rejected as strong E1 (hierarchy explicit) and moved to the P1 comparator class. E1 is chosen after a recorded audit: **East model** (C2 §6a) |
| **C2R-03** | **C2-F6:** forbidden are supplied macro-units, subject units, subsystem blocks and recovery grains. Allowed but priced are microscopic sites / spins, locality and local rules. Microscopic sites may not be relabelled as emergent architecture |
| **C2R-04** | C2-A1 is split into **A1a** (abstract: K, p_min, C_V, s) and **A1b** (family → s). s = O(η/g) with a K-independent constant is not assumed. K-dependence of every constant is reported |
| **C2R-05** | Primary metric **m_rel = min_σ max_i π(B_i Δ B̂_σ(i))/π(B_i)**, with the permutation chosen on the audit side only. It is bounded explicitly from the idempotent errors. s·K^{3/2} → 0 is a candidate sufficient condition only |
| **C2R-06** | The "gap-free ⇒ maybe K-uniform" inference is withdrawn. Gap-free concerns coefficient separation, not dimension. Whether K^{3/2} is a proof artefact remains open |
| **C2R-07** | Nesting metric ν; no forced nesting, no hierarchical-clustering optimisation; exact numerical nesting is evidence only |
| **C2R-08** | **C2-F7** disorder firewall |
| **C2R-09** | H1 and H2 get dual grades (… and ARCHITECTURE SUPPLIED BY FAMILY). H1 must account for n·λ_s/g → 0 for the 2^n slow sector to separate |
| **C2R-10** | CONJECTURE C2-B is stated only under an explicit locality / rate / temperature class, and is not inferred from the product bound spread ≥ R^d |

## C2 note (process; no repair)

- **Exploratory runs.** The C2 script was run three times while adding (i) the A2b breakdown sweep, (ii) the P1 ζ-scaling
  check, and (iii) the East Δ_HS, L = 10 and modal-configuration diagnostics. The logged run is the final one. Earlier
  outputs were identical on the shared parts.
- **Numerical-diagnostic threshold flagged.** The East depth count classifies a cut as "diverging" when the fitted
  exponent in 1/q exceeds 0.5. This is a numerical diagnostic, not a certificate. Under C0-P4 it is never a modelling
  input.
- **Upper-bound m_rel.** m_rel is computed with a Hungarian (sum-optimal) matching, which upper-bounds the min-max m_rel
  of C2R-05.
- **Improved Lemma G5-3′.** Introduced in C2. It **supersedes** the anticipated s·K^{3/2} sufficient condition in the C2
  charter §3, with s·K^{3/4} under balance. The C2 charter's §3 remains as written (preregistration text).

## C2 SCOPE REPAIR 01 (owner review of `b403594`)

**Owner ruling:**
- **C2-A PASSES:** blind recovery with K_N → ∞ is proved under the stated regularity, with the balanced sufficient
  condition s_N·K_N^{3/4} → 0 (Lemma G5-3′ accepted).
- Growing architecture from a simple rule: **NOT FOUND**.
- Unbounded partition depth without supplied scaling: **NOT FOUND**.
- Intrinsic K-growth obstruction: **OPEN**.
- TRUE COMPRESSION: **0**.
- **C3 not opened.**

| ID | correction |
|---|---|
| **C2S-01** | "No K supplied" is too broad. Correct statement: **no K is supplied to the partition-recovery map once the canonical slow space is certified. The family-level certification of that slow space may use the hidden proof architecture.** The numerical controls used the construction rank to extract the first K modes. Firewall: blind recovery, conditional rank / cut certification (RA3-03) |
| **C2S-02** | The East terminal is evidence, not a theorem: **"No endogenous growing partition architecture found; the tested spectral hierarchy requires q → 0, consistent with known East theory."** It is not claimed that fixed-q East can never produce a partition-like hierarchy for all L |

**Applied to:** `C2_GROWING_ARCHITECTURE.md`, `DA0_ZOOM_OUT_02.md`, `DA0_STATUS.md`. **C2-E2** is preregistered in the C2
charter §6b; it is **not run**.

## C2S-03 / C2S-04 and C2 CLOSURE (owner review of `c394e03`)

| ID | correction |
|---|---|
| **C2S-03** | The registered E2 sequence 2 × 2 × L (L = N/4) keeps both transverse widths fixed. It is a finite-width spin-glass tube / ladder, effectively one-dimensional in the limit, not the cubic L × L × L sequence whose T_c ≈ 1.10 (finite-size scaling up to L = 40) was cited. T = 0.5 or 0.8 therefore does **not** place this sequence inside the 3D spin-glass phase. A positive result would show only finite-width / tube structure, a negative result would not rule out 3D architecture, and an indeterminate result adds nothing. Recorded: **E2 NOT RUN — ACCESSIBLE EXACT-SPECTRAL GEOMETRY DOES NOT PRESERVE THE TARGET 3D THERMODYNAMIC LIMIT.** It is **not** recorded as "E2 NO GROWING PARTITION ARCHITECTURE FOUND", and **not** counted as a failed physical test |
| **C2S-04** | "Global Z₂ produces exact paired sectors" is withdrawn. At finite N the heat-bath Glauber chain is irreducible; global spin flip commutes with the generator and gives parity structure and configuration pairing, not two disconnected Markov sectors. A K = 2 structure explainable solely by this supplied symmetry is non-novel for C2 |

**C2 FINAL OWNER RULING:** C2 COMPLETE — GROWING BLIND RECOVERY PROVED CONDITIONALLY; ENDOGENOUS GROWING ARCHITECTURE
FROM A SIMPLE FIXED RULE NOT FOUND; INTRINSIC K-GROWTH OBSTRUCTION OPEN. C3 not opened.

**Applied to:** `C2_CHARTER.md`, `C2_GROWING_ARCHITECTURE.md`, `DA0_ZOOM_OUT_02.md`, `DA0_STATUS.md`. Created
`C2_FINAL_HANDOFF.md`.

## C3 CHARTER REPAIR 01 (owner review of `f266be2`; charter only, no computation)

**Owner ruling:** C3 CHARTER CONDITIONALLY ACCEPTED — PRIMARY MODEL ACCEPTED, REPAIR BEFORE RUN. The pre-repair boundary
is `f266be2f6432e4c8de6c73a183ea5135fe816a8e`.

| ID | correction |
|---|---|
| **C3R-01** | "n = 1 ⇒ no architecture" and "architecture needs n ≠ 1" are **withdrawn**. At T = T′ the σ-marginal is uniform, but J \| σ ~ N(s_b/μ, T/μ) keeps bond / loop correlations (plaquette ⟨ΠJ_b⟩ = 1/μ⁴). CPS mean-field has an ordered q > 0 phase for n ≤ 2, including n = 1. A1′ is renamed **EQUILIBRIUM ADAPTIVE-COUPLING COMPARATOR**; A0 / A1 remain the nulls. For T ≠ T′, the two reservoirs are supplied, and entropy production / heat flow is to be calculated at finite ε, not asserted |
| **C3R-02** | Primary analytic object: **𝓕_T(J) = (μ/2)ΣJ_b² − T log Z_T(J)**, with dJ/dt = −ε∇𝓕_T + slow noise. Stationary points depend on (T, μ), **not on T′**. Hessian μδ − (1/T)Cov_J(s_b, s_b′). Exact linear instability of J = 0 at **T = 1/μ** |
| **C3R-03** | **C3-F2:** exact local gauge symmetry σ_i → η_iσ_i, J_ij → η_iη_jJ_ij. Architecture is counted modulo all exact symmetries, using gauge-invariant observables. No gauge fixing may generate K-growth. Generated gauge-invariant frustration is distinguished from gauge copies |
| **C3R-04** | The pairwise median-PR criterion is withdrawn: it certifies independent product bits (PR ~ N/2). It is replaced by **PR_edge on elementary transitions** (independent bits → O(1)), plus a non-vanishing contrast **liminf S_ab > 0** |
| **C3R-05** | **C3-M1 resolved model-specifically.** Detector = symmetry-inequivalent stable minima of 𝓕_T joined by canonical minimum-barrier saddles, with diverging barriers / exit times. No clustering, PCCA, chosen K, lag-time partitions or fitted macrostate count. C3-A1 requires growth **and** collective edges **and** growing barriers |
| **C3R-06** | **Analytic preflight first** (landscape → symmetries → stability → loop expansion → bond-local vs collective), on separate approval. No large-lattice simulation before review |
| **C3R-07** | **CPS motivates the mechanism; it does not validate the local lattice theory.** The analysed CPS model is mean-field / infinite-range; the nearest-neighbour version is a new construction |

**C3-B stays closed. No C3 code, no C3 calculation.**

## C3 CHARTER REPAIR 02 (owner review of `38febe1`; before the analytic preflight)

**Owner ruling:** Repair 01 accepted; primary model remains accepted. **Analytic preflight APPROVED after Repair 02. No
simulation.**

| ID | correction |
|---|---|
| **C3R-08** | A plaquette (or any bounded loop) is the first gauge-invariant **inter-bond interaction**, not C3 collectivity. Its support is bounded, so its participation is O(1), and it does not satisfy C3-A1. The preflight question is whether bounded loop interactions bootstrap into extended, symmetry-inequivalent metastable structures |
| **C3R-09** | The primary dimension is **fixed at d = 3 before results**. d = 1, 2 are analytic controls only |
| **C3R-10** | Mandatory ferromagnetic-envelope test: Z_T(J) ≤ Z_T(\|J\|), so 𝓕_T(J) ≥ 𝓕_T(\|J\|), with strictness. It constrains global minima only. The key question is frustrated **local** minima with diverging barriers |
