# GRAVITY-SCOUT-1 CORRECTION LEDGER (additive)

**Rule** (carried from QFT-SCOUT-1):
- Original script and log outputs are kept as emitted.
- Where document wording is corrected, this ledger takes precedence.
- No repair creates a prediction or TRUE COMPRESSION.

## GRAVITY REPAIR 01 — G1 scope and theorem accounting (owner ruling after GRAVITY_ZOOM_OUT_01)

**Applied to:** `G1_DMIN_RESULT.md`, `GRAVITY_ZOOM_OUT_01.md`, `GRAVITY_SCOUT_1_STATUS.md`.
**Effect:** G1 is accepted **provisionally** after this repair.

| ID | correction | applied as |
|---|---|---|
| **GR1-01** | Inclusion-level terminal | The single "d_min^incl = NO CONSTRAINT" is withdrawn. It is split into three regimes:<br>**fixed-background LCQFT:** d_min^incl = 0 under the Fewster / split assumptions;<br>**perturbative gravity (Donnelly–Giddings, O(κ), arbitrary U_ε):** NO POSITIVE COLLAR SCALE DERIVED; SUBSYSTEM NOTION REPLACED / REPRICED;<br>**fine-grained gravity (Raju, scoped examples):** ORDINARY d_min^incl MAY BE NOT APPLICABLE / BLOCKED, not generalized to all quantum gravity.<br>**Overall:** no positive universal minimum splitting length; in dynamical gravity the ordinary split structure may be replaced or may fail rather than acquire a finite d_min |
| **GR1-02** | Preparation bookkeeping | The bound acts on (H_cross = 0, d, R, N, state class). It is recorded as a **PREPARATION-CLASS CONSTRAINT** on the H_cross / H_epoch side, not on Fewster's primitive. Terminal: **CONSTRAINED-NONUNIQUE — CONDITIONAL ON PREMISE L_all — HEURISTIC**. L_all: every admissible product preparation carries at least the required gravitational mass within R + O(d). The numerics establish localization only for the minimum-energy Gaussian product state |
| **GR1-03** | Necessary, not sufficient | "Every d ≥ d\* remains admissible" is replaced by **"d ≥ d\* is not excluded by this particular collapse test."** Preparability, stability and the absence of other obstructions are not proved |
| **GR1-04** | Exact lattice proposition wording | The proof is now the covariance argument (product ⇒ zero centered cross covariance ⇒ Gaussianization ⇒ displacement (K > 0) ⇒ covariance-level T-averaging ⇒ Gaussian realization). Averaging *states* is not claimed to preserve productness. The same covariance-level discipline is applied to the species-additivity and rotation-averaging arguments. The `g1/g1_dmin.py` docstring predates this and is superseded |
| **GR1-05** | Collapse coefficient scope | **Robust candidate scaling:** d\* ~ (N ℓ_P² R)^{1/3}. **Model-dependent prefactor:** (8π α₃)^{1/3}, not theorem-grade. It uses semiclassical backreaction, the flat collar estimate, L_all, and the neglect of binding and stress fluctuations. Curvature ≫ ℓ_P⁻² does not by itself guarantee the semiclassical Einstein equation if stress-tensor fluctuations are anomalously large |
| **GR1-06** | Crossed-product scope | **Witten:** a specific emergent large-N black-hole setting, III₁ → II∞. **CLPW:** a dressed de Sitter static patch, II₁. Neither is a universal statement about local regions. "Type II" is **not** taken to imply G1-P1 for geometric regions: G1-P1 applies to an exact (M, M′) pair, or geometrically once the duality / commutant identification is established |
| **GR1-07** | Semiclassical inclusion wording | "For each fixed semiclassical solution the inclusion-level infimum is unchanged" is **withdrawn**. Replacement: **no tested semiclassical-backreaction argument derives a positive universal collar scale at the inclusion level** |
| **GR1-08** | Bousso route | I(A:B)/2 ≤ A/(4ℓ_P²) for the collar's vacuum mutual information is not the Bousso conjecture itself. Grade: **HEURISTIC APPLICATION OF A CONJECTURE**. It supplies no G1 evidence |
| **GR1-09** | Novelty wording | No novelty claim for the (ℓ_P² R)^{1/3} scaling: Károlyházy / Ng–van Dam analogues exist (gr-qc/9906003). α₃ and the product-energy route are not exhaustively novelty-searched. Braunstein–Pirandola–Życzkowski is cited only as conceptual "energetic curtain" precedent, not as a derivation of α₃ or d\* |
| **GR1-10** | Source upgrades | see table below |

**Source grades after GR1-10:**

| source | grade |
|---|---|
| Donnelly & Giddings, PRD 98, 086006 (2018), arXiv:1805.11095 | **PRIMARY / SOURCE-TEXT VERIFIED** (owner): the first-order gravitational splitting; the arbitrary-U_ε setup |
| Raju, "Failure of the split property in gravity and the information paradox", arXiv:2110.05470 (2021/2022) | **PRIMARY / SOURCE-TEXT VERIFIED** (owner): the scoped failure of the split property; holography of information |
| Witten, "Gravity and the crossed product", arXiv:2112.12828, JHEP 10 (2022) 008 | **PRIMARY ARXIV ABSTRACT VERIFIED** (II∞ crossed product) |
| Chandrasekaran–Longo–Penington–Witten, arXiv:2206.10780 (2023) | **PRIMARY ARXIV ABSTRACT VERIFIED** (II₁ de Sitter static patch) |
| Hayward, "Gravitational energy in spherical symmetry", PRD 53, 1938 (1996) | **PRIMARY JOURNAL ABSTRACT VERIFIED** (spherical Misner–Sharp trapped / marginal / untrapped criterion) |
| Braunstein–Pirandola–Życzkowski, arXiv:0907.1190 | **PRIMARY ARXIV RECORD VERIFIED** (terminology / qualitative precedent only) |
| Ng & van Dam, gr-qc/9906003 | SOURCE LOCATED, TEXT NOT RE-READ (abstract: distance uncertainty (l ℓ_P²)^{1/3}); novelty comparison only |

### Owner decisions recorded with REPAIR 01

1. **Preparation-level narrowing counts**, but only as a constraint on the product-preparation class. It does not change
   d_min^incl.
2. **No QEI / localization proof before G2.** GS-P1 stays heuristic. If the bound becomes publication-relevant, a dedicated
   localization-proof gate may be opened. QEIs alone are not assumed to prove L_all.
3. **G2 is reframed:** "What replaces QFT split / type-I subsystem structure when gravitational gauge constraints are
   imposed, and does the replacement remove any frozen residual information or only relocate it into charges, dressing,
   boundary observables, or resolution?"

## G2 note (no repair; recorded with G2)

- **G2 opens QG-8 … QG-13** (asymptotic charges, perturbative order, dressing prescription, observable class /
  resolution, vacuum premises, observer / clock). See `G2_SUBSYSTEM_LEDGER.md`.
- **Toy numerics procedure** (`g2/g2_subsystem.py`). Two numerical defects were found and fixed before the log was
  emitted:
  1. an out-of-memory full SVD of the stacked commutant system;
  2. a zero / near-scalar generator stack whose purely relative threshold counted rounding noise as rank.

  The final code uses unit-norm, identity-projected generators and an absolute floor. The resolution parameter ε is
  explicit and reported. Intermediate algebra dimensions at coarse ε are **resolution-dependent by construction**: they
  are an illustration of the resolution ladder, not a physical claim.
- **The S2 ingredient list** (boundary Hamiltonian, unique vacuum, density of boundary-generated states) is summarized
  from the owner-verified source and **not re-read** in this environment.

## GRAVITY REPAIR 02 — gravitational subsystem / resolution scope (owner ruling after GRAVITY_ZOOM_OUT_02)

**Applied to:** `G2_GRAVITATIONAL_SUBSYSTEM_RESULT.md`, `G2_SUBSYSTEM_LEDGER.md`, `GRAVITY_ZOOM_OUT_02.md`, status.
**Effect:** G2 is accepted **provisionally** after this repair.
**No change to:** TRUE COMPRESSION = 0; empirical payoff = 0; zero confirmed distinctive GRUT quantitative predictions.

| ID | correction | applied as |
|---|---|---|
| **GR2-01** | Donnelly–Giddings do not literally reconstruct a replacement type-I algebra | "𝒩 is replaced by ⊕_charge 𝒩_E" is withdrawn **as a continuum claim**. Replacement: **the ordinary AQFT type-I interpolation is no longer the relevant localization structure; perturbatively, a charge-labelled gravitational splitting (of Hilbert-space subspaces, at leading order) replaces its operational role.** The toy ⊕_E B(ℋ_E) is an **algebraic caricature only**: shared charge labels, sector-wise independence, loss of arbitrary product specification. Classification: **CONSTRAINED / CHARGE-SECTOR GRAVITATIONAL SPLITTING** + **RELOCATION into charges + dressing + perturbative order** |
| **GR2-02** | Raju forbids ordinary split independence, not every conceivable type-I factor | Terminal: **ORDINARY QFT SPLIT INDEPENDENCE FORBIDDEN IN CLASS** (verified scope). A Doplicher–Longo-style 𝒩 implementing independent inside / outside specification is **BLOCKED / NOT APPLICABLE** there. No claim that no type-I factor of any kind can occur |
| **GR2-03** | Resolution vs observable-class dependence | **Source-backed:** subsystem independence is **OBSERVABLE-CLASS / COARSE-GRAINING / PERTURBATIVE-ORDER dependent** (Donnelly–Giddings: leading-order exterior observables resolve total charges only; Raju: a rich boundary algebra determines the global state, while coarse-grained observable sets and their entropy are generally state-dependent). **Illustration only:** the ε ladder → **A_resolution COUPLING CANDIDATE — ILLUSTRATION GRADE**; gravity is **not** shown to determine A_resolution. The toy 75 / 6 equality with T1 is a **dimension coincidence only** |
| **GR2-04** | Residual mapping | Σ / A_partition: **CONSTRAINED / OBSERVABLE-CLASS-PRICED**. A_interface: **RELOCATION**. A_resolution: **POTENTIALLY LOAD-BEARING; NOT YET DERIVED**. A_time: not yet tested. Gravitational charge sector: **NEW SUPPLIED LABEL**. No supplied residual information disappears. The earlier "A_resolution and Σ / A_partition are no longer independent residual entries" is superseded |
| **GR2-05** | Crossed-product controls | Witten II∞ and CLPW II₁ are kept as examples of gravity changing algebraic type. Terminal: **NO TYPE-I INTERPOLATION ESTABLISHED**. No inference of a type-I interpolation, independent product preparation, or a universal gravitational subsystem algebra |
| **GR2-06** | Revised G2 terminal | Recorded verbatim in `G2_GRAVITATIONAL_SUBSYSTEM_RESULT.md` ("Revised G2 terminal") |

### Owner ruling on the next gate (recorded with REPAIR 02)

- **Orientation is deferred.**
- **G3 = OPERATIONAL ACCESS, FINITE TIME AND GRAVITATIONAL RESOLUTION** (G3-0 … G3-10). Hard stop after
  `GRAVITY_ZOOM_OUT_03.md`.

## G3 note (no repair; recorded with G3)

- **Access vector.** G3 records resolution as an access vector (𝒪, Λ, ε_t, k, ∂, δ), not as a single ε
  (`G3_ACCESS_LEDGER.md` §1).
- **Toy numerics procedure** (`g3/g3_timeband.py`). Two defects were found and fixed before the log was emitted:
  1. mpmath's general eigensolver failed to converge; it was replaced by the Hermitian solver `eighe`;
  2. a precision artifact at the smallest ε_t for N = 8 (condition number beyond the working precision) produced a spurious
     trade-off row. Fixed with 160-digit arithmetic and a positivity guard that treats beyond-precision matrices as
     singular.
- **Mode normalizations** c_n are set to 1. The model is the standard global-AdS₄ l = 0 spectrum, not a transcription of
  the CPR paper's own truncated construction (not re-read).
- **Source caution.** A search summary attributed to arXiv:2008.01740 a coarse-grained time-band algebra with a
  non-trivial commutant in states with a macroscopic bulk observer. This may conflate it with later work, and it is not
  used as evidence. All CPR protocol facts are owner-stated, not re-read here.
