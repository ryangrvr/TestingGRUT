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
