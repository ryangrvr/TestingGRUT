# MASTER RESIDUAL LEDGER — GRUT RESIDUAL SYNTHESIS 01

**Consolidation only. NO NEW PHYSICS.** Every entry is copied from a frozen campaign at the grade it earned there.
Frozen correction ledgers take precedence over superseded wording.

## Frozen sources and hashes

| short | branch @ hash | correction / precedence source |
|---|---|---|
| **S** | `scout-2-reviewed @ 1f08ae156a346f6089ae3b3aab2be80b9ee4ef8f` | `review/REVIEW_ACCOUNTING_LEDGER.md`, `review/ERRATA_PROPOSED.md` (E-01 … E-08), precedence IR-06 |
| **B** | `grut-bridge-1-frozen @ 753b90a3e3976ed286f385ca1b455e04c92a2d26` | `bridge/BRIDGE_CORRECTION_LEDGER.md` (BR1-01 … BR5-05) |
| **Q** | `qft-scout-1-frozen @ 5b07573bf45650d332df14258554cc8428a0706c` | `qft_scout_1/QFT_CORRECTION_LEDGER.md` (Q1R-01 … Q4R-02) |
| **G** | `gravity-scout-1-frozen @ df25f29bf217931b366c0ef676cb38b7b34ad115` | `gravity_scout_1/GRAVITY_CORRECTION_LEDGER.md` (GR1-01 … GR3-06) |
| **C** | canonical GRUT `GRUT-RAI/master-w25bu9 @ b935099f61008acf111a762a6e346d947c2d0c38` | read-only; Ledger A-1 … A-16, nine supplied layers [BR5-01] |

## Vocabulary bridge (preserved, not merged)

| term | campaign | meaning (as frozen) |
|---|---|---|
| **C5 → D_dyn ⊕ [Σ ⊗ H_corr]_coupled ⊕ A_res** | S (reviewed) | the SCOUT residual |
| **H_corr\|Σ → (H_marginals, H_cross, H_epoch)** | B (B2) | a Bridge **resolution** of the same component, not a new component |
| **A_res → A_interface ⊕ A_partition ⊕ A_resolution ⊕ A_time**, with A_interface := (A_seed, A_readout) | B (B3) [BR2-01] | a Bridge resolution |
| **TRUE DERIVATION (scoped)** | S | an operational C5 structure derived from (H, ψ) inside the finite envelope. These derivations are *why* C5 reduced to the residual. They are **not** eliminations of a residual component |
| **TRUE COMPRESSION** | B, Q, G | a supplied target information item reconstructed without supplying equivalent information. Count **0** in all three |

## Ledger

**Columns:**
1. residual item;
2. original SCOUT status (S);
3. Bridge-1 result (B);
4. QFT-SCOUT refinement (Q);
5. Gravity-SCOUT refinement (G);
6. final grade;
7. supplied inputs still required;
8. frozen source hash(es);
9. controlling correction IDs.

| 1 item | 2 SCOUT (S) | 3 Bridge-1 (B) | 4 QFT-SCOUT (Q) | 5 Gravity-SCOUT (G) | 6 final grade | 7 still-required supplied inputs | 8 hash | 9 controlling IDs |
|---|---|---|---|---|---|---|---|---|
| **D_dyn** | genuine residual: the Hamiltonian / transition law, connectivity, propagation rule **not derived**. Graph / spectral / causal-order statistics, invariant measures under rigid D and structural attractors are **CONDITIONAL DERIVATIONS from D** | **RENAMING / IDENTICAL** with the canonical generator (BI-07). Layer 2 internally decomposable: D_dyn = K ⊕ generator class ⊕ class-relative orientation carrier (BI-18; B4-CUC-2). R1: the net is re-encoded in on-site drift (BI-03) | type III₁ local algebras: **injective state-independent outer modular homomorphism δ_M** (conditional specification compression; no representative, no physical time, no orientation) | not targeted (gravitational dynamics are supplied QG-1 – QG-5) | **SUPPLIED** (conditional compressions are downstream of D only) | K (static coupling), generator class, optional drift / noise law; for Q, the type class (III₁) | S, B, Q | IR-06(b); BR4-01; B4-CUC-2; Q2R-01, Q2R-02 |
| **Σ** (subsystem / local net / grouping) | coupled to H_corr; **not selected**. CRITERION-PRICED + NUMERICAL CLASSIFICATION; commutant copies are **GAUGE**. Physical non-uniqueness from objective / dimension / scale / state conflicts | **Σ REMAINS SUPPLIED**: **NO BRIDGE** from K (N = 12 isospectral inequivalent net). C1-a class: **CONDITIONAL COMPRESSION, CRITERION-PRICED** (DLS = known-result import). R1: **REDUNDANT SUPPLY / CONSISTENCY** in drift / non-uniform noise | net supplied (QP-2), so any Σ use is **RELOCATION**. HSMI: **PARTIAL CONDITIONAL COMPRESSION OF ONE-DIMENSIONAL ORDERED LOCALIZATION STRUCTURE** (a Borchers triple, not automatically a local interval net). Poincaré conditional on a supplied modular pattern | ordinary split independence: available (flat); its operational role **replaced** by charge-labelled gravitational splitting (perturbative); **forbidden in class** (fine-grained, Raju scope). **Σ / A_partition: CONSTRAINED / OBSERVABLE-CLASS-PRICED** | **SUPPLIED**; **OBSERVABLE-CLASS-PRICED** in the gravity envelope; partial conditional compression of 1D ordered localization only | a net / criterion / objective, or (Q) a modular seed + relative-commutant standardness + pattern; (G) an observable class, framework, perturbative order | S, B, Q, G | IR-05/E-06; BR1-01 … 03; Q3R-01, Q3R-04; GR2-01, GR2-02, GR2-04; GR3-06 |
| **H_marginals** | inside H_corr\|Σ; the state price is **not** H_env and **not** low entropy | H_sys **SUPPLIED** (memory compressed downstream). Bath marginal: **PARTIAL RELOCATION INTO ENVIRONMENT** (conditional on the Gibbs postulate + T_b). S_ref downstream (BI-16). Thermal class conditionally defined (BR3-01) | the state / vacuum is supplied (QP-3). With a split buffer, arbitrary normal marginals extend to normal product states | state class supplied (QP-5, QG-5) | **SUPPLIED** | T_s, T_b, Gibbs-vs-GGE postulate; (Q) the state / vacuum | S, B, Q | BR3-01; BR4-01 |
| **H_cross** | H_corr\|Σ proper: the load-bearing state price | **SUPPLIED; INDEPENDENT; LOAD-BEARING.** At equal temperature it sets the sign of X_J: X(∞) = −⟨E_int⟩_G(½ − α) (B2-P1, bridge theorem reviewed at stated scope) | **CONSTRAINED-NONUNIQUE.** H_cross = 0 **FORBIDDEN** for a sharp normal non-type-I factor / commutant pair (G1-P1); **admissible with a split buffer** | preparation level: **PREPARATION-CLASS CONSTRAINT — CONSTRAINED-NONUNIQUE — CONDITIONAL ON L_all — HEURISTIC** (d\* ~ (N ℓ_P² R)^{1/3}). Charge-sharp product restriction: **toy-level only** | **CONSTRAINED-NONUNIQUE** | the cross-correlation datum within its admissible class; (Q) sharp vs split, normality, factor type, Haag duality for the geometric reading; (G) L_all and collapse premises for d\* | B, Q, G | BR3-02, BR5-02; Q1R-01 … 06; GR1-02, GR1-03, GR1-05; GR2-01 |
| **H_epoch** | the special moment (Janus) | **SUPPLIED**: absolute coordinate gauge under global time translation; the existence of a special-form boundary event is physical | **NOT SELECTED**; relocated into (d, 𝒩, ω) | not targeted | **SUPPLIED / NOT SELECTED** | the preparation form / epoch | S, B, Q | BR4-01; Q2R-01 |
| **A_interface** := (A_seed, A_readout) | inside A_res (readout basis, effect restrictions) | A_seed **NONUNIQUE / SUPPLIED; NOT REDUCIBLE TO D + Σ**. A_readout **SUPPLIED** (completeness ≠ selection). A_closure = f(D, A_seed; R_closure) **CONDITIONAL COMPRESSION**. EA-0 **BLOCKED** at Level-0 (auxiliary branch = RELOCATION INTO D) | not targeted | **RELOCATION / LOAD-BEARING**: carries dressing, charge / boundary readout and the asymptotic-observable choice | **SUPPLIED; RELOCATION / LOAD-BEARING** | seed, readout; R_closure (convention); (G) dressing prescription, boundary-detector class, asymptotic observable | S, B, G | BR2-01, BR2-02; GR2-04; GR3-01, GR3-06 |
| **A_partition** | inside A_res (fragment grouping); coupled to Σ | **SUPPLIED**; GS1 ONE-WAY (access → geometry) | sharp vs split localization choice supplied | **CONSTRAINED / OBSERVABLE-CLASS-PRICED** (with Σ) | **SUPPLIED; OBSERVABLE-CLASS-PRICED** (gravity) | which sites / blocks / regions count | S, B, Q, G | BR2-01; GR2-04; GR3-06 |
| **A_resolution** | inside A_res (record threshold δ, coarse resolution) | **SUPPLIED** (order- and threshold-relative) | not targeted | G2: **A_resolution COUPLING CANDIDATE — ILLUSTRATION GRADE**. G3: **COARSE-GRAINING SUPPLIED**; **SUPPLIED / NOT SELECTED** (kinematic coupling to time-band robustness in a toy only) | **SUPPLIED / NOT SELECTED** | threshold / order; (G) the coarse-grained observable set, cutoff Λ, operator class | S, B, G | BR4-01; GR2-03; GR3-03, GR3-05, GR3-06 |
| **A_time** | inside A_res (history times) | **SUPPLIED** (D gives units only) | not targeted | **SUPPLIED / NOT SELECTED**. Constrained combinations: flat, finite time × Bondi-charge readout **FORBIDDEN IN CLASS**; AdS, small / infinitesimal boundary-time data determines low-energy / perturbative bulk information under the stated assumptions | **SUPPLIED / NOT SELECTED** (some (A_time, A_interface) combinations forbidden) | cadence, horizon, boundary-change time; (G) time band, window | B, G | BR4-01; GR3-01, GR3-03, GR3-06 |
| **orientation** | **NOT selected** (D0, scoped no-go) | **NO BRIDGE / CONVENTION**: GRUT supplies it tied to the preparation; class-relative orientation carrier; the reversed preparation anti-relaxes; two declared S6 forwards | **NOT SELECTED**: BW sign priced by the supplied relativistic structure; passivity orients only once imposed; the KMS sign is a convention; the hsm sign flips under antiunitaries | **SUPPLIED / NOT SELECTED** (campaign scope ruling; no gate run) | **NOT SELECTED** (supplied) | a preparation-form / sign convention; (Q) spectrum condition, wedge choice, passivity postulate; (G) retarded time, future null infinity, positive Hamiltonian, oriented time band | S, B, Q, G | E-03, E-04; BR3-02; B4-CUC-2; Q2R-03 … Q2R-06; GRAVITY REPAIR 03 orientation ruling |
| **dimension** | **not one primitive**: the local factor dimension is **SUPPLIED IN Σ**; graph / spectral / walk dimensions are CONDITIONAL DERIVATIONS (graph + probe); causal-order estimators are conditional within the calibration domain | geometry and dimension **vary across nets for the same K** (BI-06: CONDITIONAL COMPRESSION given the net) | **SUPPLIED** (encoded in the modular-position pattern, GI-09) | not targeted | **SUPPLIED** | the local factor dimension; the net / pattern | S, B, Q | Y-12, IR-06(c); BR5-03 (Lorentz wording); Q3R-01 |
| **physical scale** | scale-priced (τ, ε; stable window for δ) | D gives units only (BI-10f) | **PARAMETER NORMALIZATION CONDITIONALLY FIXED; PHYSICAL SCALE SUPPLIED** | ℓ_P assembled entirely from supplied ħ, G, c. The d\* bound uses supplied (ħ, G, c, R, N); no scale selected | **SUPPLIED** | units / scale; ħ, G, c | S, B, Q, G | Q3R-03; GR1-05 |
| **quantum lift** | finite QM envelope; composition: real / complex QM survive tested families (the field is not fixed) | **BLOCKED / NO MAPPING** (layers 6 – 8, by rule). Canonical record: NONUNIQUE-LIFT | **SUPPLIED** (QP-5), so all QFT results are auxiliary | **SUPPLIED** (QP-5 carried) | **SUPPLIED** (BLOCKED for the Bridge) | the lift, ħ, complex structure, outcome rule | S, B, Q, G | BR5-01 |
| **gravitational layer** | Lorentz structure **NOT DERIVED** | canonical layer 9 **BLOCKED / NO MAPPING**. GRUT supplies its causal cone / Lorentz structure (convergence on a residual boundary, not a shared derivation) | relativistic structure supplied where consumed (BW) | **opened and priced:** QG-1 – QG-13 **SUPPLIED**; results auxiliary to canonical GRUT | **SUPPLIED** | QG-1 – QG-13 | S, B, Q, G | BR5-01, BR5-03; Q2R-04; charter §0 (G) |
| **gravitational charge / asymptotic labels** | — | — | — | charge sector: **NEW SUPPLIED LABEL**. Asymptotic structure: **SUPPLIED and load-bearing**. Perturbative order: **SUPPLIED**. Dressing: **SUPPLIED** (RELOCATION into A_interface) | **SUPPLIED** | QG-8 – QG-10 | G | GR2-04; GR3-06 |

## Sub-items carried with their parent

| sub-item | parent | final status | sources |
|---|---|---|---|
| split buffer (d, 𝒩) | H_cross / Σ | **d:** inclusion level 0 under the split assumptions (fixed background); no positive scale derived (perturbative gravity, subsystem notion repriced); ordinary split may not apply (fine-grained, scoped). Preparation level: heuristic d\*. **𝒩:** canonical given a standard vector (Doplicher–Longo, bibliographic only) in flat AQFT; **no unique gravitational analogue**; no type-I interpolation established in the crossed-product controls | Q; G [GR1-01, GR2-01, GR2-02, GR2-05] |
| A_closure | A_interface | downstream: f(D, A_seed; R_closure), **CONDITIONAL COMPRESSION** (rule-priced; P-5 vs P-6 give 16 vs 4) | B [BR2-02] |
| R_closure | A_interface | **SUPPLIED convention / rule choice** | B [BR4-01] |
| S_ref | H_marginals (bath) | downstream of the bath state (LS-1): **CONDITIONAL COMPRESSION / RELOCATION** | B (BI-16) |
| modular-position sign / pattern | Σ / orientation / dimension | sign: **ORIENTATION RELOCATED**; pattern: **RELOCATION** (encodes the geometry type and dimension) | Q (GI-05, GI-07) [Q3R-04] |
| vacuum / state selector | H_marginals | **SUPPLIED (state-priced)** | Q (GI-03); G (QG-12) |
| Gibbs vs GGE | H_marginals | **SUPPLIED state-class postulate** | B [BR3-01, BR4-01] |

## Checks

- **No grade is stronger than the frozen campaigns earned.**
  - "CONDITIONAL COMPRESSION" appears only where a frozen ledger used it.
  - "CONSTRAINED-NONUNIQUE" appears only for H_cross (Q) and for the preparation class and access combinations (G).
- **TRUE COMPRESSION: 0 in B, Q and G.** S's scoped TRUE DERIVATIONS are recorded in `INFORMATION_FLOW_LEDGER.md` under
  their own name. They are not counted as residual compression.
