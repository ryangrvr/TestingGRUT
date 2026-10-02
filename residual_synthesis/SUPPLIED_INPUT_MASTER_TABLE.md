# SUPPLIED-INPUT MASTER TABLE

**Consolidation only. NO NEW PHYSICS.**

**Class key:**
- **PHYS:** physical datum;
- **MATH:** mathematical structure / class choice;
- **PROT:** access / protocol declaration;
- **CONV:** convention.

**Rule (BR4-01).** An accounting item is **not** a physical primitive unless a frozen result says so. The class column
records the frozen classification only.

**Columns:**
1. ID;
2. supplied datum;
3. class;
4. first explicit;
5. consumed by;
6. duplicated elsewhere?;
7. downstream / derivable given other supplied data?;
8. final status.

## A. SCOUT-2 premises (reviewed)

| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| SP-1 | finite-dimensional premise envelope ((H, ψ) on a finite Hilbert / state space) | MATH | S (charter) | all S results | superseded in scope by QP-1 (Q) | no | SUPPLIED (envelope) |
| SP-2 | D: Hamiltonian / transition law, connectivity, propagation rule | PHYS | S | E-S5; all B dynamics | = canonical layer 2 (RENAMING, BI-07) | no | SUPPLIED |
| SP-3 | locality class / CPR-type criterion (k, d, n) | MATH | S (IR-01) | E-S2 | = canonical "H is k-local" (RENAMING, BI-09) | no | SUPPLIED (criterion-priced) |
| SP-4 | objectives, thresholds, scales (record threshold δ, τ, ε) | PROT | S | Σ selection, Darwinism, histories | overlaps A_resolution (B) | no | SUPPLIED |
| SP-5 | measure (typicality / mixing) | MATH | S | S2-3 relaxation; D2 typicality | — | uniquely ergodic case: fixed by rigid D (RELOCATION) | MEASURE-PRICED |
| SP-6 | epoch / special moment | PHYS (event) / CONV (coordinate) | S | D-arrow, H2 | = H_epoch (B) | no | SUPPLIED |

## B. Bridge-1: canonical layers and hidden accounting items

| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| L1a | K static coupling (A-1) | PHYS | C / B | X_J, closure, S_ref, geometry | — | no | SUPPLIED |
| L1b | local net Σ (A-1, A-2) | PHYS | C / B | geometry, access labels | **R1**: re-encoded by on-site drift / non-uniform noise (class-scoped) | C1-a only: CONDITIONAL COMPRESSION | SUPPLIED |
| L2a | generator class (A-3) | PHYS | C / B | all dynamics | — | no | SUPPLIED |
| L2b | declared nonlinear drift (A-12) | PHYS | C / B | R1 | encodes Σ (R1) | no | SUPPLIED (optional) |
| L2c | orientation carrier (A-4) | CONV | C / B | S6 forwards | — | no | SUPPLIED; class-relative (B4-CUC-2) |
| L3a | A_seed (A-13) | PROT | B (B3) | A_closure | — | not reducible to D + Σ | SUPPLIED |
| L3b | A_readout (A-13) | PROT | B (B3) | reconstruction-given-readout | — | no | SUPPLIED |
| L3c | access sets / A_partition (A-13) | PROT | B (B3) | GS1 geometry | overlaps the S–B split (R2′, overlapping not identical) | no | SUPPLIED |
| L3d | **A_resolution** | PROT | B (B4-11, hidden) | detectability | — | no | SUPPLIED (accounting item, not a physical primitive) [BR4-01] |
| L3e | **A_time** | PROT | B (B4-11, hidden) | sampling | — | D gives units only | SUPPLIED (accounting item) [BR4-01] |
| L3g | **R_closure** | CONV | B (B4-11, hidden) | A_closure | — | no | SUPPLIED convention [BR2-02, BR4-01] |
| L3f | A_closure | — | B | — | — | **downstream** f(D, A_seed; R_closure) | NOT SUPPLIED (CONDITIONAL COMPRESSION) |
| L3h | endogenous access EA-0 | — | B | — | — | needs a lift | BLOCKED at Level-0 |
| L4a | sector label (A-9) | PHYS | C / B | — | — | no | SUPPLIED |
| L4b | statistics (A-10) | PHYS | C / B | — | — | no | SUPPLIED |
| L4c | H_marginals: system marginal (A-14) | PHYS | B (B2) | X_J | — | memory compressed downstream | SUPPLIED |
| L4d | **H_cross** (inside A-14) | PHYS (genuine boundary data) | B (B4-11, hidden) | X_J sign | — | no | SUPPLIED → CONSTRAINED-NONUNIQUE (Q) |
| L4e | **H_epoch** (inside A-14) | PHYS event / CONV coordinate | B (B4-11, hidden) | — | = SP-6 | no | SUPPLIED / NOT SELECTED (Q) |
| L4f | system–bath split (A-14) | PROT | B | S6 | **R2**: duplicate of the layer-5 "partition" name | — | SUPPLIED [list duplicate] |
| L5a | environment origin (A-11) | PHYS | C / B | noise | — | — | SUPPLIED |
| L5b | bath state (T_b + **Gibbs-vs-GGE postulate**) | PHYS + MATH (postulate) | B (B4-11 hidden for the postulate) | S_ref | **R3**: duplicate in layers 4 / 5 | S_ref downstream | SUPPLIED [BR3-01, BR4-01] |
| L5c | noise law T_i (L0-1e) | PHYS | C / B | stationary correlations | encodes Σ when non-uniform (R1) | — | SUPPLIED (optional) |
| L5e | coarse-graining map (A-15) | PROT | C / B | reduced dynamics | ⊆ A_partition for fragment restrictions | general maps NOT TESTED | SUPPLIED |
| L6–8 | physical lift, ħ, outcome rule (A-5, A-6, A-7) | PHYS / MATH | C | — | — | — | BLOCKED / NO MAPPING (B); SUPPLIED |
| L9 | gravity quadruple incl. universal causal cone (A-8); cosmological transport (A-16) | PHYS | C | — | — | — | BLOCKED (B); opened as QG-1 – QG-4 (G) |

## C. QFT-SCOUT-1 premises

| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| QP-1 | infinite-dimensional von Neumann local algebras, possibly non-type-I | MATH | Q | G1-P1, G2, G3 | extends SP-1 | — | SUPPLIED |
| QP-2 | Haag–Kastler-type net (where used) | PHYS / MATH | Q | geometric readings | = Σ (RELOCATION) | — | SUPPLIED |
| QP-3 | vacuum / standard state; positive energy; covariance | PHYS | Q | BW, HSMI, modular flows | vacuum selector reappears as QG-12 | no | SUPPLIED (state-priced) |
| QP-4 | Haag duality; split / nuclearity; type III₁ | MATH | Q | G1 geometric step (Haag duality); split admissibility; injective δ_M | — | some imported as theorems for named models | SUPPLIED or imported per use |
| QP-5 | physical quantum lift (ħ, complex structure) | PHYS | Q | all Q and G results | = canonical layers 6 – 7 | no | SUPPLIED → everything auxiliary |
| QP-6 | modular-position data (HSMI sign, pattern, relative-commutant standardness) | MATH / CONV (sign) | Q (G3) | E-Q5, E-Q6 | — | the pattern encodes the dimension | SUPPLIED (RELOCATION) |
| QX-1 | normality (relative to the chosen representation) | MATH | Q | G1-P1 | — | not identified with finite energy | SUPPLIED [Q1R-05] |
| QX-2 | split buffer d and intermediate 𝒩 | PHYS (d) / MATH (𝒩) | Q | product admissibility | — | 𝒩 canonical given a standard vector (DL, bibliographic) | SUPPLIED |
| QX-3 | KMS / modular sign convention | CONV | Q | G2 | — | — | CONVENTION [Q2R-06] |
| QX-4 | passivity / Kelvin postulate (when imposed) | PHYS (postulate) | Q | thermodynamic orientation | — | — | SUPPLIED when used [Q2R-05] |

## D. GRAVITY-SCOUT-1 premises

| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| QG-1 – QG-4 | G; spin-2 / Einstein equations; universal coupling; causal cone | PHYS | C (A-8) / G | all G | = canonical layer 9 | no | SUPPLIED |
| QG-5 | semiclassical coupling + validity regime | PHYS (hypothesis) | G | E-G1 | — | — | SUPPLIED |
| QG-6 | trapped-surface admissibility (spherical: Misner–Sharp / Hayward; general: hoop conjecture) | PHYS | G | E-G1 | — | — | SUPPLIED |
| QG-7 | perturbative gauge constraints / dressing | PHYS | G | E-G3 | — | — | SUPPLIED |
| QG-8 | asymptotic structure; Poincaré / ADM / Bondi charges | PHYS | G (G2) | E-G3, E-G4, E-G5 | — | — | **SUPPLIED and load-bearing**; charge sector a NEW SUPPLIED LABEL |
| QG-9 | perturbative order | MATH | G (G2) | E-G3, E-G6 | — | no canonical nested hierarchy | SUPPLIED |
| QG-10 | dressing prescription (non-unique) | MATH / CONV | G (G2) | E-G3 | — | — | SUPPLIED (RELOCATION into A_interface) |
| QG-11 | observable class / resolution (incl. exact boundary spectral data) | PROT | G (G2) | E-G3, E-G7 | overlaps A_resolution / A_interface | — | SUPPLIED |
| QG-12 | vacuum premises (unique vacuum; boundary density) | PHYS | G (G2) | fine-grained split failure | overlaps QP-3 | — | SUPPLIED |
| QG-13 | observer / clock bounded below (CLPW) | PHYS | G (G2) | crossed-product control | — | — | SUPPLIED |
| GX-1 | L_all (localization of product-state energy) | PHYS (premise) | G (G1) | E-G1 | — | unproved | HEURISTIC PREMISE |
| GX-2 | collar U_ε, region R, field content N | PHYS / PROT | G | E-G1, E-G3 | — | — | SUPPLIED |
| GX-3 | G3 access vector: operator class 𝒪, cutoff Λ, time band ε_t, order k, asymptotic region ∂, precision δ | PROT | G (G3) | E-G4 – E-G7 | Λ overlaps A_resolution; ε_t overlaps A_time; 𝒪 overlaps A_interface | none selected | SUPPLIED |
| GX-4 | retarded time / future null infinity / oriented time band | CONV / PHYS | G | E-G4, E-G5 | carries orientation | — | SUPPLIED (orientation input) |

## E. Careful items (owner list)

| item | consolidated reading |
|---|---|
| **R_closure** | a convention / rule choice (P-5 vs P-6 give 16 vs 4). Not a physical primitive. A_closure is downstream of it [BR2-02, BR4-01] |
| **A_resolution** | a protocol declaration (B). Gravity: coarse-graining supplied; the coupling to subsystem structure is an illustration-grade candidate only [GR2-03, GR3-05] |
| **A_time** | a protocol declaration (B). Gravity: not selected; only (A_time, A_interface) combinations constrained [GR3-01] |
| **H_cross** | genuine physical boundary data (BR4-01). CONSTRAINED-NONUNIQUE (Q). Heuristic preparation constraint (G) |
| **H_epoch** | event physical; coordinate gauge (BR4-01). NOT SELECTED (Q) |
| **(d, 𝒩)** | supplied (Q). No positive universal d_min. Heuristic d\* at preparation level. No unique gravitational 𝒩 (G) |
| **modular sign / pattern** | sign = orientation relocated; pattern = geometry type / dimension relocated (Q) |
| **gravitational charges** | new supplied label (G) |
| **asymptotic structure** | supplied and load-bearing (G) |
| **perturbative order** | supplied (G) |
| **dressing prescription** | supplied, non-unique; relocated into A_interface (G) |
| **vacuum / state selector** | supplied (QP-3; QG-12); state-priced (GI-03) |
