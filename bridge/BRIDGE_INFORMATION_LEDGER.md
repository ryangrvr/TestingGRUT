# BRIDGE INFORMATION LEDGER (Bridge-1)

**What the ledger records.** For each supplied GRUT layer and each reviewed SCOUT residual component, it records where the
information enters and what Bridge-1 has shown about it.

**Classes:**
- **TRUE COMPRESSION** — reconstructed without supplying it.
- **CONDITIONAL COMPRESSION** — reconstructed given a stated, priced condition.
- **REDUNDANT SUPPLY** — supplied twice.
- **CONSISTENCY / OVERDETERMINATION** — two supplied layers agree.
- **RELOCATION** — recoverable only from another supplied layer that encodes it.
- **RENAMING** — the same input under another name.
- **GAUGE** — not physical.
- **NO BRIDGE**
- **BLOCKED**
- **NOT YET TESTED**

**Elimination rule (owner).** A layer is eliminated only by TRUE COMPRESSION.

| ID | information item | GRUT layer | SCOUT component | status after B0 + B1 + B3 + B2 | evidence |
|---|---|---|---|---|---|
| BI-01 | local net {𝒜_x} (site frame + edges) | 1 substrate | Σ | **NO BRIDGE** from the linear generator (general and L0-1b classes) | B1 table rows 1–6; N = 12 counterexample |
| BI-02 | local net, narrow C1-a class | 1 | Σ | **CONDITIONAL COMPRESSION**, CRITERION-PRICED (unit springs + Laplacian form = a supplied net class; DLS = KNOWN-RESULT IMPORT [BR1-03]) | B1 |
| BI-03 | local net vs on-site drift | 1 ↔ 2 (declared nonlinear drift) | Σ ↔ D_dyn | **RELOCATION / REDUNDANT SUPPLY + CONSISTENCY** (odeco recovery exact; control selects the other frame) | B1-5 |
| BI-04 | local net vs site temperatures | 1 ↔ 5 (environment) | Σ ↔ H_corr / environment | **RELOCATION** if Tᵢ is non-uniform; **NO BRIDGE** if uniform | B1-4 |
| BI-05 | retained-site / bath partition | 4 (partition) | Σ grouping | **NO BRIDGE** from K (every cyclic unit vector is a retained site of some chain [BR1-01]) | B1-3 |
| BI-06 | geometry from the net | 1 → effective | S2-G2 (derived from connectivity) | **CONDITIONAL COMPRESSION** in both programs (derived *given* the net). Across nets the geometry and **dimension** vary for the same K | B0 row; B1-2, (ii) |
| BI-07 | generator / law | 2 | D_dyn | **RENAMING / IDENTICAL**: supplied in both; "not derived" in both | crosswalk |
| BI-08 | time orientation | 2 | (not selected) | **NO BRIDGE**: GRUT supplies it as a convention tied to the preparation; SCOUT shows it is not selected (D0); B2 shows the S6 dynamics are Janus-symmetric and a reversed preparation anti-relaxes | crosswalk; B2-4, B2-5, B2-6 |
| BI-09 | "H is k-local" | — | IR-01 CPR class | **RENAMING** of the same smuggled criterion (both programs mark it) | EA0 §4; IR-01 |
| BI-10 | access / readout / coarse-graining | 3, 5 | A_res | **SPLIT by B3** into BI-10a … BI-10g (see `B3_ACCESS_COMPONENT_LEDGER.md`) | B3 |
| BI-10a | access seed (what is coupled) | 3 | A_res: coupled | **NONUNIQUE / SUPPLIED**: NOT REDUCIBLE TO D + Σ | B3-2 |
| BI-10b | closure / effect algebra | 3 | A_res: effect structure | **CONDITIONAL COMPRESSION** f(D, seed; closure rule). The canonical record uses ≥ 2 rules (B3-CUC-1) | B3-1 |
| BI-10c | readout map | 3 | A_res: readout | **SUPPLIED** (completeness ≠ selection); reconstruction-given-readout = CONDITIONAL | B3-3 |
| BI-10d | access sets / fragment grouping | 3, 5 | A_res: grouping | **SUPPLIED**; GS1 ONE-WAY (A → geometry) | B3-4, B3-9 |
| BI-10e | resolution / hierarchy order | 3 | A_res: threshold | **SUPPLIED** (order- and threshold-relative) | B3-7 |
| BI-10f | sampling / horizon / boundary-change time | 3 | A_res: time | **SUPPLIED** (D gives units only) | B3-10 |
| BI-10g | endogenous access (EA-0) | 3 ↔ 6 (lift) | — | **BLOCKED at Level-0**; auxiliary non-trivial branch = **RELOCATION INTO D** | B3-5 |
| BI-11 | preparation / statistics | 4 | H_corr\|Σ | **SPLIT by B2** into BI-11a … BI-11d (see `B2_BOUNDARY_COMPONENT_LEDGER.md`). The S6 product preparation is one explicit member of the H_corr\|Σ class: **RENAMING**, independently converged | B2 |
| BI-11a | system marginal H_sys | 4 | H_corr\|Σ (marginal) | **SUPPLIED** (forgotten locally, kept in X_J) | B2-0, B2-8 |
| BI-11b | S–B correlations H_cross | 4 | H_corr\|Σ proper | **SUPPLIED; INDEPENDENT; LOAD-BEARING** (sets the sign of X_J at equal T) | B2-1, B2-2 |
| BI-11c | boundary epoch H_epoch | 4 | special moment | **SUPPLIED** (selected by the preparation form; Janus-symmetric) | B2-5 |
| BI-11d | temperature values T_s, T_b | 4, 5 | — | **SUPPLIED**; thermal class CONDITIONALLY DEFINED (stationarity also admits non-thermal GGEs) | B2-7 |
| BI-12 | environment correlations / freshness | 5 | H_corr | bath state → local asymptotic state: **PARTIAL RELOCATION INTO ENVIRONMENT** (conditional on Gibbs postulate + T_b). Relaxation mechanism: **CONDITIONAL COMPRESSION** of H_corr memory downstream (B2-P1, bridge sketch). Noise layer: **RELOCATION** of the reference state | B2-3, B2-8, B2-9, B2-10 |
| BI-13 | quantum lift, ħ, Born rule, gravity, cosmology | 6–10 | — | **BLOCKED / NO MAPPING** by rule (not reopened) | charter |
| BI-14 | commutant copies of a net; access relabelling / appended uncoupled sectors (P-5 L-E, L-A) | — | H-relative gauge | **GAUGE** (both programs) | B1 C1-a row; SCOUT handoff Q6 |

## Tally after B0 + B1 + B3 + B2

| class | count |
|---|---|
| TRUE COMPRESSION | **0** |
| CONDITIONAL COMPRESSION | 4 (BI-02, BI-06, BI-10b, BI-12 memory), plus the reconstruction-given-readout part of BI-10c; each priced by a supplied net, seed, rule or bath state |
| RELOCATION / REDUNDANT SUPPLY / CONSISTENCY | 5 (BI-03, BI-04; BI-10g auxiliary branch → D; BI-12 bath state → local limit; BI-12 noise → reference state) |
| RENAMING | 3 (BI-07, BI-09, BI-11: S6 preparation ↔ H_corr\|Σ) |
| GAUGE | 1 |
| NO BRIDGE | 3 (BI-01, BI-05, BI-08) |
| NONUNIQUE / SUPPLIED | 9 (B3: BI-10a, c, d, e, f; B2: BI-11a, b, c, d) |
| NOT YET TESTED | 0 among the mapped layers |
| BLOCKED / NO MAPPING | 1 group, plus BI-10g at Level-0 |

**Layers eliminated:** none. **TRUE COMPRESSION:** 0. **Access bookkeeping:** A_res sharpened to A_interface ⊕ A_partition ⊕ A_resolution ⊕ A_time, with A_interface := (A_seed, A_readout) and A_closure = f(D, A_seed; R_closure) [BR2-01, BR2-02]. **Preparation bookkeeping (B2):** H_corr|Σ made explicit as (H_marginals, H_cross) at the declared-form epoch; H_cross is load-bearing, temperature difference and low entropy are not the price. **Redundancies found:** the local net is supplied redundantly whenever site-local drift or
non-uniform site noise is also supplied.
