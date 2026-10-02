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

| ID | information item | GRUT layer | SCOUT component | status after B0 + B1 + B3 | evidence |
|---|---|---|---|---|---|
| BI-01 | local net {𝒜_x} (site frame + edges) | 1 substrate | Σ | **NO BRIDGE** from the linear generator (general and L0-1b classes) | B1 table rows 1–6; N = 12 counterexample |
| BI-02 | local net, narrow C1-a class | 1 | Σ | **CONDITIONAL COMPRESSION**, CRITERION-PRICED (unit springs + Laplacian form = a supplied net class; DLS = KNOWN-RESULT IMPORT [BR1-03]) | B1 |
| BI-03 | local net vs on-site drift | 1 ↔ 2 (declared nonlinear drift) | Σ ↔ D_dyn | **RELOCATION / REDUNDANT SUPPLY + CONSISTENCY** (odeco recovery exact; control selects the other frame) | B1-5 |
| BI-04 | local net vs site temperatures | 1 ↔ 5 (environment) | Σ ↔ H_corr / environment | **RELOCATION** if Tᵢ is non-uniform; **NO BRIDGE** if uniform | B1-4 |
| BI-05 | retained-site / bath partition | 4 (partition) | Σ grouping | **NO BRIDGE** from K (every cyclic unit vector is a retained site of some chain [BR1-01]) | B1-3 |
| BI-06 | geometry from the net | 1 → effective | S2-G2 (derived from connectivity) | **CONDITIONAL COMPRESSION** in both programs (derived *given* the net). Across nets the geometry and **dimension** vary for the same K | B0 row; B1-2, (ii) |
| BI-07 | generator / law | 2 | D_dyn | **RENAMING / IDENTICAL**: supplied in both; "not derived" in both | crosswalk |
| BI-08 | time orientation | 2 | (not selected) | **NO BRIDGE**: GRUT supplies it; SCOUT shows it is not selected (D0) | crosswalk |
| BI-09 | "H is k-local" | — | IR-01 CPR class | **RENAMING** of the same smuggled criterion (both programs mark it) | EA0 §4; IR-01 |
| BI-10 | access / readout / coarse-graining | 3, 5 | A_res | **SPLIT by B3** into BI-10a … BI-10g (see `B3_ACCESS_COMPONENT_LEDGER.md`) | B3 |
| BI-10a | access seed (what is coupled) | 3 | A_res: coupled | **NONUNIQUE / SUPPLIED**: NOT REDUCIBLE TO D + Σ | B3-2 |
| BI-10b | closure / effect algebra | 3 | A_res: effect structure | **CONDITIONAL COMPRESSION** f(D, seed; closure rule). The canonical record uses ≥ 2 rules (B3-CUC-1) | B3-1 |
| BI-10c | readout map | 3 | A_res: readout | **SUPPLIED** (completeness ≠ selection); reconstruction-given-readout = CONDITIONAL | B3-3 |
| BI-10d | access sets / fragment grouping | 3, 5 | A_res: grouping | **SUPPLIED**; GS1 ONE-WAY (A → geometry) | B3-4, B3-9 |
| BI-10e | resolution / hierarchy order | 3 | A_res: threshold | **SUPPLIED** (order- and threshold-relative) | B3-7 |
| BI-10f | sampling / horizon / boundary-change time | 3 | A_res: time | **SUPPLIED** (D gives units only) | B3-10 |
| BI-10g | endogenous access (EA-0) | 3 ↔ 6 (lift) | — | **BLOCKED at Level-0**; auxiliary non-trivial branch = **RELOCATION INTO D** | B3-5 |
| BI-11 | preparation / statistics | 4 | H_corr\|Σ | **NOT YET TESTED** (B2); SCOUT STRICTLY SHARPER on the identification of the price | crosswalk |
| BI-12 | environment correlations / freshness | 5 | H_corr | **NOT YET TESTED** (B2) | — |
| BI-13 | quantum lift, ħ, Born rule, gravity, cosmology | 6–10 | — | **BLOCKED / NO MAPPING** by rule (not reopened) | charter |
| BI-14 | commutant copies of a net; access relabelling / appended uncoupled sectors (P-5 L-E, L-A) | — | H-relative gauge | **GAUGE** (both programs) | B1 C1-a row; SCOUT handoff Q6 |

## Tally after B0 + B1 + B3

| class | count |
|---|---|
| TRUE COMPRESSION | **0** |
| CONDITIONAL COMPRESSION | 3 (BI-02, BI-06, BI-10b), plus the reconstruction-given-readout part of BI-10c; each priced by a supplied net, seed or rule |
| RELOCATION / REDUNDANT SUPPLY / CONSISTENCY | 3 (BI-03, BI-04; BI-10g auxiliary branch → D) |
| RENAMING | 2 |
| GAUGE | 1 |
| NO BRIDGE | 3 (BI-01, BI-05, BI-08) |
| NONUNIQUE / SUPPLIED (B3) | 5 (BI-10a, c, d, e, f) |
| NOT YET TESTED | 2 (BI-11, BI-12: B2) |
| BLOCKED / NO MAPPING | 1 group, plus BI-10g at Level-0 |

**Layers eliminated:** none. **TRUE COMPRESSION:** 0. **Access bookkeeping:** A_res sharpened to seed ⊕ partition ⊕ resolution ⊕ time, with the closure downstream (conditional). **Redundancies found:** the local net is supplied redundantly whenever site-local drift or
non-uniform site noise is also supplied.
