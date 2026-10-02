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

| ID | information item | GRUT layer | SCOUT component | status after B0 + B1 | evidence |
|---|---|---|---|---|---|
| BI-01 | local net {𝒜_x} (site frame + edges) | 1 substrate | Σ | **NO BRIDGE** from the linear generator (general and L0-1b classes) | B1 table rows 1–6; N = 12 counterexample |
| BI-02 | local net, narrow C1-a class | 1 | Σ | **CONDITIONAL COMPRESSION**, CRITERION-PRICED (unit springs + Laplacian form = a supplied net class; DLS) | B1 |
| BI-03 | local net vs on-site drift | 1 ↔ 2 (declared nonlinear drift) | Σ ↔ D_dyn | **RELOCATION / REDUNDANT SUPPLY + CONSISTENCY** (odeco recovery exact; control selects the other frame) | B1-5 |
| BI-04 | local net vs site temperatures | 1 ↔ 5 (environment) | Σ ↔ H_corr / environment | **RELOCATION** if Tᵢ is non-uniform; **NO BRIDGE** if uniform | B1-4 |
| BI-05 | retained-site / bath partition | 4 (partition) | Σ grouping | **NO BRIDGE** from K (every unit vector is a retained site of some chain) | B1-3 |
| BI-06 | geometry from the net | 1 → effective | S2-G2 (derived from connectivity) | **CONDITIONAL COMPRESSION** in both programs (derived *given* the net). Across nets the geometry and **dimension** vary for the same K | B0 row; B1-2, (ii) |
| BI-07 | generator / law | 2 | D_dyn | **RENAMING / IDENTICAL**: supplied in both; "not derived" in both | crosswalk |
| BI-08 | time orientation | 2 | (not selected) | **NO BRIDGE**: GRUT supplies it; SCOUT shows it is not selected (D0) | crosswalk |
| BI-09 | "H is k-local" | — | IR-01 CPR class | **RENAMING** of the same smuggled criterion (both programs mark it) | EA0 §4; IR-01 |
| BI-10 | access / readout / coarse-graining | 3, 5 | A_res | **NOT YET TESTED** (B3) | — |
| BI-11 | preparation / statistics | 4 | H_corr\|Σ | **NOT YET TESTED** (B2); SCOUT STRICTLY SHARPER on the identification of the price | crosswalk |
| BI-12 | environment correlations / freshness | 5 | H_corr | **NOT YET TESTED** (B2) | — |
| BI-13 | quantum lift, ħ, Born rule, gravity, cosmology | 6–10 | — | **BLOCKED / NO MAPPING** by rule (not reopened) | charter |
| BI-14 | commutant copies of a net | — | H-relative gauge | **GAUGE** (both programs) | B1 C1-a row; SCOUT handoff Q6 |

## Tally after B0 + B1

| class | count |
|---|---|
| TRUE COMPRESSION | **0** |
| CONDITIONAL COMPRESSION | 2 (BI-02, BI-06); each priced by a supplied net |
| RELOCATION / REDUNDANT SUPPLY / CONSISTENCY | 2 (BI-03, BI-04) |
| RENAMING | 2 |
| GAUGE | 1 |
| NO BRIDGE | 3 (BI-01, BI-05, BI-08) |
| NOT YET TESTED | 3 |
| BLOCKED / NO MAPPING | 1 group |

**Layers eliminated:** none. **Redundancies found:** the local net is supplied redundantly whenever site-local drift or
non-uniform site noise is also supplied.
