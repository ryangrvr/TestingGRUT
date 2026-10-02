# B3 ACCESS COMPONENT LEDGER

**Fields (owner B3-14):** SCOUT item · GRUT analogue · parent inputs · derived structure · residual choice · dependence on
Σ / H_corr / lift · classification · evidence · canonical bookkeeping compressible?

**Firewall:** a conditional result inside a supplied quantum lift is **never** Level-0 TRUE COMPRESSION.

| ID | SCOUT A_res item | GRUT analogue | parent inputs | derived structure | residual choice | Σ dep. | H_corr dep. | lift dep. | classification | evidence | bookkeeping compressible? |
|---|---|---|---|---|---|---|---|---|---|---|---|
| AC-01 | what is coupled / accessible | access seed (P-5, P-6, A-13) | — (supplied) | none | which seed | the seed is written in Σ's labels; Σ does not select it | none | none (classical) / yes (quantum seeds) | **NONUNIQUE / SUPPLIED**; NOT REDUCIBLE TO D + Σ | B3-2 (112 / 276 distinct pairs; quantum 3/3); P-6 L-B, L-C | no |
| AC-02 | effect structure | closure cl(seed; D) | D + seed + **closure rule** | the closure algebra / observable subspace (unique given the rule) | the closure rule (R_P5 vs R_P6 differ on decoupled-hidden: 16 vs 4) | none beyond the seed | none | classical: none; quantum: auxiliary class | **CONDITIONAL COMPRESSION** f(D, A_seed; rule) | B3-1 | **yes: B3-CUC-1** (book the closure as downstream; fix one rule). Not REDUNDANT SUPPLY: A-13 does not book closures |
| AC-03 | readout map | earned end-site identity h = x₁ (L0 bridge) | D + declared h | state reconstruction (globally invertible) | which h (an open dense set of complete readouts) | h uses Σ's site labels | none | none | **A_readout REMAINS SUPPLIED**; reconstruction-given-h = CONDITIONAL (identity-grade) | B3-3; L0 bridge §4 | no |
| AC-04 | fragment / locus grouping | declared access sets (A-13); GS1 boundary | D + Σ + declared sets | geometry (GS1, forward only) | which sites / blocks count | presupposes Σ; Σ does not select | the record counts depend on the preparation (B2) | none | **A_partition REMAINS SUPPLIED**; **GS1 ONE-WAY** | B3-4 (9/9, 84/84 complete; selectors trivial / tied); B3-9 | no |
| AC-05 | resolution / threshold | P-6 L-D1 order-relativity; GS1 horizon | D + seed + δ / order | detectability at a given δ / order | δ, truncation order | depth measured in Σ's graph distance | none | none | **A_resolution REMAINS SUPPLIED** | B3-7 (first differing order 32 / 35 / 40; 2 / 7 / 11 / 15 detectable) | no |
| AC-06 | time sequence / sampling / horizon | P-5 boundary-change time t* (declared) | D + schedule | when information is available | cadence, horizon, t* | none | sampling interacts with the preparation (B2) | none | **A_time REMAINS SUPPLIED** (D gives units only) | B3-10 | no |
| AC-07 | coarse-graining | layer 5 coarse-graining | D + Σ + coarse map | reduced dynamics given the map | the coarse map | yes | yes (B2) | none | **CONDITIONAL**: ⊆ A_partition where the map is a fragment restriction; general aggregating maps **NOT TESTED** | B3-9 | possibly (pending) |
| AC-08 | endogenous / state-dependent access | EA-0 (P_ρ, Γ(ρ)) | lift + ρ + D | auxiliary: identity / non-selective / block-relative | the lift itself | — | values only: (a), not (c) | **requires a lift** | **BLOCKED AT LEVEL-0 (CATEGORY MISMATCH)**; auxiliary non-trivial branch = **RELOCATION INTO D** | EA0 ruling 02; B3-5; B3-6 | no |
| AC-09 | access representation (labels, appended uncoupled sectors) | P-5 "representational" | — | invariance | none | — | — | — | **GAUGE / RENAMING** | P-5 L-E, L-A; P-6 L-F; B3-1 covariance | yes (already standard) |

## Count

| | |
|---|---|
| TRUE COMPRESSION | **0** |
| CONDITIONAL COMPRESSION | 2 (AC-02 closure; AC-03 reconstruction-given-readout), plus AC-07 conditional inclusion |
| RELOCATION INTO D | 1 (AC-08, auxiliary branch) |
| BLOCKED | 1 (AC-08 at Level-0) |
| NONUNIQUE / SUPPLIED | 5 (AC-01, AC-03 selection, AC-04, AC-05, AC-06) |
| GAUGE | 1 (AC-09) |

**Sharpened residual for A:**

    A_interface ⊕ A_partition ⊕ A_resolution ⊕ A_time,   A_interface := (A_seed, A_readout)   [BR2-01]
    A_closure = f(D, A_seed; R_closure)   [BR2-02]
