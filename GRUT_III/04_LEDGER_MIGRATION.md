# 04 — Ledger migration: register node → GRUT III status

The register's net +16 (`provenance/claims.json`, validator PASS 2026-09-30), node by node.

| Register node | Old Δ | Old tier | GRUT III status | New Δ | Where |
|---|---|---|---|---|---|
| `rung1_inin_formalism` — system/bath split | +1 | assumed | **Identification P1**: the bath is the causal complement | 0 | 01 P1 |
| `rung1_inin_formalism` — Gaussian truncation | +1 | assumed | Not assumed; linear response is an approximation derived from P0 | 0 | 01 P0 |
| `rung1_inin_formalism` — Lorentzian causal structure | +1 | assumed | Part of P0 | (in P0) | 01 P0 |
| `rung1_inin_formalism` — 4d-covariant gauge-orbit zero | +1 | assumed | Part of P0/P3 | (in P0) | 01 P0 |
| **P0 (new)** — covariant QFT framework | — | — | **POSTULATE** | **+1** | 01 P0 |
| `rung1_ontology_finite_memory` | +1 | assumed/stance | Retired: finite memory is derived where true (noise), false where not (local response) | 0 | 02 D2 |
| `rung2_kms_gate` | −1 | shown | Theorem of P0+P1 (Bisognano–Wichmann / Kay–Wald); β fixed by the horizon | 0 (no credit claimed) | 02 D1 |
| `rung3_single_pole` | 0 | derived-pending | **Retired**: replaced by closed-form kernels (D2) and the species kernel (D5) | 0 | 02 D2, D5 |
| `rung4_love_kk` | 0 | shown | Derived (Kramers–Kronig on the polynomial/threshold kernel) | 0 | 02 D5 |
| `rung5_gr_limit` — entropy ∝ area (sets G) | +1 | assumed | Area law = QFT theorem; coefficient = **datum G** | 0 | 01 P2 |
| `rung5_gr_limit` — Unruh T (sets ħ) | +1 | assumed | Theorem (Bisognano–Wichmann); ħ a datum | 0 | 01 P1 |
| **P2 (new)** — local equilibrium / Clausius | — | — | **POSTULATE** | **+1** | 01 P2, 02 D3 |
| `rung6_qm_limit` — quantization | +1 | assumed | In P0 | (in P0) | 01 P0 |
| `rung6_qm_limit` — Born measure | +1 | assumed | **P5 non-contextuality** → Born by Gleason | **+1** | 01 P5, 02 D7 |
| `rung7_wz` — ε, τ₂ ∼ 1/H₀, single-departure shape | +3 | to-derive | Retired: Λ is an integration constant; w ≡ −1 | 0 | 02 D3, 03 F1 |
| `rung8_falsifier` — κ, ω_c | +2 | to-derive | Retired: decoherence kernel derived and invisible; no tabletop claim | 0 | 02 D5 |
| `arrow_of_time` — Past Hypothesis | +1 | assumed | **P4** (same price, located) | **+1** | 01 P4, 02 D8 |
| `p_tt_ansatz` | +1 | assumed | Computed from ⟨TT⟩ under **P3** | **+1** (P3) | 01 P3, 02 D5 |
| `background_time_translation_flow` | +1 | assumed | The horizon Killing flow (P1) | 0 | 01 P1 |
| **Total** | **+16** | | | **+5** | |

Notes.
1. The reduction is by *derivation* (seven nodes now follow from P0–P3 by computations in `calc/`), by
   *re-classification as data* (two nodes: G and ħ, per the charter's Λ ruling), and by *identification* (the
   split and the flow become the causal structure). It is not by relabeling.
2. Items the record carried outside the ledger (Λ_R, I1–I11 of the consolidated paper: access seed, carrier,
   probe, universal reach, co-stretch, cross-sector cone, ħ, Born, noncommutativity, states/arrow, sector content)
   map as follows: access/seed/probe → P1 (the observer's wedge); carrier/reach/cone → P0 (Einstein causality);
   ħ → datum; Born → P5; states/arrow → P4; sector content → datum (spectrum); Λ_R → the H⁰ local term of the
   kernel, fixed by the renormalisation condition G_ab = 8πG T_ab in the IR (D3).
3. The Level-0 classical-substrate programme (L0 floor, lifts, access bridge) is not migrated: in GRUT III the
   substrate is P0 and there is no classical-to-quantum lift to select. The theorems proved there (P-06b CM
   boundary, P-17 non-identifiability, P-02 edge-velocity, P-15 Born⟺martingale, S6 net arrow) stand as
   mathematics; P-15 is used (D7).
4. This table is a proposal for the owner's audit. The counting convention (what is a datum) is the charter's
   own; its consistent application is the only accounting change made here.
