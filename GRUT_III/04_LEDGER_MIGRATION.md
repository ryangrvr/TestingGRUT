# 04 — Ledger migration (revision 2): register node → GRUT III status

Status vocabulary: POSTULATE (counted), IMPORTED (established result; zero credit), DATUM (not counted),
IDENTIFICATION (naming; no content), DERIVED (computed in `calc/` from postulates + imports), RETIRED.

| Register node | Old Δ | GRUT III status | New Δ | Where |
|---|---|---|---|---|
| `rung1_inin_formalism` — system/bath split | +1 | IDENTIFICATION (P1: commutant algebra; type III, no factorisation) | 0 | 01 P1 |
| `rung1_inin_formalism` — Gaussian truncation | +1 | not assumed (linear response is an approximation of P0) | 0 | 01 P0 |
| `rung1_inin_formalism` — Lorentzian causal structure | +1 | inside P0 | (P0) | 01 P0 |
| `rung1_inin_formalism` — gauge-orbit zero | +1 | inside P0/P3 | (P0) | 01 P0 |
| **P0** covariant QFT framework | — | **POSTULATE** | **+1** | 01 |
| `rung1_ontology_finite_memory` | +1 | RETIRED (memory is a property of (observable, state): `07`) | 0 | 07 |
| `rung2_kms_gate` | −1 | IMPORTED (Bisognano–Wichmann / Kay–Wald): no credit either way | 0 | 02 D1 |
| `rung3_single_pole` | 0 | RETIRED (replaced by DERIVED kernels D2, D5, D9, D10) | 0 | 02 |
| `rung4_love_kk` | 0 | IMPORTED (Kramers–Kronig) | 0 | 02 D5 |
| `rung5_gr_limit` — entropy ∝ area | +1 | IMPORTED (area law) with the coefficient = DATUM G, inserted in P2 | 0 | 01 P2 |
| `rung5_gr_limit` — Unruh T | +1 | IMPORTED (Bisognano–Wichmann); ħ a DATUM | 0 | 01 P1 |
| **P2** local equilibrium (Clausius) | — | **POSTULATE** | **+1** | 01, 02 D3 |
| `rung6_qm_limit` — quantization | +1 | inside P0 | (P0) | 01 P0 |
| `rung6_qm_limit` — Born measure | +1 | **P5** POSTULATE (non-contextuality) + IMPORTED Gleason | **+1** | 01 P5 |
| `rung7_wz` — ε, τ₂, shape | +3 | RETIRED (Λ is an integration constant; w ≡ −1 for the minimal framework) | 0 | 02 D3, 03 |
| `rung8_falsifier` — κ, ω_c | +2 | RETIRED (decoherence kernel DERIVED and invisible; no tabletop claim) | 0 | 02 D5 |
| `arrow_of_time` | +1 | **P4** POSTULATE (located) | **+1** | 01 P4 |
| `p_tt_ansatz` | +1 | DERIVED from ⟨TT⟩ under **P3** POSTULATE (coupling inserted) | **+1** (P3) | 01 P3 |
| `background_time_translation_flow` | +1 | IDENTIFICATION (modular/Killing flow, P1) | 0 | 01 P1 |
| **Total** | **+16** | | **+5** | |

Notes (revised).
1. **What the arithmetic means.** Eleven units left the ledger as follows: seven became IMPORTED theorems (these
   earn GRUT nothing; they were always available), two became DATA (G, ħ) by the charter's Λ ruling applied
   consistently, two became IDENTIFICATIONS. The five that remain are framework postulates. This is a
   reorganisation of the accounting, not a derivation of eleven inputs; revision 1's word "compression" is withdrawn.
2. **Where the gravitational coupling enters.** P2 (η = 1/4G) and P3 (κ). GRUT III does not derive G; `v9` shows
   the one-loop induced-gravity route fails with Standard-Model content (Σk₁ = +1/6 for ξ = 0, −1/2 for ξ = 1/6).
3. **The reviewer's charter rule.** The "Stage-3 charter §13.6" is in the owner's `vsctestinggrut` repository
   (`PROGRAM/STAGE3/STAGE3_CHARTER.md`; quotes in `08`): §13.1, NR-12(b), §13.6 RECOVERY, DEF-16 (self-referential
   modular-flow KMS earns nothing), CT-75. This table applies those rules: every thermality/FDT item is IMPORTED.
4. The Level-0 classical-substrate theorems of the record (P-06b, P-17, P-02, P-15, S6) stand as mathematics;
   P-15 is used in 07 and P-17's lesson ("reduced data identify the forcing law, not the ontology") is the ancestor
   of the Δρ framework.
