# VER0 STATUS

**State:** **PAUSED FOR OWNER PROGRAM-LEVEL DISCUSSION.** VER0-B-A accepted (with S-1 banked). V0-1-C, V0-2-C, V0-3-C (with T3-F) and V0-4-C accepted.

**Branch:** `grut-independent-verification-0`, from `scout-0 @ ab2da47` (frozen, unmodified).

**Boundaries:**
- V0-1: charter `71f0a59`; spec `063274b`; reproduction before unsealing `64b0182`; comparison `f83c929`; owner ruling
  `fdcd3a6`.
- V0-2: spec `cc4710e`; reproduction before unsealing `bde523e`; comparison `55a9aa8`; owner ruling `661cdbe`.
- V0-3: spec `335067a`; reproduction before unsealing `4ce8686`; comparison `9d673bb`; owner ruling `e250f2e`.
- V0-4: spec `4110a67`; blind reading before unsealing `f181330`; comparison `3770ce7`; owner ruling `7d09ad7`.
- VER0-B: spec `938feb6`; reproduction before unsealing `046a945`; comparison `ae56ed0`; owner ruling = this commit.

**SCOUT-0:** PROVISIONALLY SATURATED — OWED CHECKS ONLY. **Criterion 2 is not met.** Items 1 (P-17), 2 (P-02 / SF-1; TARGET 3 = F) and 3
(P-15) are complete at VER-I1 with correction. Item 5 (EDA-01) is complete at VER-I1 with correction.
Items 4 and 6 are **OUTSTANDING — PRIMARY-TEXT ACCESS DEPENDENT**. Every item is at most VER-I1.

| item | grade |
|---|---|
| V0-1 P-17 | **V0-1-C, VER-I1 (orchestrator-exposed)** (owner ruling). Correction CL-1: ontology non-identifiability, not unrestricted inverse identifiability. C-B interpretation and C3 / C4 **not reproduced** (non-blocking) |
| V0-2 P-15 | **V0-2-C, VER-I1 (orchestrator-exposed)** (owner ruling). CR-1: b = 0 a.e. CR-2: absorbing endpoints are hypotheses; SF is needed only for the general boundary-value problem and follows inside the theorem; closed interval. CR-3: ρ₁₁ vs ρ |
| V0-3 P-02 / SF-1 | **V0-3-C, VER-I1 (orchestrator-exposed)**, with subgrade **V0-3-T3-F** (owner ruling). Core results exact. **CR-e:** single-v_b theorem refuted; edge-data invariants. CR-d: z = r_e. CR-b / CR-c: scope |
| V0-4 EDA-01 | **V0-4-C, VER-I1 (orchestrator-exposed)** (owner ruling). Headline NO-FIX confirmed (22 / 22). 6 / 22 rows discordant at the CON boundary. E-13 restated against the edge-data set |
| V0-5 K1-H / K1-HS | **outstanding — primary-text access dependent** |
| V0-6 secondary primary-text checks | **outstanding — primary-text access dependent** |
| VER0-B BRI1-X1 | **VER0-B-A, VER-I1 (orchestrator-exposed / context-isolated reproducer)** (owner ruling). **S-1 banked:** X1 ∉ E₂± for every N_B ≥ 1 on a common small-time interval (internally reproduced, not externally reviewed). P-1 is editorial |

No PR. No merge.
