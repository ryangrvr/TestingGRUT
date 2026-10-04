# VER0 STATUS

**State:** **V0-3 (P-02 / HCB-SF1) DONE — V0-3-C, VER-I1 (orchestrator-exposed); AWAITING OWNER REVIEW.** V0-1-C and
V0-2-C accepted. No further item started.

**Branch:** `grut-independent-verification-0`, from `scout-0 @ ab2da47` (frozen, unmodified).

**Boundaries:**
- V0-1: charter `71f0a59`; spec `063274b`; reproduction before unsealing `64b0182`; comparison `f83c929`; owner ruling
  `fdcd3a6`.
- V0-2: spec `cc4710e`; reproduction before unsealing `bde523e`; comparison `55a9aa8`; owner ruling `661cdbe`.
- V0-3: spec `335067a`; reproduction before unsealing `4ce8686`; comparison = this commit.

**SCOUT-0:** PROVISIONALLY SATURATED — OWED CHECKS ONLY. **Criterion 2 is not met.** Items 1 (P-17), 2 (P-02 / SF-1; TARGET 3 = F) and 3
(P-15) are complete at VER-I1 with correction; items 4, 5 and 6 are outstanding. Every item is at most VER-I1.

| item | grade |
|---|---|
| V0-1 P-17 | **V0-1-C, VER-I1 (orchestrator-exposed)** (owner ruling). Correction CL-1: ontology non-identifiability, not unrestricted inverse identifiability. C-B interpretation and C3 / C4 **not reproduced** (non-blocking) |
| V0-2 P-15 | **V0-2-C, VER-I1 (orchestrator-exposed)** (owner ruling). CR-1: b = 0 a.e. CR-2: absorbing endpoints are hypotheses; SF is needed only for the general boundary-value problem and follows inside the theorem; closed interval. CR-3: ρ₁₁ vs ρ |
| V0-3 P-02 / SF-1 | **V0-3-C, VER-I1 (orchestrator-exposed)**, with subgrade **V0-3-T3-F** (owner ruling). Core results exact. **CR-e:** single-v_b theorem refuted; edge-data invariants. CR-d: z = r_e. CR-b / CR-c: scope |
| V0-4 EDA-01 | queued |
| V0-5 K1-H / K1-HS | queued |
| V0-6 secondary primary-text checks | queued |
| VER0-B BRI1-X1 | registered; not opened; statement-only target; a context-isolated reproducer is required (the orchestrator is the author) |

No PR. No merge.
