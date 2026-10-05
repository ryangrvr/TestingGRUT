# CONJECTURE MODE — FINAL SYNTHESIS 01 (`grut-conjecture-mode-1`)

**Closes** `grut-conjecture-mode-1` at head `a8dbf3f45b4502b65913dcbb709c171580ef5e34` ("CARD-01 v1R2 TERMINAL: COMPUTATIONALLY UNRESOLVED").

**Status.** This is a closing summary only. It changes no frozen result, ruling, threshold or status, and every statement below is taken from the committed record. The file inventory is in `CONJECTURE_MODE_IMPORT_MANIFEST_01.md`.

**Governing law.** `program_governance/PROGRAM_CAMPAIGN_GATE_01.md`, `program_governance/CONJECTURE_MODE_CHARTER_01.md` and `program_governance/PROGRAM_GOVERNANCE_OWNER_RULING_01.md`, frozen at `grut-program-governance-1-frozen @ 5baa8afd`.

## 1. Earned

1. **Conjecture Mode governance.** The program has two legal campaign doors:
   - Door D: derivation / selector;
   - Door C: conjecture.

   Door C cards must carry an explicit postulate, an information price, an admissible origin (C0.2), a frozen projection, a named dataset and preregistered kill-condition forms.

   Further rules:
   - **Thresholds are locked by the owner before data access.** The executor proposes, the reviewer audits, and the owner approves.
   - **Prior exposure must be disclosed.**
   - **Network access is checked in a preflight**, and an unreachable product means ACCESS-BLOCKED with no improvised substitute.
   - **Nested-null result states** are used.
   - **Anti-curve-fitting rules C-F1 to C-F10** apply.

   This was frozen by owner ruling at `5baa8afd`.
2. **Card #1 (horizon relaxor) was implemented and brought into contact with official DESI DR2 products.**
   - The exact historical toy law was transcribed from `scout-1 @ a2987fe`.
   - The sign proof (ε < 0 ⇒ w < −1) passed.
   - The registered CPL diagnostic slope `wa/(1+w0) = 1.55590467`. The local-derivative value 0.465 is unregistered, and 1.37 is discarded.
   - The model passed C1 (owner-accepted, notes N1–N6). The prescribed w(z) went directly into CAMB + Cobaya through PPF.
   - The primary data were DESI DR2 BAO plus the DESI DR2 paper's baseline CMB.
3. **Provenance was verified.**
   - All consumed likelihood data are pinned and SHA-256 hashed.
   - The DESI DR2 BAO files are traced to the DESI-designated `CobayaSampler/bao_data` tree; the pinned v2.6 tree is identical to master.
   - The official DESI chains match DESI's published checksums.
   - The implementation was conformed to DESI's own input settings **before** any evaluation (CR-1 to CR-4).
4. **The first apparent mechanical result was invalidated by convergence auditing.**
   - v1's mechanical `EXPLANATORY-SURVIVES` (I_card 4.28, F 0.61) failed the frozen 0.2 start-spread diagnostic, so it was never accepted.
   - The v1R verification then showed the v1 CPL comparator was badly under-converged: it improved by 6.36 χ² from an official-DESI start.
5. **The local-optimizer terminal was reached honestly.** The sequence was v1, then v1R, then v1R2. Each campaign was frozen before execution, values already seen were disclosed, and no criterion was adjusted toward an outcome. v1R2 was the final local-optimizer campaign, and there is no v1R3.
6. **EQ-01 equivalence banked (narrow scope).** Linear background targets in `ρ_m`, `ρ_v`, `H²` and `Ḣ` collapse to one linear-response family (`conjecture_mode/EQUIVALENCE_REGISTRY.md`).

## 2. Not earned

- **No accepted empirical verdict for Card #1.** v1, v1R and v1R2 are all unresolved. Card #1 is not EXPLANATORY-SURVIVES, EXPLANATORY-KILLED or MODEL-KILLED by any accepted run.
- **No GRUT prediction, and no validation.** A DR2 comparison is retrospective by construction.
- **No evidence that Card #1 is GRUT's generative law.** Its response form is a postulated heuristic toy law (`POSTULATED HEURISTIC TOY LAW — NOT DERIVED FROM THE GRUT EFFECTIVE STRESS TENSOR`).
- **Cards #2–#4 remain historical and unrepaired**, at their Annex-A grades with C0.2 not audited. None was rerun.
- **The supernova extensions, the ε > 0 control and the MCMC posterior were never run.**

## 3. Final status

| | |
|---|---|
| **Card #1** | **POSTULATED; COMPUTATIONALLY UNRESOLVED** (`CARD-01-v1R2 — COMPUTATIONALLY UNRESOLVED AT AVAILABLE MINIMIZATION METHOD`, commit `a8dbf3f`) |
| v1 | `826f0d1`: NUMERICALLY UNRESOLVED / NO ACCEPTED SCIENTIFIC VERDICT |
| v1R | `1c02f7d`: STILL NUMERICALLY UNRESOLVED |
| v1R2 | `a8dbf3f`: COMPUTATIONALLY UNRESOLVED. ΛCDM, CPL and free-ε were reproduced; fixed ε = −0.08 was not (+0.287 > 0.2) |
| Further optimizer campaigns | **None (no v1R3).** Any further work needs a new owner ruling, e.g. a fundamentally different, independent likelihood / minimization implementation |

**Orientation only — not a verdict, and not an S5 evaluation.** These are the reproduced v1R2 core minima.

| Quantity | Value |
|---|---|
| χ²_Λ | ≈ 10976.08 |
| χ²_CPL | ≈ 10963.81 |
| χ²_Card | ≈ 10972.72 (free-ε fit, ε ≈ −0.084) |
| I_Card | ≈ 3.36 |
| I_CPL | ≈ 12.28 |
| F | ≈ 0.27 |

## 4. Record note

The per-fit raw outputs and logs of v1, v1R and v1R2 were **not committed**. They live only in the executing container's ephemeral scratch space. The committed record holds the summary JSONs, analyzers, drivers, configurations and result documents. See the manifest §4 for status and the owner decision needed.
