# CARD #1 — OWNER THRESHOLD RULING

**Ruling on the pre-data spec:** `grut-conjecture-mode-1 @ 8c634145178eaad29153934b0d2846c6628bdd07` (`CARD_01_HORIZON_RELAXOR_SPEC.md`).

**Authority (Charter §4A, ER-2):**
- **Executor** (Claude Code): proposed S1–S8 (CLAUDE-PROPOSAL, non-binding).
- **Reviewer:** audited the proposals and recommended the modifications below, plus the O3 provenance ruling.
- **Scientific owner:** **explicitly APPROVED** the reviewer's S1–S8 recommendations and the O3 ruling **exactly as stated** ("I approve that block before Claude runs it"). The owner also accepted CARD-01-C1-PASS with notes N1–N6 and the pre-data model work.

The executor transcribes this ruling. It does not author, alter or fill any value.

**Commit order:** this file is committed **before** any network preflight or data access.

**Accepted pre-data facts:**
- The registered CPL diagnostic slope `wa/(1+w0) = 1.55590467` is correct for the frozen projection.
- `0.465` is a different, local-derivative diagnostic (unregistered).
- `1.37` remains **discarded**.

---

## S1 — `epsilon_null` — APPROVED (modified wording)

`epsilon_null = 0.03`.

This is an **operational, model-defined null neighborhood**: over `0 ≤ z ≤ 2`, it corresponds to `≲ 2.7 %` pointwise displacement in `w`. It is **not** a claimed observational-systematics floor; the executor's "systematics floor" language is deleted.

| Executor proposal | Reviewer recommendation | Owner |
|---|---|---|
| modified (wording) | approve value, delete floor language | **APPROVED** |

## S2 — Primary statistic — APPROVED

Profile likelihood:

`Delta_chi2_card(epsilon) = chi2_profile(epsilon) − chi2_profile(epsilon = 0)`

- The boundary-null significance uses the Chernoff reference `½ χ²₀ + ½ χ²₁`.
- Bayesian posterior results are **secondary only**.

| Executor proposal | Reviewer recommendation | Owner |
|---|---|---|
| adopted | approve | **APPROVED** |

## S3 — EXPLANATORY-SURVIVES — APPROVED (with precision)

EXPLANATORY-SURVIVES holds **only if all** of the following are true:
1. the best qualifying point has `epsilon < −0.03`;
2. `I_card = chi2_Lambda − chi2_card ≥ 4`;
3. `F = I_card / I_CPL ≥ 0.5`, where `I_CPL = chi2_Lambda − chi2_w0wa`.

`Delta chi2 = 3.84`, taken relative to the Card #1 branch best fit, is used **only** to define and report the descriptive 95 % profile set. It is **not** the boundary-null significance calibration.

| Executor proposal | Reviewer recommendation | Owner |
|---|---|---|
| modified (precision) | approve with precision | **APPROVED** |

## S4 — MODEL-KILLED — executor proposal REJECTED

- **MODEL-KILLED is not adjudicated in Card #1 v1.**
- The executor proposal (a relative `Delta chi2 = 25` to CPL ⇒ MODEL-KILLED) is **rejected**.
- **No relative-CPL statistic may trigger MODEL-KILLED.**
- If the whole branch performs very poorly relative to CPL, report the numerical tension as tension, without the MODEL-KILLED label.
- MODEL-KILLED would require a separately preregistered, calibrated absolute goodness-of-fit / posterior-predictive test. That test is not included in v1.

| Executor proposal | Reviewer recommendation | Owner |
|---|---|---|
| rejected | reject | **REJECTED / NOT ADJUDICATED** |

## S5 — Decision logic — APPROVED (with the S4 change)

Applied in order:

| Condition | Result |
|---|---|
| `I_CPL < 4` | `INCONCLUSIVE` |
| the 95 % Card #1 profile region lies entirely inside `−0.03 < epsilon ≤ 0` | `EXPLANATORY-KILLED` |
| `epsilon_hat < −0.03`, `I_card ≥ 4` and `I_card / I_CPL ≥ 0.5` | `EXPLANATORY-SURVIVES` |
| `I_card ≥ 4` but `I_card / I_CPL < 0.5` | `EXPLANATORY-KILLED` |
| otherwise | `INCONCLUSIVE` |

No MODEL-KILLED outcome is available in v1.

| Executor proposal | Reviewer recommendation | Owner |
|---|---|---|
| modified | approve with S4 change | **APPROVED** |

## S6 — Datasets — executor proposal MODIFIED

**Primary (verdict-bearing):** **DESI DR2 BAO + the DESI DR2 paper's baseline CMB.** The exact baseline CMB composition is:
- Planck low-ℓ TT: **Commander**;
- Planck low-ℓ EE: **simall**;
- Planck NPIPE PR4 high-ℓ TTTEEE: **CamSpec**;
- the **updated Planck + ACT DR6 lensing likelihood v1.2**.

**No SN sample is the sole primary verdict-bearing dataset.**

**Prespecified extensions** (run all three, separately):
- primary + Pantheon+;
- primary + Union3;
- primary + DESY5.

Report the Card #1 result state for each extension. **The primary DESI + CMB run controls the formal v1 verdict.** Do not switch the primary after seeing results.

| Executor proposal | Reviewer recommendation | Owner |
|---|---|---|
| modified (DESY5 removed from primary) | modify | **APPROVED AS MODIFIED** |

## S7 — Prior / range — APPROVED

- `epsilon ∈ [−1, 0]`, uniform, used only for the secondary posterior. **Profile the whole interval.**
- Positive control: `epsilon ∈ [0, 1]`, separate, labeled **CONTROL — NOT A GRUT CLAIM**.
- If the primary-branch optimum hits `epsilon = −1`, **do not accept a final verdict**; flag a prior/range-sensitivity repair.

| Executor proposal | Reviewer recommendation | Owner |
|---|---|---|
| adopted, plus boundary rule | approve | **APPROVED** |

## S8 — Diagnostic levels — APPROVED

- Report contours at 68 % and 95 %.
- Use 95 % for the Card #1 CPL-line overlay.
- This diagnostic **cannot** alter the D3 result state.

| Executor proposal | Reviewer recommendation | Owner |
|---|---|---|
| adopted | approve | **APPROVED** |

---

## O3 — Provenance ruling — APPROVED

**DESI DR2 BAO (D3).** The Cobaya-distributed DESI DR2 BAO likelihood is accepted as an **authorized equivalent distribution** only if **all** of the following hold:
1. the exact Cobaya version/commit is pinned;
2. the exact `bao_data` commit is pinned;
3. every consumed likelihood data file receives a SHA-256 hash;
4. its provenance is traceable to the DESI Collaboration DR2 BAO release/paper;
5. the implementation is the standard `bao.desi_dr2` BAO-only likelihood;
6. no arbitrary third-party recreation is substituted.

**Planck, ACT and SN likelihoods/data:**
- use collaboration-maintained or standard public likelihood distributions corresponding to the versions in the frozen S6 combination;
- pin software/data versions and hashes;
- record provenance;
- arbitrary third-party repackagings do not qualify.

**D1/D2 official `w0wa` chains:**
- only DESI Collaboration / official DESI-NERSC released chains qualify;
- if they are unreachable, report D1/D2 as unavailable;
- do not substitute reconstructed or mirror chains;
- D3 may proceed independently if all D3 likelihood components satisfy the provenance rules above.

---

## Authorized execution order (from the owner ruling)

1. Freeze the CAMB/Cobaya accuracy settings in a run-configuration commit, before inspecting any likelihood value.
2. Run the network/provenance preflight.
3. If a required D3 component fails provenance or access: **ACCESS-BLOCKED**; stop rather than improvise.
4. If the D3 components pass, acquire only the approved products.
5. Hash and record every product.
6. Run the primary DESI + CMB D3 profile.
7. Run the three prespecified SN extensions.
8. Run the CPL comparator under identical data combinations.
9. Run the positive-epsilon control separately, labeled CONTROL — NOT A GRUT CLAIM.
10. Attempt D1/D2 only if official DESI chains are reachable.
11. Apply S1–S8 mechanically.
12. Write `CARD_01_RESULT.md`.
13. Commit and push only `grut-conjecture-mode-1`.
14. Verify that all frozen, scientific and governance branches are unchanged.
15. Stop for owner review.

**Do not change any of the following after data access:**
- the postulate;
- the projection;
- the primary branch;
- the thresholds;
- the primary dataset;
- the result logic.

**Not authorized:** Card #2; any new conjecture; any PR; any merge. The charter is preserved unchanged.
