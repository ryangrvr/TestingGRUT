# V0-4 — COMPARISON PHASE (EDA-01; the original was unsealed only after the blind second reading was committed at `f181330`)

**Original:** `playground/SCOUT_0/probes/EDGE_DATA_AUDIT_01.md` (blob `a9acfd75…`); matrix in §1, tally in §1, witnesses
in §3, theorem candidate in §4.

**Second reader:** `V0_4_EDA01_SECOND_READER.md` and `code/v0_4_sr_*`. Produced by a context-isolated sub-agent and
**not** edited below.

## 1. Row-by-row comparison

| row | original | second reader | relation |
|---|---|---|---|
| E-1 | READ | READ / CON | **compatible** (dual) |
| E-2 | CON | IND (READ if inverse moments count) | **DISCORDANT** (CON vs IND) |
| E-3 | CON (definitional: makes the spectral representation exist) | IND | **DISCORDANT** (CON vs IND) |
| E-4 | READ | CON (dual READ) | **compatible** (dual) |
| E-5 | IND | IND (READ on two entries) | **compatible** |
| E-6 | IND | IND | agree |
| E-7 | IND | IND | agree |
| E-8 | CON (the one earned member selection; picks which μ_r is read) | IND ("gap 0.678" is a selection margin, not a spectral gap) | **DISCORDANT** (CON vs IND) |
| E-9 | READ (cone speed = maximal group velocity) | CON (sign J ≥ 0, floor ν ≥ ħJ/2) | **DISCORDANT**: different datum identified |
| E-10 | IND | IND / SUP | **compatible** |
| E-11 | SUP | SUP (dual READ) | **compatible** |
| E-12 | READ | READ | agree |
| E-13 | READ (datum: "edge velocity …; soft-point count (P-02 theorem)") | READ / CON, with the **full edge-data set** as the datum | **label compatible; datum description corrected** (follows V0-3) |
| E-14 | READ | READ | agree |
| E-15 | READ (moments; nonlinear class) | IND (READ if moments count) | **compatible**: the original explicitly adds low moments as an equivalent invariant; the reader flags the same fork |
| E-16 | IND | IND | agree |
| E-17 | NF | NF / IND | **compatible** |
| E-18 | IND | IND | agree |
| E-19 | READ | READ | agree |
| E-20 | IND | IND | agree |
| E-21 | SUP ("amplitude and state sector-supplied") | CON (exponent class restricted; the observed exponent depends on a supplied amplitude, witness W3) | **DISCORDANT** (SUP vs CON) |
| E-22 | IND / SUP (tally: IND) | CON (J > 0 above 4T against a zero classical floor) | **DISCORDANT** (IND vs CON) |

**Agreement on FIX status: 22 / 22** (no row is FIX for either reader).

**Label compatibility:**
- 16 / 22 compatible: 10 exact and 6 dual / justified-fork;
- 6 / 22 discordant: E-2, E-3, E-8, E-9, E-21, E-22.

**Every discordance lies at the CON boundary** (CON vs IND, or CON vs READ / SUP). None touches FIX.

## 2. Headline, theorem candidate and witnesses

**Headline verdict:** the original says "No E-entry is FIXES-EDGE-DATA"; the second reader says "no row is FIX".
**AGREEMENT.**

**The second reader adds:**
- the reasons near-misses fail the two-model test (E-13, E-4, E-21, E-12, E-1);
- explicit numerical witnesses:
  - W1: the E-13 SF-1 predicate with three parents, different edge data and classes;
  - W2: reproduces the record values r = 0.582109, X_J, κ = 1/(2√2.3) exactly, then moves them by varying pin;
  - W3: E-21 masking, giving observed exponent 3 vs 7;
  - W4: carrier cone edges.
- These are **independent of, and consistent with**, the original's W-A / W-B / W-C. The original's W-C (κ moves with
  pin) and the reader's W2 agree in substance.

**The original's "strongest obstruction" on E-8 (potential CON → FIX)** is weakened by the reader's finding that the
E-8 margin is a selection margin, not spectral edge data. Both readers agree that E-8 is not FIX at the recorded scope.

## 3. Corrections and findings

- **CR-V4-1 (row discordance, grade-bearing under the pre-fixed rule).** Six rows are materially reclassified by an
  independent reader. The disagreements are systematic: **the rubric's CON label is under-specified.**
  - Does making a spectral representation exist (E-3), or restricting support (E-2), count as "constraining edge data"?
  - Does a selection margin (E-8) count?
  - Does a sign / floor inequality (E-9, E-22) count?
  - Does a per-branch exponent restriction (E-21) count?

  **The original matrix's CON / IND / SUP / READ assignments are therefore reader-dependent at these six rows. The FIX
  column is reader-independent.**
- **CR-V4-2 (follows V0-3).** E-13's datum must be stated as the **edge-data set** (positions, orientations, local
  orders r_e, soft-point differences), not "edge velocity at the occupation edge; soft-point count (P-02 theorem)".
  - The original's audited data list should add the **edge order / jet r_e**.
  - The §4 theorem candidate's "edge data (location, exponent, velocity, visible sign, soft points)" should include the
    local order.
- **CR-V4-3 (rubric).** FIX needs an explicit statement of the class over which the primitives are held. The reader shows
  that a literal reading ("no supplied primitive varied") would make any computation a FIX. Low / inverse moments should
  be an explicit audited data class. The original found them, but outside its declared list.

## 4. Grade (by the rule fixed in the spec before the reading)

**V0-4-C — EDA-01 SECOND-READER AUDIT: HEADLINE VERDICT CONFIRMED (NO FIX); ROW MATRIX REPRODUCED WITH CORRECTION —
VER-I1 (orchestrator-exposed).**
- **Headline (no earned result FIXES edge data):** independently **confirmed**, with independent witnesses.
- **Row matrix:** 16 / 22 compatible and 6 / 22 discordant, all at the CON boundary, so CON is under-specified.
- **E-13:** re-stated against the corrected quotient; READ, with CON secondary.

**Independence:** VER-I1, orchestrator-exposed. The orchestrator had read the original matrix while extracting the
target. The second reader was context-isolated, read only allowed record files, did not open the deny-list, and gave a
file-access report. **Not VER-I2.**

**Frozen `scout-0` is not modified.** Criterion-2 item 5 ("EDA-01's row classification checked by a second reader") is
**discharged at VER-I1, with correction** (pending owner ruling).

## 5. Owner ruling (additive; review of `3770ce7`)

**V0-4 ACCEPTED.** **Final grade: V0-4-C — EDA-01 SECOND-READER AUDIT: HEADLINE NO-FIX VERDICT CONFIRMED; ROW MATRIX
REPRODUCED WITH CORRECTION — VER-I1 (orchestrator-exposed).**
- Criterion-2 item 5: **COMPLETE AT VER-I1 WITH CORRECTION.**
- CR-V4-1, CR-V4-2 and CR-V4-3 are accepted.
- **No further EDA-01 classification audit.**
- **Robust result:** FIX = 0 for all 22 rows under both independent readings. The six CON / IND / READ / SUP boundary
  disagreements remain explicitly reader-dependent and do not affect the headline.
