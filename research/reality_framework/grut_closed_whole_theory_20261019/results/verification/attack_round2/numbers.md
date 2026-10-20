# Attack round 2 — number-provenance lens

Returned by a separate agent as structured findings (the harness blocked report files). Rendered verbatim into Markdown for the record. Scratch paths refer to the session scratchpad and are not part of the repository.

Findings: 14 (2 major, 12 minor).

## 1. [major] WORKING_THEORY.md; OPEN_CLAIMS.md

**Location.** WT:216 (s08 row) and WT:332 "one-loop shift K₁/(gv) = 2.795 (statistical range ±0.11, not a certified enclosure)"; OC:76-77 "the ±0.11 allowance on K₁ is statistical, not a certified enclosure"

**Problem.** The round-1 F10 fix introduced this. The primary source says the opposite: the allowance combines analytic omitted-tail bounds with an explicitly empirical numerical allowance, and is "**not** a certified enclosure of the exact continuum saddle, a statistical confidence interval, or a bound on genuine higher-loop corrections". Round 1 truncated that quote. The packet now calls the allowance the very thing the source says it is not, both in the §7 Predictions paragraph and in OC Retained failures.

**Deciding check.** research/reality_framework/grut_candidate_1_energy_20261019/REPORT.md:20-24; SPECTRAL_NOTES.md:192 ("It is not a statistical confidence interval"); budget table REPORT.md:154-163 (0.03827594 + 0.05042782 analytic tails + 0.02 empirical + ..., "rounded upward"); ledger/ledger_s08_candidate1.md:79,81 (quoted correctly there)

**Fix proposed.** In all three places: "K₁/(gv) = 2.795 with a ±0.11 numerical allowance (analytic tail bounds plus an empirical allowance; not a certified enclosure, not a statistical confidence interval, not a bound on higher-loop error)".

## 2. [major] CALCULATIONS.md; OPEN_CLAIMS.md; WORKING_THEORY.md (results/verification/)

**Location.** CALC:253 "All records named below are in results/verification/" and CALC:261-266 (independent re-implementation); OC:111; WT:368-372

**Problem.** The shipped results/verification/independent_reimplementation_code.py (sha256 e63a7dad…, 11 992 B) is not the code the shipped record certifies. The shipped report, comparison_current.json (provenance.my_reimpl_sha256) and results JSON all name scratchpad wf2/reimpl/reimpl.py (sha256 a4740e24…, 10 450 B, 20 Oct 02:57). The shipped file is byte-identical to an earlier script, wf2/reimpl/reimpl_closed_whole.py (19 Oct 21:53). Run here, the shipped script writes a different schema and different last digits: E drift 7.105e-15 vs the record's 1.776e-14; Werner p 0.66666666733 vs 0.66666667377. The record therefore cannot be regenerated from the packet. compare.py, series.npz and xcheck_expm.py, which the report also cites, are not shipped. No headline number changes, since both scripts agree with the packet to within 1e-6.

**Deciding check.** sha256sum results/verification/independent_reimplementation_code.py = e63a7dad…; comparison_current.json provenance.my_reimpl_sha256 = a4740e24…; cmp with scratchpad wf2/reimpl/reimpl_closed_whole.py: identical; my run: scratchpad wf3/round2/numbers/reimpl_run/shipped_vs_record.txt (4 of 4 sampled values differ from the record)

**Fix proposed.** Ship reimpl.py (a4740e24), compare.py and series.npz under results/verification/ and re-run finalize.py. If the 19 Oct script is kept, label it as a separate earlier attempt.

## 3. [minor] WORKING_THEORY.md

**Location.** WT:209 (s01 row) "22–62 (GW dissipation) ... orders below their respective bounds", vs WT:252 "GW dephasing 22–62 orders below detectability"

**Problem.** The round-1 F9 fix introduced this. The calc's own caveat says the 22–62 figure "stands AS A DEPHASING STATEMENT" and that the amplitude/friction channel is "NOT SETTLED EITHER WAY". The s04 ledger lists "22–62 orders as a friction/amplitude statement" under must-not-promote. The two WT lines now disagree. "GW dissipation" is the register node's name, so the label has a source, but the scope it implies does not.

**Deciding check.** ledger/ledger_s04_grut_I_II.md:70 (calc/gw_dissipation_bounds.py:93-110), :139, :154; NO_GO_LEDGER.md:54 ("The dissipative dephasing is real but ...")

**Fix proposed.** "22–62 (GW dephasing only; the amplitude channel is not settled)" in the s01 row. Borderline major, because the label widens a computed negative beyond its scope.

## 4. [minor] WORKING_THEORY.md; CALCULATIONS.md; OPEN_CLAIMS.md (code/verify.py)

**Location.** WT:366-367 "--rerun ... compares every numeric leaf"; CALC:255-256; OC:110 "every numeric leaf identical to the committed results"

**Problem.** The statement is true today: my strict re-run finds 0 differing leaves and identical key sets. But verify.py, which the packet points to ("run it"), does not enforce it. Its comparator zips the two leaf streams and silently skips any pair whose keys differ, and zip truncates. A re-run with the whole S10 section missing passes (worst 0.0). So does a re-run with one extra leaf inserted at the front of S1 and every float doubled: worst 0.0, with only 16 of 5138 leaves compared.

**Deciding check.** code/verify.py:99-109; scratchpad wf3/round2/numbers/verify_comparator_test_output.txt; strict_compare_output.txt

**Fix proposed.** In verify.py compare dict(leaves()) key sets and FAIL on any missing or extra path. Alternatively, reword to "every leaf whose path matches".

## 5. [minor] WORKING_THEORY.md

**Location.** WT:375 "s01–s09 adversarially re-read by a second reader (8 to 15 corrections each ...)"

**Problem.** The entries' own correction lists give s01 10, s02 12, s03 8, s04 11, s05 13, s06 6, s07 12, s08 7, and s09 16 items (item 16 is "Integrity verified", so 15 corrections). The true range is 6–15, not 8–15.

**Deciding check.** Correction tables and lists: ledger_s01 §0 (10 rows), s02 §0 (12), s03 §0 (items 1-8), s04 §0 (1-11), s05 §0 (13 rows), s06 §0 (items 1-6), s07 §16 (12 rows), s08 §0 (items 1-7), s09 §0 (items 1-16)

**Fix proposed.** "6 to 15 corrections each".

## 6. [minor] CALCULATIONS.md; OPEN_CLAIMS.md (results/verification/)

**Location.** results/verification/independent_reimplementation_comparison_19oct_prereview.json, cited via CALC:266 and OC:111 ("independent_reimplementation_*")

**Problem.** The file is mislabelled. Its provenance records closed_whole_demo.py sha 6115d716…, which is commit f0e4cf0 ("apply independent code-review fixes", 20 Oct 02:57), i.e. the post-review-1 script. Only the marginal hash 1aad4b06… is pre-review. The pre-review demo is ba50f60e (da8f925). The file differs from _current.json only in these hashes.

**Deciding check.** git log hash map for code/closed_whole_demo.py: da8f925=ba50f60e, f0e4cf0=6115d716, f39cff2=4ebbd2db; diff of the two comparison JSONs (provenance lines only)

**Fix proposed.** Rename it (e.g. _20oct_postreview1.json) or state which script versions it compares.

## 7. [minor] CALCULATIONS.md

**Location.** CALC:283-284 "seven items (six prose defects and one robustness caveat), all addressed"

**Problem.** The round-1 F6 fix introduced this split. Review item 5 is assessed "accurate", so it is not a defect. The seven items are five prose defects (1, 2, 3, 4, 6), one confirmed-accurate statement (5) and one robustness caveat (7). The six fixes CALC lists as addressed map to items 1, 2, 3, 4, 6 and 7.

**Deciding check.** results/verification/code_review_calculation2.md:124-161 (item 5: "-- accurate")

**Fix proposed.** "seven items (five prose defects, one statement confirmed accurate, one robustness caveat); the six actionable ones all addressed".

## 8. [minor] CALCULATIONS.md; OPEN_CLAIMS.md; WORKING_THEORY.md

**Location.** CALC:265 "all agree to 10⁻⁶ relative or better"; OC:111 "23 comparisons agree to ≤ 10⁻⁶ relative"; WT:368-369 "re-implemented the calculation-1 quantities"

**Problem.** Row 2 (E_tot drift) is 1.776e-14 vs 7.105e-15, a relative difference of 1.5; the report itself compares it absolutely and says relative agreement is undefined for it. Row 19 is 0 vs 0. Rows 22-23 are the calculation-2 isotropic chain threshold, not calculation 1.

**Deciding check.** results/verification/independent_reimplementation_comparison_current.json rows 2, 19, 22, 23; independent_reimplementation_report.md:108

**Fix proposed.** "21 agree to ≤ 10⁻⁶ relative; the two numerical zeros agree absolutely (< 10⁻¹³)"; "calculation-1 quantities and the calculation-2 chain threshold".

## 9. [minor] WORKING_THEORY.md; CALCULATIONS.md

**Location.** WT:167 "commute with H to 10⁻¹⁴ ... ρ_A(t) − ρ_A(0) below 1.2·10⁻¹⁵"; CALC:185 "by at most 1.2·10⁻¹⁵"; CALC:171 "global purity constant to 10⁻¹⁴"

**Problem.** These bound-style statements are rounded in the wrong direction. The ground-projector ‖[F,H]‖ is 1.584e-14, which is not below 1e-14. The S9 maximum is 1.2213e-15, which is not "at most" or "below" 1.2e-15. The pure-state purity deviation is 1.343e-14. WT:153 uses the same phrase "to 10⁻¹⁴" as a strict bound (all three values there are < 1e-14), so the phrase means two different things. The magnitudes are trivial and the conclusions are unaffected.

**Deciding check.** closed_whole_demo.json S9/ground_projector_of_global_H/frobenius_norm_commutator_F_H = 1.58438e-14; S9/gibbs_beta_1_of_global_H/max_trace_distance_rhoA_t_vs_0 = 1.22128e-15; S8/pure_global_state_purity_const_dev = 1.3434e-14 (trace_numbers_r2_output.txt)

**Fix proposed.** "to 2·10⁻¹⁴" and "below 1.3·10⁻¹⁵" (or "≈ 1.2·10⁻¹⁵").

## 10. [minor] WORKING_THEORY.md

**Location.** WT:290-291 "the matched uniform-floor control retains better at every sampled time for L = 3, 4, 5"; WT:335-336 "retention probabilities that the matched uniform-floor control exceeds at every sampled time"

**Problem.** This contradicts the s10 ledger's correction 2 ("'Better at every sampled time' overstates the result"): the two agree at t = 0 (the control is −2.8e-15 lower at L = 3), and the control is strictly better only at the 200 samples with t > 0. It also contradicts WT's own s10 row (WT:218, "no worse at every sampled time ... within 10⁻¹⁰"). The wording is harsher on M1, not softer, but it does not match the record.

**Deciding check.** ledger/ledger_s10_m1_memory.md:20-23, 203

**Fix proposed.** "no worse at every sampled time (equal at t = 0) and strictly better at every sampled t > 0, for L = 3, 4, 5".

## 11. [minor] WORKING_THEORY.md

**Location.** WT:61-62 "The program's own register prices this at +16 net supplied inputs ..."; WT:220 "a 74-claim register (53 GRUT nodes and 21 vacuum-cluster nodes) ... net +16"

**Problem.** These counts are presented as current and undated. The ledger records them as "2026-08-17 sync stamps of one lineage's register", lists them under must-not-promote as a current universal, and states that the 53/21 split and the net +16 were not recomputed (only 74 was).

**Deciding check.** ledger/ledger_s12_program_nogo_postulates.md:39, 161 (POSTULATE_MAP.md:5-7)

**Fix proposed.** "+16 net supplied inputs (register sync stamp of 2026-08-17; not recomputed here)" and "53/21 (stamp; the total 74 recounted)".

## 12. [minor] OPEN_CLAIMS.md

**Location.** OC:81 "GRUT III: R1–R8, R13, R15, R19 rejections"

**Problem.** These labels occur in no ledger entry, nor in CAMPAIGN_LEDGER.md or WT, so they cannot be traced through the ledger. In the source they are numbered reviewer points that the authors marked CORRECT, not numbered rejections. Several are relabellings or restatements (R4 relabelled, R7 "compression" withdrawn, R8 P1 restated). R9–R12 were also accepted and are omitted without a stated criterion. Uncertain: the selection may be deliberate, but no criterion is given.

**Deciding check.** GRUT_III/RESPONSE_TO_REVIEWS.md:10-22, 183, 254; grep -rn 'R1–R8\|R19' ledger/ CAMPAIGN_LEDGER.md WORKING_THEORY.md: no hit

**Fix proposed.** Name them as reviewer points accepted in RESPONSE_TO_REVIEWS.md (R1–R8 and R13 retractions/relabellings; R15 v11 headline refuted; R19 v13 inference refuted) and cite the s05 ledger lines, or drop the labels.

## 13. [minor] WORKING_THEORY.md

**Location.** WT:213 (s05 row) "η = 1/4Għ"

**Problem.** Carried over from the round-1 table and not fixed. As written it reads (1/4)·G·ħ; the source prints 1/(4Għ).

**Deciding check.** ledger/ledger_s05_grut_III.md:149

**Fix proposed.** "η = 1/(4Għ)".

## 14. [minor] SOURCES.md (literature route for CALC S10 and C5)

**Location.** SOURCES.md:23-31 "Standard results used in FORMAL_CORE.md, cited there with journal references: ..."

**Problem.** The list was not updated for the round-1 additions. It omits references that FORMAL_CORE and CALC cite: Linden–Popescu–Short–Winter 2009 (CALC S10, COND 2.10), Reimann 2008, Short 2011, Short–Farrelly 2012, Popescu–Short–Winter 2006, Goldstein et al. 2006, Lieb–Robinson 1972, Wigner/Araki–Yanase/Ozawa, Araújo et al. 2013 (CALC:228), Grabert–Schramm–Ingold 1988 and Cotler–Penington–Ranard 2019. FORMAL_CORE:145,172,208-232 gives the journal references correctly.

**Deciding check.** grep -n 'Linden\|Araújo\|Lieb–Robinson\|Popescu' SOURCES.md: no hit in the list

**Fix proposed.** Add these to the SOURCES.md standard-results list.

## Confirmed by this lens

- verify.py (plain): 37 hash PASS (every MANIFEST file) + 30 check PASS = 28 headline + 2 provenance; '0 failure(s)', exit 0. This matches WT:365-366 and OC:109 "28 headline results and 2 provenance hashes" (round-1 F3 landed).
- verify.py --rerun: the above plus 4 PASS: re-run closed_whole_demo.py; re-run marginal_compatibility.py; closed_whole_demo.json worst rel dev 0.00e+00; marginal_compatibility.json worst rel dev 0.00e+00. '0 failure(s)', exit 0, wall 75.0 s. This matches CALC:257 "worst relative deviation 0".
- Strict independent re-run from scratch copies: key sets identical (5138/5138, 55/55 leaves); 0 differing leaves apart from wall_seconds; 5102 and 29 non-provenance numeric leaves.
- trace_numbers_r2.py ran 249 comparisons. 248 OK. The single 'BAD' is my own 1e-12 tolerance against 12-digit storage (0.338531123755 vs 0.338531123756), not a packet defect.
- CALC calc 1, all OK within stated rounding: S1 table (64 entries; sum rules E_A+E_Ā+E_int = E_tot to 2e-6 and p(Z=+1) = 0.5+E_A to 1e-6). Invariants 7.1e-15 / 4.0e-15 / 6.7e-16 / 1.260931 / 0.489575 / 2e-16. Reduced equation 1.9e-10 / 3e-16. S2: 13 maxima, 26 sign changes, 0.3383@2.1, 0.2353@3.6, 0.4893@37.5, 0.4042@39.5, 0.5418@37.6, 1.2035, 0.0145.
- CALC calc 1, continued: S3 3e-11 / 3e-15, condition numbers 1.1–2.3 at the 6 sample times, 15 generator norms, −0.0038 / −0.1590. S4 rows 1-3, plus the new row 4: ε = 1/512, ‖εχ‖₁ = 0.500, max D 4.16e-16, A-factor 0, [H,χ] = 0. The hand value ¼Z₀ (√2/4 = 0.353553) is re-derived. S5: 24 table values; g′ = 0.75; Bohr frequencies {0, 0.5, 1.5, 2, 2.5}; T = 4π; 7.6e-15; 0.4487. S6: all values, with joint marginals consistent. S7: series equals S1 S_A(40−t); 201/400. S8: 3e-16; 0.401@11.1; MI min 1.8e-15@0; slack 0/0.2115/0.2439.
- New S9 (rounded values OK; bound wording is N9): ‖[F,H]‖ 8.86e-15 and 1.58e-14; max D 1.22e-15; times {1, 5, 10, 40} per code:638.
- New S10, all OK: level spacing 5.898e-4; Bohr-gap separation 3.2027e-8 (1/δ = 3.1e7, so 'order 10⁷' holds); 3000 random t in [0, 1e4], seed 11 (code:666-667). d_eff 3.234/29.56; bounds 0.55605/0.18392, exactly ½√(4/d_eff); intermediate 0.298/0.095; sampled 0.287/0.088; max 0.515/0.242. The rounding in WT and FORMAL_CORE (3.2, 29.6, 3e-8) is consistent.
- CALC calc 2, all OK: C1 0 / 2.8e-17 / 3.5e-18; C2 infeasible (CLARABEL and SCS), 0.288675; C3 1.49999999899 and 0.750000000106; C4 0.666666673771 / 8.66e-4 and 0.33333333563 / 5.0e-4; C5 5.55e-17, odd cycle infeasible, even cycle optimal; C6 0.0500. 'Within 10⁻⁸ of analytic' confirmed (7.1e-9, 2.3e-9, 1.0e-9). Wall times 70.35 s and 0.7 s support 'about a minute' and 'under a second'.
- Round-1 F1 landed: CALC:267-270 and OC:112 label the re-execution as pre-review. The record is consistent: ba50f60e / 1aad4b06; 5030/24 leaves bit-exact; two runs bitwise identical (independent_reexecution_19oct_prereview.md:11-28).
- Round-1 F4 landed: OC 14 = 8 NOT FORCED (each citing calc 1 S2–S7/S9) + 3 BY DEFINITION + 3 SCOPE, matching the FORMAL_CORE §4 rows.
- Round-1 F5 landed: WT:165 reads 0.2115@3 and 0.2439@20, '0.20–0.40 over t ≥ 3'. The recomputed range is 0.2001 (t=3.5) to 0.3977 (t=37.8).
- Round-1 F7 and F8 landed. 23 comparison rows, hashes 4ebbd2db/9a633383 matched. The ε remark is dropped. Condition numbers '1.1–2.3 at the six sample times, at most 2.86 over the full 0.1 grid' match code_review_calculation1.md:35-36 (2.858 at t=17.3).
- Round-1 F2 landed only partly: the review records exist in results/verification/ and are in MANIFEST. Nine items D1–D9 (review 1) and 300 random restarts (review 2 :51-52, :175) confirmed.
- Round-1 F9 and F10 landed as text, but each introduced a new defect (N3, N1). Round-1 F6 landed the count but introduced a wrong split (N7). Round-1 F11: [[LEDGERVERIFY]], [[CRITIC]] and the SOURCES [[LITERATURE]] placeholder are filled. [[ATTACK_STATUS]] (WT:378), [[ATTACK_TABLE]] (OC:118) and the draft banner remain, as declared.
- Ledger-sourced numbers re-traced (OK-L): s01 7–47 / 22–62 / 62.7; s02 0.338, orders 6 and 8; s04 26 claims, lock −2, 62.7, 32σ→2.0σ, DESI 3.1σ at z 0.35–0.5; s06 N = 10¹⁸; s07 1.93%, 1411/640, N ≥ 2²⁸, ≥1/100 fail, ε ≤ 1/1000.
- Ledger-sourced numbers re-traced (OK-L), continued: s08 2.79528→2.795, −3553.70573→−3553.71, −0.199999152115, −3067.397336565, N = 39944; s09 ε₀ = 0.4725→0.47, Ω_K −4.540e-5 and −2.061e-9 (recomputed); s11 977 / 1.64e4 / 7.98e4 → 977 / 1.6e4 / 8.0e4; s12 74 = 53+21.
- s10 numbers: all 12 values in results/m1_matched_control_rerun.txt, plus 0.999680, 0.801481/0.751472 and 0.875748/0.764940, match ledger s10. WT commits 96a7e7e and 4e9a161 exist.
- Literature references in CALC are correct: Pollock et al. PRL 120, 040405 (2018); Chen–Ji–Kribs–Lütkenhaus–Zeng PRA 90, 032318 (2014); LPSW PRE 79, 061103 (2009) via FORMAL_CORE:208; Vorob'ev TPA 7, 147 (1962); Araújo et al. PRA 88, 022118 (2013); Pechukas 1994; Štelmachovič–Bužek 2001.
- Ledger-only note, outside the three audited files: ledger_s12:32 faithfully quotes a units slip in NO_GO_LEDGER.md:54, "~10²² (q=1) to ~10⁶² (q=2) orders below", which should read 22–62 orders.
- FINDINGS.md was not written: the harness refused the Write ('Subagents should return findings as text'). The traceability content (status of each round-1 fix, the per-number source mapping, and the findings) is returned in this structured output instead.
