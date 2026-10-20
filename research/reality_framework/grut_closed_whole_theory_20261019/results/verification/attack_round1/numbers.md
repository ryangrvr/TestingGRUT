# Number-provenance audit: packet `grut_closed_whole_theory_20261019`

Lens: NUMBER-PROVENANCE AUDITOR. Date 2026-10-20. Read-only; nothing under `/home/user/TestingGRUT` modified.
Files audited: WORKING_THEORY.md (WT), CALCULATIONS.md (CALC), OPEN_CLAIMS.md (OC). Sources: `results/closed_whole_demo.json`
(demo), `results/marginal_compatibility.json` (marg), `results/m1_matched_control_rerun.txt` (m1txt), `ledger/ledger_s*.md`,
handoff, repository, and the session scratchpad `wf2/` (where the verification records actually live).

Artifacts in this directory:
- `trace_numbers.py` / `trace_numbers_output.txt`: 201 exact comparisons of every JSON-sourced number quoted in CALC/WT/OC, all OK.
- `verify_plain.txt`: `python3 code/verify.py` -> 26 hash PASS + 27 check PASS, `0 failure(s)`, exit 0.
- `verify_rerun.txt`: `python3 code/verify.py --rerun` -> all above PASS plus `re-run closed_whole_demo.py` PASS, `re-run marginal_compatibility.py` PASS,
  `re-run closed_whole_demo.json agrees with committed results (worst rel dev 0.00e+00)`, same for marginal (0.00e+00); `0 failure(s)`, exit 0.
- `m1_control_myrun.txt` / `m1_control_time.txt`: my own run of `grut_memory_candidate_20261019_review/uniform_field_control.py`: output byte-identical
  to `control_results.txt` (diff empty); wall 10.4 s here (packet records 8.6 s; environment-dependent, not a discrepancy).

## Verdict in one line

Every number that the packet derives from its own two calculations traces exactly to the committed JSONs (201/201), both verifier modes pass,
and the M1 control re-run reproduces byte-for-byte. The defects are in the packet's **verification-record numbers**: the "5030 and 24 leaves bit-exact"
re-execution was performed on a superseded version of both scripts and both JSONs; the records behind "16 quantities", "2.86", "300 restarts",
"nine/six items" exist only in the session scratchpad; and four counts/ranges in prose (22 headline numbers; ten of fourteen; slack 0.21-0.24;
six review items) do not match what they summarize. No calculated number is wrong; no number differs from its source.

## Findings (defects)

| # | File:location | Problem | Deciding check | Severity | Minimal fix |
|---|---|---|---|---|---|
| F1 | CALC:215-216 "every numeric leaf of both JSON files bit-exact (5030 and 24 leaves)"; OC:92 "5030 + 24 numeric leaves bit-exact; statuses identical; two runs bitwise identical" | The independent re-execution was run on 19 Oct 21:38-21:43 against scripts sha `ba50f60e…` / `1aad4b06…` and JSONs sha `f7eab04b…` (100 188 B) / `1dbc40df…` (3 250 B). The committed packet files are `6115d716…` / `9a633383…` and `7f5b8df3…` (102 838 B, mtime 20 Oct 02:55) / `c3263d71…` (4 295 B, mtime 20 Oct 03:13), rewritten after the two code reviews (commits f0e4cf0, 28b3ca9, a7cdbea). The committed JSONs have 5069 and 29 non-provenance numeric leaves (5070/30 with wall_seconds). The claim therefore certifies a version that is not in the packet. | `wf2/rerun/REPORT.md:11-14` (hashes/sizes); `sha256sum results/*.json code/*.py`; `ls --time-style=full-iso results/`; `trace_numbers.py` leaf counts; `git log -- code/verify.py` | **major** | Re-run the independent comparison on the committed files and quote 5069/29 (or 5070/30), or label the record "pre-review version (19 Oct), superseded; current version verified only by `verify.py --rerun` in-session". My `--rerun` (same environment, not independent) gives worst rel dev 0.0 for both committed files. |
| F2 | CALC:213-242 "Verification record" and OC:86-95 table: "16 quantities … 10⁻⁶ relative", "worst condition number 2.86", "300 random restarts", "nine items, all addressed", "six items, all addressed", "two runs bitwise identical" | None of the four underlying records (`wf2/rerun/REPORT.md`, `wf2/reimpl/REPORT.md` + `comparison.json`, `wf2/review1/REVIEW.md`, `wf2/review2/REVIEW.md`) is in the packet, in `MANIFEST.json`, or anywhere in the repository; they exist only under `/tmp/claude-0/…/scratchpad/wf2/`. The owner cannot trace any of these numbers from the packet. | `grep -rn '5030\|2\.86\|16 quantities' /home/user/TestingGRUT` hits only the packet's own prose; `ls wf2/` shows the records | **major** | Copy the four reports (with their sha256) into `results/verification/` and add them to `MANIFEST.json`, or reword to "performed in-session; records not retained in the packet". |
| F3 | WT:338 "`code/verify.py` checks every packet hash and 22 headline numbers" | `code/verify.py` has 27 `check(` calls outside the hash loop (25 headline numeric/status checks + 2 provenance-hash checks). The count was 24 at commit 28b3ca9 and 27 since a7cdbea; the phrase "22 headline" dates from the first checkpoint da8f925 and was never updated. | `grep -c '^check(' code/verify.py` = 27; `git show 28b3ca9:…/verify.py | grep -c '^check('` = 24; `git log -S'22 headline'` → da8f925 | minor | "25 headline checks plus two provenance hashes" (or whatever the final count is at release). |
| F4 | OC:15-17 "Fourteen 'closure forces X' claims … ten of them are exhibited with numbers in calculation 1 (S2–S7)" | FORMAL_CORE §4 has 14 rows (correct). 9 rows cite calc 1; only 8 cite S2–S7 (row 1, "a pure global state", cites S1). "Ten" matches neither count. | `trace_numbers.py` row parse: 14 / 9 / 8 | minor | "nine of them are exhibited in calculation 1 (eight in S2–S7, one in S1)". |
| F5 | WT:155 "mixed trajectory: mutual information ≥ 0 and Araki–Lieb slack 0.21–0.24 away from saturation after t ≈ 3, no violation" | The JSON stores only two slack values (0.2115 at t=3, 0.2439 at t=20). Recomputed from the stored S_A, S_Ā, S_Ω series, the slack over t ≥ 3 ranges from 0.2001 (t=3.5) to 0.3977 (t=37.8); 305 of the 371 grid points after t=3 exceed 0.24 and 10 are below 0.21. Two samples were written as a range. The "no violation" conclusion is unaffected. | `trace_numbers.py` AL-slack line | minor | "slack 0.2115 at t = 3 and 0.2439 at t = 20 (0.20–0.40 over t ≥ 3), no violation". |
| F6 | CALC:223-224 "Adversarial code review of calculation 2: … six items, all addressed" | `wf2/review2/REVIEW.md` §(e) lists **seven** numbered items; the seventh ("Robustness … bare solver status") is precisely what the second solver and violation certificates that CALC lists as addressed respond to. | `sed -n 124,162p wf2/review2/REVIEW.md | grep -c '^[0-9]\+\. '` = 7 | minor | "seven items (six prose defects and one robustness caveat), all addressed". |
| F7 | CALC:218 "16 quantities (E_tot and its drift, S_A at t = 10 and 40 and its maximum, the backflow sum and D's minimum, the three record probabilities, the identical-marginal trace distance, the Werner chain threshold)" | The re-implementation report tabulates **23** comparisons ("All 23 comparisons agree"); the parenthetical names 12. "16" matches neither. Also, the report's provenance block records `marginal_compatibility.py` sha `1aad4b06…` (pre-review), so its Werner comparison was against the earlier marginal JSON (the value 0.666666673771 is unchanged, so no numeric consequence). The stated offset explanation "ε = 1/512 versus the numerically halved ε_max" is not in the report (it explains offsets as 12-digit rounding, numerical zeros compared absolutely, and solver tolerance); plausible but unrecorded. | `wf2/reimpl/REPORT.md` table (23 rows) and provenance section | minor (uncertain on the ε remark) | Quote 23 (or list the 12 named), cite the record (see F2), and drop or source the ε remark. |
| F8 | CALC:240 "the conditioning of the map inversion (worst condition number 2.86)" vs CALC:74 "condition numbers 1.1–2.3 at the sample times" | Both are right but refer to different sets: JSON sample times give 1.114–2.350; review1 computed over the full 401-point grid "max cond = 2.8582 at t=17.3". Line 240 omits "over the whole grid", so the two numbers read as contradictory. | `wf2/review1/b_out.txt:8-9`; JSON `S3/tcl/map_condition_number_thermal` | minor | "worst condition number over the full 0.1 grid 2.86 at t = 17.3 (1.1–2.3 at the six sample times)". |
| F9 | WT:193 (s01 row) "every differentiator fell 7 to 63 orders below detectability" | Ledger s01:57,70 gives three different figures with different metrics: energy-basis decoherence "7–47", GW dissipation "~22–62 (printed 21–62 elsewhere)", Γ_T "62.7"; ledger s04:74 notes these are "a different metric … not reconciled". Collapsing them into "7 to 63" (63 = 62.7 rounded) presents an unreconciled set as one range. | `ledger/ledger_s01_v7_ctp.md:57,70`; `ledger/ledger_s04_grut_I_II.md:74` | minor | "7–47 (decoherence), 22–62 (GW dissipation) and 62.7 (Γ_T) orders below their respective bounds (metrics not reconciled)". |
| F10 | WT:200 (s08 row) "one-loop shift K₁/(gv) = 2.795 ± 0.11" | Traces to ledger s08:66-71 and the Candidate-1 energy REPORT, but the ledger records the ±0.11 as "not a certified enclosure … a statistical confidence" with range [2.68528, 2.90528]; WT drops that qualifier. Number correct, boundary softened. | `ledger/ledger_s08_candidate1.md:66-68,71` | minor | "2.795 (statistical range ±0.11, not a certified enclosure)". |
| F11 | OC:97-99 "[[LEDGERVERIFY]]", "[[CRITIC]]", "[[ATTACK]]"; SOURCES.md:33 "[[LITERATURE]]"; WT:1-2 draft banner | Three cells of the OC verification table and the literature-inspected line are unfilled placeholders, so three of the nine verification rows have no result to trace. | inspection | minor (blocking only if the packet is released with them) | Fill or delete before removing the draft banner. |

Not flagged (checked and consistent, listed so the reader knows they were tested): the C2 "Frobenius violation" 0.2887 and C6 0.0500 are the
epigraph SDP value min_ρ max_X ‖Tr_X ρ − target_X‖_F (`code/marginal_compatibility.py:114-124`), which the label "frobenius" describes correctly
(0.0500 for C6 is attained by a non-product correction ½1_A ⊗ diag(−0.05, 0.05)); the "7 to 63 orders" rounding of 62.7 is arithmetic-consistent
(F9 is about mixing metrics, not rounding); the M1 wall time 8.6 s vs 10.4 s here is environment variance.

## Traceability table

Status key: OK = matches source within stated rounding; OK-L = traced to a ledger entry that cites a repository path:lines (the ledger citation
was read; the underlying repository file was not re-opened unless noted); DERIVED = recomputed by me; STALE/UNTRACED = see finding.

### CALCULATIONS.md, calculation 1 (source: `results/closed_whole_demo.json`)

| CALC line | Number (unit, context) | Source key | Source value | Status |
|---|---|---|---|---|
| 20-27 | 8 qubits, dim 256; A = qubit 0; ω₀ = 1; J = 1, J_z = 0.5, h = (0.3, −0.7, 0.5, −0.2, 0.8, −0.4, 0.6); g = 0.5; β = 1; grid [0,40] step 0.1 | `datum/*` | identical strings/values | OK |
| 37-44 | S1 table, 8 rows × 8 columns (E_tot … p(X=+1)), 6 d.p. | `S1/{E_tot,E_A,E_Abar,E_int,S_A,MI_A_Abar,pZ_plus,pX_plus}[idx(t)]` | all 64 entries match to 5·10⁻⁷ | OK |
| 46 | max\|ΔE_tot\| 7.1·10⁻¹⁵; max\|ΔS(ρ_Ω)\| 4.0·10⁻¹⁵; S(ρ_Ω) 1.260931 | `S1/invariants/{max_abs_dE_tot,max_abs_dS_tot,S_tot}` | 7.1054e-15; 3.9968e-15; 1.26093113354 | OK |
| 47 | max\|Δ Tr ρ²\| 6.7·10⁻¹⁶; Tr ρ_Ω² 0.489575 | `S1/invariants/max_abs_dpurity_tot`; `S1/purity_tot[0]` | 6.6613e-16; 0.489575386273 | OK |
| 48 | Born local vs global 2·10⁻¹⁶ | `S1/born_global_minus_local_max` | 2.2204e-16 | OK |
| 50 | S(ρ_A) < ln 2 = 0.693; first exceeds 0.4 at t = 8.8 | `S1/invariants/S_A_max_possible`; scan of `S1/S_A` | 0.693147; t = 8.8 (0.41012) | OK / DERIVED |
| 52-56 | dt = 10⁻⁴; residual 1.9·10⁻¹⁰; H_Ā commutator 3·10⁻¹⁶ | `code:251`; `S1/reduced_equation_check_t3/*` | dt = 1e-4; 1.9218e-10; 2.9623e-16 | OK |
| 61 | 13 local maxima; 26 sign changes | `S2/local_maxima_times` (len), `S2/number_of_sign_changes…` | 13; 26 | OK |
| 62-63 | 0.3383 at t=2.1; 0.2353 at t=3.6; 0.4893 at t=37.5; 0.4042 at t=39.5 | `S2/first_local_max`, `first_local_min_after_it`, `S_A_max`, `t_of_S_A_max`, `S_A_min_after_max`, `t_of_min_after_max` | 0.338264/2.1; 0.235332/3.6; 0.489339/37.5; 0.404187/39.5 | OK |
| 64-67 | D min 0.5418 at 37.6; backflow 1.2035; largest step 0.0145 | `S2/BLP/*` | 0.541801/37.6; 1.203537; 0.014545 | OK |
| 73 | min Choi 3·10⁻¹¹; TP 3·10⁻¹⁵ | `S3/cptp_check/*` | 2.9575e-11; 2.7200e-15 | OK |
| 74 | condition numbers 1.1–2.3 at sample times | `S3/tcl/map_condition_number_thermal` | 1.1137 … 2.3495 | OK (see F8) |
| 79-83 | 15 generator norms (t = 1, 2, 5, 10, 20 × 3 columns) | `S3/tcl/generator_*` | all 15 match to 5·10⁻⁵ | OK |
| 86-87 | V(t+0.5,t): −0.0038 at 0.5, −0.1590 at 5; τ = 0.5 | `S3/cp_divisibility/*` | −0.003825; −0.158979; 0.5 | OK |
| 97 | ε = 1/512 (= ½ε_max); ‖εχ‖₁ 0.500; 0.1025 at 0.7; 0.3536 | `S4/inf_T_bath_quantum_correlation_Y0X1/*` | 0.001953125 (ε_max 0.00390625); 0.5000; 0.102487/0.7; 0.353553 | OK |
| 98 | 0.0097 at 14.1; 0 | `S4/inf_T_bath_classical_correlation_Z0Z1/*` | 0.009718/14.1; 0.0 | OK |
| 99 | ε 0.25; 0.293; 0.0544 at 0.7; 0.1925 | `S4/thermal_beta1_bath_dressed_quantum_correlation/*` | 0.25 (ε_max 0.5); 0.293393; 0.054445/0.7; 0.192540 | OK |
| 103 | 2gε·2⁷ Z₀ = ¼Z₀, ‖·‖_F = √2/4 = 0.35355 | hand: [gX₀X₁,Y₀X₁] = 2igZ₀, 2·0.5·(1/512)·128 = 0.25, 0.25·√2 | 0.353553 | DERIVED, OK |
| 118-121 | S5 table 4 rows × 6 numbers | `S5/trace_distance_rhoA/*`, `S5/S_A/*` | all 24 match | OK |
| 123-125 | g′ = 0.75; T = 4π = 12.566; 7.6·10⁻¹⁵; 0.4487; Bohr multiples of 0.5 | `S5/model_b_prime_exact_recurrence/*` | 0.75; 12.56637; 7.615e-15; 0.448733; [0, 0.5, 1.5, 2.0, 2.5] | OK |
| 130-132 | t₁ = 3, t₂ = 6; 0.354910 (×2) | `S6/construction`, `born_pZ_plus_at_t1`, `register_population_0_after_write` | 3/6; 0.354909508 (×2) | OK |
| 136 | 0.2155, 0.1230, 0.1394, 0.5220 | `S6/joint_distribution_zA_t2_zR` | 0.215485, 0.123046, 0.139424, 0.522045 | OK |
| 137-140 | 0.338531; 0.449147; 0.288331; 0.338531 | `S6/pZ_plus_A_at_t2_*` | 0.338531124; 0.449147358; 0.288331120; 0.338531124 | OK |
| 143-148 | 7·10⁻¹⁵ at t=10; κ = 0.3; 0.0858/0.2273/0.2505 at 4/6/10; readout 1.000000; 0.9607/0.9672/0.9613/0.9700 at 4/6/8/10 | `S6/record_persistence/*`; `code:516` | 7.073e-15; 0.3; 0.085754/0.227307/0.250487; 1.0; 0.960744/0.967172/0.961338/0.969999 | OK |
| 153-158 | 5·10⁻¹⁵; 5·10⁻¹⁵; 0.4198→0.4005→0.4484→0.3337→0.3929→0.0000; 201 of 400; 0.4893 at 2.5 | `S7/*` | 5.051e-15; 4.930e-15; 0.419838/0.400533/0.448388/0.333707/0.392871/0.0; 201; [0.489339, 2.5] | OK |
| 163-170 | 3·10⁻¹⁶; 10⁻¹⁴; 0.401 at 11.1; S(ρ_Ω) 1.2609; min MI 1.8·10⁻¹⁵ at t=0; slack 0 / 0.2115 / 0.2439 | `S8/*`, `S8/mixed_S1_trajectory/*` | 3.331e-16; 1.343e-14; 0.401232/11.1; 1.260931; 1.776e-15 (argmin t=0); 0.0/0.211540/0.243895 | OK |

### CALCULATIONS.md, calculation 2 (source: `results/marginal_compatibility.json`)

| CALC line | Number | Source key | Source value | Status |
|---|---|---|---|---|
| 176 | 0; 3·10⁻¹⁷; 3·10⁻¹⁸ | `C1_nested_restriction/*` | 0.0; 2.776e-17; 3.469e-18 | OK |
| 177 | infeasible (CLARABEL, SCS); 0.2887 | `C2_bell_chain/*` | infeasible/infeasible; 0.288675 | OK |
| 178 | 1.5000; 0.7500 | `C3_quantitative_obstruction/*` | 1.49999999899; 0.750000000106 | OK |
| 179 | 2/3 (0.666667); 8.7·10⁻⁴ | `C4_werner_thresholds/chain_*` | 0.666666673771; 0.000866024 | OK |
| 180 | 1/3 (0.333333); 5.0·10⁻⁴ | `C4_werner_thresholds/triangle_*` | 0.33333333563; 0.000499998 | OK |
| 181-183 | 6·10⁻¹⁷; infeasible; optimal (½ on 0101, ½ on 1010) | `C5_classical_contrast/*` | 5.551e-17; infeasible; optimal | OK |
| 184 | diag(0.6,0.4) vs 1/2; infeasible (both); 0.0500 | `C6_inconsistent_overlap/*` | [0.6,0.4]/[0.5,0.5]; infeasible/infeasible; 0.0499999988 | OK |
| 186-188 | fidelity 1/4; eigenvalue 3/2; F = 3/4; state (2\|000⟩+\|110⟩+\|011⟩)/√6 | analytic; `C2/analytic_proof`; `C3` | consistent | OK (analytic, not re-derived here) |
| OC:21 | thresholds within 10⁻⁸ of analytic | computed: 7.1e-9, 2.3e-9, 1.0e-9 | — | DERIVED, OK |

### CALCULATIONS.md, verification record (lines 213-245) and OPEN_CLAIMS.md verification table (86-99)

| Line | Number | Source | Status |
|---|---|---|---|
| CALC:216 / OC:92 | 5030 and 24 leaves bit-exact | `wf2/rerun/REPORT.md:23-24` (19 Oct, pre-review scripts/JSONs) | STALE (F1): committed JSONs have 5069 / 29 non-provenance numeric leaves |
| CALC:215 / OC:92 | two runs bitwise identical | `wf2/rerun/REPORT.md:28-30` (pre-review) | UNTRACED in packet (F2); superseded version (F1) |
| CALC:218-222 / OC:93 | 16 quantities, ≤ 10⁻⁶ relative | `wf2/reimpl/REPORT.md` (23 comparisons) | UNTRACED in packet (F2); count mismatch (F7) |
| CALC:223-231 / OC:95 | six items; 300 random restarts | `wf2/review2/REVIEW.md` (7 items; 300 starts at :51-52,175) | UNTRACED in packet (F2); count (F6) |
| CALC:232-242 / OC:94 | nine items; worst condition number 2.86 | `wf2/review1/REVIEW.md` (D1–D9); `wf2/review1/b_out.txt:9` (2.8582 at t=17.3, full grid) | UNTRACED in packet (F2); nine OK; 2.86 OK but scope unstated (F8) |
| CALC:244 / m1txt:4 | 8.6 s | m1txt | OK (my run 10.4 s, same output) |
| OC:90-91 / WT:338 | verify.py all pass; --rerun every leaf within 10⁻⁶ (bit-exact) | my runs: 0 failures; worst rel dev 0.00e+00 | OK (confirmed) |
| WT:338 | 22 headline numbers | `code/verify.py` | STALE (F3): 27 checks (25 + 2 provenance) |
| OC:97-99 | [[LEDGERVERIFY]], [[CRITIC]], [[ATTACK]] | — | placeholders (F11) |

### WORKING_THEORY.md §3-§4 (source: demo / marg JSON)

| WT line | Number | Source | Status |
|---|---|---|---|
| 135-137 | eight qubits; ½Z₀; J = 1, J_z = 0.5; 0.5 X₀X₁; e^{−H_Ā}/Z; [0,40] | `datum` | OK |
| 143 | −10.228789; 1.260931; 0.489575; constant to 10⁻¹⁴ | S1 invariants (7.1e-15, 4.0e-15, 6.7e-16 all < 1e-14) | OK |
| 144 | E_A ∈ [−0.169, 0]; E_int ∈ [−0.223, 0] | `S1/invariants/range_E_A` −0.169286, `range_E_int` −0.222889 | OK |
| 145-146 | 2·10⁻¹⁶; 1.9·10⁻¹⁰; 3·10⁻¹⁶ | S1 | OK |
| 147 | 26 sign changes; backflow 1.20 | S2 | OK |
| 148 | generator changes by up to 0.97 | `S3/tcl/…K(t)-K(0.5)` max 0.968867 | OK |
| 149 | Choi down to −0.159 | S3 most_negative −0.158979 | OK |
| 150 | 0.1025; ¼Z₀ exactly | S4; hand | OK |
| 151 | up to 0.79; T = 4π | S5 a_vs_c(40) 0.793045; period 12.566 | OK |
| 152 | 0.354910; 0.3385; 0.4491; 0.2883 | S6 | OK |
| 153 | readout 1.000; 0.96–0.97 | S6 record_persistence (1.0; 0.9607–0.9700) | OK |
| 154 | 5·10⁻¹⁵; 0.4198 → 0 | S7 | OK |
| 155 | 3·10⁻¹⁶; slack 0.21–0.24 after t ≈ 3 | S8 (two samples 0.2115, 0.2439); series range 0.2001–0.3977 | WRONG RANGE (F5) |
| 28 / CALC:10 / README | reproducible in about a minute | `provenance/wall_seconds` 46.03 + 0.7 | OK |
| 162-169 | 0.289; 3/2; p ≤ 2/3; p ≤ 1/3 | marg C2/C3/C4 | OK |
| 177-180 | 0.338 (fourth cumulant) | ledger s02:85,207,327 → `P3_NC_LIFT_VERDICT_01.md:57-67` | OK-L |

### WORKING_THEORY.md §0, §5, §6, §7 and OPEN_CLAIMS.md ledger-sourced numbers

| WT/OC line | Number | Ledger entry → repository citation | Status |
|---|---|---|---|
| WT:7 | commits 96a7e7e, 4e9a161 | `git log`: both exist (handoff; M1 correction); cc24967, 1997d20 also exist | OK |
| WT:56-57 | +16 net supplied inputs; empty derived tier; empty prediction ledger | s12:39,82,113,121 → `NO_GO_LEDGER.md`, `GRUT_FINAL_SYNTHESIS_01.md:222-226,235-236`; s12:39 notes "+16 … not recomputed" (74-claim count was) | OK-L (stamp, not recomputed) |
| WT:193 (s01) | 7 to 63 orders | s01:57,70 (7–47; 22–62; 62.7) | OK-L but mixes metrics (F9) |
| WT:194 (s02) | 0.338; orders 6 and 8 | s02:85,95,207 (0.338; 17.2; 54.9; order-8 gate 67.8 "red") | OK-L; OC:68 "order-8 gate red" consistent |
| WT:196 (s04) | derived 0 across 26 V1 claims; lock −2; Γ_T 62.7 orders; 32σ → 2.0σ | s04:47,50,100 (26); s04:16,57,93 (P0s/P2 = −2); s04:85,103 (6.189e-63·H0, 62.7); s04:15 (2.0σ, 1.97 central) | OK-L |
| WT:197 (s05) | η = 1/4Għ; five postulates | s05:17,102,114 (1/(4Għ)); s05:21 (P0,P2,P3,P4,P5 = 5) | OK-L (typography "1/4Għ" should read 1/(4Għ)) |
| WT:198 (s06) | ‖ρ₀‖_∞ ≤ e^{cN}/D; O(N); αN log₂N bits | s06:25,42,44,122 | OK-L |
| WT:199 (s07) / OC:61 | 1.93%; 1411/640; N ≥ 2²⁸ | s07:178,310,374,473; s07:270-272 → `FEASIBILITY.md:49-50`; s07:71,91,273 | OK-L |
| WT:200 (s08) / WT:311-312 | Q ≈ 39 944; 2.795 ± 0.11; −3553.7; −0.199999152115; −3067.397336565; "+2.80", "−3554" | s08:66-71,80-82 → Candidate-1 energy `REPORT.md:35-42,154-163` and JSON `renormalized_K1_over_gv 2.79527765069`, `DeltaE_freeX_over_gv −3553.70573` | OK-L (qualifier dropped, F10); handoff values consistent |
| WT:201 (s09) | T_eff = H(1−ε)/2π; Ω_K = −1/(e^{2N}−1); −4.5·10⁻⁵ (N=5); −2.1·10⁻⁹ (N=10) | s09:51 → `grut_closure_fork_20261019/DERIVATION.md:146-148`, `results.json:382-416` (−4.5402e-5, −2.0612e-9); recomputed by me: −4.5402e-5, −2.0612e-9 | OK-L, DERIVED |
| WT:202 (s10) / OC:56-58 | 0.987510; 0.935065; 0.778468; ≥ 0.999680; 0.998458; 0.9 standard (L = 4, 5); flips 0.1 and 0.3 | m1txt:9-10 (and `control_results.txt`, byte-identical, my re-run identical); s10:78 (0.999680); s10:70 (0.801481, 0.751472 < 0.9); s10:72,135 (κ = 0.1 → 0.875748, κ = 0.3 → 0.764940) | OK / OK-L |
| handoff table | 12 retention numbers L = 3, 4, 5 | m1txt:7-12 | OK (all 12 match) |
| WT:203 (s11) | 977 → 1.6·10⁴ → 8.0·10⁴ sweeps at β = 1.0; one sweep; 18-qubit OTOC | s11:66 (977.37; 16431.6; 79792.2 → `tb:REPORT.md:39-43`); s11:101,111; s11:34,103 | OK-L (rounding consistent) |
| WT:204 (s12) | 74-claim register; zero banked no-gos; +16 | s12:39 (74 recomputed; "+16 not recomputed"), s12:121 | OK-L |
| WT:235 | 689 Hz; 32σ → 2.0σ; 62.7; 22–62 | s01:57,70,98; s04:15,70,103 | OK-L |
| WT:233 | memory time 1/(4H) | s05:41,115 | OK-L |
| WT:276 | flips 0.1 and 0.3; 0.9 standard | s10:72,135; s10:70 | OK-L |
| WT:291-292 | DESI DR2 + CMB, z ≈ 0.35–0.5, 3.1σ | s04:146 | OK-L |
| OC:15-17 | fourteen; ten (S2–S7) | FORMAL_CORE §4: 14 rows; 9 cite calc 1; 8 cite S2–S7 | MISCOUNT (F4) |
| OC:22 | twelve sectors | `ledger/` has 12 files; WT §5 has 12 rows | OK |
| OC:18-19 | 0.289; 3/2; 2/3; 1/3 | marg | OK |
| SOURCES:18 | V7: 332 lines, 18 380 bytes, April 2026 | handoff `NEXT_CHAT_HANDOFF.md:175-177` | OK (handoff only; file absent from repository, as stated) |

## Executed commands (for the record)

```
python3 code/verify.py                      # 0 failure(s)   -> verify_plain.txt
python3 code/verify.py --rerun              # 0 failure(s), worst rel dev 0.00e+00 both files -> verify_rerun.txt
python3 -I trace_numbers.py                 # 201 OK, BAD: none -> trace_numbers_output.txt
python3 .../grut_memory_candidate_20261019_review/uniform_field_control.py  # diff vs control_results.txt empty; 10.4 s
sha256sum results/*.json code/*.py; ls --time-style=full-iso results/; git log -- code/verify.py; git show 28b3ca9:.../verify.py
```
