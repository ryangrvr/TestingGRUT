# Attack round 2 — overclaim lens

Returned by a separate agent as structured findings (the harness blocked report files). Rendered verbatim into Markdown for the record. Scratch paths refer to the session scratchpad and are not part of the repository.

Findings: 30 (11 major, 19 minor).

## 1. [major] WORKING_THEORY.md (also CAMPAIGN_LEDGER.md:1333)

**Location.** WT:218, s10 row: 'the flying-qubit writer copies a source bit at ≥ 0.999680 with the carrier paying the register's energy change'

**Problem.** The energy flow is stated backwards. The register core loses energy and the carrier gains it. The s10 ledger's own adversarial re-read corrected exactly this wording ('The energy flows the other way'), but the correction never reached the CL entry or the WT row. The reviewer's sentence 'taken from the carrier's kinetic energy' (RF/grut_memory_candidate_20261019_review/REVIEW.md:42) contradicts the reviewer's own 'carrier +1.99'.

**Deciding check.** RF/grut_memory_candidate_20261019/REPORT.md:116-118 records carrier kinetic energy 10002.0408163 → 10004.0290991 and core 1.55 → −0.437569946. checks.py gives Δcarrier = +1.9882828 and Δcore = −1.987569946; their sum, +0.000712854, matches the recorded drift 0.000712816. See also ledger/ledger_s10_m1_memory.md:18, correction 1.

**Fix proposed.** Replace the wording with: 'the register's energy drop (1.55 → −0.44) is taken up by the carrier (kinetic energy +1.99)'. Rebuild CL:1333 from the corrected ledger file.

## 2. [major] CAMPAIGN_LEDGER.md; WORKING_THEORY.md; OPEN_CLAIMS.md

**Location.** CL:7-9 header and CL:1424 ('single-reader entry … full adversarial re-read not completed'); WT:202-205 (§5 intro); WT:375-377 (§8); OC:116; CL s10 entry CL:1303-1424

**Problem.** ledger/ledger_s10_m1_memory.md:5 says it is 'the corrected entry, written after an adversarial re-read … eight substantive corrections'; this landed in commit de88553. Four places still describe s10 as single-reader, and the CL s10 entry is the pre-correction text. Three corrections never propagated: (1) energy direction (finding above); (2) 'better at every sampled time overstates', yet CL:1331 and WT:290/335 still say it; (5) exact domain-wall conservation at ε=0 is a conditional model result, yet CL:1335 still calls it recovered known physics. A fourth was dropped: 'Radiative stability is not established' was deleted at 4e9a161 and is missing from WT's 'Not established' list. The packet says the s10–s12 rows were 'checked line by line … in the main session', but no record of that check exists, and it missed the s10 energy-direction error and the two s12 errors (+16 price, EA-0 overgeneralization).

**Deciding check.** `git diff --stat 7d5c0bc de88553 -- …/ledger/` shows ledger_s10 rewritten (404 lines changed). sed -n 14,32p ledger/ledger_s10_m1_memory.md. `grep -n 'carrier pay' CAMPAIGN_LEDGER.md` hits line 1333. ls results/ contains no row-check record. Compare `git show 1997d20:RF/grut_memory_candidate_20261019/OPEN_CLAIMS.md` (line 26, 'Radiative stability is not established') with the current OPEN_CLAIMS.md (absent).

**Fix proposed.** Regenerate the CL s10 entry from the ledger file. Restate verification status as 's01–s10 adversarially re-read (s10 in de88553, 8 corrections); s11–s12 single-reader'. Propagate corrections 1, 2, 5 and 6 to WT:218, 290 and 335. Attach the row-check record or drop 'line by line'.

## 3. [major] WORKING_THEORY.md; OPEN_CLAIMS.md

**Location.** WT:216 (s08 row) 'statistical range ±0.11, not a certified enclosure'; WT:332 (§7) same; OC:76-77 'the ±0.11 allowance on K₁ is statistical, not a certified enclosure'

**Problem.** This contradicts the source, which says the allowance is NOT a statistical confidence interval. The error came in with the round-1 fix: numbers-lens F10 misread the ledger. The ledger also requires noting that the allowance's angular term is root-derived and was not audited by the challenger; the packet omits this.

**Deciding check.** RF/grut_candidate_1_energy_20261019/REPORT.md:20-24: 'combining analytical omitted-tail bounds … with an explicitly empirical numerical allowance. It is **not** a certified enclosure of the exact continuum saddle, a statistical confidence interval, or a bound on genuine higher-loop corrections.' Also SPECTRAL_NOTES.md:192, and ledger_s08 §0 item 3 (angular term 0.05042782 root-derived; challenger bound 1.98366). checks.py prints both lines.

**Fix proposed.** In all three places write: 'numerical allowance ±0.11 (analytic tail bounds plus an empirical allowance; not a certified enclosure, not a statistical interval, not a higher-loop bound; its angular term is root-derived and was not challenger-audited)'.

## 4. [major] WORKING_THEORY.md; OPEN_CLAIMS.md

**Location.** WT:45-49 (§0: 'Every addition that a campaign proposed and tested failed … C1-as-FLT'); WT:226-232 (§6: 'Every item is a recorded withdrawal, refutation or non-implication in the repository … C1 defined as the FLT Dirichlet prescription terminates on an inverse-Gaussian coefficient'); OC:67 (Retained failures)

**Problem.** An unexecuted expectation is reported as a tested result. The inverse-Gaussian coefficient is quoted from the FLT literature, not reproduced. In the repository it exists only as source expectation (a) of the unexecuted B control. The ledger forbids treating those expectations as results.

**Deciding check.** RF/grut_closure_fork_20261019/CLOSURE_SELECTION_RESPONSE.md:120 ('The quoted inverse-Gaussian result belongs to that specified prescription') and :157-162 ('If C1 is **defined** to be (a), its inverse-Gaussian failure is a conditional termination'). RF/grut_closure_control_20261019/protocol.json:21-22 lists it under source_expectations_before_execution. STATUS.md:6-7: 'without treating its source expectations as numerical results'. Ledger: CL:1193, 1212, 1249.

**Fix proposed.** Rewrite as: 'C1 defined as the FLT prescription would terminate on the inverse-Gaussian coefficient the FLT literature reports (quoted, not reproduced; this is expectation (a) of the unexecuted B control)'. Remove it from §0's 'proposed and tested' list and move it from OC Retained failures to Unexecuted.

## 5. [major] WORKING_THEORY.md

**Location.** WT:57-63 (§0): 'The six entries of Θ, in every sector: … The program's own register prices this at +16 net supplied inputs with an empty derived tier and an empty prediction ledger (s12)'

**Problem.** A historical count is presented as the current price of this packet's Θ. +16 is the 2026-08-17 sync stamp for the 53-node responsive-vacuum stance (GRUT-as-written). That stance's Bin-1 primitive is the responsive-medium ontology, which this packet demotes to a calculational tool. The count covers none of Candidate 1, M1, CP-12, WM/B or the closure fork, so 'prices this … in every sector' quietly merges separate modules into one priced theory. GRUT III re-accounts the same register at +5. The ledger lists this exact promotion as must-not-promote, and WT:357-358 itself says ledger counts are not current.

**Deciding check.** POSTULATE_MAP.md:5-7 ('REGISTER-SYNC: 53 nodes, net +16'; 2026-08-17). GRUT_FINAL_SYNTHESIS_01.md:232-233 ('prices the full stance at +16'). GRUT_III/04_LEDGER_MIGRATION.md:27-33 (+16 → +5, 'a reorganisation'). CL:1613 ('Register counts (74 claims; 53/21; net +16) into a current universal: they are 2026-08-17 sync stamps of one lineage's register') and CL:594.

**Fix proposed.** Write: 'The 2026-08-17 register of the responsive-vacuum stance recorded net +16 underived inputs, an empty derived tier and an empty prediction ledger (GRUT III re-accounts it as five postulates). No register prices the other model instances.' Move this out of the Θ sentence.

## 6. [major] WORKING_THEORY.md; OPEN_CLAIMS.md

**Location.** WT:212 (s04 row): 'TT-only projector CHOSEN (Ward identity leaves {P⁽²⁾, P⁽⁰ˢ⁾}, GR's lock −2)'; OC:79-80 GRUT I–II failures

**Problem.** The row restates a classification that was retired as unconditional on 2026-08-14. The 'exactly two survivors' result needs the advanced-branch identity, 'a symmetry no dissipative completion sustains'. Under the diagonal identity the admissible space is strictly larger, and noise-kernel transversality became a new +1 input. This was correction 1 of the s04 adversarial re-read and is a ledger must-not-promote item. OC also omits the retirement.

**Deciding check.** provenance/claims.json:319 (WARD-SCOPE CORRECTION 2026-08-14, points 1, 2 and 5). ledger/ledger_s04_grut_I_II.md §0 item 1. CL:609: "'Exactly two Ward survivors' as an unconditional symmetry classification (RETIRED 2026-08-14)".

**Fix proposed.** Write: 'TT-only CHOSEN (GR's own kernel has P0s/P2 = −2; the two-Ward-survivor classification holds only under the advanced-branch identity and was retired as unconditional on 2026-08-14, enlarging the admissible space)'. Add the retirement and the +1 noise-transversality input to OC.

## 7. [major] WORKING_THEORY.md

**Location.** WT:210 (s02 row): 'the access seed is "a physical boundary condition, not an axiom, not a derived object"' and the answer 'Yes; no; no' to 'Are kernel and noise downstream of subsystem, state and access?'; WT:266-268 (§6 Access, same quote)

**Problem.** Two corrections from the s02 adversarial re-read were not propagated. (a) The quote is the P-5 instrument's phrase, reordered, not owner wording. The owner's binding wording is 'the access boundary is currently an unresolved structural input' and 'access seed = irreducible input, not a new law'; 'physical boundary condition' reads as a physical claim. (b) The 'Yes' drops the qualifier that the retarded kernel K was state-independent in both classes tested; the 'yes' holds for the influence functional and the noise kernel. Both are ledger must-not-promote items.

**Deciding check.** P5_ACCESS_VERDICT_01.md:67-68 (instrument) against P6_SEED_SELECTION_CHARTER_01.md:3-9 (owner). ledger/ledger_s02_influence_hierarchy.md:16 (correction 4), :22 (correction 10), :299-300. CL:250-251. DISCRIMINATOR_D1_RESULT.json:22,87 (Δ_K = 0.0). KERNEL_TRANSPORT_RESULT.json:32 (1.1e-13).

**Fix proposed.** Use the owner wording in both places. Answer (a) as: 'yes for the influence functional and the noise kernel; the retarded kernel K is state-independent in the Gaussian/free classes tested'.

## 8. [major] OPEN_CLAIMS.md; WORKING_THEORY.md

**Location.** OC:29-31 (Established 6): 'in every sector the recorded answer to the parent question is "nothing additional was forced in the stated class", with the sector's own computations as the demonstration'; WT:15-17 'asks this question in one place after another'; WT:86-88 'careful constructions of exactly such countermodels'

**Problem.** The closed-whole question is retrofitted onto sectors whose own ledger entries say they never posed or computed it. Listing this under 'Established' overstates what the ledger supports.

**Deciding check.** s01: 'no executed closed-whole computation … none of these took "bath = rest of the universe" as premise' (CL:68). s04: 'no primary source poses "bath = rest of the universe"; the bath is a declared KMS equilibrium bath' (CL:523, 613, 'Question retrofit'). s11: 'did not operate in the closed-whole class at all … cannot speak to what a closed whole with internal baths forces' (CL:1469). s03: 'none is formulated on a quantum state space' (CL:391). s12's 'computation' is the program's own ledgers (CL:1577).

**Fix proposed.** Rewrite as: 'In the sectors that computed in a closed or algebraic class (s02, s05–s10), nothing additional was forced in the stated class. s01, s03, s04, s11 and s12 did not pose the closed-whole question (open or classical classes, declared baths, or program ledgers); their verdicts are re-readings, not closure demonstrations.' Soften WT:15-17 and WT:86-88 to match.

## 9. [major] WORKING_THEORY.md

**Location.** WT:271-274 (§6 Access): 'EA-0 terminated UNFORMULABLE at the earned Level-0 scope: every seedless access candidate is global, generically trivial, a correlation readout, or generator-carried by identity'

**Problem.** (a) The stated reason is wrong. EA-0 is UNFORMULABLE because the earned substrate is classical, so the construction is undefined there; the candidate pattern comes from finite-dimensional quantum auxiliary classes, which 'do not upgrade the Level-0 terminal'. (b) 'Every' overgeneralizes. The pattern is 'a finding over the examined set, not a no-go theorem' covering ten candidates, and 'generically trivial' is verified only at d = 4 and 8. The ledger forbids promoting it to a general no-go. CL:1577 makes the same over-reach.

**Deciding check.** EA0_OWNER_RULING_02.md:10-31. EA0_ENDOGENOUS_ACCESS_FORMULATION_01.md:140-144 and 365-370. EA0_INDEPENDENT_VERIFICATION_01.md:140-163. CL:1617 (must-not-promote).

**Fix proposed.** Write: 'EA-0 terminated UNFORMULABLE at the earned Level-0 scope (the earned substrate is classical, and the construction is undefined there without an unearned lift). Separately, across the ten quantum candidates examined, each seedless candidate was global, generically trivial (proved at d = 4, 8), a correlation readout, or generator-carried by identity. This is a finding over the examined set, not a no-go theorem.'

## 10. [major] OPEN_CLAIMS.md; WORKING_THEORY.md

**Location.** OC:62-90 'Retained failures (not softened)'; WT §6 (WT:223-299)

**Problem.** Many failures that the ledger retains are missing, and the list does not say it is only a summary. s02: P-3 H3 gate failed; P-6 L-C gate red; C1-a two permanent reds; EQ-1 static gate failed and charter claim withdrawn; S2-1 run VOID. s04: x<0.358 DEMOTED; TT-auto freeze voided; κ prereg self-unblinds; mu_linear boundary 'weaker than stated'; WARD-SCOPE retirement. s06: 'FORMATION NOT ESTABLISHED'; 'Milestone A is not achieved'; G1 eight red gates. s07: A1 refutation and B1 two-event control failure (absent from OC). s09: coarse-ε criterion failure, after which only a changed criterion passes. s10: 'radiative stability not established' dropped. s12: GR-1 gate red; three negatives MUST-BE-REOPENED. The last bears directly on WT's 'did not survive' list, because the record says 'ruled-out status there is only as strong as its flag'.

**Deciding check.** CL:226-236, 576-589, 863-866, 966-967, 1234, 1603, 1605. GRUT_FINAL_SYNTHESIS_01.md:180-182. RF/cp12_unmasked_formation/REPORT.md:3. RF/cp12_source_resolved_propagation_20261015/REPORT.md:4. GRUT_II_What_Survived.md:53. `git show 1997d20:RF/grut_memory_candidate_20261019/OPEN_CLAIMS.md` line 26.

**Fix proposed.** Add one line per sector for the items above. Add a header sentence: 'Summary; the complete list is each ledger entry's Failures retained section.' Flag the three MUST-BE-REOPENED negatives wherever WT lists program negatives.

## 11. [major] WORKING_THEORY.md; OPEN_CLAIMS.md

**Location.** WT:301-312 'Still open (registered, not executed)'; OC:92-102 'Unexecuted tasks' (only 'GRUT III O5/O7 items'); WT:317-319 'Every kind of addition … was tested'

**Problem.** GRUT III pre-registered two tests that never ran, and neither appears in any unexecuted list. N1 asks whether P1 selects the massless field's horizon state without a free IR scale; the source calls either outcome 'decisive for GRUT'. N2 asks whether a horizon principle fixes Λ or another renormalisation constant. The CL s05 verdict cites both as the reason 'nothing forced' is not final. WT mentions only that 'P1 selects the horizon state remains a conjecture', so a reader of §7's 'No … state-selection rule has survived' is not told that the program's own decisive state-selection test is registered and unexecuted.

**Deciding check.** GRUT_III/13_SYNTHESIS.md:357-375. CL:690 (closure verdict) and CL:692 ('Unexecuted tasks: … N1–N4').

**Fix proposed.** Add to WT:301 and OC:98: 'GRUT III N1 (horizon-state selection by P1, pre-registered as decisive) and N2 (a horizon principle fixing a renormalisation constant), both unexecuted (s05).'

## 12. [minor] WORKING_THEORY.md

**Location.** WT:290-291 (§6) 'the matched uniform-floor control retains better at every sampled time for L = 3, 4, 5'; WT:335-336 (§7) 'retention probabilities that the matched uniform-floor control exceeds at every sampled time'

**Problem.** This is stronger than the source and contradicts s10 ledger correction 2. The source says only 'no worse … within 1e-10 … not an all-continuous-time theorem'. At t=0 the control is not better. Both places also omit that the control was specified after the review, not blinded.

**Deciding check.** RF/grut_memory_candidate_20261019/REPORT.md:28-29. checks.py on results_with_control.json: uniform−candidate at t=0 is −2.8e-15 (L=3) and −1.8e-15 (L=5); for t>0 the control is strictly better by at least 2.6e-6. control_protocol.json:6.

**Fix proposed.** Write: 'no worse at all 201 sampled times (equal at t=0, strictly better at every t>0); a post-review control, not a blind prediction'.

## 13. [minor] WORKING_THEORY.md

**Location.** WT:211 (s03 row): 'O-6 … was ruled FALSIFIED for a strict arrow from identities, the D-4 recurrence theorem and a non-member chain'

**Problem.** This misattribution came in with the round-1 fix. The O-6 ruling rests on identities I-1…I-5 and the D-1 band-edge asymptotics, not on recurrence. The review says explicitly 'it is not recurrence', and the ripples persist on the infinite chain. The row also omits that S6-1 later ran on the O-6 parent and found 'Pointwise arrow: NO. Integrated/net direction: YES', preparation-relative (ledger s03 §0 item 1).

**Deciding check.** L0_1G_PREFREEZE_REVIEW_01.md:30-33 ('This holds on the infinite chain as well. It is non-Markovian memory, and it is not recurrence.'). L0_1G_OWNER_RULING_02.md:30-37. ledger_s03_level0_floor.md §0 item 1.

**Fix proposed.** Write: 'ruled FALSIFIED from identities and band-edge asymptotics checked on a non-member chain (not from recurrence); the chartered strict-gate run was never executed; S6-1 later found a preparation-relative net arrow on the parent'.

## 14. [minor] WORKING_THEORY.md

**Location.** WT:213 (s05 row): 'in FRW two conventional temperatures of different observables exist (a detector temperature in the Wigner convention and a near-horizon Kodama–Hayward one)'

**Problem.** The Kodama–Hayward temperature is proved as a near-horizon scaling-limit KMS temperature (KPV 2021), but its cosmological case is open. Its tension with P1's definition, which assigns no temperature (WT:286), is unreconciled. 'Exist in FRW' is stronger than the ledger.

**Deciding check.** CL:672 ('cosmological extension open — the P1-definition/KPV tension is left unreconciled'). CL:702. GRUT_III/13_SYNTHESIS.md:163-167. GRUT_III/09_TWO_TEMPERATURES.md:21-33.

**Fix proposed.** Write: '… and a near-horizon Kodama–Hayward KMS temperature proved for trapping horizons (cosmological case open; its relation to P1's definition unreconciled)'.

## 15. [minor] WORKING_THEORY.md

**Location.** WT:253-255 (§6 Dynamical): 'single-pole finite memory stays a conjecture (s01, s04)'

**Problem.** This softens a retirement and a reversal. GRUT III RETIRED both rung1_ontology_finite_memory and rung3_single_pole. The in-house support suffered a 'structural reversal': exact dS forces 'infinite scale-free memory — the opposite shape of the original bet'.

**Deciding check.** GRUT_III/04_LEDGER_MIGRATION.md:13,15. GRUT_MODEL_FRAMEWORK.md:112-117. ledger_s01_v7_ctp.md:78,96-97.

**Fix proposed.** Write: 'single-pole finite memory: in-house support reversed (exact dS forces infinite scale-free memory), retired in GRUT III, UNRESOLVED in the final synthesis'.

## 16. [minor] WORKING_THEORY.md

**Location.** WT:286-287 (§6 Arrow and clock): 'The temporal generator is supplied (S5-0 GENERATOR-IRREDUCIBLE) … (s02, s05)'

**Problem.** The sector attribution is wrong. Neither the s02 nor the s05 ledger cites S5-0. It appears in s01, s03 and s12.

**Deciding check.** `grep -l 'S5-0\|GENERATOR-IRREDUCIBLE' ledger/*.md` returns ledger_s01, ledger_s03 and ledger_s12.

**Fix proposed.** Change the attribution to '(s01, s03, s12; P1 and modular time: s05)'.

## 17. [minor] WORKING_THEORY.md

**Location.** WT:50-51 and WT:241 ('the CP-12 min-entropy theorem as a selection law'); WT:58-59 and WT:239 (TT projector filed under preparation/state selection); WT:54 ('autonomous every-time retention') and WT:294 ('refuted at the permitted point')

**Problem.** Three items in the 'did not survive' lists are misstated. (a) No campaign proposed the min-entropy theorem as a selection law; the source says three times that it is not one, and it survives as a conditional obligation whose parent theorem is unaudited. (b) The TT projector is a kernel/channel choice of the dynamics (p_tt_ansatz under the P3 coupling), not a preparation. (c) A1 refutes every-time retention at one permitted point only, and other points remain free.

**Deciding check.** CL:844. RF/cp12_preparation_minentropy_audit/REPORT.md:104,111,159. GRUT_II_What_Survived.md:35-43. GRUT_III/04_LEDGER_MIGRATION.md (p_tt_ansatz row). CL:958 ('A1 refutes uniformity at one point only').

**Fix proposed.** (a) Drop it from the did-not-survive list, or say 'shown to be an obligation, not a selection law'. (b) Move the projector to the dynamics/action slot. (c) Write 'at one permitted parameter point' in both places.

## 18. [minor] WORKING_THEORY.md

**Location.** WT:62-63 ('"currently irreducible" after the EA-0 attempt') and WT:220 (s12 row hangs the same quote on EA-0); WT:17-20 ('The old parent object, an in-in action S_CTP treated as a "genotype", is the Gaussian approximation to a derived object')

**Problem.** (a) The quote predates EA-0. It comes from GRUT_PROGRAM_STATE_SYNTHESIS_01.md, dated 2026-09-28; EA-0 was ruled on 2026-09-29. (b) §0 still states V7's object as fact. The round-1 V7-attribution fix reached §2 (WT:124-127) but not §0, and per ledger_s01 §1 the S_CTP=S_IF identification is the handoff's and cannot be checked.

**Deciding check.** checks.py prints the date lines: GRUT_PROGRAM_STATE_SYNTHESIS_01.md 2026-09-28 (quote at :73); EA0_OWNER_RULING_02.md 2026-09-29. ledger_s01_v7_ctp.md §1 and CL:42.

**Fix proposed.** (a) Write: 'the record called it "currently irreducible" the day before EA-0, and EA-0 left it so'. (b) Write: 'the parent object the handoff attributes to V7 (file not in the repository)'.

## 19. [minor] WORKING_THEORY.md

**Location.** WT:209 (s01 row): '22–62 (GW dissipation) and 62.7 (Γ_T) orders below their respective bounds'; WT:251-252 (§6)

**Problem.** Each figure covers only part of its channel. The 22–62 is a dephasing statement; the amplitude channel (Γ = B·H₀/2 = 0.2 H₀ at B = 0.4, inside the few×H₀ slot) is 'NOT SETTLED EITHER WAY'. The 62.7 is horn (a) only; horn (b) is 'NOT banked'. Labelling the first 'GW dissipation' is exactly what s04 forbids.

**Deciding check.** CL:542, 550, 584, 599. calc/gw_dissipation_bounds.py:93-110.

**Fix proposed.** Write: 'GW dephasing 22–62 orders below detectability (amplitude channel unsettled); Γ_T horn (a) 62.7 orders below (horn (b) not banked)'.

## 20. [minor] WORKING_THEORY.md

**Location.** WT:309-312 (§6): 'recorded against a GRUT I–II shape exclusion … the exclusion itself is conditional on the supplied projector and arrow'

**Problem.** The status and conditionality are misstated. The no-crossing of w = −1 is not a GRUT I–II boundary. The deposit lists it as 'Held OPEN, deliberately not boundaries … rung3-gated … never ⊘', conditional on the open pole structure (rung3).

**Deciding check.** GRUT_I_What_Survived.md:26 (CORRECTED 2026-08-10c) and :34. SIGNATURE_AUDIT.md:56.

**Fix proposed.** Write: 'against the held-open, rung3-gated no-crossing of a single passive relaxor (not a banked boundary)'.

## 21. [minor] WORKING_THEORY.md

**Location.** WT:219 (s11 row): 'decides the self-correction/no-scrambling bridge claim against the bridge in that classical class'

**Problem.** The row compresses the scope. The rejection holds only in the Z-sector Metropolis 4D code at β∈{1.0,1.1}, L∈{3,4,5}, with a local cooling decoder and t_max 1e5. By the letter of the pre-registration, R2 FAILS on the face spin, and the fallback rests on O3, which was pre-registered as secondary. O2 and O3 ran in separate seeded runs, not 'the same trajectories'.

**Deciding check.** CL:1471, 1482, 1489, 1495, 1502.

**Fix proposed.** Add the scope and: 'by the letter of the pre-registration R2 fails; the fallback rests on secondary O3'.

## 22. [minor] OPEN_CLAIMS.md

**Location.** OC:51-52: 'no common action, state space, constraints or calibration exists'

**Problem.** This asserts non-existence, where WT:325-326 says correctly that none 'has been established'. The M1 packet itself carries an open conjecture relating M1 to Candidate 1 as a low-energy sector.

**Deciding check.** CL:1335 ('relation to Candidate 1 as a low-energy sector', conjecture). RF/grut_memory_candidate_20261019/COMPATIBILITY.md:87-92.

**Fix proposed.** Replace 'exists' with 'has been established'.

## 23. [minor] OPEN_CLAIMS.md

**Location.** OC:37-43 'Conditional (holds in a stated class only)'; OC:4 'Anything not listed under "Established" is not claimed'

**Problem.** Items that are unaudited, author-derived or red-gated are listed as holding. The CP-12 obligation rests on a parent theorem with no partner audit. The B1 window is author-derived, review pending and partly S9-conditional. The cumulant ladder's order-8 rung failed its frozen gate. Separately, OC:4 contradicts the existence of a Conditional section of claims.

**Deciding check.** CL:859, CL:973, CL:984, CL:226.

**Fix proposed.** Write '(stated in a class; author-derived or unaudited where marked)' and mark each item. Reword OC:4 to cover the Conditional section.

## 24. [minor] OPEN_CLAIMS.md

**Location.** OC:79-80: 'economical w(z) … the α→TT route … all retracted'

**Problem.** Settled-negative boundaries are labelled as retractions. The economical w(z) and the α→TT bridge 'failed as derivations (two survive as boundaries)', and α→TT is 'not a formal no-go'. The 689 Hz falsifier belongs to the June deposit (s01), not GRUT I–II.

**Deciding check.** GRUT_II_What_Survived.md:75. GRUT_I_What_Survived.md:30-31. GRUT_PROGRAM_CLOSURE_01.md:35-37.

**Fix proposed.** Split the list into 'retracted: …' and 'failed as derivation, retained as settled-negative boundaries: economical w(z), α→TT'.

## 25. [minor] WORKING_THEORY.md

**Location.** WT:317-319 (§7): 'Every kind of addition in FORMAL_CORE §7 except a modified global law was tested by at least one campaign'; WT:353-354: 'none has produced a condition that is not already restriction, trace matching, positivity, normalization or modular KMS'

**Problem.** 'Tested' overstates two slots. For 'measure over Θ', the CP-12 obligation is a necessary condition, and the Past Hypothesis is a supplied Bin-1 primitive that 'no calculation would discharge'. For the record slot, the M1 writer and B2 receiver are supplied apparatus, not proposed laws. The last sentence extends the C2 state-selection finding to every slot, although the ledger records supplied conditions of other kinds: P2 Clausius, the TT projector, and sequestering's global equation.

**Deciding check.** FORMAL_CORE.md §7 rows 'measure over Θ' and 'instrument / record slot'. POSTULATE_MAP.md:17-21. CLOSURE_SELECTION_RESPONSE.md:7-11. CL:657, CL:1204.

**Fix proposed.** Write: '… tested or examined, except a modified global law and a measure over Θ (supplied as P4, never tested)'. Restrict the last sentence to state-selection attempts.

## 26. [minor] WORKING_THEORY.md

**Location.** WT:328-339 (§7 'Predictions, conditional on inputs')

**Problem.** The paragraph omits that none of the model-instance predictions distinguishes its instance from known physics. GRUT III makes 'no prediction that distinguishes it from GR + QFT'. Candidate 1's lensing and Born outputs 'reproduce adopted laws', and its loop number is a consistency number. M1's retention is beaten by the control.

**Deciding check.** CL:682, CL:1125, CL:1088. RF/grut_closure_fork_20261019/REPORT.md:5-8.

**Fix proposed.** Add: 'none distinguishes its instance from the known-physics baseline'.

## 27. [minor] WORKING_THEORY.md

**Location.** WT:210 (s02 row): 'discrimination persists at orders 6 and 8'

**Problem.** The row omits that the order-8 frozen gate FAILED (ratio 67.8) and 'stays red'. Discrimination at order 8 rests on a labelled post-hoc diagnostic.

**Deciding check.** P4_HIERARCHY_VERDICT_01.md:36-41. CL:226.

**Fix proposed.** Write: '… at order 6, and at order 8 only by a labelled post-hoc diagnostic (frozen gate red)'.

## 28. [minor] WORKING_THEORY.md

**Location.** WT:213 (s05 row): 'η = 1/4Għ'

**Problem.** The formula reads as Għ/4. The source has η = 1/(4Għ). Round 1 noted this but it was not fixed.

**Deciding check.** CL:657. GRUT_III/01_PRINCIPLES.md:45-52.

**Fix proposed.** Write 'η = 1/(4Għ)'.

## 29. [minor] WORKING_THEORY.md

**Location.** WT:215 (s07 row): 'Closure raised the proof and resource cost; it supplied no principle'; B2-RA alone qualified as author-derived

**Problem.** The ledger marks this grouping as 'this ledger's synthesis, not a source statement', yet it appears as the sector's result. The row qualifies only B2-RA as author-derived, but independent mathematical review is absent for every result in the sector.

**Deciding check.** CL:951, CL:973. RF/grut_b2_feasibility_closeout_20261018/REVIEW_STATUS.md:3-7.

**Fix proposed.** Add '(ledger synthesis)' and 'all results author-derived; independent review pending'.

## 30. [minor] WORKING_THEORY.md

**Location.** WT:41-42 (§0): 'not the constants, not the action. Each negative is exhibited by a legal closed model in calculation 1'

**Problem.** The round-1 BY-DEFINITION split never reached §0. FORMAL_CORE §4 and OC:15-19 file constants and action as BY DEFINITION (elements of Θ), not as countermodels from calc 1. Uncertain: calc 1 S5 does exhibit legal closed models with different couplings.

**Deciding check.** checks.py FORMAL_CORE §4 table: 'the observed constants or a specified matter/gravity action' is BY DEFINITION. OC:15-19.

**Fix proposed.** Write: 'Each dynamical negative is exhibited by a legal closed model in calculation 1; constants and action are not forced because they are entries of Θ by definition.'

## Confirmed by this lens

- M1 numbers are the corrected 4e9a161 values: WT:218 L=4 min 0.987510 through t=10, 0.935065 at t=40, readout 0.778468; control 0.998458 at t=40. All match RF/grut_memory_candidate_20261019_review/control_results.txt and the handoff table.
- The facilitation-benefit claim is stated as withdrawn at WT:49, 218, 290 and 361, and OC:64.
- Candidate 1 is called conditional. WT:321-323 calls it a conventional Einstein+SM+scalar EFT among 'conditional model instances'. The s08 row classification (WT:216) is 'conditional one-loop constraint; supplied couplings'. WT:331-334 says 'for its supplied couplings … conditional on reduced-saddle stability'. OC:40. Tree and loop orders are kept separate: −3553.71 gv is labelled tree-plus-one-loop, and the tree gaps are attributed to the closure-fork packet.
- C2 is said to supply no new state-selection rule: WT:217 ('observer restriction selects nothing') and WT:233-235 ('adds no selection condition beyond restriction, trace matching, positivity, normalization and modular KMS').
- The B three-prescription control is said to be unexecuted at WT:217, 301 and 360-361, and OC:68 and 94.
- No combined cross-module law is claimed: WT:321-326 says Candidate 1, M1 and GRUT III live in different state spaces and that no common action, state space, constraints or calibration has been established.
- Nothing announces validated universal physics or a derivation of gravity or all reality: README:30-32, WT:65-68, OC:53-58.
- Round-1 fixes that landed cleanly: differentiator ranges split (WT:209); toric qualifiers in §0 and §6 (WT:52-54, 296-299); Araki–Lieb slack values; OC #4 counts (8 NOT FORCED / 3 BY DEFINITION / 3 SCOPE) match the FORMAL_CORE §4 table; verify.py has 28 headline checks plus 2 provenance hashes (30 check calls); V7 attribution in §2; B2 stop stated as a method and resource decision; CP-12 vacuity through N=10^18 and the missing partner audit; s01–s09 re-read correction counts 9/12/8/11/15/8/12/8/11 match '8 to 15'.
- Numbers traced or re-derived: Ω_K(N=5) = −4.5402e-5 and Ω_K(N=10) = −2.0612e-9; 977 → 1.64e4 → 7.98e4 sweeps; 0.338; 1.93%; 1411/640; N ≥ 2^28; K1 2.79528 → 2.795; −3553.70573 → −3553.71; tree gaps −0.199999152115 and −3067.397336565; ε0 = 0.4725; 74 = 53+21 and net +16 as printed at POSTULATE_MAP.md:6-7.
- python3 -I code/verify.py (hash mode) reports 0 failures, and the working tree is clean.
