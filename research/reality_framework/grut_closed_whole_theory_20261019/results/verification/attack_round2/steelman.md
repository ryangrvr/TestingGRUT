# Attack round 2 — steelman lens

Returned by a separate agent as structured findings (the harness blocked report files). Rendered verbatim into Markdown for the record. Scratch paths refer to the session scratchpad and are not part of the repository.

Findings: 20 (11 major, 9 minor).

## 1. [major] WORKING_THEORY.md; CALCULATIONS.md; OPEN_CLAIMS.md; FORMAL_CORE.md

**Location.** WT section 0:38-39 'Closure forces nothing local pointwise'. CALC:240. OPEN_CLAIMS:11 'Closure alone forces only the global theorems' and :19. WT section 0:50 and section 6:241: 'the CP-12 min-entropy theorem as a selection law' is listed under proposed additions that did not survive.

**Problem.** The headline is false as stated, and it conflicts with FC section 3:243, which lists monogamy and the entropy inequalities (pointwise local constraints) as forced. Closure alone also forces an athermality budget. For any stationary sigma, D(rho_Omega(t)||sigma) is conserved (THM 2.1). By monotonicity of relative entropy under partial trace, D(rho_A(t)||Tr_Abar sigma) <= D(rho_Omega(0)||sigma) for every subsystem A and every time t. This statement is local, pointwise and forced, and it is absent from FC sections 2-3. The CP-12 parent and min-entropy theorems are the program's own instance: ledger s06 says closure does work at two steps of the proof. Their source never proposed them as a selection law (ledger s06:68: 'not a new dynamical state-selection principle'). Listing the theorem 'as a selection law' among failed additions therefore attacks a proposal nobody made. It also contradicts OPEN_CLAIMS, whose Conditional section counts the CP-12 obligation as surviving.

**Deciding check.** checks.py A2, on the packet's own S1 datum with beta=1. The global budget D(rho_Omega||gamma_H) = 0.8911 drifts by 1.3e-13 over [0,40]. The local D(rho_A(t)||Tr gamma_H) stays below the budget at every t, with maximum 0.8017 (90%, at t=0). forbidden.py: gamma_A = diag(0.279, 0.721), and D(|0><0| || gamma_A) = 1.276 > 0.891. So closure forbids A from ever reaching |0> from this preparation, and forces p(Z_A=+1) <= 0.911 at every t (the actual maximum is 0.500).

**Fix proposed.** Replace the headline with: 'no local dynamical law, kernel, arrow direction or pointwise relaxation is forced; pointwise local inequalities are (entropy inequalities, monogamy, athermality budget).' Add an athermality-budget theorem to FC section 2. Move the CP-12 item in WT sections 0 and 6 to the surviving conditional results ('closure load-bearing; vacuous at printed sizes'), with 'not a selection law' kept as a scope note.

## 2. [major] FORMAL_CORE.md; WORKING_THEORY.md; CALCULATIONS.md

**Location.** FC section 4 arrow row: 'NOT FORCED | calc 1 S7 ... the arrow is the low-entropy preparation'. WT section 0:40 'not an arrow'. CALC S7:166. FC section 7 instrument/record row.

**Problem.** The forced half of the arrow is missing. Take rho_A (x) gamma_Abar with gamma_Abar the Gibbs state of H_Abar, which is the packet's own S1 preparation. Closure then gives the identity Sigma(t) = dS_A + beta dE_Abar = I(A:Abar)(t) + D(rho_Abar(t)||gamma_Abar) >= 0 at every t (Esposito-Lindenberg-Van den Broeck, NJP 12, 013013 (2010)). The Landauer erasure cost (Reeb-Wolf, NJP 16, 103011 (2014)) follows from the same identity. This is closure plus one named hypothesis (a passive, uncorrelated complement at the reference time), exactly the COND format of 2.10-2.13. The inequality is two-sided in t and not monotone, so it is consistent with 'direction not forced'. It is still a forced pointwise local statement, and it is the record-erasure cost that the record slot of section 7 omits. The program's own ARROW_OF_TIME.md already cites Partovi and Pusz-Woronowicz for this.

**Deciding check.** checks.py A, on the S1 datum (401 points). Minimum Sigma = -5.8e-15 (zero at t=0); Sigma(1) = 0.346, Sigma(10) = 0.739, Sigma(40) = 0.665. The identity Sigma - I - D holds to 1.4e-14. Sigma decreases on 217 of 400 steps. On the S7 reversed preparation rho(40)*, which is correlated and has a non-Gibbs bath, Sigma_rev is negative on 66% of the grid, with minimum -0.665 at t=40. The hypothesis, not closure, carries the inequality.

**Fix proposed.** Add 'COND 2.14 (second law and Landauer for an initially passive, uncorrelated complement)', stating Sigma = I + D >= 0, its two-sidedness and the erasure bound. Rewrite the arrow row as: direction NOT FORCED (S7); pointwise Sigma >= 0 FORCED under COND 2.14, not monotone. Add the Landauer cost to the FC section 7 record row.

## 3. [major] FORMAL_CORE.md; WORKING_THEORY.md

**Location.** FC section 4 rows 'pointwise or monotone relaxation', 'a Markovian or exponential memory kernel', 'a thermal complement'. WT section 0:38-41. CALC S2, S3, S5.

**Problem.** Every countermodel is finite-dimensional. In finite dimension exact relaxation is already excluded by quasi-recurrence (THM 2.9), so these countermodels cannot separate 'closure does not force relaxation' from 'finite dimension forbids it'. On the natural reading of 'the rest of the universe' the complement is infinite, and standard theorems then force these properties under named hypotheses the packet never states: (i) return to equilibrium, for an infinite KMS reservoir with mixing and spectral conditions and small coupling (Jaksic-Pillet, CMP 178, 627 (1996); Bach-Frohlich-Sigal, JMP 41, 3985 (2000)); (ii) a Markov semigroup, from Davies' weak-coupling limit (CMP 39, 91 (1974)); (iii) local thermality, from closure plus ETH plus a narrow energy window (Deutsch 1991; Srednicki 1994). The thermal-complement row also does not point to COND 2.11. None of these is closure alone, so the verdicts stand, but the conditional routes are missing.

**Deciding check.** Precise argument: finite dimension implies almost periodicity (FC THM 2.9), so lim rho_A(t) exists only for stationary rho_A. The cited theorems require an infinite reservoir with continuous spectrum. grep -i 'Davies|return to equilibrium|ETH|eigenstate' over FC and WT finds only the toric Davies-type bath of s11.

**Fix proposed.** Add COND items: an infinite KMS complement with mixing and weak coupling gives return to equilibrium; van Hove scaling gives the Davies semigroup; ETH plus an energy window gives local thermality. Annotate the three section 4 rows: 'NOT FORCED by closure alone; FORCED under COND x; the finite countermodel class excludes relaxation by recurrence and does not test the infinite reading.'

## 4. [major] FORMAL_CORE.md; WORKING_THEORY.md; CAMPAIGN_LEDGER.md

**Location.** FC section 4 records row 'NOT FORCED | records are stable for all states iff [H, Z_R] = 0'. FC section 7 instrument/record row. WT section 6:290-293, where 'exact saturation is a barrier to formation as well as a memory (s10)' is filed under 'Proposed additions that did not survive'. WT section 5 s10 classification.

**Problem.** (a) Closure alone forces a protect/form trade-off. If [H, Z_R] = 0, the Z_R statistics are constant for every state, so the closed dynamics can neither write nor erase that record. A register the dynamics can write has [H, Z_R] != 0 and is not stable for all states. M1's 'Perfect protection is not a mechanism for generating the protected state' (grut_memory_candidate_20261019/WORKED_REGIME.md:41-44) is this theorem. The corrected s10 entry reclassifies it as a conditional constraint (ledger s10 section 0, item 5), so filing it as a failure inverts its status. (b) Under non-degenerate gaps, the time-averaged distinguishability of two record preparations on R is at most D(omega_R^0, omega_R^1) plus the two LPSW terms, and omega depends only on conserved, energy-diagonal data. So records not carried by conserved quantities are forced to fade on average. (c) Quantum Darwinism (Brandao-Piani-Horodecki, Nat. Commun. 6, 7908 (2015)): closure plus a fragmented rest forces redundantly held records to be classical information about a single measurement. None of these is stated.

**Deciding check.** checks.py C: ||[H', Z_R]|| = 0.0 exactly in S6, so H' cannot write R. checks.py D, packet H with a Neel bath: |+> vs |-> on A have identical dephased states (4.7e-16) because Z-parity is conserved, so the X record must fade. The forced average bound is 0.368; the sampled average D over 1500 t in [0, 1e4] is 0.152, down from 1.000. For |0> vs |1>, D(omega_A^0, omega_A^1) = 0.169 and the sampled average is 0.170: the part that survives on average is exactly the conserved part.

**Fix proposed.** Add a protect/form theorem, a 'records need conserved data' COND and a BPH COND, each with its hypotheses. Rewrite the records row as 'persistence NOT FORCED; necessary conditions FORCED'. Move M1's formation barrier from WT section 6 'did not survive' to the surviving conditional constraints, and add 'conditional constraint' to the s10 classification.

## 5. [major] FORMAL_CORE.md; CALCULATIONS.md; WORKING_THEORY.md

**Location.** FC section 3:241-242, which lists 'the record/instrument equivalence and the multi-time structure (2.8)' as forced from closure alone. FC DEF 1.5. CALC S6:137 'At t1 = 3 apply CNOT(A->R)' and its heading 'Records and interventions are included subsystems and unitaries'. Compare FC DEF 1.1:36 'no reservoir, drive or measurement device outside Omega'.

**Problem.** This is a hidden hypothesis in the forced list. The S6 write is applied by hand at a sharp external time that, by DEF 1.1 and 1.6, no subsystem reads. H' does not generate it (H' conserves Z_R), and it changes the conserved energy. Under strict closure an intervention must be generated by H through an internal clock and controller, which brings energy bookkeeping, WAY (COND 2.13) and clock imperfection. On the last point: for H bounded below no internal clock runs monotonically with certainty (Unruh-Wald, PRD 40, 2598 (1989)), and autonomous control with a finite clock degrades (Woods-Silva-Oppenheim, AHP 20, 125 (2019)). The program's own M1 packet flags the energy point: 'An energy-conserving detector completion must receive that energy' (WORKED_REGIME.md:144-148).

**Deciding check.** checks.py C: <H'> is -10.228789 before the write and -10.021233 after. The hand-applied CNOT injects dE = +0.2076 = -E_int(3), a jump in THM 2.1's conserved quantity. S6 is not a closed whole during the write.

**Fix proposed.** Label S6 as an externally timed intervention that opens the whole at t1 and injects 0.2076 of energy. In FC section 3, move 'multi-time structure under interventions' to a line 'given interventions (external, or internal with a clock and controller)'. Add internal-apparatus costs (energy, WAY, clock imperfection) to the section 7 instrument and clock rows, or model the writer inside Omega as M1-W does.

## 6. [major] WORKING_THEORY.md; CALCULATIONS.md; FORMAL_CORE.md

**Location.** WT section 0:42 'a state rule built from H alone is shown to be dynamically empty'. CALC S9:179 'freezes every subsystem' and :183 'no subsystem has any dynamics'. CALC reading 2:242-243. FC section 4:258 'rules F(H) that depend on H alone freeze all local dynamics'. FC COND 5.3:306. WT section 7:343.

**Problem.** The equivariance theorem makes only the one-time, un-intervened marginals constant in the external t. Multi-time correlations, responses to interventions, FDT physics and relational (Page-Wootters) dynamics all remain. The packet's own THM 2.8 calls multi-time statistics the full content, and GRUT III's KMS response physics lives entirely in stationary states. Separately, WT section 0, FC section 4 and CALC reading 2 drop the 'unitarily equivariant' hypothesis: a rule that references a fixed basis is a function of H and need not be stationary.

**Deciding check.** checks.py B on the packet's H. For the Gibbs state at beta=1, ||[F,H]|| = 8.9e-15, yet the connected <Z0(t)Z0> ranges over [0.339, 0.805] for t in [0,20]. A local pi/2 Y-kick on A changes rho_A(t) by up to 0.381, and the selective outcome X=+1 (p=0.5) changes it by up to 0.921. For the ground state the corresponding values are 0.887 and 0.953.

**Fix proposed.** Reword to: 'a unitarily equivariant F(H) is stationary: its one-time marginals do not evolve in the external parameter; multi-time correlations, intervention responses and clock-relative dynamics remain; one-time local evolution requires breaking equivariance.' Restore 'unitarily equivariant' in all three places that drop it.

## 7. [major] FORMAL_CORE.md; WORKING_THEORY.md

**Location.** FC sections 2-3: COND 2.13 is the only theorem about apparatus inside Omega. FC section 7 clock row (cost 'COND 5.2; must reproduce geometric cases') and instrument/record row. FC section 4 Born row. The calc 1 instruments are sharp projective measurements at a sharp external t.

**Problem.** The theorems specific to having no exterior are missing. (i) No complete self-measurement: an apparatus inside a system cannot determine the state of the system containing it (Breuer, Philos. Sci. 62, 197 (1995)); forced by closure alone in finite dimension. (ii) Internal clocks are imperfect (Unruh-Wald 1989; Pauli); forced by closure plus H bounded below. (iii) Wigner's-friend no-gos: universal unitarity applied to an observer inside the whole, plus a superobserver with reversible control, contradicts absoluteness of observed events together with locality and no-superdeterminism (Bong et al., Nat. Phys. 16, 1199 (2020); Brukner, Entropy 20, 350 (2018); Frauchiger-Renner, Nat. Commun. 9, 3711 (2018)); forced by closure plus named access and LF assumptions. (iv) Relational access under a symmetry when no external frame exists (Bartlett-Rudolph-Spekkens, RMP 79, 555 (2007)), the finite counterpart of COND 5.4. Each constrains a slot of Theta (instrument, record, clock, access). The packet's own idealized instruments are what these theorems remove. Uncertainty: the owner may treat the Wigner's-friend items as interpretation; the clock and self-measurement items are not interpretation.

**Deciding check.** grep -i 'Breuer|Unruh|Wald|friend|Frauchiger|Brukner|reference frame|self-measure' over FC, WT and OPEN_CLAIMS finds nothing. The forcing status of each item is given in the problem text.

**Fix proposed.** Add a short subsection to FC section 2, 'closure-specific limits on internal apparatus', containing COND 2.13 and items (i)-(iv), each with its hypothesis. Put their costs into the section 7 clock, instrument/record and access rows.

## 8. [major] FORMAL_CORE.md; CALCULATIONS.md; WORKING_THEORY.md

**Location.** FC section 4 row 1, 'a pure global state | BY DEFINITION ... excluded by this definition rather than refuted'. CALC S1 (mixed Gibbs complement) and S8 (Araki-Lieb on the mixed trajectory).

**Problem.** The round-1 relabel is honest, but the packet then says nothing about the excluded reading. If the whole has no exterior, a mixed rho_Omega is either a correlation with an outside (which contradicts closure) or ignorance over pure closed wholes. Every Born prediction is linear in rho_Omega, so the second reading changes no Born statistic. It does change every nonlinear statement. S(rho_Omega) = 1.2609, I(A:Abar), the Araki-Lieb slack and Sigma (F2) become ensemble averages, not properties of one closed whole. For each actual whole, closure then forces S_A = S_Abar for every cut at every time: all local entropy is entanglement with the rest, the one place where 'bath = rest' gains content. Page-type statements follow under COND 2.11. This is the one 'not forced' row whose verdict flips on the natural reading of the parent question.

**Deciding check.** Precise argument: Tr(rho E) is linear in rho, and THM 2.9 gives pure rho_Omega => S_A = S_Abar. The packet's S8 pure-state check (deviation 3e-16) shows the forced identity this reading gains.

**Fix proposed.** Give both readings in FC section 4 row 1: BY DEFINITION under DEF 1.1; FORCED on the ontic reading of 'no exterior', where a mixed rho_Omega is an ensemble of closed wholes, Born predictions are unchanged, entropic statements refer to the ensemble, and S_A = S_Abar holds for every cut. Label the CALC S1/S8 entropic numbers accordingly.

## 9. [major] WORKING_THEORY.md; FORMAL_CORE.md

**Location.** WT section 0:30 'From closure alone, exactly the theorems of standard quantum theory'. FC section 0, where 'known physics' means QM with supplied data. FC DEF 1.1 and DEF 1.8.

**Problem.** The thesis is true by construction. DEF 1.1 defines closure as unitary evolution on a finite Hilbert space, and DEF 1.8 supplies the Born rule. In the finite class, 'nothing beyond standard QM' could not have come out otherwise; only the section 7 'modified global law' row could. The steelman goes further: read operationally, 'every state is the marginal of a pure state of system plus rest, unique up to reversible maps on the rest' is the purification postulate. Classical theory satisfies the other five principles; purification is the one that singles out finite-dimensional QT among probabilistic theories (Chiribella-D'Ariano-Perinotti, PRA 84, 012311 (2011); PRA 81, 062348 (2010)). On that reading the parent premise forces known physics, not nothing.

**Deciding check.** FC DEF 1.1:33-37 and DEF 1.8:86-93 make unitarity and the Born rule inputs. grep -i 'Chiribella|purification' over FC and WT finds nothing. CDP 2011 abstract: the five principles hold for classical theory, and purification singles out QT.

**Fix proposed.** Add one sentence to WT section 0 and FC section 0: 'In the finite class this negative is definitional, since DEF 1.1 defines closure inside QT; operationally, closure as purification (CDP 2011) selects quantum theory among probabilistic theories, so closure forces known physics, not nothing.' Optionally add a section 7 row 'theory class (GPT level)'.

## 10. [major] WORKING_THEORY.md; FORMAL_CORE.md

**Location.** WT section 7:341-354 'What extra condition would make a distinct prediction'. FC section 7 table costs (for example 'must recover known physics; introduces its own constants'). FC COND 5.2 'automatic and selects nothing'.

**Problem.** The packet under-claims a GRUT-specific result. Handoff:230-231 asks it to 'explain exactly what extra condition would be needed to make a distinct prediction'. The program already has an exact answer. GRUT III T1 (GRUT_III/07_ACCESS_RESOLVED_RESPONSE.md:20): KMS at beta implies N_f = (1/2) coth(beta omega/2) rho_f for every observable. T2 (07:34-38): with the same state and observables, a theory can differ from the conventional one only by Delta rho_f != 0, by a different state (non-KMS, different beta, or non-stationary), or by higher cumulants. The conclusion is 'In the regime its own principle defines, GRUT III can only be standard' (GRUT_III/03_PREDICTIONS_AND_FALSIFIERS.md:8). There is also a priced ten-line checklist for any sixth postulate (GRUT_III/12_BASELINE_AND_HANDOFF.md:36-105), ending 'I do not have a candidate that passes all ten lines'. Ledger s05 records these as surviving. The 0.338 fourth-order split in s02 is option (c) of T2. Also missing is the conditional theorem 'KMS => FDT', which goes beyond THM 2.5's Gaussian case; COND 5.2's 'selects nothing' omits that KMS fixes the noise from the dissipation.

**Deciding check.** The cited GRUT_III lines. grep 'Delta rho|trichotomy|12_BASELINE|ten-line|sixth postulate' over WT and FC finds nothing.

**Fix proposed.** State T2 in WT and FC section 7 as the form a distinct prediction must take: a Delta rho, a different state or a higher cumulant, for a stated observable, state and preparation. Cite 12 section 2 as the program's existing checklist for the dynamics and state rows. Add COND 'KMS => FDT (T1)' to FC section 2 or 5.

## 11. [major] WORKING_THEORY.md; CAMPAIGN_LEDGER.md; OPEN_CLAIMS.md

**Location.** WT section 5 s10 row:218 'with the carrier paying the register's energy change'. CAMPAIGN_LEDGER.md:26 and :1333 'energy conservation makes the carrier pay ... (core 1.55 -> -0.437569946 as carrier 10002.0408163 -> 10004.0290991)'. CAMPAIGN_LEDGER:1424, WT section 5:203-205, section 8:375-377 and OPEN_CLAIMS:116 'single-reader'. WT section 6:290-291 'retains better at every sampled time'. WT section 7:335-336 'exceeds at every sampled time'.

**Problem.** ledger/ledger_s10_m1_memory.md was replaced at de88553 by an adversarially re-read entry with eight corrections, and none of them propagated. Correction 1 reads 'The energy flows the other way ... the write releases register energy into the carrier', and the quoted numbers themselves show the carrier gaining energy. Correction 2: the control is equal at t=0 and strictly better only at t>0. Correction 3: the gap mechanism and 'T=0 TFIM code' are the reviewer's words, and no gap was computed, yet CAMPAIGN_LEDGER:26 asserts 'the gap of -K tau tau'. Correction 5: D-conservation is a conditional constraint. WT, OPEN_CLAIMS and the assembled CAMPAIGN_LEDGER still carry the superseded s10 entry. The packet thus under-states s10's verification and over-states its energy-flow content. This is a new defect from the round-1 cycle.

**Deciding check.** grep -c 'energy flows the other way' CAMPAIGN_LEDGER.md returns 0. See ledger/ledger_s10_m1_memory.md:1-30. git log -- ledger/ledger_s10_m1_memory.md shows de88553 after 7d5c0bc, the commit that assembled the CAMPAIGN_LEDGER s10 section.

**Fix proposed.** Rebuild the CAMPAIGN_LEDGER s10 section and summary row from the corrected file. In WT s10 write 'the write releases register energy into the carrier'. Mark s10 as adversarially re-read (8 corrections) everywhere, and check s11 and s12 the same way. In WT sections 6 and 7 write 'no worse at every sampled time; equal at t=0, better at every t>0 sample'.

## 12. [minor] WORKING_THEORY.md

**Location.** Section 0:38-42: '... not the constants, not the action. Each negative is exhibited by a legal closed model in calculation 1'

**Problem.** FC section 4:268 marks the constants and the action BY DEFINITION, but WT section 0 still says every listed negative is exhibited by a closed model. The round-1 F13 fix did not reach section 0. The constants and the action are also not 'local pointwise' items.

**Deciding check.** FC section 4:268 versus WT:41-42.

**Fix proposed.** Write: 'Each dynamical negative is exhibited by a legal closed model; the constants and the action are inputs by definition (FC section 4).' Move them out of the 'local pointwise' sentence.

## 13. [minor] WORKING_THEORY.md

**Location.** Section 0:57-61 'The six entries of Theta: field content and action, preparation, access structure, measurement rule (P5), arrow (P4), time parameter'

**Problem.** This list is not the six slots of DEF 1.7 (algebra, dynamics, preparation, access, instruments and records, clock), which WT section 0:23 itself repeats. It merges algebra with dynamics, omits instruments and records, and adds the measurement rule (DEF 1.8, not a slot) and the arrow (part of the preparation).

**Deciding check.** FC:79-80 and WT:23 against WT:57-61.

**Fix proposed.** List the six DEF 1.7 slots, then 'plus the measurement rule (DEF 1.8); the arrow sits in the preparation'.

## 14. [minor] WORKING_THEORY.md; OPEN_CLAIMS.md; FORMAL_CORE.md

**Location.** FC section 3:243 includes 'no-signalling between commuting subsystems, and monogamy constraints'; WT section 0:30-34 and OPEN_CLAIMS:11-13 omit both.

**Problem.** The three forced lists disagree. Only instantaneous no-signalling is forced by closure; no-signalling over time needs local H (COND 2.12). The Lieb-Robinson consequence for 'bath = rest' is never drawn: with local H, the effective bath of A up to time t is its light cone.

**Deciding check.** Compare the three lists. The commutator argument covers simultaneous operations only; later influence through H_AB is generic.

**Fix proposed.** Use the same list in all three files. Write 'instantaneous no-signalling (closure); finite-speed signalling (COND 2.12)' and add one sentence on the light-cone bath.

## 15. [minor] FORMAL_CORE.md; WORKING_THEORY.md

**Location.** FC COND 2.10:212 'for every subsystem small compared with d_eff'. WT section 0:35-36, which omits the purity hypothesis. FC COND 2.11:220 'a subspace H_R of large dimension'.

**Problem.** The LPSW bound (1/2)sqrt(d_A^2/d_eff) is non-trivial only when d_A^2 < d_eff. Popescu-Short-Winter need the effective bath dimension of Tr_A(1_R/d_R) to be large compared with d_A (roughly d_R >> d_A^2), not merely d_R large. WT section 0 drops COND 2.10's purity hypothesis (Short 2011 covers mixed states). Page's average-entropy theorem (PRL 71, 1291 (1993)) is the COND 2.11 plus THM 2.9 corollary and goes unmentioned.

**Deciding check.** The LPSW and PSW theorem statements.

**Fix proposed.** Write 'd_A^2 << d_eff' and 'd_R >> d_A^2', and add 'for a pure global state (mixed: Short 2011)' in WT section 0. Optionally add one line on Page.

## 16. [minor] FORMAL_CORE.md

**Location.** Section 4 chaos row: 'SCOPE | ... out-of-time-order growth at intermediate times is model specific'

**Problem.** Maldacena-Shenker-Stanford (JHEP 08 (2016) 106) give lambda_L <= 2 pi T under a thermal state plus factorization and analyticity: a conditional theorem that bounds the exponent without making it small. The M1 compatibility file already invokes the e^{2 pi T t} channel (ledger s10).

**Deciding check.** MSS theorem statement.

**Fix proposed.** Rewrite the row as: SCOPE in the finite model; the upper bound lambda_L <= 2 pi T is FORCED under the MSS hypotheses; smallness is not forced.

## 17. [minor] CALCULATIONS.md

**Location.** S4:113-114 'The inhomogeneous term vanishes for product preparations, but not only for them.'

**Problem.** This is a residual from the round-1 NZ fix. The term vanishes for rho_A (x) rho_ref. A product rho_A (x) sigma with sigma != rho_ref generically has a nonzero term; FC THM 2.3(ii) itself says so.

**Deciding check.** FC THM 2.3(ii); round-1 formal F1.

**Fix proposed.** Write 'vanishes for preparations rho_A (x) rho_ref, and not only for them'.

## 18. [minor] WORKING_THEORY.md

**Location.** Section 5 s02 row:210 and section 6:266-268: the access seed is 'a physical boundary condition, not an axiom, not a derived object', filed under 'Proposed additions that did not survive'.

**Problem.** The quoted sentence is the P-5 instrument's wording, not the owner ruling (ledger s02 section 0 item 4, :107-108, :299). The owner's binding wording is 'the access boundary is currently an unresolved structural input', with 'access is arbitrary bookkeeping' FALSIFIED in-class. That falsification is a surviving positive in-class result: moving the boundary exposes out-of-access physics (0.38). It is filed among failures.

**Deciding check.** ledger/ledger_s02_influence_hierarchy.md:16, :107-108, :299.

**Fix proposed.** Attribute the phrase to the P-5 instrument and quote the owner wording. List 'access is not arbitrary bookkeeping (P-5, in-class)' as a surviving conditional result. Keep 'access is not endogenous (EA-0)' as the failed addition.

## 19. [minor] WORKING_THEORY.md

**Location.** Section 0 overall, plus the banner at lines 1-2

**Problem.** Readability: section 0 does answer all five completion-standard questions under labelled paragraphs, with these defects. (a) The thesis sentence the owner needs first arrives only in the final Status paragraph. (b) The 'did not survive' paragraph is a single 15-item run-on list. It mixes retracted empirical claims (689 Hz, 32 sigma ISW) with mechanisms, framings, an absence ('a GRUT-specific prediction of GRUT III') and a non-proposal (CP-12 'as a selection law'), while section 6 already groups the same items by Theta slot. (c) m1 and m2 above. (d) The banner is stale: it says sections await the campaign-ledger workflow, which has finished, and it contradicts the README status line.

**Deciding check.** WT:12-68 read against handoff:244-248; WT:1-2 against README:30-32.

**Fix proposed.** Lead with the thesis. Group 'did not survive' by slot as section 6 does, and drop absences and non-proposals from it. Refresh or remove the banner after this round.

## 20. [minor] FORMAL_CORE.md; WORKING_THEORY.md

**Location.** FC section 4 Born row: 'BY DEFINITION | ... dynamical closure is not one of the assumptions that the reconstruction theorems replace it with'. FC DEF 1.8. WT section 6 measurement-rule bullets.

**Problem.** Zurek's envariance derivation (PRL 90, 120404 (2003); PRA 71, 052105 (2005)) takes 'the bath is the rest and the whole is pure' as its central premise, plus contested further assumptions (Schlosshauer-Fine, Found. Phys. 35, 197 (2005)). The program's withdrawn HJW route is a variant of it. The row is right about Gleason, Busch and MGM, but it omits the one route that uses the parent premise. The Born rule is therefore 'by definition here; conditionally derivable from closure plus contested assumptions'.

**Deciding check.** grep -i 'envarian|Zurek' over FC and WT finds nothing.

**Fix proposed.** Add a sentence with these citations to the Born row and DEF 1.8, noting that the program's HJW variant was withdrawn (s05).

## Confirmed by this lens

- The central negative survives the steelman. No candidate examined is a law beyond known physics or a state-selection rule, so 'No additional law or state-selection rule has survived' (WT sections 0 and 7) is correct.
- Status of each candidate. Forced by closure alone: athermality budget, protect/form trade-off, no complete self-measurement, instantaneous no-signalling, monogamy. Forced on the ontic reading: purity of each whole. Closure plus a named hypothesis: LPSW, records-need-conserved-data, canonical typicality and Page, ETH, return to equilibrium and Davies, Clausius/Landauer, KMS=>FDT, Lieb-Robinson, WAY, clock imperfection, quantum Darwinism, Wigner's-friend no-gos, relational access, the MSS bound, envariance. Not forced: a specific state, the constants, the action, compact topology, cosmological energy conservation.
- Under-claimed or misfiled GRUT-internal results: the CP-12 relative-entropy and min-entropy theorems (F1); M1's protect/form trade-off (F4); GRUT III T1/T2, the Delta rho instrument and the ten-line sixth-postulate sheet (F10); the s02 in-class falsification of 'access is arbitrary bookkeeping' (m7). ARROW_OF_TIME's existence/direction decomposition is reported correctly as a localization of an import.
- COND 2.12 (Lieb-Robinson) and COND 2.13 (WAY: additive conserved L, apparatus inside Omega) are stated correctly. WAY does not apply to calc 1, whose only conserved symmetry, Z-parity, is multiplicative, and the packet claims no such application.
- The S9 algebra is correct: equivariance with U = e^{iHs} gives [F(H), H] = 0. Reproduced: ||[F,H]|| = 8.9e-15 for Gibbs and 1.6e-14 for the ground state. Only the interpretation 'dynamically empty' is overstated (F6).
- The S10 method (LPSW with the exact HS^2 average) is sound. My X-record and Z-record runs use the same machinery and agree with the theorem.
- FC COND 5.3 (Fewster-Verch) now carries the primary-text hypotheses. FC section 7 is slot-by-slot with CPR scoped, and the marginal problem is attached to bottom-up families only.
- The section 4 statuses (8 NOT FORCED, 3 BY DEFINITION, 3 SCOPE) are counted correctly in OPEN_CLAIMS item 4.
- Round-1 minor fixes are all present: MGM, TCL invertibility, rho_ref, BDF, KMS sign, FV formula, clock wording, CKW, Broadbent-Grilo, THM 2.4, S4 column, Araujo, V7 attribution, toric qualifiers, O-6, differentiator ranges, K1 qualifier, counts, Araki-Lieb slack.
- WT section 0 answers all five completion-standard questions under labelled paragraphs, subject to the m1, m2 and m8 defects.
- Candidate 1, M1 and GRUT III are kept as separate model instances with no common law (WT section 7); no cross-import found.
