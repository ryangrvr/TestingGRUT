# Campaign ledger entry — S10: M1 memory (key `s10_m1_memory`)

Sector: M1 memory — finite Ising repetition register with matter feedback and programmed writer; facilitation benefit withdrawn.

Checkout note. The task names commit `96a7e7e` as the checkout; the working tree is actually at `da8f925` ("WIP checkpoint: closed-whole working-theory packet"), a direct descendant of `96a7e7e`. Both M1 directories are byte-identical at `4e9a161`, `96a7e7e` and `da8f925` (`git diff --stat` empty for both directories across those commits), so every citation below is valid at all three.

Path abbreviations used in citations:
- `pkt/` = `research/reality_framework/grut_memory_candidate_20261019/`
- `rev/` = `research/reality_framework/grut_memory_candidate_20261019_review/`
- `hand/` = `research/reality_framework/grut_closed_whole_handoff_20261019/NEXT_CHAT_HANDOFF.md`
- JSON files are cited by key path rather than line, since the exact locator is the key.

---

## 1. Question asked in this sector

The sector asked whether one closed, time-independent Hamiltonian with no external bath, drive, disorder or controller can (a) retain a readable, redundantly encoded classical record while the record alters later matter interactions and while matter information spreads; (b) write such a record autonomously from internal resources; and (c) whether a proposed extra law, "state-dependent kinetic saturation" (facilitation), is what protects the record (`pkt/protocol.json:5-8,21`; `pkt/writer_protocol.json:11`; `pkt/REPORT.md:44-47`). It was explicitly scoped away from chaos-exponent and gravitational claims (`pkt/protocol.json:41-42`; `pkt/writer_protocol.json:29`).

Under the parent closed-whole question, the test was: does closure (internal bath only, exact unitary evolution, explicit energy accounting) force a new memory principle here, or does the record persist for a reason known physics already contains?

## 2. Assumptions / model class

- State space: a ring of L≥3 cells, each a matter qubit σ_i and a record qubit τ_i; for the writer, additionally L²(R,dQ)⊗C²_R for a massive carrier with an internal source bit (`pkt/GRUT_CANDIDATE_M1.md:17-32`).
- Core law: H_M1 = Σ_i [ J(1+ατ_i^z)σ_i^zσ_{i+1}^z + h_xσ_i^x + h_zσ_i^z − Kτ_i^zτ_{i+1}^z − ΔF_iτ_i^x ], F_i = ε + (1−ε)(1−τ_{i−1}^zτ_{i+1}^z)/2 (`pkt/GRUT_CANDIDATE_M1.md:34-43`). Writer completion: H = P²/(2M) + H_M1 + Ω(1−R^z)/2 Σ_i f_i(Q)τ_i^x with fixed cos² gate profiles of half-width w (`pkt/GRUT_CANDIDATE_M1.md:53-66`).
- Fixed benchmark: ħ=J=1, (α,h_x,h_z,K,Δ,ε) = (0.5,1.05,0.5,1,1,0.05) (`pkt/GRUT_CANDIDATE_M1.md:93-96`); writer L=4, M=0.5, p₀=100, q₀=−4, σ_Q=0.35, w=0.3, Ω=π(p₀/M)/(2w) (`pkt/GRUT_CANDIDATE_M1.md:99-101`), numerically Ω=1047.1975511965977 (`pkt/writer_results_initial.json: parameters.Omega`).
- Supplied postulates, by the packet's own table: quantum mechanics, a fixed local substrate (ring, two qubits per cell), the kinetic-saturation projector, physical feedback through α, and modeled preparation resources (source qubit, massive particle, packet, blank register) (`pkt/GRUT_CANDIDATE_M1.md:75-89`). The small ratio ε and the feedback are supplied, "neither is derived from the words 'memory' or 'closure'" (`pkt/GRUT_CANDIDATE_M1.md:49-51`).
- External Hamiltonian time is fundamental; there is no gauge or gravitational constraint (`pkt/GRUT_CANDIDATE_M1.md:29-32`).
- Preparations are pure prepared states with nonzero excitation energy, "not a ground-state preparation or a demonstrated finite-temperature memory phase" (`pkt/WORKED_REGIME.md:124-125`). Sizes L=3,4,5 only; sampled time grids retention [0,40] step 0.2 (201 points), readout [0,2] step 0.02 (`pkt/protocol.json:32`; `pkt/calculate.py:110-111`).
- Success observable is the equal-prior local guessing probability P_A = 1/2 + |⟨A⟩₊−⟨A⟩₋|/4, "not the optimal global information-recovery probability" (`pkt/GRUT_CANDIDATE_M1.md:133-141`).
- The matched control replaces only F_i by the constant ε=0.05 with everything else fixed (`pkt/REPORT.md:18-20`; `pkt/control_protocol.json:7`). It was specified after the review's outcome was known: "This is not a blinded or held-out prediction" (`pkt/control_protocol.json:6`).

## 3. Executed evidence

### 3a. Core retention and readout at the M1 point (all 13 registered conditions pass)

| Run | min memory success through t=10 | success at t=40 | best matter readout |
|---|---:|---:|---:|
| M1, L=3 | 0.994999 | 0.985006 | 0.806447 (t=1.20) |
| M1, L=4 | 0.987510 | 0.935065 | 0.778468 (t=0.58) |
| M1, L=5 | 0.989132 | 0.940399 | 0.778460 (t=0.58) |

(`pkt/WORKED_REGIME.md:75-79`; `rev/control_results.txt:1,3,5`). Full precision: L=4 minimum_memory_t10 = 0.9875095925556739, memory_t40 = 0.9350645853495996, best_readout = 0.7784678210302141 at 0.58 (`pkt/summary.json: core[1]`). The L=4 readout is P_Y = 1/2 + |−0.491561535 − 0.622309749|/4 = 0.778467821 with simultaneous memory success 0.999312115 (`pkt/WORKED_REGIME.md:87-95`). The seven original runs took 11.79 s and all 13 registered acceptance conditions passed (`pkt/EXECUTION.md:11-16`; `pkt/results_initial.json: all_registered_acceptance_pass = true`, 13 `checks`). The registered acceptance thresholds were memory ≥0.9 through t=10 and readout >0.6 with contemporaneous memory ≥0.9 (`pkt/protocol.json:35-36`).

Original controls (all L=4): F=1 "unfacilitated" 0.621452 / 0.709936 / 0.741417; feedback removed (α=0) 0.996241 / 0.994449 / 0.500000; exactly frozen (ε=0) 1.000000 / 1.000000 / 0.778577 (`pkt/WORKED_REGIME.md:80-82`); half softness ε=0.025: 0.9968211359193733 / 0.9792849192829733 / 0.7785498324346294 (`pkt/summary.json: core[6]`).

### 3b. The decisive matched control (F_i ≡ ε = 0.05)

| L | M1 min through t=10 | uniform-floor min through t=10 | M1 at t=40 | uniform-floor at t=40 | M1 min over [0,40] | uniform-floor min over [0,40] | best readout M1 / control |
|---|---:|---:|---:|---:|---:|---:|---:|
| 3 | 0.994999 | 0.998738 | 0.985006 | 0.994692 | 0.984591 | 0.994266 | 0.806447 / 0.806626 |
| 4 | 0.987510 | 0.999055 | 0.935065 | 0.998458 | 0.935065 | 0.998397 | 0.778468 / 0.778483 |
| 5 | 0.989132 | 0.998479 | 0.940399 | 0.997517 | 0.940399 | 0.995717 | 0.778460 / 0.778475 |

(`pkt/REPORT.md:22-26`; `pkt/WORKED_REGIME.md:83-85`; `rev/control_results.txt:1-6`; `rev/REVIEW.md:25-32`). Full precision at L=4: uniform_floor_minimum_t10 = 0.9990547948737153, uniform_floor_t40 = 0.998457849859687, uniform_floor_best_readout = 0.7784834017608981 at 0.58 versus candidate 0.7784678210302141; readout_best_probability_difference = 1.5580730684083832e-05 (`pkt/summary.json: matched_comparison.4`). Best-readout differences at L=3,4,5 are 0.000179931, 0.000015581 and 0.000015652 (`pkt/EXECUTION.md:46`).

The control is "no worse at every sampled time through t=40, within 1e-10 numerical tolerance. This is not an all-continuous-time theorem" (`pkt/REPORT.md:28-29`). Implementation: `uniform_no_worse_on_all_samples = bool(delta.min() >= -1e-10)` on the 201-point grid (`pkt/calculate.py:110,301`). Raw minima of (uniform − candidate) are −2.7755575615628914e-15 (L=3), 8.881784197001252e-16 (L=4), −1.7763568394002505e-15 (L=5), i.e. machine-epsilon ties at t=0, which is why the protocol prescribes the 1e-10 tolerance "near t=0" (`pkt/summary.json: matched_comparison.{3,4,5}.minimum_sampled_uniform_minus_candidate`; `pkt/control_protocol.json:14`). Recorded mechanism claims: `facilitation_improves_minimum_retention_all_sizes = false`, `matched_uniform_control_retains_better_on_registered_samples = true` (`pkt/summary.json: mechanism_claims`).

Executed by two different constructions: the review's independent tensor-product script (`rev/uniform_field_control.py:10-19`, run ≈1 min, `rev/REVIEW.md:56`), re-executed with one BLAS thread and reproducing all six archived rows (`pkt/EXECUTION.md:30-32`), and the packet's bit-table `calculate.py` ten-run replay in 7.31415 s (`pkt/EXECUTION.md:34-38`), with bit-table and Kronecker matrices agreeing at every size (`pkt/summary.json: control_method_checks.matched_tensor_constructions = true`; `pkt/calculate.py:307-309`). A third execution exists at HEAD, outside this packet: the WIP checkpoint reran the unchanged review script (Python 3.13.16, numpy 2.5.3, scipy 1.18.1, 8.6 s) and reports output "IDENTICAL (diff empty)" to `control_results.txt` (`research/reality_framework/grut_closed_whole_theory_20261019/results/m1_matched_control_rerun.txt:1-12`). That file is a secondary cross-check, not a primary source for this ledger.

Mechanism reading from the review: "the memory is the ferromagnetic gap of −Kτ^zτ^z against a small transverse field; that is a zero-temperature transverse-field Ising repetition code, and it works for the same reason in both models"; facilitation is inactive on an ordered record and, after one leak, "turns both neighbours' flip amplitude from 0.05 to 1, which converts an off-resonant single-flip leakage into resonant spreading (an East-model cascade)" (`rev/REVIEW.md:36`). The packet adopts this: after a defect, opposite-neighbour flips at amplitude 1 "move domain walls without changing the ferromagnetic energy, rather than annihilating the error ... This reversible facilitation does not supply an error-correction bias" (`pkt/WORKED_REGIME.md:116-122`). "The old F=1 comparison increases the transverse field twentyfold and does not isolate neighbor dependence" (`pkt/WORKED_REGIME.md:101-102`).

### 3c. Exact regime at ε = 0

With D = ½Σ_i(1−τ_i^zτ_{i+1}^z), [H_M1(0), D] = 0; the matrix commutator is exactly zero at L=3,4,5 (`pkt/WORKED_REGIME.md:13-22`; `pkt/algebra_results.json: runs[*].hard_constraint_commutator_max = 0.0`, `checks.exact_domain_wall_conservation_at_epsilon0 = true`). In the zero-wall sector the record is perfectly retained while matter evolves under a record-dependent H_σ^{(r)} (`pkt/WORKED_REGIME.md:24-37`). Formation cost: an all-+x memory has zero-wall probability 2^{1−L} (0.25, 0.125, 0.0625 at L=3,4,5) and unitary evolution at ε=0 cannot increase it (`pkt/WORKED_REGIME.md:41-44`; `pkt/algebra_results.json: runs[*].initial_ordered_sector_probability`). At the soft point the ordered-sector probability peaks at 0.3289007412515012 / 0.3022736769655148 / 0.08277297473471425 and ends at 0.2614344170727129 / 0.29884391381968967 / 0.04729380180227547 for L=3,4,5 — "reported without claiming a persistent organized state" (`pkt/algebra_results.json: runs[*]`; `pkt/EXECUTION.md:23-26`).

Short-time analytic checks reproduced by matrix calculation: ⟨τ_0^z(t)⟩_r = r[1 − 2Δ²ε²t² + O(t³)]; readout slopes 4J(1+αr)+2h_z = 7 (r=+1), 3 (r=−1); adjacent matter squared-commutator coefficients 8J²(1+αr)² = 18 and 2 (`pkt/WORKED_REGIME.md:47-71`; `pkt/protocol.json:44-48`).

### 3d. Registered stronger attacks that failed (six conditions false, retained)

- Single-cell bit: local success at t=10 is 0.801481 (L=4) and 0.751472 (L=5), below the 0.9 standard (`pkt/WORKED_REGIME.md:155-157`; `pkt/challenge_results.json: local_bit.4.p10 = 0.801480851769999, local_bit.5.p10 = 0.751471559810829`).
- Bare matter-to-memory write: maximum success 0.553352 (L=4) and 0.526376 (L=5), "close to chance" (`pkt/WORKED_REGIME.md:158-160`; `pkt/OPEN_CLAIMS.md:22-23`; `pkt/challenge_results.json: autonomous_write.{4,5}.maximum`).
- Unconstrained flip −κΣτ_i^x: minimum success through t=10 falls to 0.875748 (κ=0.1) and 0.764940 (κ=0.3) (`pkt/WORKED_REGIME.md:161-165`; `pkt/challenge_results.json: unfacilitated_perturbations[*].minimum_t10`).
- Recorded as six false conditions alongside three true numerical-method checks (`pkt/summary.json: adversarial_conditions`; `pkt/verify.py:55-57`; `pkt/EXECUTION.md:17-22`). The challenge protocol was registered after the initial results but before the attacks ran (`pkt/challenge_protocol.json:2`).

### 3e. Writer (M1-W)

- Mechanism: carrier speed p₀/M = 200, gate crossing 2w/(p₀/M) = 0.003; conditional pulse area Ωw/(p₀/M) = π/2, "giving a controlled flip. This leading operation does not require facilitation" (`pkt/COMPATIBILITY.md:34-37`; `rev/REVIEW.md:42-44`; `pkt/REPORT.md:75-78`).
- Initial output in all four cells: "At least 0.999680" at t=0.06 (`pkt/REPORT.md:91`). Derived here from the stored per-cell means at storage t=0: minimum cell 0.999680152 (dt 0.0002, 2048 points), 0.999680089 (dt 0.0001, 2048), 0.999680260 (dt 0.0001, 4096) (`pkt/writer_results_refined.json: runs[*].{plus,minus}.storage.memory_means[*][0]`, via P = ½+|m₊−m₋|/4).
- All-cell success through sampled storage t=10: ≥ 0.991697 (`pkt/REPORT.md:92`; `pkt/summary.json: writer_minimum_memory_t10 = 0.9916972267188419`; coarse run 0.9916969879289492 in `pkt/writer_results_initial.json: runs[0]`). This storage stage "neglects residual carrier/gate interaction" (`pkt/REPORT.md:92,110-112`; `pkt/EXECUTION.md:88-94`).
- Matter readout 0.769977 at storage t=0.52 with memory 0.998568 (`pkt/REPORT.md:93`; `pkt/writer_results_refined.json: runs[2].readout.best_success = 0.769977434073726, memory_at_best = 0.9985675362234426`).
- Refinement: maximum probability change below 6.8e-7 (`pkt/REPORT.md:94`; `pkt/summary.json: writer_refinement.overall.max_memory_probability_difference = 6.705103664517154e-07`).
- Energy: carrier mean kinetic energy 10002.0408163 → 10004.0290991 in the writing branch; core 1.55 → −0.437569946; total-energy drift 0.000712816 J, "not zero; its cutoff behavior is nonmonotonic" (`pkt/REPORT.md:116-121`). Per run, writing-branch drift is 0.00010979246144415811 (dt 0.0002, 2048), 1.0308212949894369e-08 (dt 0.0001, 2048), 0.0007128160505089909 (dt 0.0001, 4096) (`pkt/writer_results_refined.json: runs[*].minus.total_energy_drift`). Plus-branch drift is at the 1e-11 level (`runs[*].plus.total_energy_drift`).
- Diagnostics: final norm error 1.33e-12, residual gate-region probability 6.87e-10, edge 1.76e-9, negative momentum 1.36e-8; "do not prove that a later carrier return is impossible" (`pkt/REPORT.md:123-128`; `pkt/writer_results_refined.json: runs[2].minus`).
- Execution history: coarse run 414.59 s; the slower dt=0.0001 local-Pauli replay was interrupted with Ctrl-C; `execution_amendment.json` replaced it by `writer_replay.py` (exact core exponentials, 2048 and 4096 clock points) — "an execution interruption, not an erased physical failure or a completed refinement" (`pkt/EXECUTION.md:64-78`; `pkt/execution_amendment.json: reason, preservation, new_execution`). Replay wall time 314.89 s (`pkt/writer_results_refined.json: wall_seconds = 314.8902410519995`).
- All six writer conditions pass at registered thresholds (`pkt/summary.json: writer_checks`; thresholds in `pkt/writer_protocol.json:21-30`). "These are internal numerical results, not independent validation" (`pkt/REPORT.md:103-106`). No writer scattering was rerun on the uniform-floor core (`pkt/EXECUTION.md:32`; `pkt/control_protocol.json:16`; `pkt/COMPATIBILITY.md:36-37`).

### 3f. Numerical integrity, replay and finalizer

- Core: largest norm deviation 1.29e-14; total-energy deviation below 3.07e-14 J; Kronecker/Krylov agreement 2.24e-13 on the retention trace and 5.2e-15 on readout means (`pkt/WORKED_REGIME.md:129-133`; `pkt/challenge_results.json: independent_full_H_max_error = 4.163336342344337e-17`).
- Ten-run replay reproduces the original seven conditional traces: `original_replay_max_conditional_mean_error = 0.0` (`pkt/MANIFEST.json`), asserted < 1e-10 (`pkt/verify.py:46-52`).
- The finalizer hashes every file except `MANIFEST.json` — 31 entries — plus the manifest it writes (`pkt/verify.py:84-95`; `pkt/MANIFEST.json: files` length 31). A read-only sha256/size comparison performed for this ledger matches all 31 entries at HEAD with no unlisted files. `verify.py` itself was not executed here because it rewrites `MANIFEST.json` (`pkt/verify.py:95`); the recorded pass is in `AGENT_COORDINATION.md:5822`.
- Review values reproduced within 5.01e-7 of the six printed rows (`pkt/verify.py:63-71`; `pkt/control_protocol.json:14`).

### 3g. Equilibrium and scrambling diagnostics (no exponent)

At β=1 the normalized connected memory autocorrelation is 0.538351 at t=10 with plateau 0.458812; the adjacent matter squared commutator reaches 1.41102 within [0,4]; the centered regulated OTOC still has initial factorization defect 0.607085; "No λ_L/(2πT) is reported" (`pkt/WORKED_REGIME.md:177-190`; `pkt/summary.json: lyapunov_exponent = null, gravitational_prediction = null`; `pkt/challenge_results.json: centered_regulated_otoc.lyapunov_exponent = null`).

## 4. Answer

What works, in the stated finite regime: a prepared ferromagnetic repetition register (L=3–5, one classical bit repeated) retains the bit at ≥0.987510 through sampled t=10 and ≥0.935065 at t=40 (L=4), while a τ-dependent matter bond gives a 0.778468 local matter readout; at ε=0 retention is exact by domain-wall conservation; a closed, time-independent flying-qubit writer copies a source bit into all four cells at ≥0.999680 with the register's energy change paid by the carrier (`pkt/REPORT.md:84-95,116-119`; `pkt/WORKED_REGIME.md:13-44`).

What was withdrawn: the attribution of that retention to the facilitation (kinetic-saturation) law. The matched F≡ε control retains better at every sampled time and size and gives nearly the same readout, so "the demonstrated retention therefore does not support the added facilitation principle. The finite calculation remains correct; the attribution of its protection to that principle is withdrawn" (`pkt/REPORT.md:3-8`; `pkt/OPEN_CLAIMS.md:16-19`). The original 1997d20 headline, "a constructed and calculated memory candidate" whose "added principle is kinetic saturation" (`git show 1997d20:pkt/REPORT.md`, lines 1-16), and its classification "a narrower working model rather than the intended universal theory" (lines 48-53 of that file) were replaced at 4e9a161 by "a finite Ising repetition register with matter feedback and a programmed flying-qubit writer; the proposed facilitation benefit is rejected in the demonstrated regime" (`pkt/OPEN_CLAIMS.md:56-61`).

What is not established: a finite-temperature memory phase (`pkt/REPORT.md:39-40`; `pkt/COMPATIBILITY.md:41-42`; the reviewer adds that a 1D Ising ring has none by Peierls, `rev/REVIEW.md:48`, which is the reviewer's statement and not a packet computation); cloning of an unknown quantum state ("the model does not clone an arbitrary unknown quantum state", `pkt/GRUT_CANDIDATE_M1.md:143-144`); gravity ("Missing. No gravitational degrees of freedom or gravitational equation have been supplied", `pkt/COMPATIBILITY.md:8`; `pkt/REPORT.md:155-167`); a chaos exponent ("no exponent or ratio to the MSS bound has been extracted", `pkt/REPORT.md:159-161`); all-times retention or a scattering certificate (`pkt/REPORT.md:29,111-112,126-128`); independent local capacity, thermodynamic localization, robustness to thermal inclusions (`pkt/OPEN_CLAIMS.md:38-40`); autonomous formation from generic resources (`pkt/OPEN_CLAIMS.md:36-37`); novelty or empirical viability (`pkt/OPEN_CLAIMS.md:50-52`).

## 5. Closure verdict

Nothing additional was forced beyond known physics in this class, and the computation that shows it is the matched-control table in §3b. The closed law was specified with no external bath, reservoir, pulse or postselection (`pkt/protocol.json:8,21`; `pkt/writer_protocol.json:11`; `pkt/COMPATIBILITY.md:14-17`). Within that closed description the record persists for the standard reason — the gap of −Kτ^zτ^z against a weak transverse field, a T=0 transverse-field Ising repetition code (`rev/REVIEW.md:36`; `pkt/REPORT.md:38-39`) — and the one candidate extra law makes retention worse (L=4, t=40: 0.935065 vs 0.998458; `pkt/REPORT.md:25`). What closure did supply is ordinary bookkeeping that known quantum mechanics already contains: energy conservation makes the carrier pay the register's energy change (core 1.55 → −0.437569946 as carrier 10002.0408163 → 10004.0290991, `pkt/REPORT.md:116-118`; `rev/REVIEW.md:42`); the nonselective readout removes −1.2010144303 and −0.8997375490 from the modeled system, which an energy-conserving detector completion must receive (`pkt/WORKED_REGIME.md:143-148`); and a conserved quantity at ε=0 forbids gaining weight in the protected sector ("Perfect protection is not a mechanism for generating the protected state", `pkt/WORKED_REGIME.md:41-44`). The packet's own statement of the gravitational implication is that closure plus structured memory "does not by itself bound the largest exponent of the other operators" and that an eikonal G_N s channel with s ∝ e^{2πTt} is untouched by a coefficient change (`pkt/COMPATIBILITY.md:60-77`). The result is "a result, not a failure" in the parent program's sense: a finite working construction in known physics, with the proposed addition rejected by its own numbers (`pkt/REPORT.md:169-176`).

## 6. Classification

Mixed, itemized:
- Recovered known physics: retention of a prepared repetition bit by a ferromagnetic gap against a weak transverse field; a π/2 pulse-area controlled flip as the writer; energy conservation in closed scattering; exact conservation of D at ε=0 (`rev/REVIEW.md:36,42`; `pkt/WORKED_REGIME.md:13-22`).
- Rejected mechanism: the facilitation-benefit interpretation of kinetic saturation, rejected by the matched F≡ε control at L=3,4,5 (`pkt/REPORT.md:137-139`; `pkt/summary.json: mechanism_claims`).
- Supplied model choice: ring topology, two qubits per cell, the projector form of F_i, ε=0.05, α=0.5, the parameter point, writer apparatus parameters, choice of local guessing observable (`pkt/GRUT_CANDIDATE_M1.md:75-89,93-101`; `pkt/PREDICTIONS_AND_TESTS.md:81-83`).
- Conditional constraint: all probabilities hold only at L=3–5, on the sampled grids, for the pure prepared states, in the benchmark; the ε=0 theorem holds only at exact saturation; the storage numbers hold only under the stated post-scattering decoupling approximation (`pkt/WORKED_REGIME.md:95-97`; `pkt/PREDICTIONS_AND_TESTS.md:72-79`; `pkt/EXECUTION.md:88-94`).
- Conjecture: an error-correction-by-cascade regime where facilitation helps against the F=ε baseline (`rev/REVIEW.md:54`; `pkt/PREDICTIONS_AND_TESTS.md:85-88`); any relation to Candidate 1 as a low-energy sector (`pkt/COMPATIBILITY.md:87-92`).
- Unexecuted task: the writer scattering on the uniform-floor core (`pkt/control_protocol.json:16`); the cascade calculation (`pkt/COMPATIBILITY.md:41-44`); thermal-inclusion/robustness attacks (`pkt/OPEN_CLAIMS.md:38-40`); a dynamical terminal detector and reset cost (`pkt/OPEN_CLAIMS.md:48-49`); the interrupted dt=0.0001 local-Pauli replay, superseded rather than completed (`pkt/EXECUTION.md:67-70`); the registered no-boundary B control, referenced from this packet as still unexecuted (`pkt/COMPATIBILITY.md:94-97`).

## 7. Remaining freedom

- The entire parameter point and ε, α, K, Δ, h_x, h_z; no calibration to seconds or to any physical system exists (`pkt/WORKED_REGIME.md:6-7`; `pkt/PREDICTIONS_AND_TESTS.md:8-14,81-83`).
- Topology (ring), dimensionality, size beyond L=5, choice of record encoding and of local versus optimal readout (`pkt/GRUT_CANDIDATE_M1.md:80,139-142`).
- Whether to keep the facilitation term at all: the packet retains the original M1 law "rather than silently replaced by the better control" (`pkt/GRUT_CANDIDATE_M1.md:151-152`); the control is itself the better-retaining model in this regime.
- A successor cascade model would need "noise, an error syndrome, correction dynamics, entropy/energy resources and an advantage over this matched baseline" (`pkt/COMPATIBILITY.md:41-44`).
- A gravitational successor would have to "specify which horizon/scattering assumption its new law changes, calculate a regulated gravitational observable, and recover the required low-energy gravitational tests" (`pkt/COMPATIBILITY.md:79-83`).
- Rejection rule for a calibrated version: disagreement with the full conditional traces beyond stated uncertainties rejects that calibrated model, "not ... every idea called GRUT" (`pkt/PREDICTIONS_AND_TESTS.md:45-49`).

## 8. Failures retained

1. Facilitation-benefit attribution withdrawn after the matched control (`pkt/REPORT.md:3-8,137-139`; `pkt/OPEN_CLAIMS.md:16-19`; `AGENT_COORDINATION.md:5820`).
2. The original F=1 control was confounded (twentyfold field) and "cannot isolate the benefit of neighbor dependence" (`pkt/REPORT.md:32-33`; `rev/REVIEW.md:21`).
3. Single-cell retention fails the 0.9 standard: 0.801481 (L=4), 0.751472 (L=5) at t=10 (`pkt/WORKED_REGIME.md:155-157`).
4. Bare autonomous write fails: 0.553352 and 0.526376 maxima (`pkt/WORKED_REGIME.md:158-160`).
5. Unconstrained-flip perturbations break the 0.9 condition: 0.875748 (κ=0.1), 0.764940 (κ=0.3) (`pkt/WORKED_REGIME.md:161-165`).
6. At ε=0 the ordered sector cannot gain weight from a disordered superposition — a formation barrier (`pkt/OPEN_CLAIMS.md:25-28`).
7. The dt=0.0001 local-Pauli writer replay was interrupted (Ctrl-C) and replaced, not completed (`pkt/EXECUTION.md:67-70`).
8. Writer total-energy drift 0.000712816 J is nonzero and nonmonotonic with cutoff; "the probability refinement is much stronger than the energy-convergence evidence" (`pkt/REPORT.md:119-121`).
9. No all-times scattering certificate; a later carrier return is not excluded (`pkt/REPORT.md:126-128`).
10. The matched control is post-hoc, "not a blind prediction" (`pkt/REPORT.md:138-139`; `pkt/control_protocol.json:6`).
11. Core-condition passing "does not establish the proposed principle" (`pkt/REPORT.md:98-100`); the explicit writer "does not retroactively pass the failed bare-writing protocol" (`pkt/REPORT.md:144-145`).
12. No finite-temperature memory phase, no Lyapunov exponent, no gravity, no priority claim (`pkt/OPEN_CLAIMS.md:30-52`).

## 9. What must not be promoted

- Do not cite "≥0.987510 retention" or any M1 number as evidence for a memory principle; the matched control does better (`pkt/REPORT.md:98-100`; `pkt/PREDICTIONS_AND_TESTS.md:50-57`).
- Do not promote "no worse at all 201 sampled times" to a continuous-time or thermodynamic theorem (`pkt/REPORT.md:29`; `pkt/WORKED_REGIME.md:111-114`).
- Do not read the ≥0.991697 storage number as a result for the complete carrier Hamiltonian; it neglects residual coupling (`pkt/REPORT.md:110-112`).
- Do not call the 13 passed core conditions or the 31-file integrity pass physical validation ("an integrity pass does not turn them into successful physics", `pkt/README.md:62-64`; `pkt/PREDICTIONS_AND_TESTS.md:78-79`).
- Do not transfer writer scattering probabilities to the uniform-floor core; they were not recomputed (`pkt/COMPATIBILITY.md:36-37`).
- Do not call the bond feedback gravity, the squared commutator 1.41102 a Lyapunov exponent, or λ_L<2πT a prediction (`pkt/COMPATIBILITY.md:8`; `pkt/WORKED_REGIME.md:182-190`; `pkt/PREDICTIONS_AND_TESTS.md:65-70`).
- Do not import the toric-branch claims quoted in `rev/REVIEW.md:48` as M1 results; the packet explicitly does not (`pkt/SOURCES.md:13-15`).
- Do not present the ε=0 exact retention as autonomous formation (`pkt/WORKED_REGIME.md:126-127`).
- Do not use the handoff as a primary source for any number above.

## 10. Sources read

Primary (all at `da8f925`, identical to `4e9a161` and `96a7e7e`):
- `pkt/REPORT.md` 1-185; `pkt/OPEN_CLAIMS.md` 1-61; `pkt/COMPATIBILITY.md` 1-97; `pkt/WORKED_REGIME.md` 1-190; `pkt/PREDICTIONS_AND_TESTS.md` 1-88; `pkt/EXECUTION.md` 1-99; `pkt/SOURCES.md` 1-30; `pkt/GRUT_CANDIDATE_M1.md` 1-167; `pkt/README.md` 1-64.
- `pkt/summary.json`, `pkt/execution_amendment.json`, `pkt/MANIFEST.json`, `pkt/results_with_control.json` (top-level keys and matched_comparison), `pkt/results_initial.json` (keys, 7 runs, 13 checks, wall time), `pkt/writer_results_refined.json`, `pkt/writer_results_initial.json`, `pkt/challenge_results.json`, `pkt/algebra_results.json` — scalar fields printed with `python3 -I`.
- `pkt/protocol.json`, `pkt/control_protocol.json`, `pkt/challenge_protocol.json`, `pkt/writer_protocol.json` in full; `pkt/verify.py` 1-103 in full; `pkt/calculate.py` via grep (lines 110-111, 154-155, 187, 229-312).
- `rev/REVIEW.md` 1-56; `rev/control_results.txt` 1-6; `rev/uniform_field_control.py` 1-42.
- `git show --stat` of `1997d20`, `cc24967`, `4e9a161`; `git show 1997d20:pkt/REPORT.md` lines 1-45 and `1997d20:pkt/OPEN_CLAIMS.md` lines 1-20, 44-53; `AGENT_COORDINATION.md` 5806-5814 (1997d20 entry) and 5816-5822 (4e9a161 entry).
- `hand/` lines 1-60, 115-200, 225-260 (row at line 123; numerical facts 138-162).
Secondary: `research/reality_framework/grut_closed_whole_theory_20261019/results/m1_matched_control_rerun.txt` 1-12.

## 11. Handoff-row discrepancies

Row text (`hand/:123`) and numerical block (`hand/:138-162`) were checked against the primary sources. No numerical contradiction was found: all twelve table entries, the 201-sample / 1e-10 statement, the L=4 readouts .778483402 vs .778467821 at t=.58, the ≥.999680 per-cell, ≥.991697 storage, .000712816 J drift, six false attack conditions, the exact seven-trace replay (max error 0.0) and the 31-file finalizer all match (§3b, §3d, §3e, §3f). Nuances and omissions, recorded rather than treated as errors:
1. "No worse at all 201 sampled times ... within 1e-10 tolerance" (`hand/:146-147`) rests on the tolerance at t=0: raw (uniform − candidate) minima are −2.78e-15 at L=3 and −1.78e-15 at L=5 (`pkt/summary.json: matched_comparison.{3,5}`), as `pkt/control_protocol.json:14` anticipates ("near t=0"). The handoff drops that qualifier.
2. The row says the programmed writer "works in the stated finite regime" but does not say that writer scattering was never recomputed for the better-retaining uniform-floor core (`pkt/control_protocol.json:16`; `pkt/COMPATIBILITY.md:36-37`).
3. "All recorded integrity checks pass" (`hand/:154`) is a recorded outcome (`AGENT_COORDINATION.md:5822`; `pkt/MANIFEST.json`); this ledger re-verified only the 31 file hashes read-only and did not execute `verify.py`, which rewrites `MANIFEST.json`.
4. The handoff says the "ten-run replay reproduces all seven original conditional traces exactly" (`hand/:151-152`); the source assertion is `< 1e-10` with a recorded maximum of 0.0 (`pkt/verify.py:52`; `pkt/MANIFEST.json`). Consistent, but "exactly" is the recorded value, not the test criterion.
5. The task text names the checkout commit as `96a7e7e`; the tree is at `da8f925`. The handoff itself says only that its commit "will be a descendant of 4e9a161" (`hand/:33-35`). Sector files are identical at all three commits, so no citation is affected.

## 12. Unread or missing

- The toric branch `claude/toric4d-memory-vs-scrambling-20261019` (commit `8e0f89f`) exists as a remote-tracking branch but was not read; the claims about it in `rev/REVIEW.md:48` and `hand/:124,189-195` are not verified here and belong to another sector.
- `pkt/verify.py` was not executed (it writes `MANIFEST.json`, a tracked file). Recorded passes are cited as recorded.
- `pkt/calculate.py`, `pkt/calculate_original.py`, `pkt/challenge.py`, `pkt/algebra.py`, `pkt/autonomous_writer.py`, `pkt/writer_replay.py`, `pkt/assemble.py` were not read in full (calculate.py only by grep); correctness of the code is taken from the recorded cross-checks, not re-derived.
- The 201-point traces in `results_initial.json`, `results_with_control.json` and `challenge_results.json` were inspected only through their recorded scalar summaries and the recorded verify assertions, not point by point.
- `pkt/core_summary.png` and `pkt/writer_summary.png` were not viewed.
- The task's "top-level REVIEW.md" was taken to be `rev/REVIEW.md`; no other `REVIEW.md` for this sector exists in the repository (a `find` to depth 3 returned none).
- The primary literature in `pkt/SOURCES.md:19-24` (MSS, Shenker–Stanford, quantum East, disorder-free localization, MBL OTOCs, avalanches) was not read; the packet's own read scope is "abstract and publisher identification" for most (`pkt/SOURCES.md:21-24`).
- No heavy computation was run; the only derived numbers are the per-cell t=0 successes computed from stored means (§3e) and a sha256 comparison (§3f).
