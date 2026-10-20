# Campaign ledger entry: S10, M1 memory (key `s10_m1_memory`)

Sector: M1 memory. A finite Ising repetition register with matter feedback and a programmed writer. The facilitation benefit is withdrawn.

Status: this is the corrected entry, written after an adversarial re-read. The earlier entry's numbers were almost all right. The re-read made eight substantive corrections, listed in §0, and added failures and dropped caveats that the earlier entry missed.

Checkout note. The branch is `codex/grut-closure-fork-20261019`. The working tree is now at `1489773` ("WIP checkpoint 7"). That commit is a descendant of the handoff commit `96a7e7e` and of the correction commit `4e9a161`. `git diff --stat` between `4e9a161` and HEAD is empty for both M1 directories, and between `96a7e7e` and HEAD it is empty for the handoff directory, so every citation below holds at all of these commits. (The earlier entry said the tree was at `da8f925`. That has since moved.)

Path abbreviations:
- `pkt/` = `research/reality_framework/grut_memory_candidate_20261019/`
- `rev/` = `research/reality_framework/grut_memory_candidate_20261019_review/`
- `hand` = `research/reality_framework/grut_closed_whole_handoff_20261019/NEXT_CHAT_HANDOFF.md`
- JSON files are cited by key path. Files at older commits are cited as `<commit>:pkt/<file>` with line numbers from `git show`.

---

## 0. What the adversarial re-read changed

1. **The energy flows the other way.** The earlier entry said the register's energy change was "paid by the carrier", and that "energy conservation makes the carrier pay." It took this from `rev/REVIEW.md:42`, which says "taken from the carrier's kinetic energy." The recorded numbers show the opposite. In the writing branch the core energy falls, 1.55 → −0.437569946, and the carrier's kinetic energy rises, 10002.0408163 → 10004.0290991 (`pkt/REPORT.md:116-118`). The carrier's mean momentum also rises, from 100 to 100.00993831588102 (`pkt/writer_results_refined.json: runs[2].minus.clock_mean_p`). So the write releases register energy into the carrier. It is consistent that the r=−1 readout preparation energy is −0.45, against 1.55 for the blank/r=+1 configuration (`pkt/WORKED_REGIME.md:142`).
2. **"Better at every sampled time" overstates the result.** The packet's statement is "no worse at every sampled time through t=40, within 1e-10" (`pkt/REPORT.md:28-29`). Recomputing from the stored traces in `pkt/results_with_control.json` shows three things:
   - The only samples where uniform − candidate is negative are at t=0: −2.7755575615628914e-15 at L=3 and −1.7763568394002505e-15 at L=5.
   - At every one of the 200 samples with t>0, the control is strictly better, by at least 2.575300059093344e-06 (L=3), 2.615226767366252e-06 (L=4) and 2.6156037454816428e-06 (L=5), each at t=0.2.
   - At t=0 the two are equal to rounding.
3. **Who said what about the mechanism.** "Zero-temperature transverse-field Ising repetition code" and "resonant spreading (an East-model cascade)" are the reviewer's words (`rev/REVIEW.md:36`), not the packet's. The packet's own wording is "a prepared ferromagnetic repetition register with weak transverse flips and matter feedback" (`pkt/REPORT.md:38-39`). It says the preparations are "pure prepared states with nonzero excitation energy, not a ground-state preparation" (`pkt/WORKED_REGIME.md:124-125`), with excitation energy 11.0731790514 above the L=4 ground state (`pkt/WORKED_REGIME.md:137-139`). It also says "an exact resonance theorem for the full coupled Hamiltonian is not claimed" (`pkt/WORKED_REGIME.md:120-121`). Neither packet nor review computes a gap or its scaling.
4. **Sign wording of the readout energy.** The nonselective readout *changes* system energy by −1.2010144303 and −0.8997375490. That means 1.2010144303 and 0.8997375490 *leave* the modeled system (`pkt/WORKED_REGIME.md:144-148`). The earlier entry said it "removes −1.20…", which is a double negative.
5. **Misclassification.** Exact conservation of D at ε=0 is an exact result about the supplied constrained law. It belongs under conditional constraint (model-specific theorem), not "recovered known physics."
6. **Dropped caveat now retained.** At `4e9a161` the sentence "Radiative stability is not established." was deleted from open-claims item 1 (`1997d20:pkt/OPEN_CLAIMS.md:26-27`, compare `pkt/OPEN_CLAIMS.md:32-35`). Nothing in the packet resolves it. "Radiative" now appears nowhere in the packet.
7. **Withdrawn original texts quoted, plus the review's narrowed overstatements.** The original texts are `1997d20:pkt/WORKED_REGIME.md:95-97`, `1997d20:pkt/GRUT_CANDIDATE_M1.md:3-8,81` and `1997d20:pkt/OPEN_CLAIMS.md:5,48-53`. The review's overstatements, which the packet narrows, are "every size and every time", "unchanged to four decimals" (false at L=3) and "writes the uniform-field control equally well" (never computed).
8. **Handoff-row checks corrected.** The handoff omits that the control was post-hoc. "Final per-cell" is a terminology mismatch. The "≥.991697" figure holds strictly only for the 4096-point refined run. Two line citations were off by one. The checkout state is updated.

---

## 1. Question asked in this sector

Can one closed, time-independent Hamiltonian with no external bath, drive, disorder or controller do three things?
- (a) Retain a readable, redundantly encoded classical record while the record changes later matter interactions and matter information spreads.
- (b) Write such a record autonomously from internal resources.
- (c) Show that a proposed extra law, "state-dependent kinetic saturation" (facilitation), is what protects the record.

Sources: `pkt/protocol.json:5-8,21`; `pkt/writer_protocol.json:11`; `pkt/REPORT.md:44-47`. Chaos-exponent and gravity claims are explicitly out of scope (`pkt/protocol.json:41-42`; `pkt/writer_protocol.json:29`).

The campaign's stated motivation was constructive: "Try a new explicit construction instead of inferring that the intuition has no possible law" (`pkt/protocol.json:5`). The law was supplied, not derived from closure: "neither is derived from the words 'memory' or 'closure'" (`pkt/GRUT_CANDIDATE_M1.md:49-51`).

Under the parent closed-whole question, the sector therefore tests something narrower than "what closure forces." It asks whether, inside a closed description, a candidate extra law is needed for (or helps) memory, or whether the record persists for a reason known physics already contains.

## 2. Assumptions / model class

- **State space.** A ring of L≥3 cells, each with a matter qubit σ_i and a record qubit τ_i. For the writer, add L²(R,dQ)⊗C²_R for a massive carrier with an internal source bit (`pkt/GRUT_CANDIDATE_M1.md:17-32`).
- **Core law.** H_M1 = Σ_i [ J(1+ατ_i^z)σ_i^zσ_{i+1}^z + h_xσ_i^x + h_zσ_i^z − Kτ_i^zτ_{i+1}^z − ΔF_iτ_i^x ], with F_i = ε + (1−ε)(1−τ_{i−1}^zτ_{i+1}^z)/2 (`pkt/GRUT_CANDIDATE_M1.md:34-43`).
- **Writer completion.** H = P²/(2M) + H_M1 + Ω(1−R^z)/2 Σ_i f_i(Q)τ_i^x, with fixed cos² gates of half-width w (`pkt/GRUT_CANDIDATE_M1.md:53-66`).
- **Benchmark.** ħ=J=1 and (α,h_x,h_z,K,Δ,ε) = (0.5,1.05,0.5,1,1,0.05) (`pkt/GRUT_CANDIDATE_M1.md:93-96`).
- **Writer parameters.** L=4, M=0.5, p₀=100, q₀=−4, σ_Q=0.35, w=0.3, Ω=π(p₀/M)/(2w) (`pkt/GRUT_CANDIDATE_M1.md:99-101`). Numerically Ω=1047.1975511965977 (`pkt/writer_results_initial.json: parameters.Omega`).
- **Supplied postulates.** Quantum mechanics, a fixed local substrate, the kinetic-saturation projector, physical feedback through α, and modeled preparation resources (`pkt/GRUT_CANDIDATE_M1.md:75-89`). The packet calls "the particular facilitation projector, spin content, topology and parameter point" construction choices (`pkt/GRUT_CANDIDATE_M1.md:86-87`).
- **Time and constraints.** External Hamiltonian time is fundamental. There is no gauge or gravitational constraint (`pkt/GRUT_CANDIDATE_M1.md:29-32`).
- **Preparations and grids.** Pure prepared states with nonzero excitation energy, not a ground state or a finite-temperature phase (`pkt/WORKED_REGIME.md:124-125`). Sizes are L=3,4,5 only. Retention is sampled on [0,40] with step 0.2 (201 points). Readout is sampled on [0,2] with step 0.02 (`pkt/protocol.json:32`; `pkt/calculate.py:110-111`).
- **Success observable.** The equal-prior local guessing probability P_A = ½ + |⟨A⟩₊−⟨A⟩₋|/4. This is "not the optimal global information-recovery probability." The record is one classical bit repeated, not L independent bits (`pkt/GRUT_CANDIDATE_M1.md:133-143`).
- **Matched control.** Only F_i is replaced, by the constant ε=0.05. It was specified after the review outcome was known: "This is not a blinded or held-out prediction" (`pkt/control_protocol.json:6-7`; `pkt/REPORT.md:18-20,137-139`).
- **Closure as imposed.**
  - The core has no quenched disorder or external controller (`pkt/protocol.json:8`).
  - The equilibrium diagnostic has no external bath in the dynamical law (`pkt/protocol.json:21`).
  - The writer has "No external pulse, reset, postselection, reservoir or time-dependent physical coupling" (`pkt/writer_protocol.json:11`).
  - "At no point is a failed transmission discarded" (`pkt/COMPATIBILITY.md:14-17`).

## 3. Evidence

### 3a. Core retention and readout at the M1 point (all 13 registered conditions pass)

| Run | min memory success through t=10 | success at t=40 | best matter readout |
|---|---:|---:|---:|
| M1, L=3 | 0.994999 | 0.985006 | 0.806447 (t=1.20) |
| M1, L=4 | 0.987510 | 0.935065 | 0.778468 (t=0.58) |
| M1, L=5 | 0.989132 | 0.940399 | 0.778460 (t=0.58) |

Sources: `pkt/WORKED_REGIME.md:75-79`; `rev/control_results.txt:1,3,5`.

- **L=4 at full precision.** 0.9875095925556739 / 0.9350645853495996 / 0.7784678210302141 (`pkt/summary.json: core[1]`).
- **L=4 readout arithmetic.** P_Y = ½ + |−0.491561535 − 0.622309749|/4 = 0.778467821, with simultaneous memory success 0.999312115 (`pkt/WORKED_REGIME.md:87-95`).
- **Execution.** The seven original runs took 11.79 s and all 13 registered conditions passed (`pkt/EXECUTION.md:11-16`; `pkt/results_initial.json: all_registered_acceptance_pass = true`, `total_wall_seconds = 11.790696546000163`).
- **Thresholds.** Memory ≥0.9 through t=10; readout >0.6 with contemporaneous memory ≥0.9 (`pkt/protocol.json:35-36`).
- **What the pass does and does not show.** "none tested superiority over the matched weak-field baseline. Their passing does not establish the proposed principle" (`pkt/REPORT.md:98-100`).

**Original controls (all at L=4).**
- F=1 "unfacilitated": 0.621452 / 0.709936 / 0.741417.
- Feedback removed (α=0): 0.996241 / 0.994449 / 0.500000.
- Exactly frozen (ε=0): 1.000000 / 1.000000 / 0.778577.

Sources: `pkt/WORKED_REGIME.md:80-82`. Half softness (ε=0.025) gives 0.9968211359193733 / 0.9792849192829733 / 0.7785498324346294 (`pkt/summary.json: core[6]`).

### 3b. The decisive matched control (F_i ≡ ε = 0.05)

| L | M1 min t≤10 | control min t≤10 | M1 t=40 | control t=40 | M1 min [0,40] | control min [0,40] | best readout M1 / control |
|---|---:|---:|---:|---:|---:|---:|---:|
| 3 | 0.994999 | 0.998738 | 0.985006 | 0.994692 | 0.984591 | 0.994266 | 0.806447 / 0.806626 |
| 4 | 0.987510 | 0.999055 | 0.935065 | 0.998458 | 0.935065 | 0.998397 | 0.778468 / 0.778483 |
| 5 | 0.989132 | 0.998479 | 0.940399 | 0.997517 | 0.940399 | 0.995717 | 0.778460 / 0.778475 |

Sources: `pkt/REPORT.md:22-26`; `pkt/WORKED_REGIME.md:83-85`; `rev/control_results.txt:1-6`; `rev/REVIEW.md:25-32`.

- **L=4 at full precision.** Control 0.9990547948737153 / 0.998457849859687 / 0.7784834017608981 at t=0.58, against M1's 0.7784678210302141. The readout difference is 1.5580730684083832e-05 (`pkt/summary.json: matched_comparison.4`).
- **Readout differences at all sizes.** 0.000179931, 0.000015581 and 0.000015652 at L=3,4,5 (`pkt/EXECUTION.md:46`). The control's best readout is slightly *higher* at all three sizes (`pkt/summary.json: matched_comparison.{3,4,5}.readout_best_probability_difference` > 0). So facilitation helps neither retention nor readout here.
- **Packet's statement.** "no worse at every sampled time through t=40, within 1e-10 numerical tolerance. This is not an all-continuous-time theorem" (`pkt/REPORT.md:28-29`).
- **How it is implemented.** `uniform_no_worse_on_all_samples = bool(delta.min() >= -1e-10)` (`pkt/calculate.py:301`). The protocol says this tolerance is for comparisons "near t=0" (`pkt/control_protocol.json:14`).
- **Re-check from stored traces (this re-read).** Negative deltas occur only at t=0, at −2.78e-15 (L=3) and −1.78e-15 (L=5). For t>0 the control is strictly better at all 200 samples, with minimum margins 2.575300059093344e-06, 2.615226767366252e-06 and 2.6156037454816428e-06, each at t=0.2. Source: `pkt/results_with_control.json: runs.{candidate,uniform_floor}_L{3,4,5}.retention.success`.
- **Recorded mechanism claims.** `facilitation_improves_minimum_retention_all_sizes = false` and `matched_uniform_control_retains_better_on_registered_samples = true` (`pkt/summary.json: mechanism_claims`).
- **Independent executions.**
  - The review's independent tensor-product script (`rev/uniform_field_control.py:10-19`) was re-executed with one BLAS thread and reproduces all six rows (`pkt/EXECUTION.md:30-32`).
  - The packet's bit-table `calculate.py` ten-run replay took 7.31415 s (`pkt/EXECUTION.md:34-38`; `pkt/results_with_control.json: total_wall_seconds = 7.314150894999784`). Bit-table and Kronecker matrices agree (`pkt/summary.json: control_method_checks`).
  - Recomputed here from stored traces: the six review rows match within 4.967923176169009e-07, inside verify's 5.01e-7 criterion (`pkt/verify.py:63-71`).
  - Secondary, not primary: a WIP-checkpoint rerun of the unchanged review script reports identical output (`research/reality_framework/grut_closed_whole_theory_20261019/results/m1_matched_control_rerun.txt:1-12`).
- **The review's mechanism reading.** "the memory is the ferromagnetic gap of −Kτ^zτ^z against a small transverse field; that is a zero-temperature transverse-field Ising repetition code." Facilitation is inactive on an ordered record and, after a leak, "converts an off-resonant single-flip leakage into resonant spreading (an East-model cascade)" (`rev/REVIEW.md:36`).
- **What the packet adopts from that reading.** It adopts the qualitative mechanism with qualifications:
  - Opposite-neighbour flips at amplitude 1 "move domain walls without changing the ferromagnetic energy, rather than annihilating the error. Matter can still change the transition energy; an exact resonance theorem for the full coupled Hamiltonian is not claimed. This reversible facilitation does not supply an error-correction bias" (`pkt/WORKED_REGIME.md:116-122`).
  - It does not adopt "zero-temperature." The states are excited prepared states (`pkt/WORKED_REGIME.md:124-125,137-139`).
- **Review overstatements the packet narrowed.** Each is recorded here as a correction.
  - The review says the control is better "at every size and every time," and the `cc24967` commit message says "identical matter readout." The packet narrows this to "no worse at every sampled time … not an all-continuous-time theorem" and "differ by 1.56e-5, rather than being exactly identical" (`pkt/REPORT.md:28-32`).
  - The review says readout is "unchanged to four decimals" (`rev/REVIEW.md:34`). This is false at L=3: 0.806447 vs 0.806626, a difference of 0.000179931 (`pkt/EXECUTION.md:46`).

### 3c. Exact regime at ε = 0 (an exact result about the supplied law)

- **Conservation.** With D = ½Σ_i(1−τ_i^zτ_{i+1}^z), [H_M1(0), D] = 0, and the computed commutator is exactly zero at L=3,4,5 (`pkt/WORKED_REGIME.md:13-22`; `pkt/algebra_results.json: runs[*].hard_constraint_commutator_max = 0.0`).
- **Perfect retention in the zero-wall sector.** The record is perfectly retained while the matter evolves under the record-dependent H_σ^{(r)} (`pkt/WORKED_REGIME.md:24-39`).
- **Formation barrier.** The zero-wall probability starts at 2^{1−L} (0.25, 0.125, 0.0625). At ε=0, unitary evolution cannot increase it: "Perfect protection is not a mechanism for generating the protected state" (`pkt/WORKED_REGIME.md:41-44`; `pkt/algebra_results.json`).
- **At the soft point (ε=0.05).** The ordered-sector probability peaks at 0.3289007412515012 / 0.3022736769655148 / 0.08277297473471425 and ends at 0.2614344170727129 / 0.29884391381968967 / 0.04729380180227547. These are "reported without claiming a persistent organized state" (`pkt/algebra_results.json: runs[*]`; `pkt/EXECUTION.md:23-26`).
- **Short-time checks reproduced by matrix calculation.**
  - ⟨τ_0^z(t)⟩_r = r[1 − 2Δ²ε²t² + O(t³)].
  - Readout slopes 7 and 3.
  - Squared-commutator coefficients 18 and 2.

  Sources: `pkt/WORKED_REGIME.md:47-71`; `pkt/protocol.json:44-48`.

### 3d. Registered stronger attacks that failed (six conditions false, retained)

The challenge protocol was registered after the initial results and before the attacks (`pkt/challenge_protocol.json:2`).

- **Single-cell bit (standard: ≥0.9 through t=10).**
  - At t=10, 0.801481 (L=4) and 0.751472 (L=5) (`pkt/WORKED_REGIME.md:155-157`).
  - The minima through t=10 are lower still, 0.6571024197046156 (L=4) and 0.6760116942976508 (L=5) (`pkt/challenge_results.json: local_bit.{4,5}.minimum_t10`).
- **Bare matter-to-memory write.**
  - Maxima 0.553352 (L=4) and 0.526376 (L=5), reached only at t=37.4 and 37.8. There are no witness times, and the minimum through t=10 is 0.5 (`pkt/WORKED_REGIME.md:158-160`; `pkt/challenge_results.json: autonomous_write.{4,5}`).
  - The packet calls this "close to chance" (`pkt/OPEN_CLAIMS.md:22-23`).
- **Unconstrained flip −κΣτ_i^x.** Minimum success through t=10 is 0.875748 at κ=0.1 and 0.764940 at κ=0.3 (`pkt/WORKED_REGIME.md:161-165`).
- **Tally.** Six conditions false, three numerical-method checks true. The checks include `independent_full_H_max_error = 4.163336342344337e-17`, and the challenge run took 7.36 s (`pkt/summary.json: adversarial_conditions`; `pkt/challenge_results.json`; `pkt/EXECUTION.md:17-22`).

### 3e. Writer (M1-W)

- **Mechanism.** The gate crossing time is 2w/(p₀/M) = 0.003. The conditional pulse area is Ωw/(p₀/M) = π/2, "giving a controlled flip. This leading operation does not require facilitation" (`pkt/COMPATIBILITY.md:34-37`; `rev/REVIEW.md:42`).
- **Only one source branch interacts.** In the R^z=+1 branch the coupling (1−R^z)/2 vanishes. Final interaction energy is 0.0 and drift is at the 1e-11 level (`pkt/writer_results_refined.json: runs[*].plus`). There the blank all-+z register *is* the record (`pkt/GRUT_CANDIDATE_M1.md:112-113`). The nontrivial operation is the R^z=−1 flip of all four cells.
- **Per-cell output immediately after scattering (t=0.06).** "At least 0.999680" (`pkt/REPORT.md:91`; `pkt/writer_protocol.json:23`). Recorded per-run minimum cells are 0.9996801522490156 (dt 0.0002, 2048 points), 0.9996800894784178 (dt 0.0001, 2048) and 0.9996802603994459 (dt 0.0001, 4096) (`pkt/writer_results_refined.json: runs[*].initial_all_cell_success`).
- **All-cell success through sampled storage t=10.** "At least 0.991697" (`pkt/REPORT.md:92`). The value is 0.9916972267188419 in the 4096-point run, which the summary uses (`pkt/summary.json: writer_minimum_memory_t10`; `pkt/verify.py:74`). The other two runs give 0.9916969879289492 (coarse) and 0.9916969501876505 (dt 0.0001, 2048 points). Both round to 0.991697 but are strictly below it (`pkt/writer_results_refined.json: runs[0..1].minimum_all_cell_success_through_t10`).
- **Storage approximation.** Storage "neglects residual carrier/gate interaction" (`pkt/REPORT.md:92,110-112`; `pkt/EXECUTION.md:88-94`).
- **Matter readout.** 0.769977 at storage t=0.52, with memory 0.998568 (`pkt/REPORT.md:93`; `runs[2].readout`: 0.769977434073726, 0.9985675362234426).
- **Refinement.** Maximum probability change is below 6.8e-7 (`pkt/REPORT.md:94`; `pkt/summary.json: writer_refinement.overall.max_memory_probability_difference = 6.705103664517154e-07`). This "Does not bound every physical approximation."
- **Energy.** In the writing branch the carrier goes 10002.0408163 → 10004.0290991 and the core goes 1.55 → −0.437569946 (`pkt/REPORT.md:116-118`). Energy therefore flows *into* the carrier.
  - Total-energy drift is 0.000712816 J (model energy units), "not zero; its cutoff behavior is nonmonotonic. The probability refinement is much stronger than the energy-convergence evidence" (`pkt/REPORT.md:119-121`).
  - Per-run writing-branch drift is 0.00010979246144415811, 1.0308212949894369e-08 and 0.0007128160505089909. The finest spatial grid has the largest drift (`runs[*].minus.total_energy_drift`).
- **Diagnostics.**
  - Norm error 1.33e-12, residual gate-region probability 6.87e-10, edge probability 1.76e-9, negative momentum 1.36e-8 (`pkt/REPORT.md:123-128`; `runs[2].minus`).
  - These do "not prove that a later carrier return is impossible."
- **Execution history.**
  - The coarse run took 414.59 s.
  - The slower dt=0.0001 local-Pauli replay "was interrupted with Ctrl-C before completion … not an erased physical failure or a completed refinement" (`pkt/EXECUTION.md:64-70`).
  - `pkt/execution_amendment.json` records replacing it with `writer_replay.py`, which took 314.89 s (`pkt/writer_results_refined.json: wall_seconds = 314.8902410519995`).
- **Status of the writer results.**
  - All six writer conditions pass (`pkt/summary.json: writer_checks`).
  - The packet calls these "internal numerical results, not independent validation" (`pkt/REPORT.md:103-106`).
  - No writer scattering was rerun or transferred to the uniform-floor core (`pkt/EXECUTION.md:32`; `pkt/control_protocol.json:16`; `pkt/COMPATIBILITY.md:36-38`). The review's line that "the same particle writes the uniform-field control equally well" (`rev/REVIEW.md:44`) was never computed.

### 3f. Numerical integrity, replay and finalizer

- **Core accuracy.** Norm deviation is at most 1.29e-14 and energy deviation is below 3.07e-14 J. Kronecker/Krylov agreement is 2.24e-13 on the retention trace and within 5.2e-15 on readout means (`pkt/WORKED_REGIME.md:129-133`). The recorded maximum readout difference is 5.10702591327572e-15 (`pkt/challenge_results.json: measurement[*].spectral_mean_difference`).
- **Replay.** `original_replay_max_conditional_mean_error = 0.0` (`pkt/MANIFEST.json`), asserted to be below 1e-10 (`pkt/verify.py:46-52`). Recomputed here from the stored traces, the result is 0.0.
- **What the finalizer does.** `verify.py` asserts source and protocol identities, recorded outcomes, syntax and links. It then *regenerates* `MANIFEST.json` with hashes of the 31 other files (`pkt/verify.py:23-95`). It does not compare those files against an earlier manifest.
- **Hash check done here.** A read-only sha256-and-size comparison against the stored manifest matches 31/31 at HEAD, with no unlisted files. `verify.py` was not executed because it rewrites a tracked file. Its pass is recorded at `AGENT_COORDINATION.md:5822`.

### 3g. Equilibrium and scrambling diagnostics (no exponent)

- **Thermal memory diagnostic.** At β=1 the connected memory autocorrelation is 0.538351 at t=10, with a finite-size plateau of 0.458812 (`pkt/WORKED_REGIME.md:177-180`; `pkt/results_initial.json: runs.candidate_L4.equilibrium.memory`).
- **Matter spreading.** The adjacent matter squared commutator reaches 1.41102 in [0,4]. The preparations and observables are kept separate from the prepared-record experiment (`pkt/WORKED_REGIME.md:182-185`).
- **No exponent.** The centered OTOC still has an initial factorization defect of 0.607085, and "No λ_L/(2πT) is reported" (`pkt/WORKED_REGIME.md:187-190`; `pkt/summary.json: lyapunov_exponent = null, gravitational_prediction = null`).

### 3h. Headline history (1997d20 → cc24967 → 4e9a161)

- **Original `1997d20`.**
  - Title "a constructed and calculated memory candidate"; "The added principle is **kinetic saturation**" (`1997d20:pkt/REPORT.md:1,13-16`).
  - "This candidate proposes state-dependent kinetic saturation … a complete nonrelativistic model of that mechanism" (`1997d20:pkt/GRUT_CANDIDATE_M1.md:3-8`).
  - Postulate row: "Protects prepared ordered patterns without freezing all transitions" (`1997d20:pkt/GRUT_CANDIDATE_M1.md:81`).
  - "Removing saturation leaves an interacting readout but substantially erases the record. Both parts of the proposed mechanism have a tested operational role" (`1997d20:pkt/WORKED_REGIME.md:95-97`).
  - Classification "a narrower working model rather than the intended universal theory … the demonstrated finite memory mechanism" (`1997d20:pkt/OPEN_CLAIMS.md:48-53`).
  - Open-claims item 1 included "Radiative stability is not established" (`1997d20:pkt/OPEN_CLAIMS.md:26-27`).
- **Review `cc24967`.** Adds the F≡ε control and reclassifies M1 as "a correct finite model of a repetition code with matter feedback and a flying-qubit writer; the claimed principle is not supported by its own numbers" (`rev/REVIEW.md:49`).
- **Correction `4e9a161`.** "a finite Ising repetition register with matter feedback and a programmed flying-qubit writer; the proposed facilitation benefit is rejected in the demonstrated regime" (`pkt/OPEN_CLAIMS.md:56-61`). The withdrawal is logged at `AGENT_COORDINATION.md:5820`. "Stable, computable regime" was narrowed to "computable finite-time regime" (`pkt/OPEN_CLAIMS.md:5`). The radiative-stability sentence was deleted without being resolved.

## 4. Answer

**What works, in the stated finite regime.**
- A prepared ferromagnetic repetition register (L=3–5, one classical bit repeated) retains the bit at ≥0.987510 through sampled t=10 and at 0.935065 at t=40 (L=4).
- A τ-dependent matter bond gives a 0.778468 local matter readout (`pkt/REPORT.md:84-90`).
- At ε=0 the zero-wall sector is exactly retained, by domain-wall conservation in the supplied law (`pkt/WORKED_REGIME.md:13-39`).
- A closed, time-independent flying-qubit writer conditionally flips all four cells from a source bit, with every cell ≥0.999680 just after scattering (`pkt/REPORT.md:91`). In the writing branch, energy released by the register (core 1.55 → −0.437569946) goes into the carrier (10002.0408163 → 10004.0290991) (`pkt/REPORT.md:116-118`).
- Subsequent storage is ≥0.991697 under the stated decoupling approximation (`pkt/REPORT.md:92,110-112`).

**What was withdrawn.** The attribution of retention to the facilitation (kinetic-saturation) law.
- The matched F≡ε control is no worse at all 201 sampled times (equal at t=0, strictly better at every t>0 sample). It has better minima and better t=40 values at every size (L=4, t=40: 0.998458 vs 0.935065), and a slightly higher best readout (0.778483 vs 0.778468 at L=4).
- Hence "The demonstrated retention therefore does not support the added facilitation principle. The finite calculation remains correct; the attribution of its protection to that principle is withdrawn" (`pkt/REPORT.md:3-8`; `pkt/OPEN_CLAIMS.md:16-19`).
- Also withdrawn: the 1997d20 inference that "Both parts of the proposed mechanism have a tested operational role" (`1997d20:pkt/WORKED_REGIME.md:95-97`, replaced by `pkt/WORKED_REGIME.md:103-105`).
- The packet's replacement explanation is "a prepared ferromagnetic repetition register with weak transverse flips and matter feedback" (`pkt/REPORT.md:38-39`). The reviewer's sharper wording ("zero-temperature … code", "resonant spreading", `rev/REVIEW.md:36`) is only partly adopted. The packet keeps the states as excited prepared states and claims no exact resonance theorem (`pkt/WORKED_REGIME.md:120-125`).

**What is not established.**
- **A finite-temperature memory phase** (`pkt/REPORT.md:39-40`; `pkt/COMPATIBILITY.md:41-42`). The reviewer's Peierls remark at `rev/REVIEW.md:48` is not a packet computation.
- **Cloning of an unknown quantum state.** "the model does not clone an arbitrary unknown quantum state" (`pkt/GRUT_CANDIDATE_M1.md:143-144`).
- **Gravity** (`pkt/COMPATIBILITY.md:8`; `pkt/REPORT.md:155-167`).
- **A chaos exponent** (`pkt/REPORT.md:159-161`; `pkt/summary.json: lyapunov_exponent = null`).
- **All-times retention, or a scattering certificate** (`pkt/REPORT.md:29,111-112,126-128`).
- **Locality, scaling and robustness.** Independent local capacity, thermodynamic localization and robustness to thermal inclusions (`pkt/OPEN_CLAIMS.md:38-40`).
- **Generic autonomous formation** (`pkt/OPEN_CLAIMS.md:36-37`).
- **Radiative stability of the facilitation law.** Stated as unestablished at 1997d20 and silently dropped (`1997d20:pkt/OPEN_CLAIMS.md:26-27`).
- **Novelty or empirical viability** (`pkt/OPEN_CLAIMS.md:50-52`).

## 5. Closure verdict

Nothing beyond known physics was found to be needed in the demonstrated finite class: L=3–5, the benchmark point, prepared pure states, sampled grids. The deciding computation is the matched-control table (§3b; `pkt/REPORT.md:22-26`):

| L | min through t=10, M1 vs control | t=40, M1 vs control |
|---|---|---|
| 3 | 0.994999 vs 0.998738 | 0.985006 vs 0.994692 |
| 4 | 0.987510 vs 0.999055 | 0.935065 vs 0.998458 |
| 5 | 0.989132 vs 0.998479 | 0.940399 vs 0.997517 |

**Scope of the verdict.** M1 did not derive any law from closure. The law, ε and feedback were supplied (`pkt/GRUT_CANDIDATE_M1.md:49-51`; `pkt/protocol.json:5`). Closure entered as an imposed constraint: no external bath, reservoir, pulse, postselection or discarding (`pkt/protocol.json:8,21`; `pkt/writer_protocol.json:11`; `pkt/COMPATIBILITY.md:14-17`). Within that closed description, the one proposed extra law is not needed and makes retention worse. The record persists in an ordinary weak-transverse-field Ising repetition register (`pkt/REPORT.md:38-39`). This rejection is itself conditional: "in this preparation and parameter regime" (`pkt/PREDICTIONS_AND_TESTS.md:55-56`). A cascade regime where facilitation might help is untested (`pkt/COMPATIBILITY.md:41-44`).

**What closed bookkeeping did supply.** All of it is standard quantum mechanics.
- **Writer energy.** It is exchanged within the closed system. In the writing branch the register releases energy to the carrier (§3e; `pkt/REPORT.md:116-118`).
- **Readout energy.** The nonselective readout changes system energy by −1.2010144303 and −0.8997375490, so that energy leaves the modeled system. "An energy-conserving detector completion must receive that energy. Detector reset is not modeled" (`pkt/WORKED_REGIME.md:144-148`).
- **A protect/form trade-off.** At ε=0 a conserved D forbids growth of the protected sector (`pkt/WORKED_REGIME.md:41-44`).

**Gravity.** The packet's algebraic remark: for H=Σ_rP_r⊗H_r, "Structured memory therefore does not by itself bound the largest exponent of the other operators." A memory sector that changes a coefficient without removing the eikonal G_N s, s ∝ e^{2πTt} channel "does not lower its exponent" (`pkt/COMPATIBILITY.md:62-77`).

**Result.** A finite working construction with the proposed addition rejected by its own numbers. "This is not a newly established memory principle, error-correction mechanism or gravitational law" (`pkt/REPORT.md:174-176,184-185`).

## 6. Classification: mixed

- **Recovered known physics.**
  - Retention of a prepared repetition bit in a weak-transverse-field Ising register. This is a qualitative reading, supported by the matched control; no gap or scaling calculation was made (`pkt/REPORT.md:38-39`; `rev/REVIEW.md:36`).
  - The π/2 pulse-area controlled flip (`pkt/COMPATIBILITY.md:34-37`; `rev/REVIEW.md:42`).
  - Unitary energy exchange in closed scattering (`pkt/COMPATIBILITY.md:14-17`).
- **Rejected mechanism.**
  - The facilitation-benefit interpretation of kinetic saturation, rejected by the F≡ε control at L=3,4,5 (`pkt/REPORT.md:137-139`; `pkt/summary.json: mechanism_claims`).
  - The original "both parts have a tested operational role" inference (`1997d20:pkt/WORKED_REGIME.md:95-97`).
- **Supplied model choice.** Ring topology, two qubits per cell, the projector form of F_i, ε, α, the parameter point, the writer apparatus parameters and gate shapes, and the local guessing observable (`pkt/GRUT_CANDIDATE_M1.md:75-89,93-101,140-142`; `pkt/PREDICTIONS_AND_TESTS.md:81-83`).
- **Conditional constraint.**
  - All probabilities hold only at L=3–5, on sampled grids, for the prepared pure states at the benchmark (`pkt/WORKED_REGIME.md:95-97`; `pkt/PREDICTIONS_AND_TESTS.md:72-79`).
  - The exact [H_M1(0),D]=0 conservation and the zero-wall retention theorem hold for the supplied law at ε=0 only (`pkt/WORKED_REGIME.md:13-44`).
  - Writer storage numbers hold only under post-scattering decoupling (`pkt/EXECUTION.md:88-94`).
  - The rejection of facilitation is itself regime-specific (`pkt/PREDICTIONS_AND_TESTS.md:55-56`).
- **Conjecture or unspecified successor.**
  - An error-correction-by-cascade regime where facilitation beats F≡ε (`rev/REVIEW.md:54`; `pkt/PREDICTIONS_AND_TESTS.md:85-88`).
  - A relation to Candidate 1 as a low-energy sector (`pkt/COMPATIBILITY.md:87-92`).
- **Unexecuted task.**
  - Writer scattering on the uniform-floor core (`pkt/control_protocol.json:16`; `pkt/EXECUTION.md:32`).
  - A cascade calculation (`pkt/COMPATIBILITY.md:41-44`).
  - Thermal-inclusion and robustness attacks, and sizes beyond L=5 (`pkt/OPEN_CLAIMS.md:38-40`).
  - A dynamical terminal detector, its reset cost, and rigorous long-time writer bounds (`pkt/OPEN_CLAIMS.md:48-49`).
  - The interrupted dt=0.0001 local-Pauli replay, which was superseded rather than completed (`pkt/EXECUTION.md:67-70`).
  - Radiative stability of the facilitation law (`1997d20:pkt/OPEN_CLAIMS.md:26-27`).
  - The registered no-boundary B control, still unexecuted (`pkt/COMPATIBILITY.md:94-97`).

## 7. Remaining freedom

- **Parameters.** The whole parameter point (ε, α, K, Δ, h_x, h_z). There is no calibration to seconds or to any physical system (`pkt/WORKED_REGIME.md:6-7`; `pkt/PREDICTIONS_AND_TESTS.md:8-14,81-83`).
- **Structure.** Topology, dimensionality, size beyond L=5, the record encoding, and local versus optimal readout (`pkt/GRUT_CANDIDATE_M1.md:80,139-142`).
- **Whether to keep facilitation at all.** The original law is preserved "rather than silently replaced by the better control" (`pkt/GRUT_CANDIDATE_M1.md:151-152`).
- **Cascade successor.** It would need "noise, an error syndrome, correction dynamics, entropy/energy resources and an advantage over this matched baseline" (`pkt/COMPATIBILITY.md:41-44`).
- **Gravitational successor.** It must "specify which horizon/scattering assumption its new law changes, calculate a regulated gravitational observable, and recover the required low-energy gravitational tests" (`pkt/COMPATIBILITY.md:79-83`).
- **Writer inputs.** Writer parameters and gate shapes are apparatus inputs (`pkt/PREDICTIONS_AND_TESTS.md:12-14`). Finite-clock errors for a uniform-floor writer remain uncalculated (`pkt/control_protocol.json:16`).
- **Rejection rule for a calibrated instance.** A disagreement rejects "that calibrated model. This does not establish that every idea called GRUT is false" (`pkt/PREDICTIONS_AND_TESTS.md:45-49`).

## 8. Failures retained

1. The facilitation-benefit attribution was withdrawn after the matched control (`pkt/REPORT.md:3-8,137-139`; `pkt/OPEN_CLAIMS.md:16-19`; `AGENT_COORDINATION.md:5820`).
2. The original inference "Removing saturation … substantially erases the record. Both parts of the proposed mechanism have a tested operational role" was withdrawn (`1997d20:pkt/WORKED_REGIME.md:95-97` → `pkt/WORKED_REGIME.md:103-105`). So was the postulate claim "Protects prepared ordered patterns" (`1997d20:pkt/GRUT_CANDIDATE_M1.md:81`).
3. The F=1 control was confounded by a twentyfold field and "cannot isolate the benefit of neighbor dependence" (`pkt/REPORT.md:32-33`; `rev/REVIEW.md:21`).
4. Single-cell retention fails the 0.9 standard: 0.801481 / 0.751472 at t=10, with minima 0.6571024197046156 / 0.6760116942976508 through t=10 (`pkt/WORKED_REGIME.md:155-157`; `pkt/challenge_results.json: local_bit`).
5. The bare autonomous write fails: maxima 0.553352 / 0.526376, reached only at t≈37–38 (`pkt/WORKED_REGIME.md:158-160`; `pkt/challenge_results.json: autonomous_write`).
6. Unconstrained-flip perturbations break the 0.9 condition: 0.875748 (κ=0.1) and 0.764940 (κ=0.3) (`pkt/WORKED_REGIME.md:161-165`).
7. At ε=0 the ordered sector cannot gain weight, which is a formation barrier (`pkt/OPEN_CLAIMS.md:25-28`).
8. The dt=0.0001 local-Pauli writer replay was interrupted with Ctrl-C and replaced, not completed (`pkt/EXECUTION.md:67-70`).
9. The writer's energy drift, 0.000712816 J, is nonzero and nonmonotonic, and is largest at the finest grid. "The probability refinement is much stronger than the energy-convergence evidence" (`pkt/REPORT.md:119-121`).
10. There is no all-times scattering certificate, and a later carrier return is not excluded (`pkt/REPORT.md:126-128`).
11. The matched control is post-hoc: "not a blind prediction" (`pkt/REPORT.md:138-139`; `pkt/control_protocol.json:6`).
12. Passing the 13 core conditions "does not establish the proposed principle." The explicit writer "does not retroactively pass the failed bare-writing protocol" (`pkt/REPORT.md:98-100,144-145`).
13. The radiative-stability caveat was deleted at 4e9a161 without being resolved (`1997d20:pkt/OPEN_CLAIMS.md:26-27`).
14. The review overstated three things, and the packet narrowed or contradicted each:
    - "every size and every time" (`rev/REVIEW.md:34`).
    - "unchanged to four decimals", false at L=3 (`pkt/EXECUTION.md:46`).
    - "writes the uniform-field control equally well", never computed (`rev/REVIEW.md:44` vs `pkt/COMPATIBILITY.md:36-38`).

    The review also has the energy direction reversed: "taken from the carrier's kinetic energy" (`rev/REVIEW.md:42` vs `pkt/REPORT.md:116-118`).
15. No finite-temperature memory phase, no Lyapunov exponent, no gravity and no priority claim (`pkt/OPEN_CLAIMS.md:30-52`).

## 9. What must not be promoted

- **M1 retention numbers.** Do not use them as evidence for a memory principle. The control does better (`pkt/REPORT.md:98-100`; `pkt/PREDICTIONS_AND_TESTS.md:50-57`).
- **The rejection of facilitation.** Do not generalize it beyond "this preparation and parameter regime" (`pkt/PREDICTIONS_AND_TESTS.md:55-56`).
- **"No worse at 201 sampled times".** Do not read it as a continuous-time or thermodynamic theorem (`pkt/REPORT.md:29`; `pkt/WORKED_REGIME.md:111-114`).
- **≥0.991697 storage.** It is not a complete-Hamiltonian result (`pkt/REPORT.md:110-112`).
- **Integrity passes.** The 13 passes and the 31-file integrity pass are not physical validation: "an integrity pass does not turn them into successful physics" (`pkt/README.md:62-64`; `pkt/PREDICTIONS_AND_TESTS.md:78-79`).
- **Writer probabilities.** Do not transfer them to the uniform-floor core (`pkt/COMPATIBILITY.md:36-38`).
- **"Memory and spreading coexist".** They coexist under one law but in separate preparations. The packet keeps "the preparations and observables … separate" (`pkt/WORKED_REGIME.md:182-185`).
- **Gravity and chaos.** The bond feedback is not gravity. The 1.41102 squared commutator is not a Lyapunov exponent, and λ_L<2πT is not a prediction (`pkt/COMPATIBILITY.md:8`; `pkt/WORKED_REGIME.md:182-190`; `pkt/PREDICTIONS_AND_TESTS.md:65-70`).
- **"Zero-temperature".** Do not attribute the reviewer's "zero-temperature" label to the packet. The states are excited prepared states (`pkt/WORKED_REGIME.md:124-125,137-139`).
- **Toric-branch claims.** Do not import them from `rev/REVIEW.md:48` as M1 results (`pkt/SOURCES.md:13-15`).
- **ε=0 exact retention.** It is not autonomous formation (`pkt/WORKED_REGIME.md:126-127`).
- **The handoff.** Do not use it as a primary source for any number.

## 10. Sources read

All primary sources were read at HEAD `1489773`, where they are identical to `4e9a161` and `96a7e7e`.

- **Packet documents.** `pkt/REPORT.md` 1-185; `pkt/OPEN_CLAIMS.md` 1-61; `pkt/COMPATIBILITY.md` 1-97; `pkt/WORKED_REGIME.md` 1-190; `pkt/PREDICTIONS_AND_TESTS.md` 1-88; `pkt/EXECUTION.md` 1-99; `pkt/SOURCES.md` 1-30; `pkt/GRUT_CANDIDATE_M1.md` 1-167; `pkt/README.md` 1-64.
- **Protocols and scripts.** `pkt/protocol.json`, `pkt/control_protocol.json`, `pkt/challenge_protocol.json`, `pkt/writer_protocol.json`, `pkt/execution_amendment.json` and `pkt/verify.py`, all in full. `pkt/calculate.py` by grep (lines 105-115, 295-312).
- **Results files (scalars and traces).** `pkt/summary.json`, `pkt/MANIFEST.json`, `pkt/results_initial.json`, `pkt/results_with_control.json`, `pkt/writer_results_initial.json`, `pkt/writer_results_refined.json`, `pkt/challenge_results.json`, `pkt/algebra_results.json`.
- **Review.** `rev/REVIEW.md` 1-56; `rev/control_results.txt` 1-6; `rev/uniform_field_control.py` 1-42.
- **Commits.** `git show --stat` of `1997d20`, `cc24967` (with its message) and `4e9a161`. `git diff 1997d20 4e9a161` (removed lines). `1997d20:pkt/REPORT.md` 1-50; `1997d20:pkt/OPEN_CLAIMS.md` 1-53; `1997d20:pkt/WORKED_REGIME.md` 70-100; `1997d20:pkt/GRUT_CANDIDATE_M1.md` 1-10, 72-82.
- **Coordination log.** `AGENT_COORDINATION.md` 5804-5824.
- **Handoff.** `hand` 20-40, 115-170, 180-200, 225-251.
- **Secondary.** `research/reality_framework/grut_closed_whole_theory_20261019/results/m1_matched_control_rerun.txt` 1-12.

## 11. Handoff-row discrepancies

These check the row (`hand:123`) and the numbers block (`hand:138-162`) against the primary sources.

1. **Numbers.** No numerical contradiction was found.
   - All twelve table entries match (`hand:141-145`).
   - The 201-sample / 1e-10 statement matches (`hand:147-148`).
   - The L=4 readouts .778483402 vs .778467821 at t=.58 match (`hand:149`).
   - The six false attack conditions match (`hand:151-152`).
   - The exact seven-trace replay, max error 0.0, matches (`hand:152-153`).
   - ≥.999680 per cell, ≥.991697 storage and .000712816 J drift match (`hand:156-160`).

   The earlier entry's line citations hand:146-147 and hand:151-152 were each off by one.
2. **Omission: the control's post-hoc status.** The handoff never says the matched control was specified after the review outcome was known ("not a blinded or held-out prediction", `pkt/control_protocol.json:6`; `pkt/REPORT.md:138-139`). It does give a post-hoc qualifier for the toric row (`hand:124`).
3. **Omission: the writer on the control core.** The row's "works in the stated finite regime" does not say that writer scattering was never recomputed for the better-retaining uniform-floor core (`pkt/control_protocol.json:16`; `pkt/COMPATIBILITY.md:36-38`).
4. **Terminology: "final per-cell recovery".** The handoff says "final per-cell recovery is at least .999680" (`hand:156`), wording that follows `AGENT_COORDINATION.md:5812`. The packet calls this "initial output in all four cells … at t=0.06", taken immediately after scattering (`pkt/REPORT.md:91`; `pkt/writer_protocol.json:23`). The later storage minimum is lower (.991697).
5. **Precision: "original M1-W writer … at least .991697".** This holds strictly only for the refined 4096-point run (0.9916972267188419). The original coarse run gives 0.9916969879289492, and the dt-only refinement gives 0.9916969501876505. Both round to .991697 (`pkt/writer_results_refined.json: runs[*]`). The .000712816 drift is likewise from the 4096-point run. The coarse writing-branch drift is 0.00010979246144415811.
6. **Consistent: 201 sampled times.** "No worse at all 201 sampled times … within 1e-10" (`hand:147-148`) is accurate. Only t=0 needs the tolerance, and for t>0 the control is strictly better at every sample (§3b). The earlier entry's "the handoff drops the 'near t=0' qualifier" is a nuance, not an error: the packet's own `pkt/REPORT.md:28-29` uses the same wording.
7. **Integrity wording.** "The corrected finalizer checks 31 files plus its manifest; all recorded integrity checks pass" (`hand:153-154`) is a recorded outcome (`AGENT_COORDINATION.md:5822`). `verify.py` regenerates the manifest rather than comparing against a prior one (`pkt/verify.py:84-95`). The read-only 31/31 hash match done here shows the files are unchanged since that run.
8. **Checkout state.** The task names `96a7e7e`. The tree is at `1489773`, a descendant. The handoff says only "a descendant of 4e9a161" (`hand:33-35`). No citation is affected.

## 12. Unread or missing

- **Toric branch.** `claude/toric4d-memory-vs-scrambling-20261019` (`8e0f89f`) was not read. The claims about it at `rev/REVIEW.md:48` and `hand:124,189-195` belong to another sector.
- **`pkt/verify.py`.** Not executed, because it rewrites the tracked `MANIFEST.json`.
- **Code not read in full.** `pkt/calculate.py` (grep only), `calculate_original.py`, `challenge.py`, `algebra.py`, `autonomous_writer.py`, `writer_replay.py` and `assemble.py` were not read in full. Code correctness is taken from the recorded cross-checks.
- **Figures.** `pkt/core_summary.png` and `pkt/writer_summary.png` were not viewed.
- **"Top-level REVIEW.md".** Taken to be `rev/REVIEW.md`. No REVIEW.md exists at the repository root, `research/` or `research/reality_framework/`. The review branch `origin/claude/m1-review-20261019` has tip `cc24967`, identical to the review directory.
- **External literature.** The papers in `pkt/SOURCES.md:19-24` were not read. The packet's own read scope is mostly "Abstract and publisher identification."
- **Computation.** No heavy computation was run. The only recomputations were from stored traces: the delta signs, the replay error, the review-row deviation, and a sha256 comparison.
