# Open claims, exact boundaries, retained failures

Packet `grut_closed_whole_theory_20261019`. This file is the claim boundary of the working
theory. Anything not listed under "Established" is not claimed.

## Established (and by what)

1. **The framework.** A closed whole is the theory datum Θ of FORMAL_CORE DEF 1.7; its rules
   are the standard theorems of FORMAL_CORE §2. Established by citation to standard results and
   by the two calculations. No novelty is claimed for any rule.
2. **Closure forces only the global theorems.** Unitary invariants, the exact reduced equation,
   the Nakajima–Zwanzig structure, CPTP reduced maps for product preparations, nested and
   necessary overlap compatibility, entropy inequalities, quasi-recurrence. Checked numerically
   in calculation 1 (S1, S8) in one supplied eight-qubit model.
3. **Closure forces nothing local.** Fourteen "closure forces X" claims are each refuted by a
   legal closed model or a standard theorem (FORMAL_CORE §4); ten of them are exhibited with
   numbers in calculation 1 (S2–S7).
4. **Compatibility has a computable edge.** Bell-pair chain infeasible; max F_AB + F_BC = 3/2;
   Werner chain p ≤ 2/3; Werner triangle p ≤ 1/3 (calculation 2). Standard quantum marginal
   problem; the numbers are solver outputs within 10⁻⁸ of the analytic values.
5. **The campaign ledger.** Twelve sectors re-read as campaigns; in every sector the recorded
   answer to the parent question is "nothing additional was forced in the stated class", with
   the sector's own computations as the demonstration (CAMPAIGN_LEDGER.md).
6. **Reproduction.** Both calculations reproduce bit-exactly on re-execution and agree with an
   independent re-implementation from the specification; the M1 matched-control script
   reproduces its recorded table line for line.

## Conditional (holds in a stated class only)

- Every sector result marked *conditional constraint* in the ledger: GRUT III D11/D13 slow-roll
  thermality; Candidate 1's one-loop free-carrier binding at its benchmark; the CP-12
  min-entropy obligation; the B1 delivery window; the Level-0 in-class non-implications; the
  influence-hierarchy cumulant ladder on constructed pairs; the GRUT I–II boundaries at their
  stated scopes. Each is conditional on its supplied Θ.
- The finite-model statements of calculation 1 are statements about that model; they illustrate
  theorems that hold generally, but the numbers are the model's.

## Not established, and not claimed

- No additional law, state-selection rule, access rule, clock rule or measurement rule beyond
  standard quantum theory. None survived any campaign (WORKING_THEORY §6).
- No combined "GRUT law" uniting Candidate 1, M1, GRUT III or the Level-0 floor. They live in
  different state spaces; no common action, state space, constraints or calibration exists.
- No derivation of gravity, of the Standard Model, of the constants, of spacetime dimension, of
  the arrow of time, or of the Born rule.
- No statement about continuum, gauge or gravitational closed wholes beyond the identifications
  in FORMAL_CORE §5; no solution to the problem of time; no total-energy statement in an
  expanding spacetime.
- No empirical validation of anything in this packet.
- No claim that a future postulate could not have content: FORMAL_CORE §7 lists exactly what it
  would have to constrain.

## Retained failures (not softened)

- M1: facilitation benefit withdrawn after the matched F = ε control; single-cell failure at the
  0.9 standard for L = 4, 5; bare matter-to-memory writing fails; unconstrained flips violate
  retention; saturation barrier. The retention probabilities survive; their interpretation does not.
- Closure fork: C1-as-FLT inverse-Gaussian termination; D2 novelty, curvature-sign falsifier and
  the universal-positive-curvature claim withdrawn; B control registered, unexecuted.
- B2: eight-event endpoint failed at 1.93%; feasibility stop at N ≥ 2²⁸ within the envelope;
  moving operation neither established nor refuted.
- GRUT I–II: ζ-numerology, info-principle items, conformalon, economical w(z), the 32σ anchor,
  the 689 Hz falsifier, the α→TT route, the "biggest contribution" framing: all retracted.
- GRUT III: R1–R8, R13, R15, R19 rejections; the HJW route to P5 withdrawn.
- Level-0: H-HERM-2, strict H-ORD, two-property reduction falsified; O-6 declared-chain run not
  executed.
- Influence hierarchy: order-8 gate red; "no untracked pair constructible" never executed.
- V7: absent from the repository; its content known only through the handoff; the rung-3
  pole-versus-cut ask drafted and never sent.
- The toric bridge: post-hoc gauge-observable correction; separate branch; not merged.

## Unexecuted tasks (registered, not run here)

- B three-prescription contour control (closure fork), with the source corrections of the C2
  response applied first.
- Candidate 1 three-charge one-loop family; lifetime; nonlinear gravitational matching.
- Level-0 O-6 declared-chain run and magnitude maps.
- GRUT III O5/O7 items; independent review of CP-04 and CP-05.
- M1 error-correction cascade against the matched baseline; thermal inclusions; dynamical
  terminal detector.
- The toric branch merge decision; no toric code was run here.
- Any continuum or algebraic version of calculation 1 (split property, regulated subregion
  density matrix): not attempted.

## Verification performed for this packet

| Check | Result |
|---|---|
| `code/verify.py` (hashes + headline numbers) | all pass |
| `code/verify.py --rerun` | both calculations re-executed, every numeric leaf within 10⁻⁶ relative (bit-exact in practice) |
| Independent re-execution by a separate agent | 5030 + 24 numeric leaves bit-exact; statuses identical; two runs bitwise identical |
| Independent re-implementation from the specification (separate agent, code unseen) | 16 quantities agree to ≤ 10⁻⁶ relative; offsets explained (rounding, ε = 1/512, solver tolerance) |
| Adversarial code review of calculation 1 | "sound with minor issues"; nine items, all addressed (CALCULATIONS.md §Verification) |
| Adversarial code review of calculation 2 | [[REVIEW2]] |
| M1 matched-control re-run | identical to `control_results.txt` |
| Ledger entries | each written by one reader and adversarially re-read by a second against the same sources; [[LEDGERVERIFY]] |
| Completeness critic over the ledger | [[CRITIC]] |
| Attack-and-revise pass on this report | [[ATTACK]] |
