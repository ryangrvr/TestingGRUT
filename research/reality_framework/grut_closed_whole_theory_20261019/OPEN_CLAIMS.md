# Open claims, exact boundaries, retained failures

Packet `grut_closed_whole_theory_20261019`. This file is the claim boundary of the working
theory. Anything not listed under "Established" is not claimed.

## Established (and by what)

1. **The framework.** A closed whole is the theory datum Θ of FORMAL_CORE DEF 1.7; its rules
   are the standard theorems of FORMAL_CORE §2. Established by citation to standard results and
   by the two calculations. No novelty is claimed for any rule.
2. **Closure alone forces only the global theorems.** Unitary invariants, the exact reduced
   equation, the Nakajima–Zwanzig structure, CPTP reduced maps for product preparations, nested
   and necessary overlap compatibility, entropy inequalities, quasi-recurrence. Checked
   numerically in calculation 1 (S1, S8) in one supplied eight-qubit model.
3. **Closure plus one named hypothesis forces four conditional local results.** Equilibration on
   average under non-degenerate gaps (computed in calculation 1 S10), canonical typicality under a
   stated measure, Lieb–Robinson light cones for local H, and conservation-law limits on records
   (FORMAL_CORE COND 2.10–2.13). The hypotheses are properties of the supplied datum.
4. **Closure forces nothing local pointwise.** Of fourteen "closure forces X" claims
   (FORMAL_CORE §4), eight are refuted by a legal closed model, each exhibited in calculation 1
   (S2–S7, S9); three are settled by how the definitions are written (pure global state,
   constants and action, Born rule); three are scope statements the finite model cannot address
   (chaos exponents, spatial topology and curvature, cosmological energy). A state rule built from
   H alone and equivariant under unitaries is stationary (S9).
5. **Compatibility has a computable edge.** Bell-pair chain infeasible (two solvers; minimum
   marginal violation 0.289); max F_AB + F_BC = 3/2; isotropic-state chain p ≤ 2/3; triangle
   p ≤ 1/3 (calculation 2). Standard quantum marginal problem; the thresholds are solver outputs
   within 10⁻⁸ of the analytic values, with positive violation certificates just above them.
6. **The campaign ledger.** Twelve sectors re-read as campaigns; in every sector the recorded
   answer to the parent question is "nothing additional was forced in the stated class", with
   the sector's own computations as the demonstration (CAMPAIGN_LEDGER.md). Verification status
   of each entry is printed in the ledger.
7. **Reproduction.** Both calculations re-execute to identical results (`verify.py --rerun`) and
   agree with a separate agent's re-implementation from the specification (23 comparisons); the M1
   matched-control script reproduces its recorded table line for line.

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
- B2: eight-event endpoint failed at 1.93%; its contact-cancel control fails the ≥ 1/100 gate;
  the right-boundary parent target ε ≤ 1/1000 unachieved; B2-RA author-derived with independent
  review pending; stop at N ≥ 2²⁸ is a method and resource decision; moving operation neither
  established nor refuted.
- CP-12: the min-entropy bound is vacuous at every size its checker prints (through N = 10¹⁸);
  the parent relative-entropy theorem has no partner audit; its constants are "too loose for
  useful small-N probability bounds".
- Candidate 1: acceptance thresholds consolidated after execution, not preregistered; the
  ±0.11 allowance on K₁ is statistical, not a certified enclosure; small targets below the found
  minimum unresolved.
- GRUT I–II: ζ-numerology, info-principle items, conformalon, economical w(z), the 32σ anchor,
  the 689 Hz falsifier, the α→TT route, the "biggest contribution" framing: all retracted.
- GRUT III: R1–R8, R13, R15, R19 rejections; the HJW route to P5 withdrawn; "GRUT derives the
  initial state" withdrawn in the synthesis's own adversarial review; two unreconciled graviton
  verdicts in the secular-infrared survey.
- Level-0: H-HERM-2, strict H-ORD, two-property reduction falsified; O-6 declared-chain run not
  executed.
- Influence hierarchy: order-8 gate red; "no untracked pair constructible" never executed.
- V7: absent from the repository; its content known only through the handoff; the rung-3
  pole-versus-cut ask drafted and never sent; the bath behind the in-in action never specified
  (NO-PHYSICAL-BATH-CANDIDATE).
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

| Check | Result | Record |
|---|---|---|
| `code/verify.py` (hashes + headline checks) | [[NCHECKS_TABLE]] | run it |
| `code/verify.py --rerun` | both calculations re-executed; every numeric leaf identical to the committed results | run it |
| Independent re-implementation from the specification (separate agent, code unseen) | 23 comparisons agree to ≤ 10⁻⁶ relative; re-compared against the committed results after the final edits | `results/verification/independent_reimplementation_*` |
| Independent re-execution (separate agent, 19 Oct) | bit-exact, but of the pre-review version (superseded) | `results/verification/independent_reexecution_19oct_prereview.md` |
| Adversarial code review of calculation 1 | "sound with minor issues"; nine items, all applied | `results/verification/code_review_calculation1.md` |
| Adversarial code review of calculation 2 | "sound with minor issues"; seven items, all applied | `results/verification/code_review_calculation2.md` |
| M1 matched-control re-run | identical to `control_results.txt` | `results/m1_matched_control_rerun.txt` |
| Ledger entries | [[LEDGER_VERIFY_TABLE]] | `ledger/`, CAMPAIGN_LEDGER.md |
| Primary literature | [[LIT_TABLE]] | SOURCES.md |
| Attack on the packet | [[ATTACK_TABLE]] | `results/verification/attack_round*/` |
