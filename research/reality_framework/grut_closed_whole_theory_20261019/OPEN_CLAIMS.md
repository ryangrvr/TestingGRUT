# Open claims, exact boundaries, retained failures

Packet `grut_closed_whole_theory_20261019`. This file is the claim boundary of the working
theory. Only what is listed under "Established" and "Conditional" is claimed, each in the scope
stated there.

## Established (and by what)

1. **The framework.** A closed whole is the theory datum Θ of FORMAL_CORE DEF 1.7; its rules
   are the standard theorems of FORMAL_CORE §2. Established by citation to standard results and
   by the two calculations. No novelty is claimed for any rule. In the finite class this is
   definitional (closure is defined inside quantum theory); read as purification, closure
   singles out quantum theory among probabilistic theories (FORMAL_CORE §0).
2. **What closure alone forces.** Unitary invariants, the exact reduced equation, the
   Nakajima–Zwanzig structure, CPTP reduced maps for product preparations, nested and necessary
   overlap compatibility, quasi-recurrence, and pointwise local inequalities: the entropy
   inequalities, monogamy, instantaneous no-signalling, the athermality budget (THM 2.14) and
   protect-or-form for records (THM 2.15). Checked numerically in calculation 1 (S1, S8, S11) in
   one supplied eight-qubit model.
3. **What closure plus one named hypothesis forces.** Equilibration on average for a pure
   global state under non-degenerate gaps (computed in S10); the second-law identity
   Σ = I + D ≥ 0 for a complement prepared in a Gibbs state uncorrelated with A (computed in S11; general
   passivity is not sufficient); relaxation for an
   infinite complement under stated mixing or weak-coupling conditions; canonical typicality
   under a stated measure; Lieb–Robinson light cones for local H; fluctuation–dissipation and
   the three-way form of any difference for KMS states; the chaos bound; conservation and clock
   limits on apparatus inside the whole (FORMAL_CORE COND 2.10–2.20). The hypotheses are
   properties of the supplied datum.
4. **What closure does not force.** Of fourteen "closure forces X" claims (FORMAL_CORE §4), six
   are refuted by a legal closed model exhibited in calculation 1 (a thermal complement,
   pointwise relaxation, a Markovian kernel, a unique local law, determination of the future by
   ρ_A(0), record persistence); one is split (the arrow: its direction is not forced, Σ ≥ 0 is
   forced given an uncorrelated Gibbs reference preparation); four are settled by how the definitions are
   written (a pure global state, a preferred initial state, the constants and action, the Born
   rule); three are scope statements the finite model cannot address (chaos exponents, spatial
   topology and curvature, cosmological energy). A unitarily equivariant state rule built from H
   alone is stationary, so its one-time marginals do not evolve (S9); multi-time and response
   physics remain.
5. **Compatibility has a computable edge.** Bell-pair chain infeasible (two solvers; minimum
   marginal violation 0.289); max F_AB + F_BC = 3/2; isotropic-state chain p ≤ 2/3; triangle
   p ≤ 1/3 (calculation 2). Standard quantum marginal problem; the thresholds are solver outputs
   within 10⁻⁸ of the analytic values, with positive violation certificates just above them. This
   is the one-time necessary condition; the process marginal problem is not computed.
6. **The campaign ledger.** Twelve sectors re-read as campaigns. In the sectors that computed in
   a closed or algebraic class (s02, s05–s10), nothing additional was forced in the stated class.
   s01, s03, s04, s11 and s12 did not pose the closed-whole question (open or classical classes,
   declared baths, or program ledgers); their verdicts are re-readings, not closure
   demonstrations. Verification status of each entry is printed in the ledger.
7. **Reproduction.** Both calculations re-execute to identical results (`verify.py --rerun`) and
   agree with a separate agent's re-implementation from the specification (23 rows: 21 to 10⁻⁶
   relative, two numerical zeros absolutely); the M1 matched-control script reproduces its
   recorded table line for line.

## Conditional (stated in a class; author-derived or unaudited where marked)

- Sector results marked *conditional constraint* in the ledger: GRUT III D11/D13 slow-roll
  thermality; Candidate 1's one-loop free-carrier margin at its benchmark (conditional on
  reduced-saddle stability); the CP-12 min-entropy obligation (parent theorem unaudited; vacuous
  at printed sizes); M1's protect/form constraint at ε = 0; access is not arbitrary bookkeeping
  in the P-5 class; the B1 delivery window (author-derived, review pending, partly conditional on
  the S9 control); the Level-0 in-class non-implications; the influence-hierarchy cumulant
  ladder on constructed pairs (order-8 gate red); the GRUT I–II boundaries at their stated
  scopes. Each is conditional on its supplied Θ.
- The finite-model statements of calculation 1 are statements about that model; they illustrate
  theorems that hold generally, but the numbers are the model's. Its interventions are
  externally timed.

## Not established, and not claimed

- No additional law, state-selection rule, access rule, clock rule or measurement rule beyond
  standard quantum theory. None survived any campaign (WORKING_THEORY §6). A modified global law,
  a measure over Θ and relational constraints between slots were never tested.
- No combined "GRUT law" uniting Candidate 1, M1, GRUT III or the Level-0 floor. They live in
  different state spaces; no common action, state space, constraints or calibration has been
  established.
- No derivation of gravity, of the Standard Model, of the constants, of spacetime dimension, of
  the direction of the arrow of time, or of the Born rule.
- No statement about continuum, gauge or gravitational closed wholes beyond the identifications
  in FORMAL_CORE §5; no solution to the problem of time; no total-energy statement in an
  expanding spacetime.
- No prediction that distinguishes any model instance from the known-physics baseline, and no
  empirical validation of anything in this packet.
- No claim that a future postulate could not have content: FORMAL_CORE §7 lists the kinds the
  program considered and names known kinds outside that list.

## Retained failures (not softened)

Summary; the complete list is each ledger entry's "Failures retained" section.

- M1: facilitation benefit withdrawn after the matched F = ε control, which was specified after
  the review and is no worse at every sampled time (equal at t = 0, strictly better at every
  t > 0); single-cell failure at the 0.9 standard for L = 4, 5; bare matter-to-memory writing
  fails; unconstrained flips violate retention; radiative stability not established. The
  retention probabilities survive; their interpretation does not.
- Closure fork: D2 novelty, the curvature-sign falsifier and the universal-positive-curvature
  claim withdrawn; the coarse-ε criterion failed and only a changed criterion passes; the B
  control is registered and unexecuted, and a C1-as-FLT termination would rest on an
  inverse-Gaussian coefficient quoted from the FLT literature, not reproduced.
- B2 and WM: the eight-event endpoint failed at 1.93%; its contact-cancel control fails the
  ≥ 1/100 gate; the right-boundary parent target ε ≤ 1/1000 unachieved; B2-RA author-derived with
  independent review pending; the stop at N ≥ 2²⁸ is a method and resource decision; moving
  operation neither established nor refuted; A1 refutes autonomous every-time retention at one
  permitted point; B1's two-event control fails.
- CP-12: the min-entropy bound is vacuous at every size its checker prints (through N = 10¹⁸); the
  parent relative-entropy theorem has no partner audit; its constants are "too loose for useful
  small-N probability bounds"; formation not established; Milestone A not achieved.
- Candidate 1: acceptance thresholds consolidated after execution, not pre-registered; the ±0.11
  allowance on K₁ is a numerical allowance (analytic tail bounds plus an empirical allowance),
  not a certified enclosure, not a statistical interval and not a higher-loop bound, and its
  angular term was not challenger-audited; small targets below the found minimum unresolved.
- GRUT I–II. Retracted: ζ-numerology, info-principle items, the conformalon, the 32σ anchor,
  the "biggest contribution" framing. Failed as derivations but retained as settled-negative
  boundaries: the economical w(z) and the α→TT route ("not a formal no-go"). The two-survivor
  Ward classification retired as unconditional (2026-08-14), adding noise-kernel transversality as
  an input; x < 0.358 demoted; the κ pre-registration self-unblinds; the mu_linear boundary
  weaker than stated.
- June deposit and V7: the 689 Hz falsifier withdrawn; V7 absent from the repository, its content
  known only through the handoff; the rung-3 pole-versus-cut ask drafted and never sent; the
  bath behind the in-in action never specified (NO-PHYSICAL-BATH-CANDIDATE).
- GRUT III: reviewer points accepted in RESPONSE_TO_REVIEWS.md (R1–R8 and R13 retractions or
  relabellings; R15, the v11 headline coefficient, refuted; R19, the v13 inference, refuted); the
  HJW route to P5 withdrawn; "GRUT derives the initial state" withdrawn in the synthesis's own
  review; two unreconciled graviton verdicts in the secular-infrared survey.
- Level-0: H-HERM-2, strict H-ORD and the two-property reduction falsified; O-6 declared-chain
  run not executed.
- Influence hierarchy: order-8 gate red; "no untracked pair constructible" never executed; P-6's
  L-C gate red; S2-1's first execution void.
- Program ledgers: three negatives flagged MUST-BE-REOPENED (the ROOT-1 battery's hard-coded
  gate, the composite no-pole negative, the tt_worldline decay).
- The toric bridge: post-hoc gauge-observable correction; by the letter of the pre-registration
  R2 fails; separate branch; not merged; independent review outstanding.

## Unexecuted tasks (registered, not run here)

- B three-prescription contour control (closure fork), with the source corrections of the C2
  response applied first; it would also reproduce, or not, the FLT inverse-Gaussian coefficient.
- GRUT III N1 (horizon-state selection by P1, pre-registered as decisive) and N2 (a horizon
  principle fixing a renormalisation constant); GRUT III O5/O7 items; independent review of
  CP-04 and CP-05.
- Candidate 1 three-charge one-loop family; lifetime; nonlinear gravitational matching.
- Candidate 1 Z₂ domain walls and the hidden sector's initial thermal history; an external
  calibration of the total vacuum energy for its cosmological application (s08 addendum).
- Level-0 O-6 declared-chain run and magnitude maps.
- M1 error-correction cascade against the matched baseline; thermal inclusions; dynamical
  terminal detector.
- The toric branch merge decision and its quantum OTOC; no toric code was run here.
- The partner audit of the CP-12 parent theorem; independent review of B2-RA and B1.
- The completeness critic over the ledger and full adversarial re-reads of s11 and s12.
- An internal writer for calculation 1 (a clock and controller inside Ω); any continuum or
  algebraic version of calculation 1 (split property, regulated subregion density matrix); the
  process marginal problem for the sectors.

## Verification performed for this packet

| Check | Result | Record |
|---|---|---|
| `code/verify.py` (hashes + headline checks) | every file in MANIFEST.json hashed; 31 headline results and 2 provenance hashes; all pass | run it |
| `code/verify.py --rerun` | both calculations re-executed; identical leaf paths; every numeric leaf identical to the committed results | run it |
| Independent re-implementation from the specification (separate agent, code unseen) | 23 rows: 21 agree to ≤ 10⁻⁶ relative, two numerical zeros agree absolutely; re-compared against the committed results after the final edits | `results/verification/independent_reimplementation_*` |
| Independent re-execution (separate agent, 19 Oct) | bit-exact, but of the pre-review version (superseded) | `results/verification/independent_reexecution_19oct_prereview.md` |
| Adversarial code review of calculation 1 | "sound with minor issues"; nine items, all applied | `results/verification/code_review_calculation1.md` |
| Adversarial code review of calculation 2 | "sound with minor issues"; seven items (six actionable), all applied | `results/verification/code_review_calculation2.md` |
| M1 matched-control re-run | identical to `control_results.txt` | `results/m1_matched_control_rerun.txt` |
| Ledger entries | s01–s10 adversarially re-read by a second reader, corrections folded in; s11–s12 single-reader; claims drawn from s10–s12 checked in the main session; completeness critic not run | `ledger/`, CAMPAIGN_LEDGER.md, `results/verification/s10_s12_row_checks.md` |
| Primary literature | Fewster–Verch, BFV, Masanes–Galley–Müller, Witten, Pollock et al. (PRL, PRA), Chen et al., Cotler–Penington–Ranard: full text; Feynman–Vernon: not accessed | SOURCES.md |
| Attack on the packet | round 1: 2 blocking, 11 major, 22 minor; round 2: 0 blocking, 30 major, 50 minor; all addressed; no third round on the final version | `results/verification/attack_round1/`, `attack_round2/` |
