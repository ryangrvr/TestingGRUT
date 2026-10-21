# Brief: can a principle fix the total renormalized vacuum energy? (issue #13, Codex 02b4515)

21 October 2026. Candidate 1 dynamics is paused. Codex's question (issue #13 comment 6324148604; packet
`research/reality_framework/grut_scale_feasibility_20261021/NEXT_DIRECTION.md` at commit 02b4515, read with
`git show 02b4515:<path>`):

**Energy-origin obstruction.** On a fixed background, H -> H + cI leaves density-matrix evolution, reduced states,
measurement probabilities, energy gaps and normalized Gibbs states unchanged, so a rule built only from those closure
properties cannot select c. With a dynamical metric, adding -C ∫ sqrt(-g) d^4x to the matter action shifts
T_mu_nu by -C g_mu_nu and changes curvature. A compensated counterterm relabelling (rho_m,vac + C, rho_Lambda,ct - C)
changes nothing physical. The obstruction is elementary and scoped to fixed-background closure arguments; it is not a
no-go theorem for quantum gravity, closed gravitational constraint systems, nonlocal proposals or new dynamics.

**Task.** Try to refute the obstruction and specify the minimal additional principle/equation that escapes it, if any.
Compare its actual assumptions with existing adjustment and sequestering proposals before calling it new. Write the
independently testable consequence before building any candidate.

**Success gate for any proposal** (Codex's four questions). Target: the TOTAL renormalized rho_vac / M_Pl^4 with the
reduced-Planck convention explicit. Write the principle as equations and list every parameter, boundary condition,
state measure and cosmological input. Then:
1. Which extra structure makes it constrain the gravitational vacuum term that the shift-invariant matter argument
   leaves free?
2. Does it fix a residual (sign and scale) or only remove a large contribution while keeping an arbitrary constant?
3. Is the value stable when a heavy threshold/loop contribution is included at the stated order, without retuning to
   data?
4. What independent dynamical or observable consequence follows beyond fitting the cosmological constant, and which
   already-exposed data cannot be called held out? (H0 is not a parameter-free target: matter, radiation and curvature
   inputs remain.)

**GRUT context.** The parent question of the program is what a closed system, whose every subsystem's bath is the rest
of the universe, forces beyond known physics (`research/reality_framework/grut_closed_whole_theory_20261019/`:
WORKING_THEORY.md, FORMAL_CORE.md, CAMPAIGN_LEDGER.md; earlier closure work in
`research/reality_framework/grut_closure_fork_20261019/` including CLOSURE_SELECTION_RESPONSE.md, and GRUT III
material referenced in the ledger). Closure alone did not force new physics in earlier work; no state-selection rule
survived. Do not repeat those negatives as new.

**Rules for every track agent.**
- Write ONLY inside `research/reality_framework/grut_vacuum_principle_20261021/tracks/<key>/`: NOTES.md (equations,
  derivations, verdicts, read levels, limits), code/ if any (numpy/scipy/sympy), results/*.json if any.
- No git commands other than read-only `git show`/`git log`. Do not edit other packets.
- Literature: give arXiv/DOI identifiers and the read level actually achieved (abstract / section / full). Do not copy
  paper text or PDFs into the repository; keep downloads in the scratchpad
  /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/.
- Distinguish: derived [D], computed [C], quoted [Q], estimate [E], argued [A]. Never claim a new theory, an empirical
  validation or a confirmatory test. Label any already-observed data as retrospective.
- Novelty: anything matching an existing proposal must be credited to it; "GRUT-native" only if it uses an assumption
  that the cited literature does not.
