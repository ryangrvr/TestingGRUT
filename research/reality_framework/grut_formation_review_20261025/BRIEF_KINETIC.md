# Brief: kinetic case-L threshold at the heavy-mediator spot point (issue #13)

26 October 2026. Claude owns this calculation at Codex's request (issue #13 comment 6413510190); Codex retains the
independent source, collision-inventory and conservation review. Read BRIEF.md first; all of its rules apply (write
only inside `tracks/kin_<key>/`; read Codex packets only via `git show`/`git archive` into the scratchpad; no commits;
own code; labels [D]/[C]/[Q]/[E]/[A]; retrospective vs held out; failures reported, not repaired).

## The question

The heavy-mediator window (REVIEW_HEAVYWINDOW.md; tracks/hw_critic/CONSOLIDATED.md) passes the mediator-debris test by
a wide margin. Its only failing gate is the decay of primordial longitudinal vectors into hard X X̄ pairs, whose
verdict depends on which case-L threshold is right. The two equilibrium proxies differ by 10×:
- isentropic proxy: R_L ≈ 80–88 (passes at e = 1e-7);
- energy-conserving proxy: R_L ≈ 7–8 (≈ 15–20 with the on-shell s route) (fails by 6–9 at e = 1e-7; passes for
  e ≳ 3e-6).

Replace both proxies by a kinetic treatment and report **f_sym(e, source, release spectrum, formation kernel)**: the
final charge-symmetric hidden mass relative to the net-charge mass, against the 1% engineering allowance (f = 0.01).
Give a conditional R crossing only where monotonicity and closure hold. Energy conservation is required. The task is
not to choose whichever threshold keeps the point alive.

## Fixed inputs

- **The point.** Spot point of 6406685915: H_i = 1.3335e11 GeV, M_Y = 30 H_i, h = 40 H_i, μ = 7.15945616 H_i,
  κ = 0.003458684721, r = 0.45, λ_Y = λ_Σ = 1, e = 1e-7 (and 1e-6 as a second case); light sector g_X = 0.1,
  λ_X = 1e-3, λ_φ = 0.03, v = 0.924089637 GeV, Q = 200000, Y_Q = 4.8294e-9. Conditional heavy preparation and
  normalised light sector: Codex's `research/reality_framework/grut_heavy_vector_budget_20261025/INPUTS.json` at
  commit 225028b. Release: ε = −0.14572 (primary window a = 30–40), N_X = −2788.18, E_X = 38386.25 (H_i units);
  mediator debris R = 1.14 (N_YH 0.01725, N_YL 0.33447, N_σ 3e-9), to be carried as a small extra source.
- **The vector source.** Codex's `research/reality_framework/grut_kinetic_contract_20261026/DECAY_SOURCE.json` at
  commit 9877518: time profile and comoving momentum distribution of the daughters per sign, with normalisation and
  moments, for e = 1e-7 (pairs/net = 7.406624e-4, R = 40.5626, unprocessed rest-mass f = 0.0015025, mean p/m_X at the
  case-L density coordinate 77200.9 without collisions) and e = 1e-6. This is the finite band k/k* = 0.2–5 (a
  conservative partial source). The broader-grid source of our tracks/hw_gates (R ≈ 48.5, N_A = 2.47 H_i⁻³) is the
  second input; do not mix one spectrum's normalisation with the other's decay mass or Higgs profile.
- **Scalar reservoir.** Codex's THERMODYNAMIC_CONTROL.json (same packet): one neutral scalar s (m_s = 0.22635 GeV,
  trilinear 0.018482 GeV, vacuum width 1.7333e-5 GeV) and two charged degrees; at T = 12 m_X the scalar holds 33.24% of
  the gas energy in chemical equilibrium; the matched-energy T/m shifts from 8.118 to 7.344 at the central R.
- **Our existing code (reuse with changes stated):** tracks/ls_threshold/code (lsrelic.py, lsthermo.py, run_ke*.py),
  tracks/sp_yledger/code (thermal_relic.py, rates.py), tracks/gp_leftovers/code (amp42.py, channels.py,
  equilibrium.py, catalysis.py), tracks/sp_critic/code/rmax_relic.py, tracks/hw_gates/code (gate1_vector.py).
  Codex's KINETIC_SPEC.md (9877518) states the reaction and conservation requirements; follow it.

## Requirements (from KINETIC_SPEC.md, condensed)

1. **State variables.** At minimum X, X̄, the light scalar s (number and energy), the remaining parent source, and the
   classical/condensed (object) sector. Either momentum distributions f_±(q), f_s(q) plus an object sector, or
   per-species number and energy moments after checking that elastic collisions enforce the assumed spectral shape
   fast enough. A single temperature and a single pair chemical potential cannot be imposed before the relaxation that
   creates them: either evolve the pre-onset relaxation or bracket it. Match the releasing condensate's energy and
   charge as well as the vector source.
2. **In-medium masses and phase.** Specify the scalar expectation value and thermal masses; the old rates use vacuum
   masses in the broken phase and cannot establish that phase at every temperature. Evaluate or bound the real
   scalar's cubic and quartic self-reactions.
3. **Collision network.** A → X X̄ (source); elastic X/X̄/s; s ↔ X X̄ forward and reverse with s number and energy
   storage; ss ↔ X X̄ and its reverse; charge-conserving 4 ↔ 2 with energy-dependent amplitudes (non-relativistic
   threshold constants are not the kernel for p/m_X ≈ 77000 daughters; they may be an adequate late limit after
   relaxation, with a validity check and a matched high-energy treatment or bound); scalar self-reactions; formation,
   evaporation, capture and merger with an explicit kernel or a labelled bracket. Remove overlap between real
   sequential s processes and their resonant-scattering description, and state the subtraction.
4. **Formation.** Gas kinetics alone cannot select the object population. Report under labelled alternatives: a
   gas-only control; the old instantaneous-condensation control; and any supported nucleation/fragmentation kernel.
   Keep finite-Q energy and ω = dE/dQ distinct; include the object-number work term or an explicitly linear
   comparator; track the hard X emitted by X̄ capture (bb211d1) rather than calling capture erasure.
5. **Checks.** Vanishing internal collision energy moment; conserved total charge; stationarity of equilibrium
   distributions under forward/reverse kernels; nonnegative distributions; momentum and time grid refinement with
   residuals reported separately from pass/fail; parent and mean-field work included; source number and energy
   moments verified against Codex's analytic and ODE values; a forced-thermal-initialisation comparison with the
   relaxation processes that justify or invalidate it; late s quanta that decay into pairs tracked; anti-object charge
   reported; a no-source control so the release's own leftover population is not attributed to the vectors; an R
   rescaling scan at fixed spectral shape (R = 1–90, including 40.5626 and 48.5) labelled as a diagnostic, not a new
   cosmological point.

## Clarifications from Codex (issue #13 comment 6414910614, added 26 October)

6. **Initial spectrum.** The release's energy and charge alone do not determine its antiparticle population. The
   solver must derive the free X̄ population of the releasing condensate from the release spectrum (our certified
   free-mode extraction in tracks/spectrum, scaled as stated), or bracket it with labelled alternatives, or demonstrate
   that the relaxation stage erases its influence on f_sym. Report which of the three applies and the sensitivity.
   **Matching caveat (Codex, comment 6415201060):** the certified spectrum in tracks/spectrum was extracted from the
   prepared spherical collapse of the *light-target* condensate (ff49e5d family), not from this heavy-field release.
   Before it is used here, state and justify the matching (which quantities are carried over, in what units, and why
   the heavy release's free X̄ fraction should follow it); otherwise treat the release X̄ population as a labelled
   bracket (for example 0, the light-collapse fraction, and a thermal-onset value) and report the spread.
7. **Benchmark normalisation.** Use one and the same denominator for the 1% allowance in every run (the net-charge
   mass of the final object population at the fixed yield Y_Q), with the finite-object masses (E(Q), not Q m_X)
   accounted for explicitly, so that runs with different formation rules are compared on the same footing.

## Deliverable

f_sym for: e = 1e-7 with the finite-band source; e = 1e-7 with the broader source; e = 1e-6 with its own source; each
under the formation alternatives; plus the diagnostic R scan and the conditional R_L crossing. State which of the two
old proxies the kinetic result is closer to and why. Every approximation with a stated regime or error. No
cosmological viability, observational evidence or GRUT prediction is claimed.
