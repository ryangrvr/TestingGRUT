# Brief: one calculable vacuum mechanism with a solved example (issue #13)

21 October 2026. This follows the vacuum-energy principle pass
(`research/reality_framework/grut_vacuum_principle_20261021/`, tracks/*/NOTES.md and SYNTHESIS.md). Candidate 1
remains paused. Its status, as Codex corrected it: "Candidate 1 is not independently test-ready; universal equivalence
to cold dark matter remains unproved."

## Mandate

The owner (21 October): "a coherent, approximate, unproven working hypothesis is a legitimate goal now ... The
immediate deliverable should contain equations and a solved example ... Allow transparent refinement. A failed
prediction can motivate a new version. Keep the original result and identify what changed. A borrowed mechanism can
be a useful foundation when credited accurately. GRUT's contribution must then be stated precisely. The next milestone
is one mechanism, one consistent set of equations, and one worked scenario showing how its claimed behavior emerges."
Matter, gravity and vacuum behaviour must follow from compatible equations, with supplied assumptions identified. The
distinction that matters is between approximations we can quantify and contradictions in the construction.

Codex (issue #13 comment 6324495415): write the action or equivalent closed equations, the degrees of freedom, the
state and boundary assumptions, and the provenance of every parameter. Demonstrate at least one complete scenario,
including the gravitational response of the vacuum mechanism. The result is either a **candidate construction with a
solved example**, or the strongest partially specified construction together with its **single next mathematical
dependency**. "At most one calibrated parameter" is a design goal, not a law. List every input. If a calibration fixes
the residual, say so. Internal inconsistencies cannot be labelled refinement. A parameter update is a new version, and
the old predictions are kept.

## What the vacuum pass found (short)

- The energy-origin obstruction holds on its scope (fixed-background closure is blind to H -> H + cI; gravity sees
  -C g_mu_nu). No route surveyed fixes the total renormalized rho_vac/M^4 with sign and scale, stays stable without
  retuning, and gives a held-out consequence. Nothing is GRUT-native.
- Sequestering: global KP removes every constant vacuum term but its residual is a functional of a supplied history.
  Its linear-trigger realization is retrospectively inconsistent (<R> = -115.6 s H0^2, no root). Local and Omnia
  Sequestra keep a free residual. Unimodular/HT relocate the value to the state and protect only spacetime-constant
  terms; phase transitions gravitate as in GR. Everpresent Lambda fixes a spread, not the mean, and our history is
  improbable in its simplest model.
- The reviewers flagged one uncredited, closest closure-type argument: **Klinkhamer-Volovik q-theory**. A
  self-sustained vacuum with no external environment has zero gravitating vacuum energy in equilibrium, through a
  conserved vacuum variable q. GRUT's "the whole has no bath, so its vacuum level should not gravitate" is almost
  verbatim this argument.

## The chosen mechanism (borrowed, credited)

Klinkhamer and Volovik, q-theory: PRD 77, 085015 (2008), arXiv:0711.3170; PRD 78, 063528 (2008), arXiv:0806.2805;
PRD 80, 083001 (2009), arXiv:0905.1919; JETP Lett. 91, 259 (2010), arXiv:0907.4887; Volovik, IJMPD 15, 1987 (2006),
gr-qc/0604062; and their later q-theory papers. Every element taken from them is credited to them.

Starting transcription (ours; every agent must check it against the sources or derive it independently).
Conventions: signature (-,+,+,+), hbar = c = 1, M = reduced Planck mass = 2.435e18 GeV, K = M^2 = 1/(8 pi G).
Lambda symbols are energy densities unless stated. rho_Lambda,obs = 2.52e-47 GeV^4 = (2.24 meV)^4, rho_Lambda/M^4 =
7.16e-121 (Planck 2018 inputs, exposed data).

    S = ∫ d^4x sqrt(-g) [ (1/2) K(q) R - eps(q) ] + S_m[psi, g]
    F = dA (A a 3-form),  F_{kappa lambda mu nu} = q sqrt(-g) [kappa lambda mu nu]   =>  q = F~ / sqrt(-g)

with F~ metric-independent. Our derivation (to be verified) gives:

    (q)  eps'(q) - (1/2) K'(q) R = mu                      mu: a spacetime constant (integration constant)
    (g)  K G_mu_nu = T^m_mu_nu - (eps - mu q) g_mu_nu - (g_mu_nu Box - nabla_mu nabla_nu) K
    (m)  nabla^mu T^m_mu_nu = 0
    gravitating vacuum energy  rho_V = eps(q) - mu q   (in flat equilibrium: eps - q eps')

Equilibrium (Minkowski): eps'(q0) = mu0 and eps(q0) - mu0 q0 = 0.
**Proposed GRUT formalization (to be tested, credited to KV):** "the closed whole has no exterior, so its equilibrium
vacuum carries no external pressure", i.e. the state of the whole has mu = mu0[eps]. Which eps (pre- or
post-transition, bare or fully renormalized) the condition refers to is a central open question, not a detail.

## What "done" means for the construction

1. Action or closed equations; degrees of freedom; state and boundary assumptions; provenance of every parameter
   (derived / supplied [P] / calibrated / literature).
2. Coupled evolution and conservation equations, with the Bianchi identity verified.
3. One complete solved scenario including the gravitational response of the vacuum mechanism (FRW with matter and
   radiation, a vacuum-energy change at a transition, the relaxation, and the expansion history it produces).
4. One relation between measurable quantities beyond the calibration target, with its uncertainty and every input
   listed; or an explicit statement that none follows, plus the single next mathematical dependency.
5. GRUT's contribution stated precisely (it may be small; say so if it is).
6. The four-question gate (Codex): (1) which extra structure constrains the gravitational vacuum term; (2) does it fix
   a residual (sign and scale) or only remove a large term; (3) stable under heavy thresholds/loops at the stated order
   without retuning; (4) what independent consequence follows, and which data are already exposed.

## Rules for every agent

- Write ONLY inside `research/reality_framework/grut_q_vacuum_construction_20261021/tracks/<key>/`: NOTES.md
  (equations, derivations, numbers, verdicts, read levels, limits), code/ (numpy/scipy/sympy), results/*.json.
- No git commands other than read-only `git show`/`git log`. Do not edit other packets.
- Literature: give arXiv/DOI identifiers and the read level achieved (abstract / section / full). Do not copy paper
  text or PDFs into the repository; keep downloads in the scratchpad
  /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/qtheory/.
- Label claims: derived [D], computed [C], quoted [Q], estimate [E], argued [A], recalled [R]. Never claim a new
  theory, an empirical validation or a confirmatory test. Label already-observed data as retrospective.
- Approximations must be quantified (state the error or the regime). An inconsistency is reported as an
  inconsistency, never as a refinement.
- Every calibrated quantity is named as calibrated, with the datum used. Do not tune a parameter to the observed
  vacuum energy and then call the match a prediction.
- Credit: anything matching published work is credited to it. "GRUT-native" only for an assumption the cited
  literature does not use.
- Return value: a digest of at most about 700 words with the key equations, numbers and verdicts, and pointers to
  NOTES.md sections. The NOTES.md is the record.
