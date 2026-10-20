# GRUT working theory: the closed whole

D. Ryan Grover's program, synthesis of 19–20 October 2026. Packet `grut_closed_whole_theory_20261019`,
branch `codex/grut-closure-fork-20261019`, written from the handoff at commit `96a7e7e` and the
corrected M1 record at `4e9a161`. Detailed derivations are in [FORMAL_CORE.md](FORMAL_CORE.md),
every number in [CALCULATIONS.md](CALCULATIONS.md), every sector in
[CAMPAIGN_LEDGER.md](CAMPAIGN_LEDGER.md) and `ledger/`, and claim boundaries in
[OPEN_CLAIMS.md](OPEN_CLAIMS.md).

## 0. What GRUT now says, in one page

**Status.** GRUT is a working closed-system framework of known physics with conditional model
instances and computed boundaries. No additional law or state-selection rule has survived any
campaign. Nothing here is validated universal physics, and neither gravity nor "all reality"
has been derived.

**The parent object is a question.** What does a closed system, whose every subsystem's bath is
the rest of the universe, force that known physics does not already contain? GRUT asks this
question in one place after another and reports the answer, including "nothing additional in
this class, and here is the computation". The old parent object, which the handoff describes as
V7's in-in action S_CTP treated as a "genotype" (the V7 file is not in the repository), is the
Gaussian approximation to a derived object, the Feynman–Vernon influence functional, and sits
downstream of a supplied state, dynamics and access structure. It is a calculational tool.

**What calculates.** A model instance is a theory datum
Θ = (algebra, dynamics, preparation, access structure, instruments and records, clock), all
supplied. Given Θ, the global law ρ_Ω(t) = U(t)ρ_Ω(0)U(t)†, the restriction ρ_A = Tr_Ā ρ_Ω and the
Born rule produce every prediction, conditional on Θ. Calculation 1 carries one fully specified
eight-qubit Θ from preparation through evolution to measurement with exact energy and
correlation accounting; calculation 2 computes the compatibility conditions that overlapping
descriptions must satisfy. Both run in about a minute.

**What closure constrains.** From closure alone, the theorems of standard quantum theory:
unitary invariants, the exact reduced equation dρ_A/dt = −i[H_A, ρ_A] − i Tr_Ā[H_{AĀ}, ρ_Ω], the
Nakajima–Zwanzig structure with its initial-correlation term, complete positivity of reduced
maps for product preparations, nested and necessary overlap compatibility, quasi-recurrence in
finite dimension, and some local inequalities that hold at every instant: the entropy
inequalities, monogamy, a relative-entropy "athermality budget" that no subsystem can exceed,
and protect-or-form for records (a record the closed dynamics can write is not protected for
every state). From closure plus one named hypothesis on the supplied datum: equilibration on
average when energy gaps are non-degenerate (for a pure global state; its size is set by the
preparation); the second-law identity Σ = I + D ≥ 0 when the complement starts in a Gibbs state
uncorrelated with A (passivity alone is not enough); relaxation when the complement is infinite and satisfies stated mixing or
weak-coupling conditions; canonical typicality under a stated measure; light cones for local H;
fluctuation–dissipation for KMS states; and conservation-law and clock limits on apparatus
inside the whole. Closure forces **no local dynamical law**: not monotone relaxation, not
Markovianity, not a kernel, not a unique local law, not thermality of the complement, not the
direction of the arrow, not record persistence. Each of these is exhibited by a legal closed
model in calculation 1. A preferred state, the constants and the action are not forced because
they are entries of Θ by definition. A unitarily equivariant state rule built from H alone is
stationary, so its one-time marginals do not evolve.

**Which proposed additions did not survive** (details by slot in §6):
- *State selection*: C2 as a selection rule; the novelty of slow-roll thermality; the
  curvature-sign falsifier.
- *Dynamics*: the CTP action as a parent principle; in-in selection of the Einstein–Hilbert
  action; the TT projector as forced; any GRUT-internal kernel-selection principle; the founding
  ζ-hypothesis, the conformalon and the information principle; the facilitation benefit of the
  M1 memory.
- *Access*: an endogenous access structure (EA-0).
- *Measurement*: decoherence as a derivation of the Born weights; the HJW route to the outcome
  rule; influence-data positivity as a new principle.
- *Retracted empirical claims*: the 689 Hz falsifier; the 32σ ISW anchor; a GRUT-specific
  prediction of GRUT III.
- *Mechanisms*: autonomous every-time retention (refuted at one permitted parameter point); the
  self-correction/no-scrambling bridge in the classical 4D toric code (as reported on an unmerged
  branch, not re-run here, independent review outstanding).

**What remains a supplied choice.** The six slots of Θ in every model instance: the algebra
(GRUT III P0 and data), the dynamics (Candidate 1's couplings; GRUT I–II's projector), the
preparation (CP-12 preparation classes; the no-boundary state that C1/C2 did not select; the
arrow, which GRUT III locates as postulate P4), the access structure (EA-0; M1's register), the
instruments and records (M1's writer), and the clock; plus the measurement rule (P5). The
2026-08-17 register of the responsive-vacuum stance recorded net +16 underived inputs, an empty
derived tier and an empty prediction ledger (GRUT III re-accounts the same stance as five
postulates); no register prices the other model instances.

**Two caveats about the answer.** In the finite class closure is defined inside quantum theory,
so "nothing beyond standard quantum theory" could not have come out otherwise there; read
operationally, closure as purification singles out quantum theory among probabilistic theories,
so in that reading closure forces known quantum physics rather than nothing (FORMAL_CORE §0).
And several sectors did not themselves pose the closed-whole question (s01, s03, s04, s11, s12:
open or classical classes, declared baths, or program ledgers); their verdicts here are
re-readings, not closure demonstrations.

## 1. The parent question and what counts as an answer

A closed whole has no exterior. Its subsystems are whatever the access structure names, and
each subsystem's bath is, by definition, the commutant of that subsystem inside the whole
(FORMAL_CORE DEF 1.1–1.3; the twisted commutant for fermions). The phrase "bath = rest"
therefore adds a definition, not a property: the rest is not thereby thermal, Markovian,
infinite, weakly coupled, uncorrelated or resettable. Every one of those is a property a model
instance may or may not have.

An answer to the parent question in a sector has the form

> in model class M, closure together with the stated supplied data implies X,

with either a derivation or a countermodel. Three kinds of answer count:

1. **Forced by closure alone.** All such X found are standard theorems (FORMAL_CORE §3).
2. **Forced given a named hypothesis.** Closure plus a checkable property of the supplied datum
   implies X (FORMAL_CORE COND 2.10–2.20). The hypothesis, not closure, carries the content.
3. **Not forced.** A legal closed model in M lacks X. Most program campaigns built exactly such
   countermodels, or tested proposed extra postulates that did not survive.

"Nothing additional was forced in this class, and this computation establishes it" is a
complete result. It is not proof that no future postulate could have content; FORMAL_CORE §7
lists the kinds of postulate the program considered, slot by slot of Θ, with their costs, and
names known kinds outside that list.

## 2. The formal core in equations

Finite regulated case, ħ = 1, time independent H, supplied factorization ℋ_Ω = ℋ_A ⊗ ℋ_Ā and
H = H_A + H_Ā + H_{AĀ}:

    ρ_Ω(t) = e^{−iHt} ρ_Ω(0) e^{iHt},                       (global law; unitary)
    ρ_A(t) = Tr_Ā ρ_Ω(t),                                     (restriction; exists for every state)
    dρ_A/dt = −i[H_A, ρ_A] − i Tr_Ā [H_{AĀ}, ρ_Ω(t)].          (exact; the last term is global)

With 𝒫ρ = (Tr_Ā ρ) ⊗ ρ_ref and ℒ = −i[H, ·] the same equation reads

    𝒫ρ̇(t) = 𝒫ℒ𝒫ρ(t) + ∫₀ᵗ 𝒦(t−s) 𝒫ρ(s) ds + 𝒫ℒ e^{𝒬ℒt} 𝒬ρ(0),     𝒦(τ) = 𝒫ℒ e^{𝒬ℒτ} 𝒬ℒ𝒫,

so the memory kernel is a function of all of H, of ρ_ref and of the factorization. The
inhomogeneous term vanishes when ρ_Ω(0) = ρ_A(0) ⊗ ρ_ref, and also for some correlated
preparations, for instance a correlation that commutes with H (calc 1 S4 exhibits one). For a
product preparation the map ρ_A(0) ↦ ρ_A(t) is linear and completely positive but in general
not CP-divisible; for a correlated preparation there is in general no such map at all.

Records are registers R in the access structure. Interventions are of two kinds. Externally
timed ones are supplied controls that open the whole at a chosen time (calc 1 S6's write
changes the energy by 0.2076); their complete description is the process tensor, a multilinear
completely positive map from the controls to the outcomes, determined by Θ. Internal ones need a
clock and controller inside Ω, with energy, conservation-law and clock-imperfection costs. Two
standard identities hold in either case: an instrument on A and a unitary record into an idle,
traced-out register give the same later statistics; and the statistics after an intervention
are not functions of the one-time map.

For a factorized preparation the Feynman–Vernon influence functional
F[x, x′] = Tr_Ā[U_Ā[x] ρ_Ā U_Ā[x′]†] is the path-integral form of the same reduced dynamics;
correlated preparations add a preparation functional. If the complement is linear, linearly
coupled and prepared thermal, F is Gaussian and fixed by a dissipation kernel and a noise kernel
built from a spectral density J(ω). The closed-time-path generating functional is that object in
field-theory notation. A quadratic S_CTP therefore presupposes J(ω), the temperature, linearity
and the preparation, and earns nothing by itself (ledger s01).

Compatibility: restrictions of one global state agree on nested regions and on overlaps
(necessary, automatic). The converse fails: a family of local states that agrees on every
overlap need not come from any positive global state. This quantum marginal problem is the
one-time necessary condition for sectors to be compatible with one global model and
preparation. Full compatibility also requires the sectors' reduced processes to be marginals of
one global process, a stronger problem this packet does not compute.

Beyond the finite case the rules must be read algebraically: bath = commutant 𝔄(O)′, equal to
the causal-complement algebra only under Haag duality; restriction ω|_{𝔄(O)}; no partial trace
for type III algebras; modular KMS automatic and not a temperature in a physical clock unless the
flow is geometric. GRUT III's P1 is this identification and supplies neither dynamics nor a
state (ledger s05, s09). External time is an input throughout; nothing here touches the problem
of time, the gravitational constraints, or the absence of a conserved total energy in an
expanding spacetime.

## 3. One closed whole, from preparation to measurement

The supplied datum (CALCULATIONS §1): eight qubits; A = qubit 0 with H_A = ½Z₀; the rest a
disordered XXZ chain (J = 1, J_z = 0.5, fixed fields); coupling 0.5 X₀X₁; preparation
|+⟩⟨+|_A ⊗ e^{−H_Ā}/Z (a mixed global state); projective Z and X instruments on A; t ∈ [0, 40].

What the numbers show (all from `results/closed_whole_demo.json`):

| Question | Result |
|---|---|
| Global invariants | E_tot = −10.228789, S(ρ_Ω) = 1.260931, Tr ρ_Ω² = 0.489575, each constant to 10⁻¹⁴ |
| Local energy | E_A ranges over [−0.169, 0], E_int over [−0.223, 0]: not conserved separately |
| Born rule | local p(Z=+1), p(X=+1) equal the global expressions to 2·10⁻¹⁶ |
| Exact reduced equation | residual 1.9·10⁻¹⁰; the H_Ā commutator traces to 3·10⁻¹⁶ |
| Monotone relaxation? | No: S(ρ_A) has 26 sign changes of its derivative; information backflow 1.20 on the 0.1 grid (Breuer–Laine–Piilo) |
| Fixed Lindbladian? | No: the exact time-convolutionless generator changes by up to 0.97 along the trajectory and differs between thermal, ground and infinite-temperature baths |
| CP-divisible? | No: intermediate maps have Choi eigenvalues down to −0.159, which certifies memory |
| Does ρ_A(0) fix the future? | No: two preparations with identical ρ_A(0) and identical ρ_Ā(0) separate to trace distance 0.1025; the A-factor of the initial-correlation term is ¼Z₀ exactly |
| Is that term zero only for products? | No: a correlation that commutes with H (the conserved Z-parity) leaves the futures identical to 4·10⁻¹⁶ although the preparation is correlated |
| Unique local law? | No: same H_A and ρ_A(0) with three different complements, trace distances up to 0.79; a commensurate complement returns A exactly at T = 4π |
| Records | an externally timed write changes the energy by 0.2076; the register population then equals the Born probability 0.354910 (an identity); p(Z=+1) at t₂ is 0.3385 with the record, 0.4491 without it, and 0.2883 from composing one-time maps with the instrument |
| Record persistence | exact for [H, Z_R] = 0, but then H cannot write the register; a register rotation loses nothing (optimal readout 1.000); coupling the register to the bath lowers the optimal readout to 0.96–0.97 |
| Arrow | forward then backward evolution returns ρ_Ω(0) to 5·10⁻¹⁵; the conjugated final state is a legal preparation whose S(ρ_A) returns non-monotonically from 0.4198 to 0 |
| Second law given a Gibbs complement uncorrelated with A | Σ = ΔS_A + βΔE_Ā equals I + D to 1.1·10⁻¹⁴, is ≥ 0 at every time and decreases on 22 of 40 unit steps |
| Athermality budget | D(ρ_A(t)‖σ_A) never exceeds the conserved D(ρ_Ω‖σ) = 0.891 for the stationary global Gibbs state σ |
| Entropy theorems | pure global state: S(ρ_A) = S(ρ_Ā) to 3·10⁻¹⁶; mixed trajectory: mutual information ≥ 0 at all 401 times, Araki–Lieb slack 0.2115 at t = 3 and 0.2439 at t = 20 (0.20–0.40 over t ≥ 3), no violation |
| A state rule built from H alone? | the Gibbs and ground states of H commute with H to 2·10⁻¹⁴ and their one-time marginals stay fixed to below 1.3·10⁻¹⁵; multi-time and response physics remain |
| Closure plus generic gaps | equilibration on average: ⟨D(ρ_A(t), ω_A)⟩ ≤ ½√(d_A²/d_eff); bound 0.184 vs sampled 0.088 for d_eff = 29.6 (the d_eff = 3.2 case is outside the non-trivial regime); the gap condition holds only marginally (closest Bohr-gap coincidence 3·10⁻⁸) |

Every line is standard quantum mechanics. Computing them in one supplied model turns the parent
question's answers, negative and conditional, into exhibits rather than assertions.

## 4. Compatibility and underdetermination

**Compatibility** (CALCULATIONS §2, `results/marginal_compatibility.json`). Two Bell pairs Φ⁺ on
AB and on BC agree exactly on the shared qubit (ρ_B = 1/2 from both sides) and have no common
three-qubit state: the semidefinite program is infeasible with two solvers, and analytically a
pure ρ_AB forces ρ_ABC = ρ_AB ⊗ ρ_C. The obstruction is quantitative: the smallest achievable
Frobenius violation of the two marginals is 0.289; the largest achievable F_AB + F_BC is 3/2 (it
would be 2 if compatible); isotropic states W(p) = pΦ⁺ + (1−p)1/4 on AB and BC are jointly
realizable only for p ≤ 2/3, and on all three pairs of a triangle only for p ≤ 1/3, the
separability threshold. Classically, two marginals sharing one variable always extend, and an
odd cycle of perfect anticorrelations does not (an even one does); the quantum failure on a chain
is the specifically quantum constraint (monogamy). Consequence for GRUT: effective environments
asserted sector by sector are not automatically restrictions of anything; the family must pass
this one-time test, and the test is standard.

**Underdetermination.** The local data underdetermine the global datum in three computed ways.
(i) ρ_A(0), even with ρ_Ā(0), does not determine ρ_A(t) (calc 1, S4). (ii) The same local
Hamiltonian and preparation admit different complements with different local laws (S5).
(iii) Two-point influence data do not determine higher influence data: the program's own
dependence campaigns found exactly matched two-point baths that split a probe coherence by 0.338
through the fourth cumulant (ledger s02). So no amount of local response data selects the
closed whole; the global datum is supplied, and the "response" vocabulary adds no content beyond
the restriction rules (ledger s01, s02).

## 5. The campaign ledger

Each row is one campaign; the full line-cited entry, with its assumptions, executed evidence,
retained failures and sources, is in [CAMPAIGN_LEDGER.md](CAMPAIGN_LEDGER.md) and `ledger/`.
Every entry was written from the primary sources by one reader. Entries s01–s10 were then re-read
adversarially by a second reader (6 to 15 corrections each, folded into the entries). Entries s11
and s12 are single-reader; the claims below drawn from s10–s12 were also checked against the
sources in the main session, and the round-2 attack corrected further errors in them
(`results/verification/s10_s12_row_checks.md`).

| # | Campaign | Question it asked | Answer | Classification |
|---|---|---|---|---|
| s01 | V7 / CTP organization | Does CTP doubling, normalization or the response vocabulary derive dynamics, Born probabilities or the asserted predictions? | No. The Feynman–Vernon / Keldysh form is generic ("any local causal open system yields it"); in-in does no work at the gravity rung (the Ward identity allows EH, R², f(R), Lovelock); the declared action "fixes nothing" in its own family; the Born measure is inherited; the differentiators fell 7–47 (decoherence), 22–62 (GW dephasing only; the amplitude channel is not settled) and 62.7 (Γ_T horn (a); horn (b) not banked) orders below their respective bounds; the bath behind the in-in action was never specified. The V7 file is absent from the repository; its content is known only through the handoff. | mixed: recovered physics on supplied choices; rejected mechanisms retained |
| s02 | Influence hierarchy and dependence | Are kernel and noise downstream of subsystem, state and access? Do two-point data fix higher data? Is positivity a new principle? | Yes for the influence functional and the noise kernel (the retarded kernel was state-independent in the Gaussian and free classes tested); no; no. Matched two-point baths split a probe coherence by 0.338 at fourth order; discrimination persists at order 6, and at order 8 only by a labelled post-hoc diagnostic (frozen gate red); every positivity constraint reduced to state positivity (NULL-AS-NEW-PRINCIPLE). In-class, access is not arbitrary bookkeeping (P-5); the owner ruling calls the access boundary "an unresolved structural input". | mixed: dominantly recovered physics, bounded by in-class conditional constraints |
| s03 | Level-0 floor and earned theorems | Do the floor theorems constrain a closed whole beyond standard quantum theory? | No. Every theorem is about a supplied classical class (rings, OU chains, declared readouts); the lift to quantum theory is "IRREDUCIBLE/SUPPLIED". The one closed-system item, O-6 (conservative chain with internal bath), was ruled FALSIFIED for a strict arrow from identities and band-edge asymptotics checked on a non-member chain, not from recurrence; the chartered run on the declared chain was never executed; S6-1 later found a preparation-relative net arrow on the parent. | mixed: recovered mathematics; conditional constraints; three hypotheses falsified |
| s04 | GRUT I–II vacuum response | What was recovered with imports, what was a projector choice, which exclusions are conditional, is the arrow supplied? | Derived count 0 across 26 V1 claims and V2. TT-only projector CHOSEN: GR's own kernel has P0s/P2 = −2, and the two-Ward-survivor classification holds only under the advanced-branch identity; it was retired as unconditional on 2026-08-14, enlarging the admissible space and adding noise-kernel transversality as a +1 input. The four boundaries hold at stated scopes; the arrow's direction is imported through the KMS gate; Γ_T horn (a) is 62.7 orders below its bound; the ISW anchor retired 32σ → 2.0σ. | mixed: recovered physics; supplied choices; conditional exclusions; many retractions retained |
| s05 | GRUT III observer bath | Is P1 a law? Is modular KMS a temperature? What is imported in GR recovery and universal coupling? | P1 is an identification with zero credit, no dynamics, no state; modular KMS is a physical temperature only where the flow is geometric (imports); GR recovery imports the area law, Jacobson and η = 1/(4Għ); universal coupling inserts κ; "GRUT III makes no prediction that distinguishes it from GR + QFT in curved spacetime". Slow-roll thermality (D11, D13) is conditional. In FRW, besides a detector temperature in the Wigner convention, a near-horizon Kodama–Hayward KMS temperature is proved for trapping horizons (cosmological case open; its relation to P1 unreconciled); the secular-infrared survey holds two unreconciled graviton verdicts. The pre-registered decisive tests N1 and N2 are unexecuted. | mixed: dominated by recovered physics; five postulates supplied |
| s06 | CP-12 geometry and preparation | Which geometry-selection obstructions are model- and preparation-specific? What does the min-entropy theorem assume? | All obstructions are model- and preparation-specific. The min-entropy theorem needs ‖ρ₀‖_∞ ≤ e^{cN}/D and O(N) excess energy; closure is load-bearing through unitary invariance of the spectrum and energy conservation, converting the question into an initial-state obligation of αN log₂N bits, "not a state-selection law". The bound is vacuous at every size its checker prints, through N = 10¹⁸, and its parent relative-entropy theorem has no partner audit. Pure coherent preparations are outside it, escaping "without any proof that it assembles geometry"; formation is not established; no theorem excludes every emergent geometry. | mixed: a surviving conditional constraint (unaudited); rejected local-geometry mechanisms |
| s07 | WM / B1 / B2 matter and moving channels | Keep acceptance, recoverability, delivery and retention distinct; is B2 an impossibility proof? | Conditional externally powered preparation (P1); autonomous every-time retention refuted at one permitted point (A1); a conditional delivery window (B1, whose two-event control fails); a fixed-path two-source pass that does not transfer to the moving law; B2's eight-event endpoint failed (1.93%). The closeout is a method and resource decision: no phase-sensitive estimate is available and every certificate is vacuous or out of envelope at the required scale (1411/640 majorant, N ≥ 2²⁸); moving operation is "neither established nor refuted". All results are author-derived with independent review pending. Ledger synthesis: closure raised the proof and resource cost and supplied no principle. | mixed: recovered physics; conditional windows; one refutation; one feasibility stop |
| s08 | Candidate 1 scalar EFT | What is the action, what is supplied, what does the one-loop calculation show? | Einstein–Hilbert + SM + three singlet scalars; everything structural supplied; a classically supported charged profile at Q ≈ 39 944 (with an adverse branch and a second higher state). One-loop shift K₁/(gv) = 2.795 with a ±0.11 numerical allowance (analytic tail bounds plus an empirical allowance; not a certified enclosure, not a statistical interval, not a higher-loop bound; its angular term is root-derived and was not challenger-audited); tree-plus-one-loop margin −3553.71 gv for the soliton against N = 39 944 free carriers, conditional on reduced-saddle stability. The tree-level single-carrier and equal-split gaps −0.199999152115 and −3067.397336565 gv come from the closure-fork packet. No lifetime, no quantum fission, no gravity matching, no novelty (FLS-type soliton, Graham's spectral method). | mixed: conventional EFT; conditional one-loop constraint; supplied couplings |
| s09 | Closure / no-boundary | Does C2 supply a state-selection rule? Is the slow-roll coefficient new? Is curvature sign a discriminator? | No. The restriction map exists for every state and modular KMS is automatic, so observer restriction selects nothing; the Ivo–Li–Maldacena trace-matching conditions are already in the published no-boundary density-matrix prescription. T_eff = H(1−ε)/2π is the first-order term of Obadia's 2008 coefficient (D2 withdrawn). The fork's own report: closure "alone fixes neither the observer's restriction nor its thermal spectrum", "not a theorem that every closed universe is thermal"; the present universe is not slow-roll (ε₀ ≈ 0.47). Compactness does not imply k = +1 (flat torus); closed de Sitter gives Ω_K = −1/(e^{2N} − 1), −4.5·10⁻⁵ at N = 5 and −2.1·10⁻⁹ at N = 10, so the sign alone bounds nothing, although a model can be disfavored by its curvature distribution. C1 defined as the FLT prescription would end on the inverse-Gaussian coefficient the FLT literature reports, quoted and not reproduced here. The B three-prescription control is registered and unexecuted (its files entered the tree in the M1 commit, so the registration order is self-declared). Program decision: option C, preserve the map. | mixed: recovered physics; supplied prescription choices; rejected falsifiers; one unexecuted control |
| s10 | M1 memory | What works in the stated finite regime, what was withdrawn, what is not established? | In the stated finite regime the repetition register retains the bit (L = 4: min 0.987510 through sampled t = 10, 0.935065 at t = 40; matter readout 0.778468), and the flying-qubit writer copies a source bit at ≥ 0.999680, the write releasing register energy (1.55 → −0.44) into the carrier (kinetic energy +1.99). The matched F = ε control, specified after the review, is no worse at every sampled time (equal at t = 0, strictly better at every t > 0; L = 4, t = 40: 0.998458), so the facilitation-benefit interpretation is withdrawn. The independent review reads the retention as an ordinary gapped repetition register, a reading no gap calculation tested. Exact domain-wall conservation at ε = 0 is a model-specific protect/form constraint. Not established: a thermal memory phase, cloning of an unknown state, gravity, a chaos exponent, all-times retention, radiative stability. | mixed: recovered physics; rejected mechanism (facilitation benefit); conditional constraints at L = 3–5 on sampled grids |
| s11 | Toric memory / scrambling bridge (separate branch) | What was computed, what was corrected post hoc, what is it not? | An open-system computation: single-face Metropolis dynamics of the Z-sector 4D toric code at β ∈ {1.0, 1.1}, L ∈ {3, 4, 5}, local cooling decoder, t_max = 10⁵ sweeps. The logical half-life grows with L (977 → 1.6·10⁴ → 8.0·10⁴ sweeps at β = 1.0) while stabiliser observables decorrelate in about one sweep; the face spin was reclassified post hoc as gauge-variant, so by the letter of the pre-registration R2 fails and the fallback rests on the secondary damage observable; the half-life estimator was changed after an unpreserved first attempt; the 18-qubit OTOC illustration "does not bear on B". As reported it recovers the known thermal stability of the 4D toric code and decides the bridge claim against the bridge in that classical class. It is not a closed-unitary gravitational scrambling computation; its de Sitter and holography statements are uncomputed literature claims. Separate branch, not merged, independent review outstanding. | mixed: recovered physics; conditional rejected mechanism; conjectures not computed; unexecuted quantum OTOC |
| s12 | Program no-go ledger, postulate map, endogenous access | Which recorded no-gos bear on "closure forces X"? Is access ever derived? | The program's ledgers: a register of 74 nodes (53 GRUT and 21 vacuum-cluster at the 2026-08-17 sync stamp; total recounted here) banks no structural-impossibility (FORBIDDEN) no-go, its strongest exports being settled-negative, empirically excluded, suppressed or borrowed; at that stamp net +16 underived inputs, an empty derived tier and an empty prediction ledger; every kernel alternative "excluded by an input, never by a principle"; the arrow's direction, the Born measure, the system/bath split, ħ, c, G, the action, the projector, the generator, the lift and the access seed recorded as supplied. The record called the access seed "currently irreducible" the day before EA-0 was ruled, and EA-0 left it so: UNFORMULABLE at the earned Level-0 scope because the earned substrate is classical, with the bridge gate ending NONUNIQUE-LIFT. Three negatives are flagged MUST-BE-REOPENED. Pure global state, compact topology and finite four-volume are not addressed in these sources. | mixed: dominantly supplied choices; rejected mechanisms; unexecuted tasks; conditional in-class constraints |

## 6. Proposed additions that did not survive, surviving conditional results, and what remains open

Grouped by the slot of Θ each addition would have filled (FORMAL_CORE §7). This is a summary;
each ledger entry's "Failures retained" section is the complete list.

**State selection.**
- C2 identifies the observer-state object and adds no selection condition beyond restriction,
  trace matching, positivity, normalization and modular KMS; the no-boundary density-matrix
  construction already contains the trace-matching conditions (s09).
- D2 "slow-roll thermality" novelty withdrawn (Obadia 2008 is the prior source for the
  first-order coefficient); the sign-only curvature falsifier and "all sequestering
  formulations force positive curvature" withdrawn; the coarse-ε criterion failed, and only a
  changed criterion passes (s09).
- The sharp Sorkin–Johnston state is rejected; "GRUT derives the initial state" was withdrawn
  in the GRUT III synthesis's own review; "P1 selects the horizon state" remains a conjecture
  whose pre-registered decisive test (N1) never ran (s05).

**Dynamics.**
- The Schwinger–Keldysh formalism does not select the Einstein–Hilbert action; the in-in
  machinery does no work at the gravity rung, which reverts to Jacobson with G inserted
  (s01, s04, s05). One-loop induced gravity with Standard-Model content does not supply G (s05).
- The TT-only projector of GRUT I–II is a choice, and its two-survivor classification was
  retired as unconditional (s04). The x < 0.358 bound is demoted; the κ pre-registration
  self-unblinds; the mu_linear boundary as banked is weaker than stated (s04).
- a₀ via Jacobson, a memory time 1/(4H), a bulk Planckian fluid, c_eff at LIGO frequencies and
  w ≡ −1 as a GRUT discriminator are rejected (s05).
- The 689 Hz falsifier is withdrawn (s01); the ISW anchor retired from 32σ to 2.0σ (s04).
- The founding ζ-hypothesis, the conformalon, the information principle and the
  single-particle interference memory mechanism are rejected. Single-pole finite memory: its
  in-house support was reversed (exact de Sitter forces infinite scale-free memory), it was
  retired in GRUT III, and it is UNRESOLVED in the final synthesis (s01, s04, s05).
- Candidate 1 is a quartic extension of Friedberg–Lee–Sirlin with supplied couplings; its
  "coefficients vanish at every scale" rule and a pointwise Coleman–Weinberg shortcut are
  rejected; its acceptance thresholds were consolidated after execution (s08).
- The classical relaxational carrier is eliminated inside the supplied quantum class (CA-1);
  H-HERM-2, strict H-ORD and the two-property reduction of the floor are falsified; P-6's L-C
  gate stays red and S2-1's first execution is void (s02, s03).
- The bath behind the in-in action was never specified: the noise-origin campaign closed with
  NO-PHYSICAL-BATH-CANDIDATE, and GRUT III retired the finite-memory ontology ("memory is a
  property of (observable, state)") (s01, s02).
- The facilitation benefit of M1 is withdrawn: the matched control is no worse at every sampled
  time (equal at t = 0, strictly better at every t > 0) for L = 3, 4, 5, a post-review control;
  a single cell fails the 0.9 standard at L = 4, 5; bare matter-to-memory writing fails;
  unconstrained flips at 0.1 and 0.3 violate retention (s10).

**Access structure.**
- An endogenous access structure: EA-0 terminated UNFORMULABLE at the earned Level-0 scope,
  because the earned substrate is classical and the construction is undefined there without an
  unearned lift. Separately, across the ten quantum candidates examined, each seedless candidate
  was global, generically trivial (proved at d = 4 and 8), a correlation readout, or
  generator-carried by identity; this is a finding over the examined set, not a no-go. The
  bridge gate ended NONUNIQUE-LIFT; no earned Level-0 structure selects the physical lift (s03,
  s12). Locality-from-spectrum was cited and ruled CRITERION-SMUGGLED, not tested (s12).

**Measurement rule.**
- Decoherence does not yield outcomes or Born weights in the tested branch-preserving unitary
  class; P5 remains a postulate; the HJW route was withdrawn; a positivity "principle" reduces to
  state positivity (NULL-AS-NEW-PRINCIPLE) (s01, s02, s05).

**Clock and arrow.**
- The arrow's direction is imported with the passive, uncorrelated past-boundary state; the
  exactly solvable model has Γ(t) = Γ(−t) and full recurrence (s01, s04). The temporal generator
  is supplied, S5-0 GENERATOR-IRREDUCIBLE (s01, s03, s12). P1 assigns no temperature in FRW with
  ε ≠ 0, and modular time equals physical time only in the imported geometric cases (s05).

**Mechanisms.**
- Autonomous every-time retention is refuted at one permitted parameter point (A1); B1's
  two-event control fails; B2's eight-event endpoint failed; the fixed-path pass does not
  transfer to the moving law (s07). CP-12 formation is not established and its Milestone A is
  not achieved (s06).
- The toric report compares thermal Metropolis trajectories and decoded memory with a post-hoc
  gauge-observable correction; it is not a closed-unitary gravitational scrambling-exponent
  computation, and its rejection of the bridge claim is as reported on the unmerged branch, not
  re-run here (s11).

**Surviving conditional results** (they hold in their stated classes; none is a new law):
- The CP-12 min-entropy obligation, in which closure is load-bearing (unitary invariance of the
  spectrum and energy conservation): an initial-state requirement for loop-rich geometry in the
  frozen model; vacuous at the sizes printed and resting on an unaudited parent theorem (s06).
- M1's protect/form constraint at ε = 0: exact domain-wall conservation protects a record and
  forbids its formation, an instance of FORMAL_CORE THM 2.15 (s10).
- Access is not arbitrary bookkeeping in the P-5 class: moving the access boundary reveals
  physics that identical in-access data hide (s02).
- The influence-hierarchy cumulant ladder on the constructed pairs, with its order-8 gate red
  (s02); the Level-0 in-class non-implications (s03); the GRUT I–II boundaries at their stated
  scopes (s04); slow-roll thermality D11/D13 in its stated class (s05); Candidate 1's one-loop
  margin at its benchmark (s08).

**Still open (registered, not executed).** The B three-prescription control (FLT, whole-line
Hartle–Hawking, Di Tucci–Lehners Robin; Yamada as a distinct prescription), including the
reproduction of FLT's inverse-Gaussian coefficient on which a C1-as-FLT termination would rest
(s09); GRUT III N1 (horizon-state selection by P1, pre-registered as decisive) and N2 (a horizon
principle fixing a renormalisation constant) (s05); Candidate 1's three-charge one-loop family,
lifetime and nonlinear gravity (s08); the O-6 declared-chain run (s03); the rung-3
pole-versus-cut ask, drafted and never sent (s01); independent review of CP-04 and CP-05 (s05);
the partner audit of the CP-12 parent theorem (s06); independent review of the B2-RA pass and of
B1's formation/use composition (s07); M1's error-correction cascade against the matched baseline
(s10); the toric branch's merge decision and its quantum OTOC (s11); the universality probe, the
criterion-unification fork, EA-1 and S-5 (s12); the three MUST-BE-REOPENED negatives (the ROOT-1
battery's hard-coded gate, the composite no-pole negative, the tt_worldline decay) (s12); the two
unreconciled graviton verdicts of the GRUT III secular-infrared survey (s05). One external signal
bears on a GRUT I–II item: the DESI DR2 + CMB preference (3.1σ) for a w(z) crossing the phantom
divide at z ≈ 0.35–0.5, against the held-open, rung-3-gated no-crossing of a single passive
relaxor, which is not a banked boundary (s04).

## 7. Physical scope, predictions and open claims

**Has any additional law or state-selection rule survived?** No. State selection, dynamical
constraints, access, the clock and the measurement rule were tested by campaigns; the algebra
and instrument slots were examined; a modified global law, a measure over Θ and relational
constraints between slots were never tested (FORMAL_CORE §7). The surviving GRUT objects are the
question; the formal core, which is standard; the compatibility and underdetermination
computations, which are standard; the conditional results of §6; and the conditional model
instances, each with its own Θ. They do not combine into one law. Candidate 1 (a conventional
Einstein plus Standard Model plus scalar effective field theory), M1 (a finite Ising repetition
register with matter feedback and a programmed writer) and GRUT III (five postulates with
imports) live in different state spaces, and no common action, state space, constraint set or
calibration has been established; the handoff's no-cross-import rule stands.

**Predictions, conditional on inputs.** Within any supplied Θ the framework predicts Born
probabilities and multi-time statistics; calculation 1 is such a prediction set, calculation 2
is a set of compatibility constraints, and COND 2.10 bounds time-averaged local deviations once
d_eff is known. Candidate 1 predicts, for its supplied couplings, a one-loop shift K₁/(gv) = 2.795
with a ±0.11 numerical allowance (not a certified enclosure, not a statistical interval) and a
tree-plus-one-loop margin of −3553.71 gv for its soliton against free carriers in the benchmark
channel, conditional on reduced-saddle stability (s08). M1 predicts, for its supplied sizes and
fields, retention probabilities that the matched uniform-floor control equals at t = 0 and
exceeds at every sampled t > 0 (s10; the control script re-run in this packet, line for line).
GRUT III predicts geometric-horizon temperatures where the modular flow is geometric, all
imported (s05). None of these distinguishes its instance from the known-physics baseline, and
none is a prediction of the closed-whole framework itself.

**What extra condition would make a distinct prediction.** The program already states the form
it must take (GRUT III T1–T2, FORMAL_CORE COND 2.18): for a stated observable, state and
preparation, either a different dissipative response, a different state, or different higher
cumulants than the conventional theory. Which slot of Θ supplies that difference is then the
question. For the state: a unitarily equivariant rule built from H alone is stationary (calc 1
S9), so one-time local evolution needs a rule that breaks equivariance; in locally covariant QFT
an embedding-covariant rule is excluded under the hypotheses of Fewster–Verch Theorem 6.13
(FORMAL_CORE COND 5.3); a bottom-up family of local states must pass the marginal problem and,
beyond it, a process marginal problem. For the dynamics: a constraint on H with recovered known
physics and its own constants. For access: an endogenous structure that does not presuppose the
factorization it derives. For the clock: an identification of modular time with physical time
beyond the geometric cases. Or a measure over Θ, or a kind outside the program's list (two-time
boundary conditions, self-location measures, changed composition rules, state-dependent
dynamics). The program's record of attempts is the ledger; among the state-selection attempts,
none produced a condition that is not already restriction, trace matching, positivity,
normalization or modular KMS. The pre-registered N1 test is the program's own open attempt.

**Claim boundaries.** The finite statements are the regulated model's; the algebraic statements
are identifications; the gravitational and cosmological statements are scope limits, not
results. Historical bounds and ledger counts are not current predictions. The toric memory
report on its separate branch is not merged and is not a gravitational scrambling computation
(s11). The B three-prescription control remains registered and unexecuted (s09). The
facilitation-benefit claim of M1 is withdrawn (s10). Interventions in calculation 1 are
externally timed; an internal writer was not modelled here.

## 8. Verification record

- `python3 code/verify.py` checks the hash of every file listed in MANIFEST.json, 31 headline
  results and 2 provenance hashes; `--rerun` re-executes both calculations and requires identical
  leaf paths and agreement on every numeric leaf.
- A separate agent re-implemented the calculation-1 quantities and the calculation-2 chain
  threshold from the specification alone: 23 rows, 21 agreeing to 10⁻⁶ relative and the two
  numerical zeros agreeing absolutely; it was re-compared against the committed results after
  the final edits. A separate re-execution verified the pre-review version. Two adversarial code
  reviews found nine and seven items, all applied. Records: `results/verification/`, summarized
  in [CALCULATIONS.md](CALCULATIONS.md) §Verification.
- The M1 matched uniform-field control script was re-run unchanged and reproduced the recorded
  table line for line (`results/m1_matched_control_rerun.txt`).
- Ledger: s01–s10 adversarially re-read by a second reader (6 to 15 corrections each, folded into
  the entries); s11–s12 single-reader; claims drawn from s10–s12 checked in the main session
  (`results/verification/s10_s12_row_checks.md`). The completeness critic over the ledger did not
  run (session limits).
- Primary literature: see SOURCES.md (Fewster–Verch, BFV, Masanes–Galley–Müller, Witten, Pollock
  et al. PRL and PRA, Chen et al. read in full; Cotler–Penington–Ranard in full by an attack
  agent; Feynman–Vernon not accessed).
- Attack on the packet by separate agents in two rounds: round 1 (equations, handoff compliance,
  number provenance) found 2 blocking, 11 major and 22 minor items; round 2 (overclaim, steelman,
  equations, number provenance) found 0 blocking, 30 major and 50 minor items. All were
  addressed in this version; records in `results/verification/attack_round1/` and
  `attack_round2/`. No third round was run on this final version.
- Not done: no old solver reopened beyond the M1 control re-run; no cosmology contour executed;
  no toric code run; no determinant recomputed; no field added to M1; no new branch; no external
  contact.
