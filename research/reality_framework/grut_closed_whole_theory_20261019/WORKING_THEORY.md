> **DRAFT CHECKPOINT (19 Oct 2026, in progress).** Sections marked `[[...]]` await the campaign-ledger
> workflow; the calculations and formal core are complete. Do not cite this file until the banner is removed.

# GRUT working theory: the closed whole

D. Ryan Grover's program, synthesis of 19 October 2026. Packet `grut_closed_whole_theory_20261019`,
branch `codex/grut-closure-fork-20261019`, written from the handoff at commit `96a7e7e` and the
corrected M1 record at `4e9a161`. Detailed derivations are in [FORMAL_CORE.md](FORMAL_CORE.md),
every number in [CALCULATIONS.md](CALCULATIONS.md), every sector in
[CAMPAIGN_LEDGER.md](CAMPAIGN_LEDGER.md) and `ledger/`.

## 0. What GRUT now says, in one page

**The parent object is a question.** What does a closed system, whose every subsystem's bath is
the rest of the universe, force that known physics does not already contain? GRUT is the
program that asks this question in one place after another and reports the answer, including
the answer "nothing additional in this class, and here is the computation". The old parent
object, an in-in action S_CTP treated as a "genotype", is the Gaussian approximation to a
derived object (the Feynman–Vernon influence functional) and sits downstream of a supplied
state, dynamics and access structure. It is a calculational tool, not a principle.

**What calculates.** A closed whole is a theory datum
Θ = (algebra, dynamics, preparation, access structure, instruments and records, clock), all
supplied. Given Θ, the global law ρ_Ω(t) = U(t)ρ_Ω(0)U(t)†, the restriction ρ_A = Tr_Ā ρ_Ω and the
Born rule produce every prediction, conditional on Θ. Packet calculation 1 carries one fully
specified eight-qubit Θ from preparation through evolution to measurement with exact energy and
correlation accounting; calculation 2 computes the compatibility conditions that overlapping
descriptions must satisfy. Both are reproducible in about a minute.

**What closure constrains.** From closure alone, exactly the theorems of standard quantum
theory: unitary invariants, the exact reduced equation dρ_A/dt = −i[H_A, ρ_A] − i Tr_Ā[H_{AĀ}, ρ_Ω],
the Nakajima–Zwanzig structure with its initial-correlation term, complete positivity of
reduced maps for product preparations, nested compatibility, necessary overlap compatibility,
the entropy inequalities, quasi-recurrence in finite dimension. From closure plus one named
hypothesis on the supplied datum, four conditional theorems: equilibration **on average** when
the energy gaps are non-degenerate, with its size set by the preparation; canonical typicality
for a stated measure over preparations; effective light cones for a local H; and
conservation-law limits on records made inside the whole. Closure forces **nothing local
pointwise**: not monotone relaxation, not Markovianity, not a kernel, not a unique local law,
not thermality of the complement, not an arrow, not record persistence, not a preferred state,
not the constants, not the action. Each negative is exhibited by a legal closed model in
calculation 1 (§3 below), a state rule built from H alone is shown to be dynamically empty, and
the compatibility requirement has a quantitative edge in calculation 2 (§4).

**Which proposed additions did not survive.** Every addition that a campaign proposed and tested failed to survive as a law. The
CTP action as a parent principle; the in-in selection of the Einstein–Hilbert action; the TT
projector as forced; decoherence as a derivation of the Born weights and the HJW route to the
outcome rule; influence-data positivity as a new principle; any GRUT-internal kernel-selection
principle; the facilitation benefit of the M1 memory; C1-as-FLT and C2 as state selection; the
novelty of slow-roll thermality; the curvature-sign falsifier; the CP-12 min-entropy theorem
as a selection law; the 689 Hz falsifier and the 32σ ISW anchor; the founding ζ-hypothesis,
the conformalon and the information principle; the self-correction/no-scrambling bridge in the
classical 4D toric code (as reported on an unmerged branch, not re-run here, independent review
outstanding); autonomous every-time retention; a GRUT-specific prediction of GRUT
III. Section 6 lists each with its sector; the ledger retains every failure.

**What remains a supplied choice.** The six entries of Θ, in every sector: the field content
and action (GRUT III P0 and data; Candidate 1's couplings), the preparation (GRUT I–II
projector choice; CP-12 preparation classes; the no-boundary state that C1/C2 did not select),
the access structure (EA0; the register and writer in M1), the measurement rule (P5), the
arrow (P4), and the time parameter. The program's own register prices this at +16 net
supplied inputs with an empty derived tier and an empty prediction ledger (s12); the access
structure in particular is "currently irreducible" after the EA-0 attempt to make it endogenous.

**Status.** A working closed-system framework of known physics with conditional model
instances and computed boundaries. No additional law or state-selection rule has survived.
Nothing here is validated universal physics, and neither gravity nor "all reality" has been
derived.

## 1. The parent question and what counts as an answer

A closed whole has no exterior. Its subsystems are whatever the access structure names, and
each subsystem's bath is, by definition, the commutant of that subsystem inside the whole
(FORMAL_CORE DEF 1.1–1.3). The phrase "bath = rest" therefore adds a definition, not a
property: the rest is not thereby thermal, Markovian, infinite, weakly coupled, uncorrelated
or resettable. Every one of those is a property a model instance may or may not have.

An answer to the parent question in a sector has the form

> in model class M, closure together with the stated supplied data implies X,

with either a derivation or a countermodel. Two kinds of answer count:

1. **Forced.** X follows from closure and the supplied data in class M. All such X found so far
   are standard theorems (FORMAL_CORE §3).
2. **Not forced.** A legal closed model in M lacks X. The program's campaigns are, for the most
   part, careful constructions of exactly such countermodels, plus the record of proposed
   extra postulates that were tested and did not survive.

"Nothing additional was forced in this class, and this computation establishes it" is a
complete result. It is not proof that no future postulate could have content; FORMAL_CORE §7
lists, slot by slot of the theory datum, what such a postulate would have to constrain and what
it would cost.

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
preparations, for instance a correlation that commutes with H (calc 1 S4 exhibits one); it is
generically nonzero otherwise. For a product preparation the map
ρ_A(0) ↦ ρ_A(t) is linear and completely positive but in general not CP-divisible; for a
correlated preparation there is in general no such map at all. The complete multi-time object
is the process tensor, a multilinear completely positive map from interventions to outcomes,
itself determined by Θ. Records are registers R in the access structure written by unitaries;
an instrument on A and a unitary record into an idle, traced-out register give the same later
statistics (FORMAL_CORE THM 2.8).

For a factorized preparation the Feynman–Vernon influence functional
F[x, x′] = Tr_Ā[U_Ā[x] ρ_Ā U_Ā[x′]†] is the path-integral form of the same reduced dynamics;
correlated preparations add a preparation functional. F is Gaussian, and fixed by a
dissipation kernel and a noise kernel built from a spectral density J(ω), only when the
complement is linear, linearly coupled and Gaussian-prepared. The closed-time-path generating
functional is that object in field-theory notation. A quadratic S_CTP therefore presupposes
J(ω), the temperature, linearity and the preparation. That is where the parent object that the
handoff attributes to GRUT V7 lives (the V7 file is not in the repository), and it earns
nothing by itself (ledger s01).

Compatibility: restrictions of one global state agree on nested regions and on overlaps
(necessary, automatic). The converse fails: a family of local states that agrees on every
overlap need not come from any positive global state. That is the quantum marginal problem,
and it is the exact form of the requirement that every sector's effective environment be
compatible with one global model and preparation.

Beyond the finite case the rules must be read algebraically: bath = commutant 𝔄(O)′, equal to
the causal-complement algebra only under Haag duality; restriction ω|_{𝔄(O)}; no partial trace
for type III algebras; modular KMS automatic and not a temperature in a physical clock unless
the flow is geometric. GRUT III's P1 is this identification and supplies neither dynamics nor
a state (ledger s05, s09). External time is an input throughout; nothing here touches the
problem of time, the gravitational constraints, or the absence of a conserved total energy in
an expanding spacetime.

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
| Monotone relaxation? | No: S(ρ_A) has 26 sign changes of its derivative; information backflow 1.20 (Breuer–Laine–Piilo) |
| Fixed Lindbladian? | No: the exact time-convolutionless generator changes by up to 0.97 along the trajectory and differs between thermal, ground and infinite-temperature baths |
| CP-divisible? | No: intermediate maps have Choi eigenvalues down to −0.159 |
| Does ρ_A(0) fix the future? | No: two preparations with identical ρ_A(0) and identical ρ_Ā(0) separate to trace distance 0.1025; the initial-correlation term is ¼Z₀ exactly |
| Unique local law? | No: same H_A and ρ_A(0) with three different complements, trace distances up to 0.79; a commensurate complement returns A exactly at T = 4π |
| Records | The register population after a unitary write equals the Born probability 0.354910 (an identity); p(Z=+1) at t₂ is 0.3385 exactly with the record, 0.4491 without it, and 0.2883 from composing one-time maps with the instrument |
| Record persistence | exact for [H, Z_R] = 0; a register rotation loses nothing (optimal readout 1.000); coupling the register to the bath lowers the optimal readout to 0.96–0.97 |
| Arrow | forward then backward evolution returns ρ_Ω(0) to 5·10⁻¹⁵; the conjugated final state is a legal preparation whose S(ρ_A) returns non-monotonically from 0.4198 to 0 |
| Entropy theorems | pure global state: S(ρ_A) = S(ρ_Ā) to 3·10⁻¹⁶; mixed trajectory: mutual information ≥ 0 at all 401 times, Araki–Lieb slack 0.2115 at t = 3 and 0.2439 at t = 20 (0.20–0.40 over t ≥ 3), no violation |
| Is the initial-correlation term zero only for product states? | No: a correlation that commutes with H (the conserved Z-parity) leaves the futures identical to 4·10⁻¹⁶ although the preparation is correlated |
| Can a state rule built from H alone select local physics? | No: the Gibbs state and the ground state of H commute with H to 10⁻¹⁴ and freeze every subsystem (ρ_A(t) − ρ_A(0) below 1.2·10⁻¹⁵); any such equivariant rule is stationary |
| What does closure plus generic gaps force? | Equilibration on average: ⟨D(ρ_A(t), ω_A)⟩ ≤ ½√(d_A²/d_eff); bound 0.556 vs sampled 0.287 for d_eff = 3.2, bound 0.184 vs sampled 0.088 for d_eff = 29.6; the gap condition holds only marginally (closest Bohr-gap coincidence 3·10⁻⁸) |

Every line is standard quantum mechanics. The point of computing them in one supplied model is
that the parent question's answers, negative and conditional, become exhibits rather than
assertions.

## 4. Compatibility and underdetermination

**Compatibility** (CALCULATIONS §2, `results/marginal_compatibility.json`). Two Bell pairs
Φ⁺ on AB and on BC agree exactly on the shared qubit (ρ_B = 1/2 from both sides) and have no
common three-qubit state: the semidefinite program is infeasible, and analytically a pure
ρ_AB forces ρ_ABC = ρ_AB ⊗ ρ_C. The obstruction is quantitative: the smallest achievable
Frobenius violation of the two marginals is 0.289; the largest achievable F_AB + F_BC is 3/2
(it would be 2 if compatible); isotropic states W(p) = pΦ⁺ + (1−p)1/4 on AB and BC are jointly
realizable only for p ≤ 2/3, and on all three pairs of a triangle only for p ≤ 1/3, the
separability threshold. Classically, two marginals sharing one variable always extend, and
an odd cycle of perfect anticorrelations does not (an even one does); the quantum failure on a
chain is the specifically quantum constraint (monogamy). Consequence for GRUT: an "effective environment"
asserted sector by sector is not automatically a restriction of anything; the family must pass
this test, and the test is standard.

**Underdetermination.** The local data underdetermine the global datum in three computed ways.
(i) ρ_A(0), even with ρ_Ā(0), does not determine ρ_A(t) (calc 1, S4). (ii) The same local
Hamiltonian and preparation admit different complements with different local laws (S5).
(iii) Two-point influence data do not determine higher influence data: the program's own
dependence campaigns found exactly matched two-point baths that split a probe coherence by
0.338 through the fourth cumulant (ledger s02). So no amount of local response data selects
the closed whole; the global datum is supplied, and the "response" vocabulary adds no
content beyond the restriction rules (ledger s01, s02).

## 5. The campaign ledger

Each row is one campaign; the full line-cited entry, with its assumptions, executed evidence,
retained failures and sources, is in [CAMPAIGN_LEDGER.md](CAMPAIGN_LEDGER.md) and `ledger/`.
Every entry was written from the primary sources by one reader. Entries s01–s09 were then re-read
adversarially by a second reader, and corrections are recorded inside them; for s10–s12 the full
re-reads did not complete, and every claim used in the rows below was instead checked line by line
against the cited sources in the main session.

| # | Campaign | Question it asked | Answer | Classification |
|---|---|---|---|---|
| s01 | V7 / CTP organization | Does CTP doubling, normalization or the response vocabulary derive dynamics, Born probabilities or the asserted predictions? | No. The Feynman–Vernon / Keldysh form is generic ("any local causal open system yields it"); in-in does no work at the gravity rung (Ward identity leaves EH, R², f(R), Lovelock); the declared action "fixes nothing" in its own family; the Born measure is inherited; the differentiators fell 7–47 (decoherence), 22–62 (GW dissipation) and 62.7 (Γ_T) orders below their respective bounds, metrics not reconciled; the bath behind the in-in action was never specified. The V7 file itself is absent from the repository; its content is known only through the handoff. | mixed: recovered physics on supplied choices; rejected mechanisms retained |
| s02 | Influence hierarchy and dependence | Are kernel and noise downstream of subsystem, state and access? Do two-point data fix higher data? Is positivity a new principle? | Yes; no; no. Matched two-point baths split a probe coherence by 0.338 at fourth order and discrimination persists at orders 6 and 8; every positivity constraint reduced to state positivity on the generated algebra (NULL-AS-NEW-PRINCIPLE); the access seed is "a physical boundary condition, not an axiom, not a derived object". | mixed: dominantly recovered physics, bounded by in-class conditional constraints |
| s03 | Level-0 floor and earned theorems | Do the floor theorems constrain a closed whole beyond standard quantum theory? | No. Every theorem is about a supplied classical class (rings, OU chains, declared readouts); the lift to quantum theory is "IRREDUCIBLE/SUPPLIED". The one closed-system item, O-6 (conservative chain with internal bath), was ruled FALSIFIED for a strict arrow from identities, the D-4 recurrence theorem and a non-member chain; the chartered run on the declared chain was never executed. | mixed: recovered mathematics; conditional constraints; three hypotheses falsified |
| s04 | GRUT I–II vacuum response | What was recovered with imports, what was a projector choice, which exclusions are conditional, is the arrow supplied? | Derived count 0 across 26 V1 claims and V2; TT-only projector CHOSEN (Ward identity leaves {P⁽²⁾, P⁽⁰ˢ⁾}, GR's lock −2); the four boundaries hold at stated scopes; arrow direction imported through the KMS gate; Γ_T 62.7 orders below the bound; ISW anchor retired 32σ → 2.0σ. | mixed: recovered physics; supplied choices; conditional exclusions; many retractions retained |
| s05 | GRUT III observer bath | Is P1 a law? Is modular KMS a temperature? What is imported in GR recovery and universal coupling? | P1 is an identification with zero credit, no dynamics, no state; modular KMS is a physical temperature only where the flow is geometric (imports); GR recovery imports the area law, Jacobson and η = 1/4Għ; universal coupling inserts κ; "GRUT III makes no prediction that distinguishes it from GR + QFT in curved spacetime". Slow-roll thermality (D11, D13) is conditional; in FRW two conventional temperatures of different observables exist (a detector temperature in the Wigner convention and a near-horizon Kodama–Hayward one); the secular-infrared survey holds two unreconciled graviton verdicts. | mixed: dominated by recovered physics; five postulates supplied |
| s06 | CP-12 geometry and preparation | Which geometry-selection obstructions are model- and preparation-specific? What does the min-entropy theorem assume? | All obstructions are model- and preparation-specific. The min-entropy theorem needs ‖ρ₀‖_∞ ≤ e^{cN}/D and O(N) excess energy; closure enters only as unitary invariance of the spectrum, together with energy conservation, converting the question into an initial-state obligation of αN log₂N bits, "not a state-selection law". The bound is vacuous at every size its checker prints, through N = 10¹⁸, and the parent relative-entropy theorem it sharpens has no partner audit. Pure coherent preparations are outside it, escaping "without any proof that it assembles geometry"; no theorem excludes every emergent geometry. | mixed: conditional constraints; rejected local-geometry mechanisms |
| s07 | WM / B1 / B2 matter and moving channels | Keep acceptance, recoverability, delivery and retention distinct; is B2 an impossibility proof? | Conditional externally powered preparation (P1); autonomous every-time retention refuted at one permitted point (A1); a conditional delivery window (B1); a fixed-path two-source pass, author-derived with independent review pending, that does not transfer to the moving law; B2's eight-event endpoint failed (1.93%). The closeout is a method and resource decision: no phase-sensitive estimate is available and every certificate is vacuous or out of envelope at the required scale (1411/640 majorant, N ≥ 2²⁸); moving operation is "neither established nor refuted". Closure raised the proof and resource cost; it supplied no principle. | mixed: recovered physics; conditional windows; one refutation; one feasibility stop |
| s08 | Candidate 1 scalar EFT | What is the action, what is supplied, what does the one-loop calculation show? | Einstein–Hilbert + SM + three singlet scalars; everything structural supplied; classically supported charged profile at Q ≈ 39 944 (with an adverse branch and a second higher state); one-loop shift K₁/(gv) = 2.795 (statistical range ±0.11, not a certified enclosure) and a tree-plus-one-loop margin of −3553.71 gv for the soliton against N = 39 944 free carriers, conditional on reduced-saddle stability; the tree-level single-carrier and equal-split gaps −0.199999152115 and −3067.397336565 gv come from the closure-fork packet; no lifetime, no quantum fission, no gravity matching, no novelty (FLS-type soliton, Graham's spectral method). | mixed: conventional EFT; conditional one-loop constraint; supplied couplings |
| s09 | Closure / no-boundary | Does C2 supply a state-selection rule? Is the slow-roll coefficient new? Is curvature sign a discriminator? | No. The restriction map exists for every state and modular KMS is automatic, so observer restriction selects nothing; the Ivo–Li–Maldacena trace-matching conditions are already in the published no-boundary density-matrix prescription. T_eff = H(1−ε)/2π is the first-order term of Obadia's 2008 coefficient (D2 withdrawn). Compactness does not imply k = +1 (flat torus), and closed de Sitter gives Ω_K = −1/(e^{2N} − 1), −4.5·10⁻⁵ at N = 5 and −2.1·10⁻⁹ at N = 10, so the sign is no discriminator. The fork's own report states that closure "alone fixes neither the observer's restriction nor its thermal spectrum" and that the result is "not a theorem that every closed universe is thermal"; the present universe is not slow-roll (ε₀ ≈ 0.47), so the first-order thermal form does not describe it. A sign alone bounds nothing, although a model can be disfavored by its curvature distribution. The B three-prescription control is registered with fixed grids and acceptance rules and unexecuted (its files entered the tree in the M1 commit, so the registration order is self-declared); its target widths are "established physics, not a GRUT prediction". Program decision: option C, preserve the map. | mixed: recovered physics; supplied prescription choices (action, state, cycle, measure, N-distribution); rejected falsifiers; one unexecuted control |
| s10 | M1 memory | What works in the stated finite regime, what was withdrawn, what is not established? | In the stated finite regime the repetition register retains the bit (L = 4: min 0.987510 through sampled t = 10, 0.935065 at t = 40; matter readout 0.778468) and the flying-qubit writer copies a source bit at ≥ 0.999680 with the carrier paying the register's energy change. The matched F = ε control is no worse at every sampled time and size, within 10⁻¹⁰, and better at the tabulated points (L = 4, t = 40: 0.998458), so the facilitation-benefit interpretation is withdrawn. The independent review reads the retention as the standard mechanism, a ferromagnetic gap against a weak transverse field (a T = 0 transverse-field Ising repetition code), with facilitation "a liability" for a prepared record. Not established: a thermal memory phase, cloning of an unknown state, gravity, a chaos exponent, all-times retention. | mixed: recovered physics; rejected mechanism (facilitation benefit); conditional numbers at L = 3–5 on sampled grids |
| s11 | Toric memory / scrambling bridge (separate branch) | What was computed, what was corrected post hoc, what is it not? | An open-system computation: single-face Metropolis dynamics of the 4D toric code at fixed β. The logical half-life grows with L (977 → 1.6·10⁴ → 8.0·10⁴ sweeps at β = 1.0) while gauge-invariant local observables decorrelate in about one sweep; the face spin was reclassified post hoc as gauge-variant; the half-life estimator was changed after an unpreserved first attempt; the 18-qubit OTOC illustration "does not bear on B". It recovers the known thermal stability of the 4D toric code and, as reported, decides the self-correction/no-scrambling bridge claim against the bridge in that classical class, on the logical half-life and damage healing alone if the post-hoc observable correction is rejected; it is not a closed-unitary gravitational scrambling-exponent computation, and its de Sitter and holography statements are uncomputed literature claims. Separate branch, not merged. | mixed: recovered physics; conditional rejected mechanism; conjectures not computed; unexecuted quantum OTOC |
| s12 | Program no-go ledger, postulate map, endogenous access | Which recorded no-gos bear on "closure forces X"? Is access ever derived? | The program's own ledgers answer: a 74-claim register (53 GRUT nodes and 21 vacuum-cluster nodes) that banks no structural-impossibility (FORBIDDEN) no-go, its strongest exports being settled-negative, empirically excluded, suppressed or borrowed; net +16 supplied inputs, an empty derived tier and an empty prediction ledger; every kernel alternative "excluded by an input, never by a principle"; the arrow's direction, the Born measure, the system/bath split, ħ, c, G, the action, the projector, the generator, the lift and the access seed all recorded as supplied. EA-0, the gate that attacked endogenous access directly, terminated UNFORMULABLE at the earned scope and the bridge gate NONUNIQUE-LIFT: "nothing selects the access seed; it is currently irreducible". Pure global state, compact topology and finite four-volume are not addressed anywhere in these sources. | mixed: dominantly supplied choices; rejected mechanisms; unexecuted tasks; conditional in-class constraints |


## 6. Proposed additions that did not survive, and what remains open

Grouped by the kind of addition each would have been (FORMAL_CORE §7). Every item is a
recorded withdrawal, refutation or non-implication in the repository; the ledger entries give
the lines.

**State selection.**
- C1 defined as the Feynman–Lehners–Turok Dirichlet prescription terminates on an
  inverse-Gaussian coefficient; that is a conditional termination of that definition, not of
  every meaning of closure (s09).
- C2 identifies the observer-state object and adds no selection condition beyond restriction,
  trace matching, positivity, normalization and modular KMS; the no-boundary density-matrix
  construction already contains the trace-matching conditions (s09).
- D2 "slow-roll thermality" novelty withdrawn (Obadia 2008 is the prior source for the
  first-order coefficient); the sign-only curvature falsifier and "all sequestering
  formulations force positive curvature" withdrawn (s09).
- The TT-only projector of GRUT I–II is a choice, not a consequence (s04); the sharp
  Sorkin–Johnston state is rejected and "P1 selects the horizon state" remains a conjecture (s05).
- The CP-12 min-entropy theorem is an initial-state obligation, not a selection law (s06).
- The HJW route from P0 + P1 to the outcome rule is withdrawn (s05).

**Dynamical constraint.**
- The Schwinger–Keldysh formalism does not select the Einstein–Hilbert action; the in-in
  machinery does no work at the gravity rung, which reverts to Jacobson with G inserted
  (s01, s04, s05).
- One-loop induced gravity with Standard-Model content does not supply G (s05).
- a₀ via Jacobson, a memory time 1/(4H), a bulk Planckian fluid, c_eff at LIGO frequencies, and
  w ≡ −1 as a GRUT discriminator are rejected (s05).
- The 689 Hz falsifier is withdrawn, the ISW anchor retired from 32σ to 2.0σ, Γ_T is 62.7 orders
  below its bound, GW dephasing 22–62 orders below detectability (s01, s04).
- The founding ζ-hypothesis, the conformalon, the information principle, the single-particle
  interference memory mechanism and the symmetry route to forced TT are rejected; single-pole
  finite memory stays a conjecture (s01, s04).
- Candidate 1 is a quartic extension of Friedberg–Lee–Sirlin with supplied couplings; its
  "coefficients vanish at every scale" rule and a pointwise Coleman–Weinberg shortcut are
  rejected (s08).
- The classical relaxational carrier is eliminated inside the supplied quantum class (CA-1);
  H-HERM-2, strict H-ORD and the two-property reduction of the floor are falsified (s02, s03).
- The bath behind the in-in action was never specified: "what bath Hilbert space was integrated
  out to make the influence functional" is recorded as the deepest open item, the noise-origin
  campaign closed with NO-PHYSICAL-BATH-CANDIDATE, and GRUT III retired the finite-memory
  ontology ("memory is a property of (observable, state)") (s01, s02).

**Access structure.**
- The access seed is load-bearing under boundary change and nothing in the record fixes it:
  "a physical boundary condition, not an axiom, not a derived object" (s02); no earned Level-0
  structure selects the physical lift, and observability at the earned readout is an identity
  (s03). 
- EA-0 terminated UNFORMULABLE at the earned Level-0 scope: every seedless access candidate is
  global, generically trivial, a correlation readout, or generator-carried by identity; the
  successor bridge gate terminated NONUNIQUE-LIFT; "nothing selects the access seed; it is
  currently irreducible" (s12).

**Measurement rule.**
- Decoherence does not yield outcomes or Born weights in the tested branch-preserving unitary
  class (EXPERIMENT_P_FRONTIER closed at that class); P5 remains a postulate; a positivity
  "principle" reduces to state positivity (NULL-AS-NEW-PRINCIPLE) (s01, s02, s05).

**Arrow and clock.**
- The arrow's direction is imported with the passive, uncorrelated past-boundary state; the
  exactly solvable model has Γ(t) = Γ(−t) and full recurrence; closure plus partition plus
  product start does not force a strict arrow (band-edge memory, anti-arrow at switch-on;
  recurrence obstruction D-4) (s01, s03, s04).
- The temporal generator is supplied (S5-0 GENERATOR-IRREDUCIBLE); P1 assigns no temperature in
  FRW with ε ≠ 0; modular time equals physical time only in the imported geometric cases (s02, s05).

**Memory and structure mechanisms.**
- M1's facilitation benefit is withdrawn: the matched uniform-floor control retains better at
  every sampled time for L = 3, 4, 5; a single cell fails the 0.9 standard at L = 4, 5; bare
  matter-to-memory writing fails; unconstrained flips at 0.1 and 0.3 violate retention; exact
  saturation is a barrier to formation as well as a memory (s10).
- Autonomous every-time retention is refuted at the permitted point (A1); B2's eight-event
  endpoint failed; the fixed-path pass does not transfer to the moving law (s07).
- The toric report compares thermal Metropolis trajectories and decoded memory with a post-hoc
  gauge-observable correction; it is not a closed-unitary gravitational scrambling-exponent
  computation. Its rejection of the bridge claim is as reported on the unmerged branch, not
  re-run here, independent review outstanding (s11).

**Still open (registered, not executed).** The B three-prescription control (FLT, whole-line
Hartle–Hawking, Di Tucci–Lehners Robin; Yamada as a distinct prescription) (s09); Candidate 1's
three-charge one-loop family, lifetime and nonlinear gravity (s08); the O-6 declared-chain run
(s03); the rung-3 pole-versus-cut ask never sent (s01); independent review of CP-04 and CP-05
(s05); M1's error-correction cascade against the matched baseline (s10); the toric branch's
merge decision (s11); the partner audit of the CP-12 parent relative-entropy theorem (s06);
independent review of the B2-RA fixed-path pass and of B1's formation/use composition (s07); the
two unreconciled graviton verdicts of the GRUT III secular-infrared survey (s05). The universality probe deferred, the criterion-unification fork never run, EA-1 closed and S-5
held (s12); the quantum OTOC and Lyapunov question in the 4D code (s11). One external signal is
recorded against a GRUT I–II shape exclusion: the DESI DR2 + CMB preference for a w(z) that
crosses the phantom divide at z ≈ 0.35–0.5 (3.1σ), "pending refutation if the signal
consolidates" (s04); the exclusion itself is conditional on the supplied projector and arrow.


## 7. Physical scope, predictions and open claims

**Has any additional law or state-selection rule survived?** No. Every kind of addition in
FORMAL_CORE §7 except a modified global law was tested by at least one campaign; the ledger
reports the result of each. The surviving
GRUT objects are: the question; the formal core, which is standard; the compatibility and
underdetermination computations above, which are standard; and the conditional model
instances, each with its own Θ. They do not combine into one law. Candidate 1 (a
conventional Einstein plus Standard Model plus scalar effective field theory), M1 (a finite
Ising repetition register with matter feedback and a programmed writer), and GRUT III (five
postulates with imports) live in different state spaces, and no common action, state space,
constraint set or calibration has been established; the handoff's no-cross-import rule stands.

**Predictions, conditional on inputs.** Within any supplied Θ the framework predicts Born
probabilities and multi-time statistics; calculation 1 is such a prediction set and
calculation 2 is a set of compatibility constraints, and COND 2.10 bounds time-averaged
local deviations once d_eff is known. Candidate 1 predicts, for its supplied couplings, a
one-loop energy shift K₁/(gv) = 2.795 (statistical range ±0.11, not a certified enclosure) and a
tree-plus-one-loop margin of −3553.71 gv for its soliton against free carriers in the benchmark
channel, conditional on reduced-saddle stability (ledger s08). M1 predicts, for its supplied
sizes and fields, retention probabilities that the matched uniform-floor control exceeds at
every sampled time (ledger s10; re-run in this packet, line for line). GRUT III predicts
geometric-horizon temperatures where the modular flow is geometric, all imported (ledger s05).
None of these is a prediction of the closed-whole framework itself; each is a prediction of a
model instance.

**What extra condition would make a distinct prediction.** A rule that fixes a slot of Θ
beyond consistency. For the state: in the finite class a rule built from H alone and
equivariant under unitaries is stationary and freezes every subsystem (calc 1 S9), so a state
rule with content must use the access structure, a clock or some non-equivariant datum; in
locally covariant QFT an embedding-covariant rule is excluded under the hypotheses of
Fewster–Verch Theorem 6.13 (FORMAL_CORE COND 5.3); and a bottom-up family of local states must
pass the marginal problem. For the dynamics: a constraint on H with recovered known physics and
its own constants. For access: an endogenous structure that does not presuppose the
factorization it derives, or a stated relation such as locality of H, which generically picks
at most one factorization. For the clock: an identification of modular time with physical time
beyond the geometric cases. Or a measure over Θ, whose consequences COND 2.11 illustrates. Each
must be compared to this packet's baseline, in which the prediction is conditional on a fully
supplied Θ. The program's record of attempts is the ledger; none has produced a condition that
is not already restriction, trace matching, positivity, normalization or modular KMS.

**Claim boundaries.** The finite statements are the regulated model's; the algebraic
statements are identifications; the gravitational and cosmological statements are scope
limits, not results. Historical bounds and ledger counts are not current predictions. The
toric memory report on its separate branch is not merged and is not a gravitational
scrambling computation (ledger s11). The B three-prescription control remains registered and
unexecuted (ledger s09). The facilitation-benefit claim of M1 is withdrawn (ledger s10).

## 8. Verification record

- `python3 code/verify.py` checks the hash of every file listed in MANIFEST.json, 28 headline results
  and 2 provenance hashes; `--rerun`
  re-executes both calculations and compares every numeric leaf with the committed results.
- A separate agent re-implemented the calculation-1 quantities from the specification alone and
  matched 23 comparisons; it was re-compared against the committed results after the final
  edits. A separate re-execution verified the pre-review version. Two adversarial code reviews
  found nine and seven items, all applied. Records: `results/verification/`, summarized in
  [CALCULATIONS.md](CALCULATIONS.md) §Verification.
- The M1 matched uniform-field control script was re-run unchanged and reproduced the recorded
  table line for line (`results/m1_matched_control_rerun.txt`).
- Ledger: s01–s09 adversarially re-read by a second reader (8 to 15 corrections each, all
  folded into the entries); s10–s12 single-reader entries whose report-row claims were checked
  line by line in the main session. The completeness critic over the ledger did not run.
- Attack on the packet: [[ATTACK_STATUS]]
- Not done: no old solver reopened; no cosmology contour executed; no toric code run; no
  determinant recomputed; no field added to M1; no new branch; no external contact.
