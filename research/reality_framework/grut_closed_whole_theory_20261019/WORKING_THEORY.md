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

**What closure constrains.** Exactly the theorems of standard quantum theory: unitary
invariants, the exact reduced equation dρ_A/dt = −i[H_A, ρ_A] − i Tr_Ā[H_{AĀ}, ρ_Ω], the
Nakajima–Zwanzig structure with its initial-correlation term, complete positivity of reduced
maps for product preparations, nested compatibility, necessary overlap compatibility, the
entropy inequalities, quasi-recurrence in finite dimension. Closure constrains **nothing
local**: not relaxation, not Markovianity, not a kernel, not a unique local law, not
thermality of the complement, not an arrow, not record persistence, not a preferred state, not
the constants, not the action. Each negative is exhibited by a legal closed model in
calculation 1 (§3 below), and the compatibility requirement has a quantitative edge in
calculation 2 (§4).

**Which proposed additions did not survive.** [[LEDGER_ADDITIONS]]

**What remains a supplied choice.** The six entries of Θ, in every sector: the field content
and action (GRUT III P0 and data; Candidate 1's couplings), the preparation (GRUT I–II
projector choice; CP-12 preparation classes; the no-boundary state that C1/C2 did not select),
the access structure (EA0; the register and writer in M1), the measurement rule (P5), the
arrow (P4), and the time parameter. [[LEDGER_SUPPLIED]]

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
lists exhaustively what such a postulate would have to constrain and what it would cost.

## 2. The formal core in equations

Finite regulated case, ħ = 1, time independent H, supplied factorization ℋ_Ω = ℋ_A ⊗ ℋ_Ā and
H = H_A + H_Ā + H_{AĀ}:

    ρ_Ω(t) = e^{−iHt} ρ_Ω(0) e^{iHt},                       (global law; unitary)
    ρ_A(t) = Tr_Ā ρ_Ω(t),                                     (restriction; exists for every state)
    dρ_A/dt = −i[H_A, ρ_A] − i Tr_Ā [H_{AĀ}, ρ_Ω(t)].          (exact; the last term is global)

With 𝒫ρ = (Tr_Ā ρ) ⊗ ρ_ref and ℒ = −i[H, ·] the same equation reads

    𝒫ρ̇(t) = 𝒫ℒ𝒫ρ(t) + ∫₀ᵗ 𝒦(t−s) 𝒫ρ(s) ds + 𝒫ℒ e^{𝒬ℒt} 𝒬ρ(0),     𝒦(τ) = 𝒫ℒ e^{𝒬ℒτ} 𝒬ℒ𝒫,

so the memory kernel is a function of all of H, of ρ_ref and of the factorization, and the
inhomogeneous term vanishes only for product preparations. For a product preparation the map
ρ_A(0) ↦ ρ_A(t) is linear and completely positive but in general not CP-divisible; for a
correlated preparation there is in general no such map at all. The complete multi-time object
is the process tensor, a multilinear completely positive map from interventions to outcomes,
itself determined by Θ. Records are registers R in the access structure written by unitaries;
an instrument on A and a unitary record into an idle, traced-out register give the same later
statistics (FORMAL_CORE THM 2.8).

The Feynman–Vernon influence functional F[x, x′] = Tr_Ā[U_Ā[x] ρ_Ā U_Ā[x′]†] is the
path-integral form of the same reduced dynamics. It is Gaussian, and fixed by a dissipation
kernel and a noise kernel built from a spectral density J(ω), only when the complement is
linear, linearly coupled and Gaussian-prepared. The closed-time-path generating functional is
that object in field-theory notation. A quadratic S_CTP therefore presupposes J(ω), the
temperature, linearity and the preparation; it is where GRUT V7's parent object actually
lives, and it earns nothing by itself (ledger s01).

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
| Entropy theorems | pure global state: S(ρ_A) = S(ρ_Ā) to 3·10⁻¹⁶; mixed trajectory: mutual information ≥ 0 and Araki–Lieb slack 0.21–0.24 away from saturation after t ≈ 3, no violation |

Every line is standard quantum mechanics. The point of computing them in one supplied model is
that the parent question's negative answers become exhibits rather than assertions.

## 4. Compatibility and underdetermination

**Compatibility** (CALCULATIONS §2, `results/marginal_compatibility.json`). Two Bell pairs
Φ⁺ on AB and on BC agree exactly on the shared qubit (ρ_B = 1/2 from both sides) and have no
common three-qubit state: the semidefinite program is infeasible, and analytically a pure
ρ_AB forces ρ_ABC = ρ_AB ⊗ ρ_C. The obstruction is quantitative: the largest achievable
F_AB + F_BC is 3/2 (it would be 2 if compatible); Werner states W(p) on AB and BC are jointly
realizable only for p ≤ 2/3, and on all three pairs of a triangle only for p ≤ 1/3, the
separability threshold. Classically, two marginals sharing one variable always extend, and
an odd cycle of perfect anticorrelations does not; the quantum failure on a chain is the
specifically quantum constraint (monogamy). Consequence for GRUT: an "effective environment"
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

[[LEDGER_TABLE]]

## 6. Proposed additions that did not survive, and what remains open

[[LEDGER_FAILURES]]

## 7. Physical scope, predictions and open claims

**Has any additional law or state-selection rule survived?** No. Each row of FORMAL_CORE §7
was tested by at least one campaign; the ledger reports the result of each. The surviving
GRUT objects are: the question; the formal core, which is standard; the compatibility and
underdetermination computations above, which are standard; and the conditional model
instances, each with its own Θ. They do not combine into one law. Candidate 1 (a
conventional Einstein plus Standard Model plus scalar effective field theory), M1 (a finite
Ising repetition register with matter feedback and a programmed writer), and GRUT III (five
postulates with imports) live in different state spaces, and no common action, state space,
constraint set or calibration has been established; the handoff's no-cross-import rule stands.

**Predictions, conditional on inputs.** Within any supplied Θ the framework predicts Born
probabilities and multi-time statistics; calculation 1 is such a prediction set and
calculation 2 is a set of compatibility constraints. Candidate 1 predicts, for its supplied
couplings, a one-loop energy shift of about +2.80 gv against a classical free-carrier binding
margin of about −3554 gv in its benchmark channel (ledger s08). M1 predicts, for its supplied
sizes and fields, retention probabilities that the matched uniform-floor control exceeds at
every sampled time (ledger s10; re-run in this packet, line for line). GRUT III predicts
geometric-horizon temperatures where the modular flow is geometric, all imported (ledger s05).
None of these is a prediction of the closed-whole framework itself; each is a prediction of a
model instance.

**What extra condition would make a distinct prediction.** A rule that fixes an element of Θ
beyond consistency: a positive, normalizable state functional ρ_Ω(0) = F(H, access) that
survives the natural-state obstruction (FORMAL_CORE COND 5.3) and the marginal problem; or a
constraint on H with recovered known physics and its own constants; or an endogenous access
structure that does not presuppose the factorization it derives; or an identification of
modular time with physical time beyond the geometric cases. Each must be compared to this
packet's baseline, in which the prediction is conditional on a fully supplied Θ. The program's
record of attempts at each is the ledger; none has produced a condition that is not already
restriction, trace matching, positivity, normalization or modular KMS.

**Claim boundaries.** The finite statements are the regulated model's; the algebraic
statements are identifications; the gravitational and cosmological statements are scope
limits, not results. Historical bounds and ledger counts are not current predictions. The
toric memory report on its separate branch is not merged and is not a gravitational
scrambling computation (ledger s11). The B three-prescription control remains registered and
unexecuted (ledger s09). The facilitation-benefit claim of M1 is withdrawn (ledger s10).

## 8. Verification record

- `python3 code/verify.py` checks every packet hash and 22 headline numbers; `--rerun` re-executes
  both calculations and compares every numeric leaf (identical on two independent re-runs).
- An independent re-execution and an independent re-implementation from the specification
  alone were performed by separate agents and are summarized in [OPEN_CLAIMS.md](OPEN_CLAIMS.md)
  §Verification, with any discrepancies and their resolution.
- The M1 matched uniform-field control script was re-run unchanged and reproduced the recorded
  table line for line (`results/m1_matched_control_rerun.txt`).
- Every ledger entry was written by one reader and re-read adversarially by a second against the
  same sources; corrections are recorded inside each entry.
- Not done: no old solver reopened; no cosmology contour executed; no toric code run; no
  determinant recomputed; no field added to M1; no new branch; no external contact.
