# GRUT handoff: write the working theory around the closed-whole question

19 October 2026. Owner: D. Ryan Grover. Status: **handoff before the new
working-theory synthesis is written**. This packet adds no physical result.

## The owner's current direction

> Our old V7's parent object was S_CTP as a "genotype," which is standard
> open-system physics wearing a new name. The honest parent object is the
> question the whole program has been asking since v1: what does a closed
> system whose every subsystem's bath is the rest of the universe force,
> that known physics does not already contain? Every sector becomes a
> campaign that asked that question in one place and reports the answer.
> Most answers are "nothing beyond known physics, and here is the
> computation that shows it." That is a result, not a failure. Work out
> our working theory.

The owner then requested a fresh chat because the conversation was too long.
The next chat should **produce the working theory**, with equations,
calculable examples and a sourced campaign ledger. Do not consume its turn
writing another handoff or merely proposing to continue. A question can
organize the research; each physical model still needs actual dynamics,
state/preparation and measurement rules.

## Repository and checkpoint

- Repository: `ryangrvr/TestingGRUT`.
- Working branch: `codex/grut-closure-fork-20261019`.
- Last scientific amendment before this handoff:
  `4e9a161a041ed23ca04052cbbe83a5c878e8b447` (M1 matched-control correction).
- Its parent is the independent M1 review
  `cc249673086b984f090d170f394d4d635e9c451d`; original M1 packet is `1997d20`.
- The handoff commit will be a descendant of `4e9a161`; its exact hash is
  supplied in the completion message. Read the branch before assuming a
  newer review or campaign has not appeared.
- Session checkout: `/workspace/scratch/9a57d63bee22/TestingGRUT-review`.
  This transient path is not required in a new environment.
- Read `AGENT_COORDINATION.md` before editing, claim the intended area,
  preserve other workers' claims, and log verification actually performed.
- Research, bounded local calculations, additive reports and posting to
  this repository are authorized. Do not ask for routine confirmation.
  No spending, account access, external messaging or indefinite jobs.
- Stay on the existing branch when practical. Earlier B-control work used
  a one-branch-per-step / at-most-two-branches-per-day budget. This handoff
  does not authorize a proliferation of new branches or architectures.

The checkout contains unrelated untracked temporary artifacts outside this
packet. Leave them untouched and stage only the intended files. No physical
calculation from this session remains running. No subagents were spawned
for the matched-control amendment or this handoff.

## What has already been established about the framing

The proposed common core is a **closed global description with mutually
compatible subsystem descriptions**. "Bath = rest" means an internal
complement relative to a specified subsystem and access structure. It does
not mean that the complement is automatically thermal, Markovian, infinite,
weakly coupled, uncorrelated, or a resettable reservoir.

For a finite regulated model with a supplied factorization, the starting
equations can be written explicitly:

\[
\rho_\Omega(t)=U(t)\rho_0U(t)^\dagger,\qquad
\rho_A(t)=\operatorname{Tr}_{\bar A}\rho_\Omega(t),
\qquad H=H_A+H_{\bar A}+H_{A\bar A}.
\]

For time-independent H, the exact reduced equation is

\[
\dot\rho_A=-\frac{i}{\hbar}[H_A,\rho_A]
-\frac{i}{\hbar}\operatorname{Tr}_{\bar A}
[H_{A\bar A},\rho_\Omega].
\]

This exposes the physical question: what fixes the global law, state,
interaction and access structure? Renaming the last term "response" does
not supply their values. Eliminating the complement can yield memory and
noise; its CTP influence functional is a derived calculational object.
Keep all required assumptions for Gaussian, stationary, KMS, Markov or
overdamped approximations visible.

For nested regions, restrictions must agree. For overlapping regions,
agreement on the overlap is necessary but need not be sufficient for a
common positive global extension. Every sector's effective environment
must be compatible with one global model and preparation. These are
standard quantum consistency requirements, not a newly discovered law.

In continuum/gauge/gravitational settings do not assume this tensor
factorization. Use a declared physical algebra and state restriction.
GRUT III's P1 identifies the bath with a commutant, equating it to the
causal-complement algebra only where Haag duality holds. Local type-III
algebras do not have the elementary reduced-density-matrix construction.
The finite Hamiltonian's external time is an input; this is not a supplied
solution to quantum gravity or the problem of time.

**Not forced by closure alone:** a pure global state, compact spatial
topology, finite four-volume, positive curvature, a preferred initial
state, a thermal bath, universal relaxation, an exponential memory kernel,
a chosen arrow, a small gravitational Lyapunov exponent, observed constants,
or a specified matter/gravity action. A closed system can have any of
several dynamics and preparations. Global energy conservation in a fixed,
time-independent Hamiltonian model must not be confused with a universally
defined total energy in an expanding GR spacetime.

These paragraphs are a construction brief, not a newly executed general
theorem. The next chat should derive the exact statements and test proposed
extra implications with explicit countermodels.

## Campaign evidence: what survives and what must not be promoted

| Campaign | Supported answer and boundary | Source to read |
|---|---|---|
| Original V7 / CTP organization | V7 supplies an action-based organization and openly imports the SM. CTP doubling, normalization and a response vocabulary do not independently derive dynamics, Born probabilities or all its asserted predictions. Historical claims require their later grades. | Original `GRUT_V7.md`; details below |
| Influence hierarchy / dependence campaigns | Kernel and noise are downstream of subsystem, state and access choices. Matching two-point data need not match higher influence data. Positivity constraints reducing to quantum state positivity are not a new principle. | [September synthesis](../../../GRUT_WORKING_THEORY_SYNTHESIS_01.md), [record index](../../../GRUT_RECORD_RECONCILIATION_INDEX_01.md) |
| GRUT I–II vacuum response | Preserve the recorded recovery-with-imports, projector-choice results, conditional exclusions and supplied arrow. Do not turn old bounds or historical ledger counts into current universal predictions. | [I](../../../GRUT_I_What_Survived.md), [II](../../../GRUT_II_What_Survived.md) |
| GRUT III / observer bath | P1 is an identification, not a new state-selection law. Modular KMS is not automatically a temperature for a physical clock. GR recovery and universal coupling import physical inputs. | [principles](../../../GRUT_III/01_PRINCIPLES.md), [synthesis](../../../GRUT_III/13_SYNTHESIS.md), [C2 response](../grut_closure_fork_20261019/CLOSURE_SELECTION_RESPONSE.md) |
| CP-12 geometry / preparation | There are model- and preparation-specific geometry-selection obstructions. The min-entropy extension concerns a particular incoherent high-entropy class and a loop-rich graph projector; pure coherent preparations are outside that assumption. This is not a theorem excluding every emergent geometry. | [min-entropy report](../cp12_preparation_minentropy_audit/REPORT.md), [gravity adjudication](../../../GRUT_III/14_CP12_GRAVITY_ADJUDICATION.md) |
| WM/B1/B2 matter and moving channels | Keep preparation acceptance, recoverability, delivery and retention distinct. The B2 closeout is a feasibility stop within a specified computational envelope, not a physical impossibility proof. Results from incompatible state spaces are not composable. | [closeout disclosure](../grut_b2_feasibility_closeout_20261018/REVIEW_STATUS.md), original WM/B reports via the continuation/index documents |
| Candidate 1 scalar EFT | A conventional Einstein/SM-plus-scalar candidate has specified action, classically supported charged objects and a model-specific one-loop calculation. The result is conditional on supplied couplings; no novelty or validated universe follows. | [candidate](../grut_candidate_1_20261018/REPORT.md), [energy](../grut_candidate_1_energy_20261019/REPORT.md), [closure fork](../grut_closure_fork_20261019/REPORT.md) |
| Closure / no-boundary | C2 as stated adds observer interpretation, not the missing state-selection condition. The slow-roll coefficient is not a new general thermality result; curvature sign without a magnitude is not a discriminator. The registered B three-prescription control remains unexecuted. | [C2 response](../grut_closure_fork_20261019/CLOSURE_SELECTION_RESPONSE.md), [B status](../grut_closure_control_20261019/STATUS.md) |
| M1 memory | A finite Ising repetition register with matter feedback and a programmed writer works in the stated finite regime. The matched F=epsilon control retains better; the facilitation-benefit interpretation is withdrawn. | [corrected report](../grut_memory_candidate_20261019/REPORT.md), [open claims](../grut_memory_candidate_20261019/OPEN_CLAIMS.md) |
| Toric memory / scrambling bridge | Separate report on another branch; read before importing. It compares thermal Metropolis trajectories and decoded memory, with a post-hoc gauge-observable correction and a separate quantum illustration. It is not a closed-unitary gravitational MSS-exponent computation. | Branch and commit below |

Negative results are part of GRUT's output: specify the hypothesis, model
class, calculation, counterexample or obstruction, and what remains free.
"Nothing additional was forced in this campaign" is a result. It is not
proof that no possible future extra postulate could have content.

Do not import two successful modules into one theory merely because both
have GRUT in their names. If the working theory is a shared framework with
different model instances, say so; if claiming one combined physical law,
actually establish a common action/state space, constraints and calibration.

## Recent numerical facts that a new chat must retain

M1's decisive matched-control comparison (same 0.05 floor, all other inputs
fixed) was executed twice by different Hamiltonian constructions:

| L | M1 min P through t=10 | Uniform-floor min P | M1 P at t=40 | Uniform-floor P at t=40 |
|---|---:|---:|---:|---:|
| 3 | 0.994999 | 0.998738 | 0.985006 | 0.994692 |
| 4 | 0.987510 | 0.999055 | 0.935065 | 0.998458 |
| 5 | 0.989132 | 0.998479 | 0.940399 | 0.997517 |

Uniform-floor retention is no worse at all 201 sampled times at each size,
within 1e-10 tolerance. This is not an all-continuous-time theorem.
At L4 the best matter readout is .778483402 versus .778467821 at t=.58;
it is nearly, not exactly, identical. The old F=1 control raised the field
twentyfold and failed to isolate facilitation. Six original stronger
scientific attack conditions remain false. The ten-run replay reproduces
all seven original conditional traces exactly. The corrected finalizer
checks 31 files plus its manifest; all recorded integrity checks pass.

The original M1-W writer's final per-cell recovery is at least .999680.
Its subsequent storage value at least .991697 through sampled t10 uses
the unconditional spin state with residual carrier coupling neglected.
The total-energy drift is .000712816 J in its model units, nonmonotonic
with cutoff. No all-times scattering certificate, unknown-state cloning,
thermal memory phase, gravity or chaos exponent follows. The impulsive
writer is a designed controlled flip, not evidence for facilitation.

Candidate 1's approximately +2.80 gv one-loop energy shift versus a
classical free-carrier binding margin about -3554 gv supports the specified
channel in that benchmark, not lifetime or all quantum fragmentation.
The cheap tree-level channel calculations give single-carrier gap
-.199999152115 gv and equal-split gap -3067.397336565 gv at Q=39944.
Keep tree and loop orders separate; do not recompute expensive determinants
merely to repeat that weak-coupling conclusion.

## V7 and other sources actually read in this interrupted synthesis

The original `GRUT_V7.md` was resolved in the owner's files and read in
full, lines 1–332, 18,380 bytes, headed D. Ryan Grover, April 2026. It is
not present at a repository path in this checkout. The parent-object
"genotype" wording is in section 2; the SM-import admission is section 15;
sections 34–43 contain the old claims and limits. No PDF variant was audited.
The file itself ends with 191+ tests; do not substitute a later PDF's
250+ count or treat either as physical validation. Use the Library skill
to resolve/read the exact file if needed. Do not publish private file IDs.

Read current M1 documents and raw correction, GRUT III principles,
selected synthesis sections, closure/C2 reports, CP-12 min-entropy report
and B2 execution disclosure. Some earlier multi-document reads were output
truncated: this is not a complete re-audit of the entire repository.
Locate CP-04/CP-05 primary reports through the GRUT III source references
before attaching precise termination labels or a fixed count of five.

Toric branch fetched, not merged:
`claude/toric4d-memory-vs-scrambling-20261019`, commit
`8e0f89f5d675f13fe8aa62e5ee17196b406dca12`.
The first 105 lines of its REPORT were read, including its method, estimator
change and gauge-observable correction. No toric code or trajectories were
rerun here. Preserve that distinction from the reviewer's broader claims.

Primary literature located for the next derivation: Feynman–Vernon,
Ann. Phys. 24 (1963) 118–173, DOI 10.1016/0003-4916(63)90068-X;
Pollock et al., arXiv:1801.09811 (operational Markov condition);
Fewster–Verch, arXiv:1106.4785 (natural-state limitations under stated
axioms). Only search/abstract metadata was checked for this synthesis;
inspect the relevant primary text before applying its theorems. No fresh
global closure theorem or finite-model demonstration was executed here.

## Next-chat deliverable and bounded execution order

Write a readable `WORKING_THEORY.md` under a separately named closed-whole
packet, with detailed evidence behind it. The core deliverables are:

1. **Parent question and formal core.** State what closure means, the exact
   global-to-local restriction/evolution rules, what is supplied, how records
   and interventions are represented, and the requirement of a common global
   realization. Give an explicit theory datum. Keep CTP downstream.
2. **Calculable demonstration.** Use a small specified closed Hamiltonian
   to carry preparation through evolution to measurements and energy/correlation
   accounting. Exhibit why closure does not force monotone relaxation or a
   unique local law. An initial-correlation example can show why a reduced
   state alone does not specify an arbitrary preparation's future.
3. **Compatibility check.** Derive nested/overlap consistency and its limits.
   A useful test is two proposed Bell-pair marginals sharing one qubit:
   identical one-qubit overlaps do not guarantee a common global state.
   This is a suggested standard-physics example, not a new computation already run.
4. **Campaign ledger.** For each sector give its question, assumptions,
   executed evidence, answer, remaining freedom and source. Classify
   recovered known physics, conditional constraint, supplied model choice,
   rejected mechanism, conjecture or unexecuted task. Retain every failure.
5. **Physical scope and open claims.** State whether any additional law or
   state-selection rule has survived. If none has, call the result a working
   closed-system framework with conditional model instances and tested
   boundaries, not a new empirically validated universal dynamics. Explain
   exactly what extra condition would be needed to make a distinct prediction.

Use three passes: formal specification, a small discriminating calculation,
then synthesis and attack. Do not reopen every old solver, add fields to
rescue M1, resume the cosmology dispute automatically, or inflate the task
into unlimited architectures. A fully specified new postulate is allowed
only with explicit costs and a separate comparison to the known-physics
baseline; it is not required to manufacture novelty.

The prior two-calculation idea and Bell-compatibility example are pending.
The original theory-writing plan was interrupted by this handoff request;
there is no unpublished completed working theory to retrieve.

## Completion standard

The next result should answer: **what does GRUT now say, what calculates,
what does closure actually constrain, which proposed additions did not
survive, and what remains a supplied choice?** It can honestly be a theory
framework of known physics plus computed boundaries. The owner explicitly
values that result. It must be more than another reconciliation document,
and it must not announce that gravity or all reality has been derived.
