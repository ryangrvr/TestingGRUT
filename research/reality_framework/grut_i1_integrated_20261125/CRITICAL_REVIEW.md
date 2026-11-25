# GRUT-I1 critical review

A separately tasked internal AI reviewer read the frozen protocol, the parent
D7 recovery proof, the D3–D6 access assumptions and the primary comparison
literature. The reviewer independently implemented density evolution and
formed the Bell scores. This is a separate implementation inside Work.
External RAI and hardware replication remain open. [External handoff](RAI_HANDOFF.md)

## A useful task value with a declared information boundary

The value difference is monotone for nested admissible policy sets. Blind
policies must remain implementable when the record is exposed. The same
task, initial data and resource caps apply in both cases. If the intervention
changes the Hamiltonian, preparation, time cost or correlations with the bath,
its measured difference is a configured mechanism contrast. Operational
evaluation of that contrast and the optimal information value answers two
separate questions.

The difference lies in \([0,1-V_{\mathrm{blind}}]\) for a bounded unit score.
Adding records in succession gives a telescoping identity. Complementary
records can be individually worthless, as with two fair shares whose XOR is
the target. Duplicates have zero incremental value once either copy is
available. A single-task value cannot determine what an observer knows
generally. Simulating the exposed record from blind-accessible data within
the same budget makes the gain zero. This links the definition to quantum
decision theory and statistical sufficiency. Characterizing a whole
comparison order can require a larger family of tests with side systems.
[Buscemi, quantum sufficiency](https://arxiv.org/abs/1004.3794) ·
[Chiribella, D'Ariano and Perinotti, quantum networks](https://arxiv.org/abs/0904.4483) ·
[Takagi and Regula, operational resource comparisons](https://arxiv.org/abs/1901.08127)

| Earlier case | Reviewer scope check |
|---|---|
| D3 | The fault register holding the past is excluded. The surviving active states without M have identical historical marginals immediately after erasure. If the writer vacuum weight is \(q\), a free M-only discriminator gets \(1-q/2\). Coherent access to M and the writer raises it to \((1+\sqrt{1-q^2})/2\). These discriminators are separate from the measured fixed-H joint repair and retention. |
| D4 | The collision of distinct histories with identical triples and reservoirs caps any restoration across that class. A blind score on a different single-fault ensemble needs its fault-site prior and accessible bath data. With a fair historical bit, no ancillary evidence and S-only access, that guess optimum is \(\max(p,1-p)\), where \(p\) is the probability of corrupting S. |
| D5 | The seed h remains available in the head with G. A controller allowed to use both can regenerate the target with M withheld. The original route's response to scrambling or disabling M diagnoses its own writer. The direct-seed null prevents assigning its full intervention contrast to unrestricted information value. |
| D6 | The one-query 0.75 weak-information score and 0.50 blind score apply at the declared response cut. The direct source and forecast witness are excluded. A constant witness itself permits a perfect blind forecast. The complete clock before the cut may carry different information. |
| D7 | The trace-norm optimum is exact at the isolated recovery interface. Blind optimization gives 0.853553390593 shallow and 0.50 deep. The accessible programme changes those into 0.895284707521 and 0.75 for the weak assay. A source bypass or access to the reference changes the interface. |

## Physical channel and thermodynamics

The source sign must persist through the field steps until the specified
local flips. Averaging an independent sign at each step instead destroys the
assayed correlation and defines a different model. Memory contrast relative
to the initial sign can remain different from contrast relative to the
current sign. The phase in the coherence block carries the entire source
history. Acting on B and M with independent fresh cells respects the stated
access pattern. The isolated R is used exclusively to score. The backward
effect of a later trace-preserving operation on any earlier probability is
zero by trace preservation.

The common hot and cold contacts leave a population map independent of
source and memory. Write \(u=P(1\mid0)\) and \(v=P(0\mid1)\);
then \(0\le u+v\le1\). In the Bell test, the summed correct-population
contribution before the factor \(1/4\) is
\(u+v+(1-u-v)\sum_i\operatorname{Tr}(V_{ii}\bar\rho V_{ii}^\dagger)\)
for a Stinespring recovery. Each stay block is a contraction. The upper
bound \(2-u-v\) and the coherence upper bound \(2\|C\|_1\) can be
attained by the same basis-preserving conditional phase. This checks the
conventional comparator for ordinary damping plus the phase family.
Population transfers conditional on the source, access to correlated bath
cells, early encoding and arbitrary nonseparable payload–programme
preparation would require a new optimum. The restricted jitter comparison
optimizes among the exposed conditional phase and exact idle at the declared
interface. It supplies no universal result about multi-pulse hardware.

Every fresh contact has a Stinespring dilation. A bath cell can be traced
after its only contact for the system predictions. A reused cell remains in
the joint state. The reduced maps from the initial state can remain CPTP
while failing semigroup composition. A finite-window generator uses fresh
product inputs and specified pulse-angle scaling. It is an effective
description under that preparation and timing, with no claimed stationary
field-theory derivation or renewal after depletion. [Attal and Pautrat](https://arxiv.org/abs/math-ph/0311002)

For initially factorized Gibbs contacts the audit checks
\(\Delta S_{\rm sys}+\beta\Delta E_{\rm bath}=I_{\rm out}+D(\rho_{\rm bath,out}\|\tau_{\rm bath})\).
The system side includes the source, record, payload and its reference for
this ledger. The identity follows from conservation of the full joint
entropy and the Gibbs logarithm. It concerns the finite contact, and does not
license a resetting rule. Bare energy balances because the resonant
exchanges preserve total excitation energy and B and M are degenerate.
Supplied scheduling, storage preparation, reservoir refresh and laboratory
controller work are separate. Heat flow from the hot supply to the cold
supply is a conventional conduction result. A reduction of payload
excitations and an improvement in transmitting an unknown state are
separate performance criteria. [Esposito, Lindenberg and Van den Broeck](https://arxiv.org/abs/0908.1125) ·
[Barsky-Giles et al., alternating hot and cold contacts](https://arxiv.org/abs/2602.18300)

## Separate reproduction and limits

The checker constructs the four-register \(B,M,Q,R\) density matrix for each
source preparation. It derives Kraus operators by exponentiating each
two-site interaction, keeping all initial source weights, and extracts each
Bell programme block after tracing B. It also constructs the matrix
Lindblad generator separately. There is no call to the production solver.
All 36 registered fresh coarsest-grid cases run. Anchor error is at most
\(5.45\times10^{-15}\). The maximum per-contact bare-energy residual is
\(2.42\times10^{-15}\), the entropy identity residual is
\(2.78\times10^{-15}\), and the minimum terminal eigenvalue is above
\(-3.0\times10^{-15}\). These are floating-point checks alongside the
analytic dilation, with no evidence for a positivity failure. The
reviewer's fresh solver does not test reuse or the geometry grids.
[Checker](reviewer_check.py) · [Output](REVIEW_CHECK.json)

For example, at sixteen slices, weak full-cold shallow fidelity is
0.35347704 versus 0.35238848 blind, and the weak deep pair is 0.29504815
versus 0.28807792. Both figures permit improvement in Bell score below
the 0.50 entanglement-witness threshold. A score below 0.50 by itself
leaves the actual separability question open. At angular jitter 0.4 the
exact idle removes the weak-memory advantage in these thermal cases.
The shallow and deep ideal anchors remain the strongest first laboratory
contrasts. A noisy contrast below the registered combined uncertainty
cannot adjudicate record usefulness on hardware.

## Organization, locality and verdict

Independent damage has an exact product probability for restoring all
cells and a binomial tail for a density predicate. A sharp change of
that tail with size comes from the chosen predicate applied to
independent cells. It supplies no collective singularity. The matched
two-cell test exposes bright and dark behaviour under a common bath
and states which throughput is matched in each arm. Neither match
equates every detailed interaction or proves a thermodynamic scaling
transition. The vacuum states are supplied resources.

The static XY test inserts one local excitation and measures where
its information is available later. Both the source perturbation and
the readout region are supplied. It tests retention and transfer,
without demonstrating unsupervised discovery of a repair rule. Matching
the incident interaction sum leaves edge structure and local energy
variance distinct. Relabelling the system and instruments gives an
exact counterexample to recovering spatial embedding from a task
score. Bounds in an interaction-graph metric already assume that
metric. CP-12's failure to derive locality therefore remains decisive
for its original claim. [Lieb and Robinson](https://www-thphys.physics.ox.ac.uk/talks/CMTjournalclub/sources/LR.pdf)

The integrated construction merits an effective-theory reading if
the complete production gates pass. Its gain, dissipation and negative
controls follow from the conventional null. The adopted causal dynamics
and the information access rule do not derive a matter EFT, dynamical
spacetime, a cosmological state or subjective experience. Candidate 1
remains an independent matter-sector investigation.
