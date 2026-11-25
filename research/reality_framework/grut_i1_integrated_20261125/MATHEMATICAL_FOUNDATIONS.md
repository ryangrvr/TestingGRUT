# Mathematical foundations of GRUT-I1

## 1. Finite causal physical model

Use the specification
$$
\mathfrak Q=(\mathcal H,\mathcal A,\rho_0,H_0,H(t),
              \mathcal I,\mathcal K,\mathcal T,\mathfrak B).
$$
Here $\mathcal A=\mathcal B(\mathcal H)$ and $\rho_0$ is a
positive unit-trace state. $H_0$ is the reference bare energy
used for the resource ledger. $H(t)$ generates the ordered
microscopic unitary. A time-independent Hamiltonian suffices
in the existing closed constructions and the network test.
The collision experiment uses explicitly timed local interactions.
The clock, gates and any trajectory conditioning are physical
inputs to each model.

$\mathcal I$ declares the accessible subsystem embeddings and
instruments, $\mathcal K$ is the permitted time order and ports,
$\mathcal T$ is the scored task and $\mathfrak B$ is the budget.
In the finite hardware below, subsystems are tensor factors.
More generally one can specify a concrete unital subalgebra
and restrict the state to it. A complement or commutant does
not automatically constitute a product thermal reservoir.
Gauge constraints can obstruct naive tensor factorization.

An instrument with branch maps $\{\mathcal J_y\}$ consists of CP,
trace-nonincreasing maps with trace-preserving sum. Subsequent
instruments and microscopic intervals compose chronologically.
Each branch probability is the trace of its composed map on
$\rho_0$. Equivalently the experiment pulls the final effect
back through the adjoints. Summing over later trace-preserving
choices leaves the earlier probabilities unchanged. The
restriction is operational causality in the declared process.
No general space-like locality follows without a geometry,
an interaction bound and a limit that preserves the bound.

The explicit I1 implementation has factors B (physical phase
source), P (probe), ready memory M, probe dump, payload Q, isolated
reference R, used and unused local reservoir cells, plus an optional
control-error bit. The initial states, preparation and all local
interactions are specified in the frozen protocol. The scoring
instrument acts on Q,R only. Applying a CP reduced model to this
episode is justified by the bath preparation and the stated trace,
so arbitrary initial correlations with those same baths would change
the derivation. [Protocol](REALITY_TEST_PROTOCOL.md) · [Effective theory](EFFECTIVE_THEORY.md)

## 2. Records, generation and observation

- A **historical record** is a physical subsystem with a conditioned
  family of states indexed by an earlier alternative $x$, with a
  stated access map and distinguishing effect. The history label
  can be retained in a reference for calibration. It is not part
  of the controller unless the declared port admits it. Retention
  compares accessible ensembles after a stated evolution and time.
- **Generative information** is carried by an initial parameter
  or state whose intervention changes the statistics of a later
  configuration under a fixed process. D5 changes the seed or rule;
  a compact generator can grow a long output while supplying the
  grammar and route. Exact configuration and coarse configuration
  correspond to different effects. A chosen success projector is
  an operational definition, not a dynamical explanation.
- **Predictive information** is tested against an independent
  target intervention or output with a known chronologically
  allowed calibration. In D6 the isolated witness supplies a
  separate physical readout. In D7 and I1 the scored Bell effect
  is independent of the hidden label. Prediction can be trivial
  if the target is constant, as the D6 sealed witness shows.
- A **physical response** is a causal coupling from the accessible
  state into a downstream interaction. Subjecting different
  memory states to the same actuator tests the channel it uses.
  A side port can provide equivalent or better information.
- **Internal observation** means that the transformed state of
  an exposed physical register predicts an operational statistic,
  and an allowed internal instrument can use that register.
  This is the finite measurement construction. It does not select
  a fundamental factorization, physical clock or outcome law.

Energy follows the declared H and control schedule. A local resource
budget can include preparation, ancilla spectra and counts, time,
support, work, entropy, readout and available communication. Finiteness
and unitarity preserve the global spectrum and von Neumann entropy.
Reduced entropy need not be conserved or monotone. Gibbs states and
low-entropy preparations require separate physical choices.
Neither closure nor recurrence forces an arrow in every state.

## 3. Operational value of memory

Let $T$ denote a task and its duration, prior, allowed operations
and normalized payoff effect $0\le E_T\le I$. Write $A$ for
all baseline accessible data, M for the additional port, and X for
any inaccessible reference labelling the experiment. At fixed
physical preparation and resource budget put
$$
 V_T(A)=\sup_{\pi\in{\cal P}(A;\mathfrak B)}
             \operatorname{Tr} E_T\,\sigma_\pi,\qquad
 R_T(M\mid A)=V_T(AM)-V_T(A).
$$
In a multistage experiment $\pi$ is a causal admissible policy,
not an arbitrary dependence on a later label. Blind policies
ignore M or use an admissible independent dummy. An intervention
that destroys M, deposits heat or changes the permitted hardware
defines a related experimental contrast. It need not equal this
difference of optima.

**Order and bounds.** For nested feasible sets, $0\le R_T\le
1-V_T(A)$. A score difference for a fixed actuator can be negative.
If the blind interface can freely prepare a dummy state $\nu_M$
and emulate the same downstream operations on it, then for the
joint input $\omega_{XAM}$
$$
 R_T\le \tfrac12\|\omega_{XAM}-\omega_{XA}\otimes\nu_M\|_1.
$$
The difference for each fixed policy is bounded by the trace
distance before optimizing, by the variational bound for effects
and channel contraction. With no additional A and an admissible
$\nu_M=\omega_M$, this is at most
$\sqrt{I(X:M)_\omega/2}$ in nats by quantum Pinsker.
There is no universal tight conversion from memory entropy
to task performance or from $R_T$ to a work rate.

For binary guessing with equal priors and arbitrary allowed POVM,
$$
 R_T=\frac14\big(\|\rho_0^{AM}-\rho_1^{AM}\|_1
                   -\|\rho_0^{A}-\rho_1^{A}\|_1\big).
$$
It follows by maximizing the effect for one guess on the positive
part of the state difference. In general tasks the complete
quantum process, possible reference systems and physical operation
budgets determine the effects. Discrimination is one special case.

**Sufficiency.** If an admissible CPTP map from A simulates the
complete needed family $AM$, including the correlated reference
for the registered experiments, then $R_T=0$. Simulation composes
with the optimal assisted policy. Data processing contracts
accessible task performance under the same feasibility convention.
A single zero $R_T$ does not imply sufficiency for all tests.
Comparisons over full decision families and quantum ancillary
systems are the standard statistical-comparison problem.
[Buscemi](https://arxiv.org/abs/1004.3794) ·
[Chiribella, D'Ariano and Perinotti](https://arxiv.org/abs/0904.4483)

**Composition.** With nested access to M then N at fixed
feasibility, $R_T(MN|A)=R_T(M|A)+R_T(N|AM)$ by telescoping.
There is no general addition of the marginal gains. Duplicate
bits are redundant. Two independently fair bits whose XOR is
the target have zero individual guessing value but joint gain
$1/2$. For independent single-shot discrimination tasks and a
payoff requiring both labels correct, joint POVMs give product
optimal success. Product POVMs attain it, and if $Y_i$ are
optimal SDP duals with $Y_i\ge p_{x_i}\rho_{x_i}\ge0$,
then $Y_1\otimes Y_2$ is dual feasible. The difference of
products is not a sum of gains. Correlated tasks, new access,
shared budgets or interference invalidate the product argument.

## 4. Apply the definition to D3–D7

These rows fix the decision cut, the payoff and the competitor.
Their actual protocol scores have the additional conditions in
the original packets. The cut-level tasks are distinct from joint physical arrival,
archive retention and finite-time execution.

| Original packet | Declared optimal task and result | Configured physics and limits |
|---|---|---|
| [D3](../grut_d3_record_repair_20261124/REPORT.md) | Fair earlier label and common erasure. With M alone at the repair cut the write-erasure probability is $\lvert a_0\rvert^2$. Blind success is $1/2$ and $R=(1-\lvert a_0\rvert^2)/2$ for optimal label restoration. Adding unrestricted coherent access to the writer W changes it to $R=\sqrt{1-\lvert a_0\rvert^4}/2$. The W marginal alone is identical for both labels. | The finite Hamiltonian's joint restoration and retention bound is about 0.991568. The matched shuffle halves its configured score. That physical result includes the actual arrival and resource constraints. |
| [D4](../grut_d4_internal_recognition_20261124/REPORT.md) | For one unflagged flip in S, M1 or M2 and a fair historical bit, the three current bits identify it by majority. S alone has success $\max(p,1-p)$ with $p=\Pr(\text{S is flipped})$. Thus $R=\min(p,1-p)$ for this cut and unrestricted label correction, giving $1/3$ for equal fault locations. | Bounds for up to four separated faults use the supplied writer and finite baths. Simultaneous faults collide in accessible state and can defeat history restoration. Knowing S always flipped can make the blind optimum perfect. |
| [D5](../grut_d5_prospective_20261124/REPORT.md) | For whole-word completion with supplied G and accessible source pilot h, the direct writer reaches the intact optimum even without M. The pilot-enabled access admits $R=0$ for the optimal word task. Erasing the pilot as well would change the cut and the task's information problem. | The original fixed writer uses the copied M. Its intervention effects are large and preserved, including finite early-fault resumption, yet the direct pilot is an explicit sufficient alternative. |
| [D6](../grut_d6_adaptive_20261125/REPORT.md) | At its exact main response cut the weak pointer has contrast $1/2$. With isolated source and witness and the same decision controlling both layouts, $V=3/4$, $V_\mathrm{blind}=1/2$, $R=1/4$. Strong gives $R=1/2$. | The direct port and pre-probe clock coherences are separate access regimes. With sealed source the layout advantage is zero, even though a constant witness can be predicted perfectly. |
| [D7](../grut_d7_functional_recovery_20261125/REPORT.md) | Recover a phase-damaged qubit and its Bell reference through the specified interface. For the weak record, $R=0.041731316928$ shallow and $R=0.25$ deep. Both blind and assisted scores optimize downstream channels. | The projector is independent of the hidden sign. No extra query or pre-channel encoding is allowed. A coherent probe has the same optimum as this pointer for the registered family. |

For the D3 coherent entry the conditional pure cut states have the
form $a_0|\mathrm{vac}\rangle_M|0\rangle_W+
\lvert b\rangle_M\lvert\chi\rangle_W$, with identical W marginals and
overlap magnitude $|a_0|^2$. The binary pure-state trace
distance gives the formula. The original record-plus-arrival
score and possible readout restrictions remain in the D3 packet.

The shared object is a normalized physical decision problem on
an accessible quantum process. It connects the preparations, records
and responses when their cuts can be composed. It supplies neither a
single scalar notion of organization nor a selection rule over
physical laws. [Evidence map](EVIDENCE_MAP.md) · [Novelty audit](NOVELTY_AUDIT.md)
