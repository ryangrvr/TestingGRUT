# GRUT-D7 · Physical recovery using an acquired internal record

25 November 2026. [D6](../grut_d6_adaptive_20261125/REPORT.md) measured
the agreement of two layouts and a forecast. Here a learned record
sets a phase compensator before a new quantum payload is loaded.
The outcome is the physical state transferred to a receiver. Both
the disturbance and the compensation are inside one finite
time-independent Hamiltonian.

## Preparation and geometry

Use one occupied head on a folded lattice path of length $L\ge8$.
A path site has a vacuum state and 48 occupied states. The head
carries probe P, memory M, alternative selector A, payload Q and
a three-state port μ. Start at site 0 with P along +X, M=A=Q=0,
and μ in the selected main, passive or direct mode. The reference R
at position (−1,−2) has no dynamics and is initially in
$|\Phi^+\rangle=(|00\rangle+|11\rangle)/\sqrt2$ with the input I.
An evaluator retains R. The receiver O starts in zero.

| Clock site | Position | Additional stationary registers |
|---|---|---|
| 0 | (0,0) | — |
| 1 | (1,0) | Source B, probe dump D=0, fair passive coin, 3-valued assay strength and 2-valued phase severity |
| 2 | (2,0) | Fair independent mask coin, disable sink E=0, 4-valued memory intervention |
| 3 | (3,0) | — |
| 4 | (3,1) | — |
| 5 | (2,1) | Ready compensator bit G=0, input I, a 4-valued compensator angle setting |
| 6 | (1,1) | — |
| 7 | (1,2) | — |
| 8 | (2,2) | Receiver O=0 |
| j≥9 | (j−6,2) | — |

The source B is a conserved qubit prepared as $|b\rangle$, with
$b=0,1$ equally likely. Only B is initially correlated with b.
All settings carry the same values in the two source conditions
of a trial. Both mixed coins are independent. The assay angle θ is
0, π/12 or π/4 for none, weak and strong. The payload sees a
separate rotation of size φ=π/8 or π/4, set independently of b.
The balanced results average both source conditions.

## Time-independent exchanges

Use ℏ=J=1 for the measurements. Define
$\omega=4J/(L+1)$ and $a_j=(\omega/2)\sqrt{(j+1)(L-j)}$.
Each forward exchange $T_j(U_j)$ moves the head from j into
the vacant site j+1 and acts by the listed unitary on the
data registers. It vanishes outside that occupancy pattern.
The Hamiltonian is

$$
H=\sum_{j=0}^{L-1}a_j\big(T_j(U_j)+T_j(U_j)^\dagger\big).
$$

The reverse exchange applies the adjoint. Each term is supported
on sites of Manhattan diameter at most 2. Links 3, 5 and 6 access
an extra adjacent register site, as shown by the geometry above.
Each controlled unitary is a local multiregister coupling.
All arms at a fixed L are sectors of the same H with stationary
or conserved control settings. No timed pulses act during a trial.

| Link j | Forward data unitary |
|---|---|
| 0 | Apply $e^{-i\theta Z_BZ_P}$ to assay P. |
| 1 | Rotate P by $e^{-i\pi X_P/4}$, XOR its pointer into M, then swap P with D. The passive port XORs its separate coin into A. The direct port XORs B into A. The main port leaves A zero. |
| 2 | Intact leaves M alone. Scrambled XORs the mask coin into M. Inverted applies X. Disabled swaps M with E. |
| 3 | XOR the selected bit into G at site 5. The selector is M in the main port and A otherwise. This fixes the forward compensator before any operation on the payload. |
| 4 | Swap the input I at site 5 into the mobile payload Q. R remains isolated. |
| 5 | Apply the physical phase disturbance $e^{-i\phi Z_BZ_Q}$ via B at site 1. This link touches Q, not any decision bit. |
| 6 | Apply $e^{+i\beta Z_GZ_Q}$ via the previously set G at site 5. This link acts on G and Q and has identity action on B, M, A, and R. |
| 7 | Swap Q into receiver O at site 8. |
| j≥8 | Identity data links. |

The actuator angle β is set at preparation independently of b.
The blind choice is β=0, the full correction uses β=φ, and the
calibrated correction uses
$\beta=\tfrac12\arg(\cos(2\phi)+i\sin(2\theta)\sin(2\phi))$.
Use 0 when both components vanish. The resulting four physical
angles are 0, $\tfrac12\arctan(1/2)$, π/8, π/4. The design supplies
the channel family, magnitudes and phase response law. Only the
sign is acquired. The compensator keeps its state on the ordered
forward path. The reversible clock can still backtrack in
laboratory time. Evidence about causal order comes from the
endpoint propagator and sector analysis.

## Readouts and access boundary

Measure the projector onto $|\Phi^+\rangle$ for O and R. The
entanglement fidelity is its probability. The projector and
target are identical for both b values. There is no source label
in the scoring rule or feedback from the reference. Measuring
$X_O X_R$, $Y_O Y_R$ and $Z_O Z_R$ on separate preparations
implements the same statistic through
$(1+\langle XX\rangle-\langle YY\rangle+\langle ZZ\rangle)/4$.

The attempted communication input is independent of the assay and
its controls. Up to commitment, the primary actuator can use M
only. The passive coin and direct source port are explicit
alternative selectors. During recovery the registered coupling
uses the received Q and G only. The independent reference, source
and discarded storage are excluded from any recovery optimization.
The bound may grant the recovery the entire available programme
state and arbitrary quantum processing jointly with Q. Opening
the source port or accessing R would change the benchmark.

The perfect-transfer time is $T=\pi/\omega$. No clock reading
is used to filter receiver data. At arbitrary time the released
receiver score includes all unfinished branches, whose receiver
is still zero and whose reference has marginal $I/2$. Record
the head's ready probability $j\ge8$ separately. The final readout
window has $|\omega t-\pi|\le0.2$. At $2T$ the state recurs and the
receiver no longer holds the payload.
