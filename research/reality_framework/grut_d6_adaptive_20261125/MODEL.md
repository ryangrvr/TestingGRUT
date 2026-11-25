# GRUT-D6 · Adaptive development from an acquired internal record

25 November 2026. This finite construction extends the local moving-head
Hamiltonian of [D5](../grut_d5_prospective_20261124/MODEL.md). The head starts
without the information needed to select a successful layout. It probes a
conserved internal coupling sign, records the outcome, leaves the source
behind, and builds two different layouts. A second spin at the source evolves
independently of the head and serves as a held-out prediction target.

## Physical registers

Set ℏ = J = 1 for production, with J setting the maximum exchange scale.
Use n ≥ 2, m = 3n + 4 and L = 2m. The open chain has sites 0 through L.
Each site has a head degree of freedom with a vacuum and 48 occupied states.
Exactly one site is occupied in the prepared state. The head carries four
qubits P, M, A and Q and one three-state port μ.

| Location | Stationary degrees of freedom | Preparation |
|---|---|---|
| 1 | Source B, witness W, probe dump D, passive coin Cₚ, strength S | B = b, W = +, D = 0, Cₚ fair mixed, S sealed/weak/strong |
| 2 | Intervention mode F, scramble coin Cₘ, disable sink E | F as registered, Cₘ fair mixed, E = 0 |
| 4 through n+3 | n bank A display qutrits | All blank |
| 2n+4 | Stage challenge Z | Bit z independent of b |
| 2n+5 through 3n+4 | n bank B display qutrits | All blank |

Inputs b=0 and b=1 have equal prior probability. Every conditional run is
also tested separately, with z=0 and z=1. The head starts at site 0 with
P in |+⟩, M=A=Q=0 and μ in its registered basis state. The port μ takes
values main, passive and direct. Only the source starts correlated with
b. In particular the witness, the probe, memory, clock, geometry, modes,
challenge and tape have the same initial state for both b values. Coins
are independent fair mixtures in every arm.

## Local Hamiltonian

Let `U_j` be the local unitary listed below, acting on an occupied head
and the registers at sites j and j+1. For each basis matrix element of
`U_j`, put the corresponding transition between
`|head at j, vacuum at j+1⟩` and
`|vacuum at j, head at j+1⟩` in `T_j(U_j)`. The operator vanishes
outside that local occupancy pattern. Define

$$
\omega=\frac{4J}{L+1},\qquad \kappa=\frac{\omega}{4},\qquad
H=K+\sum_{j=0}^{L-1}\frac{\omega}{2}\sqrt{(j+1)(L-j)}
        \left[T_j(U_j)+T_j(U_j)^\dagger\right].
$$

Here K is the witness term at site 1. In the weak and strong sectors
it is `κ Z_B Z_W`. In the sealed sector it is `κ I_B Z_W`. The phase
assay and direct port also use S. B and S are stationary conserved
labels. No head transition acts on W. The same H includes all the
registered arms, whose settings are prepared independently of b.
Interactions are nearest-neighbor on the lattice of supersites, with
bounded local dimensions. An exchange can involve several internal
registers within its two sites.

| Link j | Local unitary, listed in the order applied on a forward hop |
|---|---|
| 0 | `exp(−i θ(S) Z_B Z_P)`, with θ=0, π/12, π/4 for sealed, weak, strong respectively |
| 1 | Rotate P by `R_x(π/2)=exp(−i π X_P/4)`, XOR P into M, then swap P with D. If μ is passive, XOR Cₚ into A. If μ is direct and S is weak or strong, XOR B into A. Otherwise leave A alone. |
| 2 | Apply the F-controlled action to M. Intact is identity, scrambled XORs Cₘ, inverted applies X, and disabled swaps M with E. |
| 3 through n+2 | One A cell per link. Exchange blank with color d in the cell on the arrival site, where d=M for μ main and d=A otherwise. |
| n+3 through 2n+2 | Identity. |
| 2n+3 | XOR the stationary challenge z into Q. |
| 2n+4 through 3n+3 | One B cell per link. Exchange blank with color `d XOR Q XOR P` in the cell on the arrival site, then apply X to P. P starts this stage in 0, having left its assay value in D. |
| m through L−1 | Identity. |

A display has distinct symbols blank, color 0 and color 1. Its exchange
leaves the third symbol alone, so it has a well-defined inverse on the
entire qutrit. The reverse hop applies each complete unitary's adjoint.
The two banks are physically separate. No B-cell write consults any
A-cell value. The stage changes when the head traverses the challenge
link, which supplies a new input to the second task without another
interrogation of B.

## Success and observables

The structural targets are A all color b, and B at index i (zero-based)
color `b XOR z XOR (i modulo 2)`. These encode two uses of the inferred
polarity, a uniform deposit and a challenge-dependent alternating
layout. The unknown sign governs an actual phase response. The designer
supplies its relationship to the layouts and to the second spin.

The physical forecasting target is the Y-basis outcome of W at
terminal time `T=π/ω`. Decode Y=+1 as 0 and Y=−1 as 1. W has not
interacted with any tape or head register. In an open port (weak or
strong) this terminal outcome is b. With S sealed it is 0 for both
source states. Report both forecasts and structural scores, including
the sealed mismatch between predicting a constant and learning b.

Define the recorded decision as the basis value of the bit used to
control writes. Define record accuracy separately as the chance that
M equals b. For j≥2 the pre-intervention value can be read on the
ordered history through the stationary D. After the subsequent
intervention, M remains unchanged under all forward write links.
For an unconditional duration-dependent score, include the head
projector onto the required stage. The main late scores count
unfinished banks as failures. (A) reaching its structural target
requires a complete A bank; (B) reaching its own target requires a
complete B bank. The joint late score counts structural completion of
both and can additionally include W agreeing with the decision. No
conclusion is drawn from discarding missed arrivals. Readout of D,
B and W is permitted to an evaluator and supplies no feedback.

At the stage diagnostic, record the probability of head index in
[n+3,2n+4], where A is complete and B is still blank. Its target
subscore also requires the A word to match b. At the terminal time
check the exact clock endpoint, both words, record retention,
forecast, and resource/energy measurements. Check the whole-layout
and whole-layout-plus-forecast probabilities at the flanking times
specified in [PROTOCOL.json](PROTOCOL.json).

## Channel boundary and timing

On the ordered path, after arriving at site 3 and before the first write,
the only mobile state that depends on b in a main arm is M. P was
swapped with the prepared blank D. A is still zero, Q is still zero,
and the port setting is fixed. The source, the witness, both coins and
the discarded storage remain at sites 1 and 2. No forward response
link addresses them. Later challenge z is independent of b. The passive
arm retains M but chooses the separate Cₚ-derived bit. The direct arm
adds an explicit source-to-A bypass at the last source link. It is a
diagnostic alternative channel, gated off in every sealed run.

A time-independent reversible walk can backtrack. This is a claim
about the terminal map and the projected complete response sector
of the declared unperturbed evolution. It does not assign a classical
monotone path to a superposition or forbid a new laboratory
intervention that re-enters the probe chamber. In the sealed sector,
the Hamiltonian itself factors as identity on B for all allowed
internal head states, including the direct-port selector. Operations
on this isolated B or changing the locked strength mode would be
additional access outside that control.

The initial localized head is the prepared physical clock. Laboratory
Schrödinger time and the readout times remain specified. No gates are
clocked or switched externally after preparation. The only
interventions in the protocol occur as stationary control links.
