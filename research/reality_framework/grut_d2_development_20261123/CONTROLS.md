# Controls, resources and retained failures

All scores use the frozen model and `PROTOCOL.json`. The correction to one
control forecast is explicit in [CONTROL_CORRECTION_01.md](CONTROL_CORRECTION_01.md).
The mathematical proofs concern the unmodified model; none promotes it into
target-free development or a new law of physics.

The uniform-window guarantees concern a single read at any chosen time in
that window. They are not a joint probability that a continuously monitored
trajectory never leaves the aligned family. Repeated measurements require
their backaction and apparatus resources to be included separately.

| Comparison | Held fixed | Changed or spent | What it establishes |
|---|---|---|---|
| Connected model / interface removed | Hilbert dimension, full initial state, physical chain, reservoir length, onsite energy, all other hoppings | Remove one exchange edge; initial mean energy identical, initial variance lower by g^2/2 | The ordering and count-transfer pathway requires access to this reservoir; this does not isolate memory from cooling |
| Prepared vacuum / maximally mixed active reservoir | H and dimensions | Bath entropy, mean energy, initial resource; reference correlations remain if its tag is retained | No ordering from a maximally mixed active marginal; tagged reference information can still be redistributed |
| Baseline / staggered bath phase at tau | Same H before and after the kick; same present physical-system state; same joint occupation and count readouts | O(L) prescribed coherent phase control and external timing; bath coherences change; variance generally changes, though mean work is zero here | Present system state and the passive count record do not suffice to predict its future |
| Baseline / single physical endpoint flip at tau | Same H and reservoir, no bath reset | Timed physical X kick, approximately Delta of positive mean work | Autonomous recovery after a declared disturbance, using remaining reservoir capacity; no cost-free regulation claim |
| Short bath | Same local rule at n=2,L=1 | Far shorter reservoir and storage lifetime | Exact Rabi recurrence defeats a permanent-attractor interpretation |

## Removing the interface

For the initially maximally mixed chain, internal wall hopping preserves the
identity state. With the reservoir disconnected, `P_order=2^(-m)` at all times,
`N_B=0` and `I(K:N_B)=0`. Its initial mean energy is unchanged. For fair walls,

\[
 \operatorname{Var}H_{\rm on}=\Delta^2m/4+g^2m/2,
 \qquad
 \operatorname{Var}H_{\rm off}=\operatorname{Var}H_{\rm on}-g^2/2.
\]

Thus the complete energy distribution is not matched. The fixed decoder
`K_guess=N_B` succeeds with probability `2^(-m)` initially or disconnected.
This fixed-decoder number is not the optimal no-information guessing baseline
for a nonuniform binomial count distribution. The information baseline is
exactly zero. Merely dephasing a classical number pointer in its own basis
would not remove that pointer's record and is not used as a disabling test.

## Bath phase intervention

With wall sites indexed `0,...,m-1`, full-path staggering is
`D_j=(-1)^j`. Apply its restriction to the bath only; on the system it is the
identity. This unitary commutes with every occupation/count readout, including
the reference tag. It leaves the complete current physical-system marginal
unchanged, and changes some off-diagonal bath coherences. There is no claim
that the complete bath marginal is unchanged.

If `A=U_SS(tau)` and `D_S` is the system part of full staggering, the system
block at `2*tau` after this intervention is exactly

\[
 A_{\rm echo}=D_S+A(I-D_S)A.
\]

The almost-empty source refills, approaching its initial order probability
`2^(-m)`. The count was still present immediately after the kick. Its values
alone therefore fail to describe the reservoir influence. This is coherent
information return, not a decision made by a classical developmental memory.
No readout-controlled actuator or learned rule was introduced.

The initial covariance has purely imaginary coherence on nearest-neighbor
links, so every local hopping energy has zero mean before and after these
diagonal phase flips. Mean work is zero for this preparation. The intervention
can change the energy variance; `RESULTS.json` records it for every size.
Zero mean work does not supply a free physical controller or clock.

## Disturbance, work and record

At `tau=8*S2/g`, flip the last physical spin. It flips the last domain-wall
occupation, without resetting the bath. Its mean work is
`Delta*(1-2*<w_last>(tau))`, with vanishing hopping mean before and after the
kick in this preparation. The exact evolved expected number of walls gives
`P_order >= 1-<N_wall>`. A separate uniform analytic bound gives recovery
above `0.9922` at `2*tau`; the direct samples give stronger lower bounds.

The post-kick state can have anomalous fermionic correlations, so the
unperturbed Gaussian determinant is **not** applied to it. The derivation
uses only the exact normal covariance for the expected number, then a union
bound. The old count record is not assumed unchanged: the kick usually adds
one total excitation. Guessing the shifted count `K+1` from the reservoir
has lower bound `1-<N_wall>_final-<w_last>_before`. This weaker record bound is
reported separately, with an analytic value above `0.9891` for the family.

## D1 corrections retained

The three-CNOT D1 endpoint is exact. For its specified coherent pure input,
the regulated qubit has von Neumann entropy zero both initially and finally.
Its designated computational-basis outcome entropy decreases; that is a
different quantity. The final memory-actuator pair is pure, with entanglement
entropy `h2(p)` between its members. For a diagonal mixed source it instead
has classical correlations and no entanglement. Gate scheduling and blank
registers were supplied resources. An unused reference register was not a
clock. D2 supplies autonomous baseline dynamics but keeps Schrodinger time
external and explicitly budgets its reservoir preparation.

## Negative claims and limits

The preferred alignment family, conserved orientation seed, local interface,
initially pure bath, coarse readout, dynamics and preparation are supplied.
No microscopic target-free emergence theorem, active record-dependent
regulation, emergence of an observer, distinctive GRUT prediction, cosmic
selection rule or consciousness test has been earned. No asymptotic attractor
or exact new plateau is possible in this finite closed system. Candidate 1's
quantum-binding calculations, parameters and conclusions are unchanged.
