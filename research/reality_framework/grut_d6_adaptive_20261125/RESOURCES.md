# D6 resource accounting

## Hardware and preparation

| Resource | Count or setting | Role |
|---|---:|---|
| Chain sites | L+1 = 6n+9 | One occupied site carries a moving head |
| Local head states per site | 49 | Vacuum plus 48 occupied labels |
| Mobile registers | P, M, A, Q and port μ | Four qubits and one qutrit |
| Stationary site 1 | Dimension 48 | Source, witness, dump, passive coin and strength |
| Stationary site 2 | Dimension 16 | Four-valued control, scramble coin and blank sink |
| Stationary stage challenge | One qubit | Independent second input |
| Display media | 2n qutrits | Initially blank substrate |
| Dissipative reservoir | None | Correlations remain in the closed system |

The full tensor-product dimension is
$49^{L+1}	imes1536	imes 3^{2n}$. The one-head sector is
$(L+1)	imes48	imes1536	imes3^{2n}$. Each link
uses at most two neighboring supersites and bounded local dimensions.
The assay, copying and port controls are explicit multiregister
unitaries on these two sites.

Preparation supplies the head at the endpoint, the two known coherent
spin states, blank storage and an unoccupied route. The chosen modes,
geometry, coupling strengths, Y decoder, grammar, and relationship
between the assay and witness signs are engineering inputs. The mode
and challenge carry no initial information about b. The source carries
one hidden bit under the fair prior. The mixed comparison coins supply
two bits of entropy in every conditional run. If b and z are also
averaged independently and fairly, the full ensemble entropy is four
bits. Fixed settings and blanks add no ensemble entropy.

## Energy and work

At fixed n all arms share the exact endpoint-preparation energy law

$$
E_{k,w}=\omega(k-L/2)+w\kappa,\qquad
\Pr(k,w)=2^{-(L+1)}\binom{L}{k},\quad
k=0,\ldots,L,\quad w=\pm 1.
$$

The mean is zero with the chosen energy reference and the variance
is $\omega^2 L/4+\kappa^2$. The witness energy is initially symmetric
because it starts along +X. The law does not depend on the sign of
its coupling. Every exchange term has norm at most J and the one-head
Hamiltonian has norm less than 2J. The episode requires no externally
applied work after preparation. Choosing the initial state, building
the inhomogeneous couplings, maintaining their precision and
reading or reinitializing the device have separate costs.
The model leaves their laboratory work budgets unspecified.

## Information storage and resets

For a fixed source the phase decoder and dump retain a two-branch
pure state across the record and stationary dump. Tracing the dump
gives a classical pointer with conditional entropy
$h_2((1+\sin(2\theta))/2)$. The fair-source marginal of M is one bit
in the intact arm. Its mutual information with b is
$1-h_2(p)$, or 0.1887218755 bits at weak strength, one bit at strong
strength and zero at sealed strength. Retaining this pointer through
the change of stage requires no new source query or new memory bit.
The two layouts carry the same inferred source bit, so the two
successes are correlated. The stage challenge also determines the
second word. Its conditional contribution is supplied independently.

The dump D and optional sink E are each one qubit. They are blank
resources used in the reversible implementation. The scramble coin
retains the mask needed to undo the scramble. The second spin W is
an inaccessible verification resource and has a full bit of
information about b at T in the open arms. It is independent of
the response apparatus conditional on b. No transfer couples it
to the head. Separating these resources follows from the Hamiltonian
and the readout algebra, including coherences in the late clock sector.

The total von Neumann entropy is conserved. No finite bath
thermalizes, and no new syndrome or error reservoir is needed for
this one acquisition. The record has no guaranteed long-term
storage under the same freely evolving dynamics. At 2T the path
reverses all its data transformations and the witness differs from
its initial state. The entire joint density operator returns at 4T.
None of these recurrences gives accumulated knowledge across cycles.

At bounded J, the clock length, time $T=(L+1)\pi/(4J)$ and blank
substrate grow linearly in n. One acquired bit controls one polarity
through both banks. The apparatus supplies route, two grammars and
stage change, so output length does not measure the number of
learned parameters. Coupling perturbation tolerance and equilibration
of the head have no certification in this packet.

[Exact proof](PROOF.md) · [Computed moments and dimensions](RESULTS.json)
