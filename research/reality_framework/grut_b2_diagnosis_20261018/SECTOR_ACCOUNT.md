# Complete sector and source-flow account

Unchanged B2, fixed eight-event host and original failed theta=6 endpoint.
Full 773-state coherent propagation; no sector measurement or classical
rate replacement. Raw 193-point data are in `RAW_0.json` through
`RAW_6.json`; `SECTORS.json` includes integer-time summaries and balances.
Quantum identities are in `IDENTITIES.md`. Numerical values below are
evaluations under the pinned B2 arithmetic contract, not external validation.

## Final sectors and the fuel correction

| Mutually exclusive sector | Probability at theta=6 |
|---|---:|
| No matter, both fuels excited | 48.889814% |
| One matter-pair quantum, left fuel excited | 21.571333% |
| One matter-pair quantum, right fuel excited | 21.571333% |
| Two matter-pair quanta, both fuels empty | 7.967519% |

Each fuel's marginal occupation is 70.461148%. **The joint both-full
probability is 48.889814%, not 70.46%.** This corrects the marginal/joint
ambiguity in the recent discussion without editing the archive.

q=N_+=N_- counts matter quanta, not necessarily fully co-located composites.
At this endpoint, 43.127186% is in a one-co-pair/one-fuel configuration;
7.961763% contains two co-located pairs. Other q=1,2 configurations retain
separated charged constituents. The prescribed outgoing event is only
1.932705% and is nested inside two-co-pair configurations. None of these
numbers is multiplied by another success or readability estimate.

The mean remaining matter-pair count is 0.590777. Total resource is exactly
F_L+F_R+q=2; energy is conserved. Bare excited configurations have probability
0.00021238, while the initial high SPECTRAL band weight is much smaller.
Those are different projectors: virtual charged configurations are part of
the low-energy eigenstates and have not been removed from propagation.

## What refills the fuels?

For either endpoint use L_e=F_e+(n_+,e+n_-,e)/2. Conversion cannot change
L_e; only hopping through its boundary bond changes it. The exact balance
from theta=2 to theta=6 is:

| Change per endpoint | Left | Right |
|---|---:|---:|
| Fuel increase | 0.70348587 | 0.70348587 |
| Net resource returned across endpoint bond | 0.69012938 | 0.69012938 |
| Matter change already at endpoint | -0.01335649 | -0.01335649 |

Delta F=Delta L-Delta M closes to below 2e-16 in evaluated arithmetic.
Thus local conversion of matter already sitting at the endpoint cannot
account for most of its net fuel refill. The dynamics necessarily includes
substantial net spatial return to each source region. This is a balance of
expectations, not a probability assigning an individual path to every particle.
It does not split coherent collision and reflection mechanisms additively.

The exact zero-hopping isolated-cell reference retains the original fuel
energy, counterterm and lambda. It gives p(q=2)=0.201974 at theta=2 and
0.041904 at theta=6, compared with physical values 0.997750 and 0.079675.
It cannot reproduce the spatial resource return because L_e stays exactly
one in that reference. A resemblance in late populations is not evidence
that the physical loss was simply an isolated-cell Rabi oscillation.

At the failed endpoint, I01=0.794544 and I12=0.329483 in theta units:
both net matter-creation currents are already positive again. Resource
flux across each source bond is outward again, about 0.044724. This is a
recurrent conversion/transport cycle passing through an unfavorable endpoint,
not evidence that the fuels have become permanently trapped full.

## Coherence and causal scope

The full complex z01 and z12 are evaluated before taking their imaginary
parts for sector current. All within-sector wavefunctions continue under
the complete Hamiltonian; fields, source occupations, collisions and return
amplitudes remain coherent. `RAW_*.json` records these complex coherences
and source return amplitudes, not merely a sequence of population updates.
The pinned law, input and code reconstruct the full states; the files archive
observables rather than every coefficient of every wavefunction.

This source-flow result distinguishes required spatial backflow from
conversion-only motion. It does not uniquely attribute the backflow to
collision, finite-boundary reflection or source hybridization. B2's contact-
cancelled reference also fails the same endpoint gate, so contact repulsion
alone was never established as the sole cause. No Hamiltonian was revised.
