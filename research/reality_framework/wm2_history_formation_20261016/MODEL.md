# Physical law, preparation and numerical representation

The physical Hamiltonian and preparation are exactly S4's, with common
g in [0.9,1.1] and independent storage detunings
delta_i=bar_delta+epsilon_i, |bar_delta|<=1/2, |epsilon_i|<=1/4.
See the immutable S4 MODEL.md at c6ac0c68480d4164ad620bbf2a4b3d8c65ac85b1.
S5 changes the computation, not the law. CP-12 remains a decoupled frozen
benchmark; the original charged CP-12 matter does not inherit this metric.
Units set hbar=1 and the internal receiver hopping energy to one.
The reported times are model units, with no derived laboratory conversion.

Every ordered pair has the same single-end bond append/pop grammar.
No event coordinates, supplied event adjacency, lattice, or pair
affinity is added. A numerical ramp of detunings assigns on-event
energies; it does not restrict their attachment order. Every ordering
and reverse operation is retained.

## Preparation, receiver and energy

The seed is the original neutral product preparation with a uniform
root, bondless one-head auxiliary seed, root storage bit one, other
storage bits zero, receiver vacuum and neutral-probe vacuum. The
receiver has L_N=64N^2 supplied bosonic modes, positive one-particle
energies in [2,6], and exact capacity N-1 through the conserved charge.
The g/N and 1/(4sqrt(N)) normalizations are supplied microscopic inputs.
Root averaging evaluates a block-diagonal observable in the original
coherent root preparation; it does not replace that state by a selected
root or a newly prepared dephased mixture.

In each root component r the initial energy mean is 16N-4+delta_r,
with conditional variance g^2(N-1)/N^2. The uniform-root state adds
Var_r(delta_r), bounded by 1/16. This is not an extensive hidden energy
source. Binding, storage, receiver and interaction energies are part
of one reversible closed Hamiltonian. Numerical projection of receiver
modes is an error-accounted approximation, not cooling or a physical
energy-transfer operation. The exact physical evolution conserves energy.

During formation the probe remains vacuum. Hence the probe-dependent
bare geometry interaction vanishes exactly and storage consistency is
invariant. The full formation component consists of all words and all
receiver occupations with k photons at word length k+1. Omitting
one-probe mismatch configurations here does not omit reachable states.

## What is exact in the receiver calculation

All L_N physical frequencies remain in the moving one-particle vectors.
Free receiver evolution is the complete finite-chain unitary in its
discrete sine basis. All type words are exact. The equal-detuning class
isometry is exact, including original root weights. No photon-number
cutoff beyond conserved charge is used.

The approximation consists of product-formula time steps, a bounded
finite exponential approximation, and projections of receiver orbitals.
MATHEMATICS.md bounds their complete accumulated state error, including
the future effect of every omitted component. The active mode count is
a numerical resource; it is not the physical receiver's mode count.
An unsuccessful resource-limited run is not a physical formation failure.

## Recording and the joint response protocol

At T_N=4N^2, accept every word of length at least max(2,ceil(N/2)).
Record its root and complete ordering. The receiver is not measured,
refreshed or selected. Inject S4's normalized root-pair probe, then
apply a local number phase at the root and measure the adjacent
relational receiver under the complete one-probe Hamiltonian.
Rejected trials contribute zero signal. Original root/formation
probabilities are retained in the response per original trial.

Recording and probe injection are supplied interventions. The word
measurement dephases only the emission interaction in the formation
component. Its mean work is -<H_emission>, with magnitude <=g. Its
full outcome entropy is

    log2(N) + (1/N) sum_roots H(word | root),

between log2(N) and log2[N sum_(k=0)^(N-1) (N-1)!/(N-1-k)!]. It is
not zero device cost, and no reset-cost model is derived. Probe
preparation supplies number variance 1/4 and mean hopping energy -J;
its apparatus must account for that exchange.

S4's post-recording theorem applies to every conditional receiver state
and every accepted word. A valid S5 formation lower bound can therefore
be joined to it without a favorable field reset or word postselection.
Such a demonstration is finite-instance unless formation is controlled
uniformly with increasing N.

The single-end grammar still supplies the path geometry class and its
dimension. No common geometry for original charged matter, selected
dimension, gravitational continuum, equilibrium phase, observer law or
distinctive empirical prediction is established here.
