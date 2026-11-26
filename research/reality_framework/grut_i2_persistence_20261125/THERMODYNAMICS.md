# Operational thermodynamic ledger

All heat and work are expressed in units of $k_B T_b$.
Entropy and information use natural logarithms. The solvent
is a single fixed-temperature bath. The stock of excited
chemical molecules is finite. The source and receptor are
degenerate, so their purely thermal flips exchange no bare
energy. Driven receptor flips spend or recover the same
physical stock as the rebinding channel.

## Exact energy and entropy balance

For the count probabilities p, the physical microscopic
entropy is
$$
 S_{\rm mic}(p)=H(p)+\sum_z p_z\log g(k_z).
$$
The second term is the entropy of the uniform labelled-fuel
states inside a count class. Omitting it would misprice the
medium entropy and a charged cell's reset. The mean solvent
heat is positive into the bath and obeys
$$
 Q_b=\langle E\rangle_0-\langle E\rangle_\tau
 =J(\langle B\rangle_\tau-\langle B\rangle_0)
   +A(\overline{N_{\rm repair}}+\overline{N_{\rm sens}}).
$$
Here the signed net conversion counts sum to
$K-\langle k\rangle_\tau$. Positive conversions use F.
The reverse jumps can synthesize F, so a temporarily empty
stock is not an absorbing state. At most K net conversions
can be supplied by the initial stock. At fixed rates there
is no steady dissipative flux after equilibration.

For any of the fixed control settings detailed balance gives
$$
 \Sigma(\tau)=S_{\rm mic}(p_\tau)-S_{\rm mic}(p_0)+Q_b
 =D(p_0\Vert\pi)-D(p_\tau\Vert\pi)\ge 0.
$$
The identity follows directly from the fixed Gibbs
law. Differentiating D gives the signed current expression
on each reverse pair. Its loss rate is a sum of terms
$(u-v)\log(u/v)\ge0$ for directional probability fluxes
u,v. The generator proof covers all values of the
registered rates, before numerical propagation. It does
not assert that each channel's heat or marginal entropy
change is positive. Count-level medium entropy equals
$Q_b+\Delta\langle\log g(k)\rangle$.

## Preparation, acquisition and storing the record

The starting gas law for persistence is $\pi_X$ restricted
to O. Its excess free energy is $-\log\pi_X(O)$. A
uniform gas start has cost $D(\mathrm{unif}\Vert\pi_X)$.
The K initially all-F molecules cost
$K\log(1+e^A)$ relative to their equilibrium. Both source
and receptor begin in the fair independent Gibbs
mixture. It supplies zero information about e and no
extra nonequilibrium bit preparation. Prepared gas
cohesion is a physical binding energy, not a labelled
target pattern.

During acquisition the receptor tracks e by reversible
fuel-coupled flips. The signed acquisitions are in
$N_{\rm sens}$. The receptor and gas share a finite
stock, so this use competes with direct rebinding.
Its thermal flipping and environment switching set
record retention and appear in the same joint law.
With the sensing channel isolated, the receptor's
distinguishability in this degenerate two-state thermal
model decays as $e^{-2\nu t}$. Changing the repair
barrier through m raises or lowers forward and reverse
attempt rates together. Its energy transfers are charged
through $N_{\rm repair}$. There is no separately measured
binary outcome written into an unpriced register.

## Reset and minimum work

An ideal isothermal reset restores the **whole input
distribution**, including gas and fuel. The original
Hamiltonian is recovered. For the joint reset the
work bound in these units is
$$
 W_{\rm joint}=D(p_0\Vert\pi)-D(p_\tau\Vert\pi)=\Sigma(\tau).
$$
It is a thermodynamic infimum with arbitrary
distribution control and quasistatic time after the
scored episode. To attain the infimum, for a microscopic
distribution r first quench to Hamiltonian
$-\log r$ up to a constant. Connect that equilibrium
slowly to the Hamiltonian $-\log r_0$ and quench
back. The energy of zero-probability states is handled
by a large-gap limit. This protocol realizes the
free-energy difference, including the possibility of
recovering remaining stock free energy. It does not
construct a physical actuator that achieves it at finite
time. Paying for a new initial stock with no recovery
instead incurs the full initial preparation price.

For resets that act separately on X, e, m and the entire
fuel block, let $p_\tau^i$ be their final marginals and
$\pi_i$ their Gibbs marginals. The ideal marginal work
sum is
$$
 W_{\rm sep}=\sum_i\big[D(p_0^i\Vert\pi_i)-D(p_\tau^i\Vert\pi_i)\big].
$$
Energy is additive over these blocks, including all bonds
inside the gas block. The starting law is a product over
the blocks. Expanding relative free energy into marginal
terms and total correlation gives
$$
 W_{\rm sep}-W_{\rm joint}=
  I_{\rm tot}(p_\tau)=\sum_i H(p_\tau^i)-H(p_\tau).
$$
Fuel multiplicity cancels in this difference. Each
marginal reset can realize its own reversible infimum
without access to the correlations, and the joint
correlation work potential is lost. This accounts for
stored information and its erasure. It is not justified
to add $k_B T_b\log 2$ per receptor flip to this ledger.
Acquisition exchanges, final receptor free energy and
conditional correlations are distinct quantities.

Both comparisons restore the same state preparation
and stock per trial. Arbitrary joint reset may use all
blocks after the scored episode, including e. It
supplies no source port during the episode. Keeping a
trial's preset control label allows a reset conditioned
on that setting. Fixed barriers, kinetic calibration
and materials are shared in the declared access
comparison. Extra apparatus to implement a direct
cue changes the allowed coupling. The infima price
state preparation, storage, response and reset under
the model. They do not price creating the material
itself, finite-rate reset hardware, or an external
instrument used for verification.

## Registered example

At $L=8,K=8,J=1.2$ the cyclic persistence preparation
costs 24.500639. The charged stock accounts for 24.388699
and the shape condition for .111940. Passive operation
retains its fuel and costs .108592 for ideal joint
reset after $\tau=4$. It relaxes the shape distribution,
despite using no fuel. The best survival among the
settings that actually use an acquired record is
.752614. It uses 1.027347 net conversions for the
receptor and .314163 for assisted rebinding. Its
finite-episode heat is 3.972826, its microscopic
entropy increase is 3.570846 nats and its joint
reset costs 7.543672. The separate reset costs
7.704701, including .161029 nats of unrecovered
total correlation. The comparator with $a_0=a_1=1$
and sensor disabled reaches .754112 for 2.195984
joint or 2.448584 separate. The scores and total
works refer to unconditional trials.

The least average joint work for .75 survival in
the frozen whole menu is 1.010066. The least average
separate-reset work is 1.121588. The inputs are
calculated for each setting and then mixed as specified.
A full-episode minimum that also prices finite reset
time or operating different hardware remains open.
[All ledgers](RESULTS.json) · [Cost calculation](DERIVATIONS.md)
