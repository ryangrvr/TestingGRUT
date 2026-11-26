# GRUT-I2 microscopic model

25 November 2026 · Parent draft [PR #21](https://github.com/ryangrvr/TestingGRUT/pull/21).
All energies are in units of the fixed solvent temperature k_B T_b.
Use τ for episode duration. The state describes a small overdamped
ring of adhesive particles and one local catalytic patch. It is a
finite stochastic microscopic model of particle occupations and
chemical states. A fast solvent and local molecular mixing are
the physical assumptions behind the activated jump description.
No microscopic Hamiltonian for that solvent is inferred from this
coarse graining.

## Physical state

There are L labelled sites on a ring with L/2 identical particles.
Only single occupancy is allowed. The cohesive bond count is
$B(x)=\sum_i x_i x_{i+1}$ and the particle energy is $-J B(x)$.
Neighbors can exchange a particle and vacancy anywhere. The
patch on link (0,1) has a thermally switching, energy-degenerate
conformation e. Its independent forward and reverse flip rates
are κ. It changes the activation barrier of hopping on its
own link. It does not add a field to the energy.

At the same patch is an energy-degenerate, two-state receptor m.
It flips at symmetric thermal rate ν. It can sample any one
of K locally mixed fuel molecules in that pocket. A molecule
has two states F and W with energies A and 0. Its conversion
can supply work for a local rebinding reaction or change the
receptor. Both reverse reactions are included. At all times
the same equation applies. There is no externally timed
actuator, measurement feedback, target pattern or orientation flag.
A readout of the receptor is a chemical dependence of a rate.

Let k be the number of F molecules. The full microscopic state
has the K individual fuel bits and total dimensionless energy
$E=-JB+A k$. For the exchangeable fuel preparation with k=K,
the count is an exact Markov lumping. Its phase volume is
$g(k)=\binom Kk$. The gas, source, receptor and fuel identities
are retained as physical degrees of freedom. The equilibrium
measure of the count process is

$$
\pi(x,e,m,k)=Z^{-1} e^{JB(x)-A k}\binom Kk .
$$

## Elementary rates

When a site pair on link i has different occupations let the
exchanged configuration be y and $d=B(y)-B(x)\in\{-1,0,1\}$.
The background hopping rate is $r_{\rm bg}e^{Jd/2}$.
At the patch replace $r_{\rm bg}$ by $r_e$. A catalytic
rebind at the patch adds a separate channel, only for
$\lvert d\rvert=1$, ending at k−d. With coefficient a its rate is

$$
w_{\rm a}(x,k\to y,k-d)=a e^{(J+A)d/2}
\begin{cases}
 k/K,&d=1,\\
 (K-k)/K,&d=-1.
\end{cases}
$$

At K=0 both fuel channels are absent. A concerted exchange and
conversion at a physical contact supplies these rate ratios.
The geometric mean prefactor is an activated-barrier model.
Sampling a particular available fuel molecule has coefficient
a/K, which gives the count dependence without an undepleted
chemical potential. There is no fuel action on neutral hops.

The receptor's additional fuel-coupled flip m→1−m has
$d_s=+1$ if the new bit equals e and $-1$ otherwise.
Its rate is $s e^{A d_s/2}$ times the same available
fuel or waste fraction, with k→k−d_s. The two conformations
select the direction of the reversible chemical coupling.
They do not change receptor energy. This idealized local
transduction assumes physical kinetic recognition of e.
E, m and k are not directly available to other hops.

Holding e and m fixed, the repair coefficient a is either a
constant, $a_m$ from the local receptor, or $a_e$ in the
explicit direct-access control. It is the same coefficient
for a transition and its reverse. All other couplings,
energies and temperature are the same across these arms.
Every rate obeys microscopic thermal detailed balance
against the stated energy including fuel. Thus the
finite model can persist away from equilibrium for a time
while spending prepared fuel. It cannot sustain a steady
stationary current with finite equilibrium resources.

## Preparation and persistence

For persistence, take the canonical attractive gas conditioned
on $C(x)\le\lfloor(L/2)/2\rfloor$, where $C=L/2-B$ is the number
of nonempty connected clusters. This is a shape-class preparation
with no preference among translations or exact configurations.
Both e and m are initially fair and independent. All K fuel
molecules initially occupy F. Use the same preparation for
all controls at the same L,J,K. The shape preparation and
fuel charge are priced.

A first-exit survival is the probability of staying in this
shape class throughout [0,τ]. Endpoint membership can include
reentry and is reported separately. Other declared observables
are bonds per particle and longest connected component per
particle. For formation use a uniformly random fixed-occupancy
gas rather than a prepared aggregate. Passive dynamics can
then select bonds by cohesion. Cooling, holding or repairing
an exact microstate is a different task.

## Resources

All arms at fixed capacity start with the same finite materials
and high-energy stock. The bath temperature, patch and background
passive rates, principal clock τ and catalyst rate menu are
matched. The no-cue arms can spend the entire stock on rebinding.
Receptor-driven arms must share it with acquisition. The direct
cue comparator opens an additional local coupling to e and is
identified as such. Program settings are fixed reaction barriers,
and can be measured independently rather than inferred from score.

The physical ledger counts free energy of preparation, net
chemical conversions for both channels, bath heat, entropy and
correlations. Reset costs are thermodynamic infima with known
end distributions and ideal quasistatic operations. Compare an
ideal joint reset and separate resets of gas, e, receptor and
fuel. These operations act after the scored duration. Restricted
or finite-rate apparatus can cost more. Apparatus fabrication
and external diagnostic measurements are outside the operational
episode and cannot be estimated by the model. The rates and
material resources specify the episode without those measurements.
