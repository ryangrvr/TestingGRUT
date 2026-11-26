# Conventional comparisons and experimental discriminator

The [preproduction feasibility gate](FEASIBILITY.md) closes the
claim that a memory-dependent barrier adds new physics to this
model. The production tests determine which physical regimes
exist under a finite resource constraint. They do not treat
recovery on an engineered device as evidence of natural formation.

| Subject | Primary comparator and outcome |
|---|---|
| Passive cohesion and collective motion | Conserved neighbor exchanges in an Ising/lattice-gas model are the Kawasaki class. The exact run statistics here explain the equilibrium trend, sharp shape readout and absence of a macroscopic aggregate in the finite-J limit. [Kawasaki](https://doi.org/10.1103/PhysRev.145.224) |
| Driven aggregation and active matter | Chemical consumption can drive transient assembly. Nonequilibrium interacting particles can also cluster through activity-dependent mobility rather than attractive bonds. I2 uses a local chemical reaction and explicit attraction. It supplies no motility-induced phase-separation calculation or universality claim. [Boekhoven et al.](https://pubmed.ncbi.nlm.nih.gov/26339025/) · [Tailleur and Cates](https://arxiv.org/abs/0803.1069) |
| Autonomous sensing and conventional feedback control | A switching environment, receptor state and chemical bias are standard autonomous sensory-network ingredients. Optimizing admissible policy mixtures and comparing an instantaneous cue are ordinary restricted stochastic control. The receptor information statistics are operational, and require a work comparison to establish advantage. [Barato, Hartich and Seifert](https://doi.org/10.1103/PhysRevE.87.042104) · [Lan et al.](https://doi.org/10.1038/nphys2276) |
| Full work, memory and depletion | The microscopic reverse transitions, chemical counts, stochastic entropy production and reset from relative free energy obey ordinary stochastic thermodynamics. A finite reservoir changes activity and ratios as it depletes. The marginal versus joint work difference accounts for erased correlation, with no universal fixed acquisition cost per bit. [Rao and Esposito](https://arxiv.org/abs/1805.12077) · [Sagawa and Ueda](https://doi.org/10.1103/PhysRevLett.102.250602) |

## Measurable relationship

Observe one allowed bond-increasing hop at the patch. Select
identical local x, y, e and m and hold fixed the affinity
and state-specific activity convention. Distinguish the
assisted channel by one F→W conversion. For K=8 define
$$
 L(k,m)=\log \frac{w_{\rm a}(x,k\to y,k-1\mid e,m)}
                      {w_{\rm a}(y,k-1\to x,k\mid e,m)} .
$$
The conventional result is
$$
 L(k,m)=J+A+\log[k/(K-k+1)].
$$
Only compare nonzero gate settings. At J=1.2, A=3
the primal rate ratios from k=8 and k=4 are about
533.491 and 53.349. Their log difference is
$\log 10=2.302585$, despite the same single-molecule
gap and contact energy. The expected difference in log
affinity across receptor labels at the same k and e
is zero. Both contrasts follow ordinary finite chemical
availability and microreversibility. A model that treats
the eight-cell stock as an undepleted drive misses the
fill contrast.

One possible implementation would observe adhesive particles
or small supramolecular units in an annular trap with a local
reactive patch. This is a proposed test geometry. Existing
fuel-driven assemblies establish the interaction motivation,
not that these particular receptor kinetics or state-resolved
measurements have been achieved. Calibrate the bond energy
and molecular gap independently of persistence outcomes.
At fixed preparations of x and y, measure conditional dwell
times and separable forward and reverse channel counts.
Prepare the corresponding reverse stock k−1. Use a waiting-time
likelihood with right censoring for the rates, with distinct
data for channel misclassification, temperature, patch
configuration and drift. Long measurements cannot silently
replace the finite stock by a constant reservoir.

The frozen target resolution for the fill contrast is .25
in log-rate units. At a fixed k the m contrast can be
estimated by matching both nonzero gates without fitting
J or A to the difference. A confirmed confidence interval
excluding the conventional prediction, after resolving
hidden reactions, reverse-channel stoichiometry and
calibration, would falsify this kinetic model. An apparent
violation from an unobserved path, wrong coarse states
or an external work source would instead change the
matched conventional model. No I2 mechanism predicts a
nonzero residual at fixed microscopic conditions.
This is a quantitative null discriminator, with no
GRUT-specific departure to claim.
