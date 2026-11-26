# GRUT-I2 · Natural organization and thermodynamic persistence

25 November 2026. A finite, attractive lattice gas forms more bonds
from an unorganized start through thermal hopping. A local chemical
patch raises the bond fraction further while spending prepared fuel.
An autonomous receptor carries information about the patch's past
state, but using it does not improve the registered persistence–work
frontier once sensing and reset are priced.

The model, menus, diagnostic material and gates were frozen at
4dcc33c4720969b2898d7a0320a4a9ec3de9ea39 before production.
No named arrangement or future event enters the rates.
The solvent, bond energy and catalytic patch are specified physical
inputs. [Model](MODEL.md) · [Protocol](PROTOCOL.md) · [Freeze](FREEZE.json)

## Cohesion, drive and record

All primary runs use a ring of eight sites with four particles.
Formation starts with uniformly random occupancy at that particle
number. The dimensionless episode lasts four time units.

| Formation arm | Initial bonds per particle | Final bonds per particle | Final probability of declared shape class |
|---|---:|---:|---:|
| Flat passive diffusion | .428571 | .428571 | .628571 |
| Thermal cohesion | .428571 | .518685 | .844046 |
| Cohesion and constant fuel | .428571 | .536343 | .883409 |
| Cohesion and registered receptor gate | .428571 | .528048 | .864722 |

Cohesion adds .090114 bonds per particle, clearing the frozen .02
gain. The constant active patch adds .017657 over the cohesive
passive arm, clearing .01. The phase space is translation invariant
in energy. The patch acts only on its local link and gives no
particular layout extra binding energy.

For persistence the material is prepared in the shape class,
defined by at most two occupied clusters at this size. Thermal
cohesion raises first-exit survival from .493785 to .746719.
The J=0 and J=1.2 conditional Gibbs starts weight shapes
differently within the class. That comparison tests preparation
and passive dynamics together. The formation test starts with
the same uniform population in all arms.

| Persistence arm | First-exit survival | Endpoint in shape class | Joint reset work |
|---|---:|---:|---:|
| Cohesive passive | .746719 | .907562 | .108592 |
| Constant fuel at a=.5 | .753174 | .932807 | 1.883085 |
| Strongest sensor-off setting | .754112 | .937022 | 2.195984 |
| Best survival using acquired information | .752614 | 0.934105 | 7.543672 |

Endpoint membership can count repairs after exit. The
work comparison uses the full ensemble in every arm.

The receptor gate stores real physical information. The best
record-using setting has .085320 nats about the terminal patch
state and .032450 nats about its state at time 2 after conditioning
on the terminal state. It spends 1.027347 net fuel conversions on
sensing and .314163 on the gas. This is an acquired history signal
rather than merely an initial correlation. The top sensor-off
setting uses a constant coefficient with no history dependence.
The joint and separately reset costs are derived and recorded.
[Resource ledger](THERMODYNAMICS.md) · [All rows](RESULTS.json)

## Minimum cost within the frozen menu

The initial charged stock and physical gas preparation are identical
at fixed J and K. The joint ideal reset credits unused free energy
and retained correlations. The separate ideal reset loses the
correlation credit. Both restore the full ensemble. Work is measured
in units of $k_B T_b$.

| Allowed settings at first-exit reliability .75 | Minimum pure-setting joint | Minimum mixture joint | Minimum mixture separate |
|---|---:|---:|---:|
| Passive | Infeasible | Infeasible | Infeasible |
| Full receptor-capable menu | 1.084608 | 1.010066 | 1.121588 |
| Entire no-feedback control | 1.084608 | 1.010066 | 1.121588 |
| Direct local cue, with extra access to e | .850788 | .674274 | .757763 |

The table's mixture minimum is exact for the declared finite
menu in the ideal thermodynamic model, evaluated numerically.
None of the menus attains .80 at the primary duration.
Extra access through an instantaneous local cue can lower
the cost. That does not imply a benefit from the receptor.

Every common feasible registered budget gives a tie between the
full record-capable and no-feedback family. The frozen record
benefit gate therefore fails. A derived dominance check on the
saved rows gives a stronger finite-menu statement. The constant
setting with a=1 beats all 32 settings with sensing, including
unused-sensor controls, by at least .001498 in survival and
saves at least 4.035881 joint or 3.875468 separate reset work.
It can replace any sensing component of a mixture. This
establishes equality of the two cost fronts across this menu,
including work caps between the registered values. The check
does not change the registered gate. It extrapolates neither to
continuous rate choices nor to untested sensing architectures.
[Derivation](DERIVATIONS.md) · [Numerical frontier](RESULTS.json)

## Collective and finite-resource limits

The exact half-filled equilibrium multiplicity gives a typical
cluster fraction tending smoothly to .354344 at J=1.2.
At 512 sites the declared coarse shape class has probability
greater than .99999999999. The probability of a single aggregate
is \(3.11\times10^{-94}\). These are correlated equilibrium
clusters, and the exact counting proof also bounds largest runs
on a logarithmic scale. The coarse predicate does not establish
a phase transition or a thermodynamic-limit single aggregate.
The one-patch driven size runs cannot make a bulk claim.

At two starting fuel molecules the representative sensing gate
finishes with zero F in .492517 of trials. The constant-active
comparator does so in .049256. Zero stock removes every assisted
and acquisition transition and reproduces passive survival.
Reverse fuel moves are still possible at nonzero stock capacity.
The controls therefore distinguish finite preparation-driven
activity from indefinite chemical drive. [Negative results](NEGATIVE_RESULTS.md)

## Reproduction and reality test

All 64 primary policies, formation, stock and size controls, and
the large canonical grid are recorded. The separate labelled-fuel
solver uses positive uniformization rather than the count generator's
exponential action. It agrees on 89 scalars over six prescribed
cases within \(9.77\times10^{-15}\). The largest detailed-balance
error is \(7.11\times10^{-15}\) and the largest free-energy
identity error is \(8.23\times10^{-13}\), below the frozen tolerance.
Numerical agreement concerns the specified model. A solvent and
material-specific microscopic realization has yet to be calibrated.

The single-hop forward/reverse rate ratio changes by a factor of
ten as available F decreases from eight to four of eight molecules
in matched comparisons. Its coefficient is independent of whether
the nonzero gate uses a receptor. That is a concrete conventional
laboratory prediction. An extra measured memory-dependent affinity
under matched microscopic conditions would reject this model.
No such deviation is predicted or observed here.
[Novelty audit](NOVELTY_AUDIT.md) · [Reproduction](REPRODUCTION.md) ·
[External review packet](REVIEW_PACKET.md)
