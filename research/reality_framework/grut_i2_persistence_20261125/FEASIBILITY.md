# Literature and analytical feasibility gate

25 November 2026. This gate precedes production. Accept I1 verdict B
and test an independent lattice-gas mechanism for natural aggregate
formation and finite thermal persistence. The question is whether
a chemically acquired local record improves the material's survival
at matched physical cost. The task can yield a negative result.

## Physical motivation

Nearest-neighbor particle exchange is the particle-conserving
lattice-gas dynamics associated with Kawasaki. Cohesive bonds have
a physical energy and position-independent attraction. Chemical
fuel can transiently maintain assemblies in laboratory systems.
The I2 patch models a local analogue of fuel-aided association.
It retains reversible waste branches and a finite fuel inventory.
Chemical feedback by a receptor is physically conventional.
These sources motivate the interaction classes, not the chosen
numerical coefficients or a claim that the entire I2 device
has been realized. [Kawasaki](https://doi.org/10.1103/PhysRev.145.224) ·
[Boekhoven et al.](https://pubmed.ncbi.nlm.nih.gov/26339025/) ·
[Barato et al.](https://doi.org/10.1103/PhysRevE.87.042104)

The local patch avoids requiring a fast global sensing or control
broadcast. It samples local chemical conformers with simple
activated rates. The solvent supplies thermal noise without
encoding a shape, and one patch has no extensive fuel throughput
at fixed episode time. Large aggregates under driven bulk
chemistry would require a different registered model.

## Analytical checks before choosing the test

1. A hop changes at most one bond. Comparing two reversed
   elementary transitions including the identity of the fuel
   molecule gives the ordinary Boltzmann ratio. The finite count
   reduction adds the multiplicity ratio and retains reversibility.
2. At fixed particle number, the multiplicity with c clusters is
   $L\binom{L/2-1}{c-1}^2/c$. Attraction favors low cluster count
   without selecting an arrangement. The half-filled bulk limit has
   typical cluster fraction $1/(1+e^{J/2})$ with finite J. It is
   smooth. A fixed density threshold can still become sharp.
3. A top score defined as exactly one cluster has little room
   for internal rebinding before first exit. The registered class
   at the primary size includes precursors with two clusters,
   allowing a difference between first-exit persistence and return.
4. With one temperature and finite fuel, the expanded chain
   is detailed balanced. Its loss of relative free energy is
   entropy production. Recharging to a specified initial
   distribution has a reversible work infimum. Allowing separated
   resets can waste correlations and raises that cost. These
   results enable an exact cost frontier in a finite policy menu.
5. The state at L=8 and K=8 fits a sparse generator, and its
   persistence can be solved by restricting that generator with
   its full original escape rates. An independently enumerated
   individual-fuel calculation is small at L=6,K=2. No trajectory
   Monte Carlo is needed for acceptance.

## Novelty disposition before computation

Causal local detailed balance, fuel accounting, information
dependent rates and free-energy/reset identities are established.
This model predicts a conditional rate ratio and a cost frontier
from those laws. It supplies no additional thermodynamic affinity
associated with a record. A positive gain or failure of the record
can answer the scoped physical question, and neither outcome is
a GRUT-specific deviation. Conventional equivalence closes
that novelty claim at this gate. The forthcoming computation
tests the effective description and retains any negative result.
[Rao and Esposito](https://arxiv.org/abs/1805.12077) ·
[Sagawa and Ueda](https://doi.org/10.1103/PhysRevLett.102.250602) ·
[Lan et al.](https://doi.org/10.1038/nphys2276)
