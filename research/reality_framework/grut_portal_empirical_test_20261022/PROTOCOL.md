# Portal empirical test: declared calculations and assumptions

22 October 2026 (America/Chicago). Parent fa3bc97b52f44bdd625da8b40fdc7015837517a7.
Owner instruction: "run it ainstein" following the proposed detector and
portal-history test. This is exploratory model construction after known
data, with no blind prediction or empirical discovery claim.

Use the parent's proposed subtracted Higgs portal, the existing classical
Q=19972,39944,79988 profiles, g_X=.1, lambda_phi=.01, lambda_X=.001,
h_D=.01, lambda_D=.001, m_D0=g_X*v. The scales .5,.75,1,1.25 GeV
are explicitly hypothetical variants. They put the radial mediator below
the two-muon threshold, avoiding uncontrolled GeV hadronic decay widths.
These choices preserve neither the original scale nor its zero portal;
they are the already-proposed extension being assessed.

Detector calculation:

- Use the official XENON1T approximate likelihood release v1.0.0,
  commit 8254684105193d4a84cb8b5752f5e450d5ab70b3. Fold true spectra through
  its migration matrix and use its recommended non-asymptotic threshold.
- Reproduce the official supplied 50 GeV point-WIMP example before using
  the new profiles. Independently implement the likelihood combination and
  compare selected results with the unmodified public implementation.
- Include natural xenon isotopes and compare with the parent's mono-isotope
  calculation. Check energy/form-factor integration and profile resolution.
- Recast all 1000 supplied XENONnT 20 tonne-year background-only toy
  likelihoods if feasible. This is the named historical detector projection,
  not an updated 2026 performance forecast or observed data.
- Obtain the official LZ 2024 efficiency and SI limits as an additional
  dated screen. A count-rescaling estimate must be labelled approximate;
  it is not the unavailable full LZ likelihood for a changed spectrum.
- Report limits on kappa at supplied halo fraction 1, and the dependence
  on that fraction. No global allowed region is inferred from one release.

History calculation at detector-relevant couplings:

- Derive the h-to-ss coupling including first-order mixing, s-to-electron
  decay and low-energy visible production through the same portal.
- Evaluate a hot, initially empty hidden sector with a broken-phase SM
  bath starting at 80 GeV. Follow conserved energy transfer from Higgs
  decays, a collisionless second-moment diagnostic and a rapid-internal-
  equilibrium estimate. They are distinct approximations. Test phase
  restoration in each; do not replace either with a complete Boltzmann
  solution of all species or a domain-wall simulation.
- Evaluate an explicitly cold history whose maximum SM bath temperature
  after the supplied preparation is 20 MeV. Include electron, muon and
  leading low-energy hadronic/photon production where relevant, and an
  energy-budget bound allowing all injected energy to survive as matter.
  Check the effect of the scalar's long lifetime using its produced energy,
  not a lifetime bound that assumes a thermal abundance.
- The cold +v phase and prepared charged-object population are initial
  data. Formation, baryogenesis and the earlier reheating action are not
  derived. A bound on maximum temperature is stronger than specifying a
  low final reheating temperature while allowing an earlier hot bath.
- Radiation-dominated cosmology is conditional on a supplied total vacuum
  matching compatible with the observed late universe. The original zero
  counterterm does not achieve this. This input is explicitly calibrated;
  the calculation does not claim to select the vacuum density.

Report an overlap only with the datasets, projection and history actually
computed, propagating approximation limits. Counterexamples can exclude a
history without excluding all portal completions. Conversely, passing a
low-production history does not prove formation or full cosmological
viability. One author; no independent replication or subagents.
