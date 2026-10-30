# Stabilized M1 holder and interacting transport protocol

Parent 94f660d0e04e7b40b668ffb316793a575f3e3813. The attached Argon proposal
is a set of hypotheses, not an independent numerical reproduction.

1. Use the cubic curvature operator as the primary construction. Check its
canonical dimension, the full quartic form including kappa, and an explicit
global lower bound at fixed finite R. Distinguish scalar boundedness from
positivity, positive effective Planck mass, quantum consistency and a UV
completion. Audit both suggested higher-degree completions analytically.
Choose mu_eff/H = mu/H + (10 i - mu/H) R/(12 H_i^2), explicitly encoding
the requested initial relative phase. Keep the old light-sector parameters.
Solve the de Sitter stationary state and evolve through H_i t=800. Calculate
a fixed-dimensionless-family abundance refit, specifying which heavy masses
co-scale. Record gauge and energy constraints, timestep/tolerance controls,
window sensitivity, and nonminimal gravitational correction along the run.

2. Inspect the actual inherited 12x12 engine before reuse. It omits metric,
inflaton and radiation perturbations and does not initialize Bunch-Davies
modes. Extend the scalar/longitudinal calculation where controlled; propagate
a declared de Sitter vacuum covariance through reheating and compare the
longwave physical subspace with finite-k, finite-e controls and homogeneous
finite differences. Do not label spectator charge power as physical CDI.
The inflation duration and the source of primordial curvature are additional
inputs. Do not claim the scale-invariant uncorrelated observational bound
without its transfer and statistical assumptions.

3. Evolve the inherited 64/96 source measures under nonrelativistic identical
XX elastic scattering in comoving momentum so expansion does not deposit
particles into a new fixed physical-energy floor. Test number/energy,
positivity, angular/grid convergence, and thermal stationarity. Quantify
relativistic and high-energy-tail errors against the inherited kernel.
Use electron-positron entropy transfer for the temperature clock; integrate
cumulative scattering rather than infer thermalization from Gamma/H.
Measure the final shape and fit a declared smooth moment-matched ansatz,
retaining fit error. Attempt a controlled perturbation/linear-power
comparison; expose collision-model, primordial/source and initial-condition
assumptions, k-range, solver convergence and failures. A homogeneous
Boltzmann solution cannot by itself determine acoustic perturbations.

4. Keep the frozen M1 and previous numerical packets unchanged. Save code,
data, derivations, source pins, reproduction commands and limitations.
Publish an additive branch after verification; no independent review claimed.

## Implementation addendum (recorded after exploratory calculations)

The density calculation was implemented as an internal conserving BGK test
cosmology: supplied formation matching at 2 MeV, a fixed comoving Maxwell or
generalized-Maxwell shape, all other matter cold, one barotropic radiation
fluid, Phi=Psi, and radiation-oscillation averaging beyond a tested k/aH
threshold. These were implementation choices, not earlier preregistered
physical predictions. Controls vary q nodes, multipole truncation, tolerance,
radiation-averaging threshold, temperature seed and collision multiplier.
The exact linearized collision kernel, pre-checkpoint perturbation sources,
and the full photon/neutrino/baryon hierarchy remain missing. Negative CLASS
solver outcomes are retained. See REVIEW_RESPONSE.md for the replay-threshold
revision and its separate tighter control.
