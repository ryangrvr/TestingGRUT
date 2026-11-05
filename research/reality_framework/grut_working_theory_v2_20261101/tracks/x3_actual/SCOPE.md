# Actual-history continuation: scope committed before new results

Base: 98ccd20d136bbb42b65c3fa85c98bb003e22e8cc. Existing numerical methods and frozen parameters remain read-only. This is a new, explicitly conditional study, not a revision of the original two-node M2-X registration.

1. Use the actually evolved H = 5.5e10 GeV scalar history in x2_strip/results/lower_H_lattice.json. The old hLF55_mid history reused the H84 lattice: it is a proxy control, not the new evolved point.
2. Reconstruct the kinetic state with a documented moment/spectrum rule, preserve energy at the splice, evolve to a converged dilution integral, and compare the evolved history with the old proxy. Vary momentum resolution and handover assumptions. Record all numerical failures.
3. Test the fermion closure with controlled spatial/time resolution and explicit quantum-mode checks. A 1+1D ray with a 3D weight is not a full-3D Dirac solution. Distinguish tests of the reduced algorithm from certification of the cosmological calculation; do not label an unresolved closure as accepted.
4. For nuggets, carry separate persistent-core, free-particle, antiparticle, coherent-field and neutral-reservoir ledgers. Enforce charge and energy accounting across conversion, reverse conversion and decay. Do not infer a population from a stationary root or count the same scalar and its decay daughters simultaneously. A population output requires supplied/evolved formation and transport data, with their limitations stated.
5. The original blind result remains 0/135 within its registered band [1.03e-6,4.12e-6]. New roots may be reported only with explicitly simulated scale support. Extrapolated roots remain rejected. No new cosmological point is accepted unless spectrum, closure and numerical controls support it.
6. Numerical controls: energy continuity at the kinetic splice <= 1e-6 relative; assess timestep and momentum-resolution changes directly and report >5% changes as unresolved for affected observables. Preserve the original scalar budget error (~1.33%) rather than renormalizing it away. For quantum tests report unitarity, vacuum subtraction, spatial and temporal errors separately.
7. Outcomes: conditional pass, numerical/closure failure, or unresolved. A failure closes only the implemented instance/approximation in its tested domain, never the field content or GRUT as a whole.

This scope supplies no new state-selection principle, vacuum-density derivation, or consciousness/intentionality claim. Report actual executions and provenance, not promised calculations.
