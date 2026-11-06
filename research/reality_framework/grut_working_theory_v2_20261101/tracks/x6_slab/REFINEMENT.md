# Additional source diagnostics

The initial cutoff control changes final integrated bare work by
0.55%, while the full instantaneous bare work differs by 7.62%
relative to its peak. Both comparisons are archived. Independently
halve the timestep on the wide cylinder and vary its longitudinal
and transverse cutoffs separately.

The finite action places derivative subtractions in scalar
kinetic coefficients. To compare the fermion exchange with canonical
scalar stress, also move those explicitly known stress and work
terms into the fermion ledger. Check that the reallocation
preserves the combined budget and that its own fermion pressure
work agrees with the independent time integral of the matched
force. Apply the same numerical comparison to that work.
The original raw-work verdict is recorded alongside it.

## Initial-state and edge tests

The matched ledger reduces the wide-cutoff instantaneous-work
difference from 7.62% to 6.41%. A separate half-step on the wide
cylinder moves the matched work by 0.376%. Test the sources of
the remaining cutoff difference by extending the first-order
case to longitudinal integer 12 and transverse cutoff 16 at the
fine timestep, and by doubling the reference transverse resolution
again. Evaluate the exact decomposition of the spatial matching
coefficient into interior pairs and missing states at the
longitudinal edge.

As a distinct preparation control, add the second-order adiabatic
graph term using the initial scalar acceleration and metric rate.
Compare its reference and wide cylinders at the fine timestep.
Check that it improves the approximation to the exact stationary
rotating vacuum in a uniform flat control. Keep the original
first-order production protocol and verdict separate from this
preparation test.
