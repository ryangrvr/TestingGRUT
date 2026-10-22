# After-result strengthening: an exact finite-charge sum rule

The first prescribed run finished with 45 internal conditions and no failures.
It found the planned first-Born mass/readout ratio equal to one at all three
main charges: 1.00000000000218, 1.00000000000167 and 1.00000000000119.
The protocol had only proposed its joint large-charge/weak-probe limit.

Interpreting that unexpected pattern led to a stronger analytic result:
for every localized stationary fixed-charge solution of this *classical
flat-space action*, not just a thin-wall limit,

E = lambda_phi*v^2 * integral(v^2-phi^2) d^3x,

and therefore a_D^Born = -h_D^2 E/(4*pi*lambda_phi*v^2).
Two derivations are given in MATHEMATICS.md: the fixed-charge envelope/scaling
identity, and the integrated field equations plus virial identity. The
methods are standard; literature novelty has not been established.

This amendment adds numerical checks of that identity and its independent
integrated equations. It also adds a controlled finite-coupling interval
for the exact zero-energy scattering length when I=integral r|V| dr<1.
No physical coupling, charge target, momentum range or success criterion
is changed. The strengthening was motivated by already-seen results and
must not be described as a held-out or preregistered discovery.

The word "Born" matters: the exact finite-h_D scattering lengths differ
from the leading formula by 4.14%, 7.21% and 11.97% in the first run. The
sum rule does not remove those effects. Quantum corrections, new matched
mass scales, gravity beyond the flat saddle and extra operators can also
modify the tree identity; no all-orders protection is asserted.

Numerical and scope follow-up: double the packet momentum quadrature from
61 to 121 points to assess the integer readout counts. Explicitly identify
the short-range flat-space scattering length, rather than attaching it to
the full gravitational Coulomb problem. Estimate the finite-momentum
gravitational hierarchy without claiming a computed metric phase correction.
These are verification/scope additions; physical inputs remain unchanged.
