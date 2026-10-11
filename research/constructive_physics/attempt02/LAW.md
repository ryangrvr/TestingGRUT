# CP-02 / D003: shared induced geometry

**One selected construction; author hypothesis; not established physics of nature.**
Freeze occurs after the analytic design disclosed in DESIGN_LOG.json, before
executable checks. No modified potential or repair of CP-01 is used.

## Domain and primitives

Use units hbar=c=1. A timelike four-dimensional hypersurface is embedded in flat
five-dimensional spacetime with metric eta_AB=diag(-1,1,1,1,1). Embedding fields
X^A(x), A=0,...,4, and one massive real matter scalar chi are supplied. Worldvolume
coordinates x^mu have mu=0,...,3. Ordinary local quantum field quantization,
preparations, the Born rule and tensor composition are supplied.

Postulate ambient ISO(4,1), worldvolume reparametrization invariance, normal
reflection, and that the matter field and the embedding motion share the same
physical induced metric

\[
\gamma_{\mu\nu}=\eta_{AB}\partial_\mu X^A\partial_\nu X^B.
\]

This is an active physical symmetry and coupling hypothesis, not passive
coordinate covariance alone. No evidence for an extra dimension is asserted.
In static gauge X^A=(x^mu,Y), choose canonical pi=sqrt(T_R)Y, where the measured
renormalized tension T_R>0 has mass dimension four:

\[
\gamma_{\mu\nu}=\eta_{\mu\nu}+T_R^{-1}\partial_\mu\pi\partial_\nu\pi.
\]

Require a nondegenerate Lorentzian induced metric, a regular graph, weak momenta
and curvature below an EFT cutoff Lambda, and a positive perturbative kinetic
spectrum. Initial data are smooth small perturbations of Y=chi=0; scattering has
asymptotically flat vacuum boundary conditions. Local classical solutions are
claimed only until graph degeneracy, large curvature or loss of hyperbolicity.
No global existence, UV completion or finite-dimensional exact quantum model
is claimed.

## Action and full parameter inventory

\[
S_0=-T\int d^4x\sqrt{-\gamma}
-\frac12\int d^4x\sqrt{-\gamma}
(\gamma^{\mu\nu}\partial_\mu\chi\partial_\nu\chi+\mu^2\chi^2).
\]

mu>0 is independently calibrated. Add **every** local counterterm allowed by
the stated symmetries, not just the ones convenient for the prediction:

\[
S_{\rm EFT}=S_0+\sum_a c_a\int d^4x\sqrt{-\gamma}\,
 {\cal I}_a(\gamma,K,\nabla_\gamma,\chi).
\]

Here n_A is a unit spacelike normal with n_A partial_mu X^A=0;
K_mu nu=n_A partial_mu partial_nu X^A, and nabla_gamma is the Levi-Civita
derivative of gamma. The index a runs over local polynomial scalar contractions
of K and its covariant derivatives and of chi and its derivatives, with even
total normal parity and chi parity. Redundant terms may be removed by integration
by parts and local field redefinitions. Pure zero-curvature, zero-matter density
is the volume term; matter mass, wavefunction and self-couplings are independently
matched. R=K^2-K_mu nu K^mu nu, with Riemann fixed by the Gauss equation. Examples
include sqrt(-gamma)R, R^2, R_mu nu R^mu nu, tr(K^4), chi^2 R and chi^4.
There is no artificial finite operator closure. Wilson coefficients at the
matching scale are independent finite data. Higher derivatives are treated
perturbatively; extra poles from resumming them lie outside the claimed EFT.

At any specified derivative order this is a finite EFT truncation. Beyond that
order an infinite coefficient inventory remains. T_R, mu, Lambda, matching
coefficients, initial conditions and detector operations are supplied.
No thermostat, state-dependent Hamiltonian or consciousness variable is used.

Static heavy probes of calibrated masses M_1,M_2 use the same induced geometry:
S_probe=-sum_i M_i int ds_gamma. Finite-size/curvature probe operators are also
allowed and separately matched. The controlled force calculation uses the
weak-source limit, r much larger than source sizes, inverse matter gaps and
inverse cutoff. Probe back-reaction is retained only to order M_1 M_2/T_R^2.
No new independent graviton is included. Claims are in the flat-ambient,
gravity-decoupling domain; actual finite-G experiments require conventional
gravitational effects and background bounds absent here.

## Quantum prescription and predeclared claims

Use geometrically invariant background-field regularization and dimensional
regularization d=4-2epsilon with renormalized canonical pi. Check the nonlinear
ambient-boost Ward identity, ordinary mixture consistency and leading
perturbative unitarity. Do not infer a complete anomaly-free UV theory from a
one-loop calculation. Heavy-matter determinants and allowed geometric
counterterms must be included in the protection argument.

The proposed restricted observables are:

1. The coefficient kappa of the leading soft, on-shell four-pi amplitude,
   M_4=kappa(s^2+t^2+u^2)+higher derivative terms, massless s+t+u=0.
2. C_7=-lim_{r->infinity} r^7 V_pi(r)/(M_1 M_2), defined for the scalar-induced
   connected interaction in the controlled flat vacuum, at leading weak-source
   order. This is not the full gravitational potential.

The analytic design proposes kappa=1/(4T_R) and
C_7=3/(128 pi^3 T_R^2), hence C_7=3 kappa^2/(8 pi^3).
These are hypotheses to derive/check, not installed independent equations.
The coefficient relation is an asymptotic EFT statement. Finite-energy and
finite-distance errors are not assumed bounded without matching information.

Check a pure-pi one-loop amplitude, including its UV counterterm and RG running,
and the geometric matter determinant. No counterterm's finite matching
coefficient may be set to zero to claim an exact finite-energy prediction.

## Hostile comparisons and stopping tests

* Full ordinary 4D Lorentz/shift/reflection EFT, with the same fields, states,
  masses and measured kappa, permits an independent matter disformal coefficient
  b. Use metric eta+(b/T_R)partial pi partial pi. b=0,1,2 are predeclared
  comparator assignments; b not equal to 1 abandons the extra physical 5D
  symmetry and must be identified as such.
* The conventional EFT sharing **all** CP-02 symmetries includes all geometric
  counterterms. Test whether it already implies the proposed leading relation.
  If it does, no new GRUT exclusion over that complete baseline is established.
* Compare zeta=0 and zeta=1 multiplying int sqrt(-gamma)R^2, with other matching
  data equal. Both share the physical symmetry. Derive their differing finite
  four-point amplitudes; do not call finite-energy rigidity protected.
* Test whether induced geometry and an allowed induced Einstein-Hilbert term
  can produce the observed Newtonian 1/r force and two spin-2 polarizations.
  If not, record the obstruction; do not introduce a new bulk graviton, change
  codimension or add a new law to rescue the construction.

Predeclared software parameters are T_R=1 in cutoff-normalized units,
angular cosines {-1,-1/2,0,1/2,1}, scattering energies {1/20,1/10},
source separation ratios {1,2,5} for rescaled integral controls, kinematics
(s,t,u)={(1,-1/2,-1/2),(2,-1/2,-3/2)}, scale factors {1/2,1,2},
and determinant gradients {(0,0,0,0),(1/5,1/10,0,0),(0,1/2,1/3,0)}.
These are software controls, not measured calibration or laboratory forecasts.

## Success boundary

Derive a quantum-stable **leading** cross-observable restriction, with explicit
premises and the full conventional comparison. A surviving conditional EFT is
not a fundamental GRUT discovery. Failure of the gravity limit or novelty is
recorded without rewriting this law. No Stage-3 score, L_0 compression claim,
experimental lock, canonical integration, CP-03 or Card 3 is authorized.
