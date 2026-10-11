# CP-01 / D002: frozen reciprocal curvature-response construction

Freeze date: 2026-10-11 UTC. Author hypothesis, not validated physics.
This text precedes executable checks. Analytic selection work is disclosed
in ../DRAFT_LOG.json.

## 1. Physical domain and supplied foundation

Use c = hbar = 1, signature (-,+,+,+), and M_P^2 = 1/(8 pi G).
M_P > 0 and m > 0 are independent parameters. Spacetime is a smooth,
globally hyperbolic four-manifold with a Lorentzian metric g and positive
dimensionless field F. The connection is the metric Levi-Civita connection;
R is its Ricci scalar.

For a fully declared matter inventory take one real scalar chi with
mass mu_chi > 0 and self-coupling lambda_chi >= 0:

S_m = - integral sqrt(-g) [ (partial chi)^2/2
                          + mu_chi^2 chi^2/2 + lambda_chi chi^4/4! ].

The classical field equations apply to smooth constraint-compatible initial
data. The flat-vacuum sector uses asymptotically flat boundary conditions,
F -> 1 and chi -> 0. A timelike Dirichlet boundary, if used, includes the
usual generalized boundary term M_P^2 integral sqrt(|h|) F K, fixing F
and the induced metric there. No global regular solution is assumed for
arbitrary large data.

Matter inventory, geometry, initial data, preparations and quantum measurement
postulates are supplied. No observer/interface chain is claimed to emerge.

## 2. Complete proposed action

S_J = (M_P^2/2) integral sqrt(-g) [ F R - U(F) ] + S_m[g,chi],

U(F) = (3/2) m^2 F(F-1)^2,        F > 0.

There is no explicit Jordan-frame kinetic term for F. Its propagating kinetic
term comes from metric-scalar mixing, not from an undefined operator.

The action, absence of additional operators, universal metric coupling and
the function U are POSTULATED. They are not consequences of the old GRUT
substrate. Their physical motivation is a positive reciprocal susceptibility
with a stable flat vacuum and the lowest nonconstant reciprocal harmonic.
Reciprocity alone does NOT fix the function: removing higher harmonics is
a separate load-bearing model hypothesis.

## 3. Explicit equivalent variables and proposed restriction

For F > 0 set g_E = F g, a = sqrt(2/3),
phi = (M_P/a) log F, x = a phi/M_P, A = exp(-x/2).

S_E = integral sqrt(-g_E) [
       M_P^2 R_E/2 - (partial phi)^2/2 - V(phi)]
       + S_m[A^2 g_E,chi],

V(phi) = (m^2 M_P^2/a^2) [cosh(x)-1].

This equivalence and the resulting scalar coupling are to be DERIVED.
The vacuum potential has phi -> -phi symmetry; the matter coupling A does
not. This is not invariance under arbitrary physically performed local
transformations, nor a protection theorem for the complete quantum action.

Define lambda_n = d^n V/dphi^n at phi = 0, in the canonical variable.
The proposed tree-level restriction is:

lambda_2 = m^2, lambda_3 = 0,
lambda_4 M_P^2/m^2 = 2/3,
M_P d(log A)/dphi = -1/sqrt(6).

Only a derived operational response or amplitude makes these useful
predictions; a definition of lambda_n is not itself evidence.

## 4. Baseline and hostile comparison, fixed before execution

The comparison class is local scalar-tensor effective physics with the
same metric/matter inventory, canonical Einstein scalar, positive kinetic
terms, flat vacuum, scalar mass, universal A coupling and operational
preparations. Potential coefficients beyond the mass are independent.
It has no stipulated equilibrium state: KMS/FDT/Onsager assumptions are
not invoked for vacuum scattering or coherent classical evolution.

An explicit comparator uses the same A and:

V_epsilon = (m^2 M_P^2/a^2) {
 cosh(x)-1 + epsilon [cosh(2x)-4 cosh(x)+3] }, epsilon >= 0.

It preserves vacuum energy, tadpole, mass, parity and the leading force
coupling. It must be tested for stability and for violating the proposed
quartic restriction without changing those measured quantities.

## 5. Predeclared failure tests and correspondence obligations

- Derive the Jordan/Einstein equations and total stress conservation.
- Check positivity, the vacuum spectrum and the local Cauchy domain.
- Derive the long-distance Einstein limit and the weak-field Yukawa force.
- Exhibit a nonlinear response or scattering observable separating the
  candidate from V_epsilon using independently fixed m and M_P.
- Test the global Legendre branch; do not claim a global single-valued
  metric f(R) action without an invertibility certificate.
- Test operational mixture consistency and local causality under the supplied
  ordinary quantum interpretation. Do not use H[rho].
- Compute the scalar one-loop local potential divergence. Test whether its
  field dependence is absorbed by the declared family, allowing vacuum
  energy, parameter and field-normalization renormalization.
- Disclose omitted graviton/matter loops, finite matching constants and
  measurement assumptions. A scalar-sector closure failure is not a
  theorem excluding every possible ultraviolet completion.

No exact all-orders prediction is preclaimed. If quantization requires an
independent operator affecting the same observable, the finite-parameter
fundamental-closure claim fails. An effective model with an unspecified
matching coefficient is not a rescued parameter-free theory.

## 6. Information and evidence boundary

Physical parameters: M_P, m, mu_chi and lambda_chi, plus initial and boundary
data and any experimental control history. The novel supplied content is
the exact curvature potential/minimal-harmonic closure; universal metric
coupling and quantum/relativistic structure are existing supplied physics.

No frozen L_0 bit score, joint compression credit or Stage-3 passing verdict
is asserted. The Appendix-B evaluation kit remains independently unchecked.
The complete function contains more physical specificity than the single
number m; free parameters alone are not its information price.

This attempt may produce a conditional classical restriction. It establishes
a fundamental theory only if the stronger quantum, physical and empirical
obligations are actually met.
