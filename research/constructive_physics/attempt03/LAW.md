# CP-03 / D004: scalar determinant mechanism

**Selected once, after analytic design and literature inspection; before executable
scientific checks. Published mechanism, not an originality claim.**

## 1. Objects and action

In natural units, take a smooth oriented four-dimensional manifold M, N real
scalar fields phi^a (N >= 4), and positive field-space metric delta_ab. Set

    Q_mu_nu = sum_a (partial_mu phi^a)(partial_nu phi^a),
    V(phi) = V0 + (m^2/2) sum_a (phi^a)^2,
    S[phi] = integral_M d^4x sqrt(abs(det Q)) / V(phi).

V0 > 0 and m > 0 are supplied. The real classical differentiable branch has
rank Q = 4 and V != 0. Boundary variations have compact support; no assertion
of global existence is made. Degenerate configurations are boundary points,
not silently excluded solutions of a proven complete quantum theory.

No spacetime metric, clock fields, graviton or Einstein-Hilbert term is in this
action. Dimension, differentiability, scalar inventory, positive field-space
metric, the determinant interaction and potential ARE inputs. The action is the
clockless variant of Carone, arXiv:1911.13166, Sec. IV, in D=4. It is not CP-02's
single transverse scalar embedding action.

## 2. Auxiliary description and signature test

For the classical algebraic comparison, introduce a nondynamical symmetric
invertible g and

    S_aux = integral sqrt(abs(det g)) [ (1/2) g^mu_nu Q_mu_nu - V ].

Its metric equation is

    T_full_mu_nu = Q_mu_nu - g_mu_nu [(1/2) g^rho_sigma Q_rho_sigma - V] = 0.

In D>2 the solution is g_mu_nu = (D/2-1) Q_mu_nu / V. Substitution yields
the determinant action with factor [(D/2-1)/V]^(D/2-1); in D=4 it gives the
action above. This is a classical elimination, not a proof that constrained
and unconstrained auxiliary-metric quantum measures coincide.

The test must determine whether a real Lorentzian branch exists. An analytic
continuation is an additional hypothesis, not a derived physical time.

## 3. Quantum calculation being tested

Use ordinary perturbative scalar quantization, a large-N expansion, and the
covariant dimensional prescription D=4-epsilon. Finite epsilon is a formal UV
proxy, not a demonstrated physical fractional spacetime dimension or a unitary
UV completion. A complete nonperturbative measure, reflection positivity and
Lorentzian reconstruction are obligations, not supplied proofs.

The published organization expands by adding and subtracting the same tensor:

    Q = [V0/(D/2-1)] eta + [Q - c1 eta],
    c1 = V0/(D/2-1),

and uses c2 to normal order the scalar mass operator. Eta in this organization
is an expansion/background assumption; it is not credited as an emergent
Lorentzian metric. Normal-ordering removes leading tadpoles. The consequent
leading interaction is to be independently checked:

    L_int = -[T_mu_nu Pi^mu_nu|rho_sigma T_rho_sigma]/(4 V0),
    T_mu_nu = sum_a [partial_mu phi^a partial_nu phi^a
               - eta_mu_nu ((partial phi^a)^2/2 - m^2 (phi^a)^2/2)],
    Pi = (D/2-1) I_sym - (eta tensor eta)/2,
    I_sym^mu_nu|rho_sigma = (eta^mu_rho eta^nu_sigma + eta^mu_sigma eta^nu_rho)/2.

Here T is the FREE MATTER stress tensor, not T_full. The distinction must be
preserved in the Weinberg-Witten and constraint arguments.

For incoming equal-mass on-shell momenta, define

    E_mu_nu(p1,p2) = -(p1_mu p2_nu + p2_mu p1_nu)
                     + eta_mu_nu (p1 dot p2 + m^2),
    q = p1+p2,       q^mu E_mu_nu = 0,
    A0 = -Pi/(2 V0),
    A = A0 + K A.

The one-loop kernel is to be computed from two E insertions and two scalar
propagators. In the conserved-source quotient its corrected result is

    K = -lambda [1 - D q^2/(12 m^2)] I_sym + O(q^4),
    lambda = N(D/2-1) Gamma(-D/2) m^D / [2 V0 (4 pi)^(D/2)].

The massless critical value is lambda=-1, NOT the lambda=+1 in the original
2016 paper's uncorrected signs. The later papers explicitly correct this.

This calculation sums the bubble chain at leading 1/N; it is not an all-orders
or finite-N graviton existence theorem. Small q^2/m^2 and a consistent
covariant regulator are required. The leading Laurent controls in D=4-epsilon:

    V0_critical = -N m^4 / (32 pi^2 epsilon),
    M_P^2 = -N m^2 / (48 pi^2 epsilon),
    A_pole = -[I_sym - (eta tensor eta)/2]/[M_P^2 q^2].

Use epsilon < 0 for positive M_P^2 in this formal prescription. Terms O(epsilon)
and finite matching contributions are not declared zero physical corrections.

## 4. Selected hypotheses and failure rules

POSTULATED: manifold/dimension, real scalar inventory and positive target metric,
the determinant action and potential, general reparametrization invariance,
ordinary quantum operations, large N, the formal regulator/measure organization,
flat-background/tadpole conditions and a possible analytic continuation.

TO BE DERIVED: classical metric elimination, constraint identities, determinant
expansion, kernel normalization, critical pole, conserved-source residue,
conditional universal scalar coupling and Newtonian response.

NOT DERIVED BY NAMING A POLE: a real Lorentzian quantum theory, complete first-
class quantum constraint algebra, absence of every ghost, finite-N stability,
full nonlinear Einstein equations, observed matter inventory or new GRUT physics.

Do not modify the action, field-space signature, matter inventory or regulator to
rescue a failure. A conventional EFT with the same low-energy pole and scalar
inventory is the novelty comparator. Full allowed matching freedom must be
included; an arbitrary cutoff-to-Planck relation is not credited as a new
prediction without an independently defined physical cutoff.

## 5. Fixed computation contract

1. Exact four-dimensional trace/stress algebra and quartic expansion using
   rational matrices. Three fixed symmetric kinetic matrices and two potential
   values (0, 2/5), with V0=7.
2. Classical Gram/signature and Hessian controls at N=4,6,8 using fixed rational
   gradient matrices. Analytic sign and rank proofs govern the whole branch.
3. Independently derive the leading dimensional pole of the scalar bubble using
   Feynman parameters, radial integrals and angular moments. Do not copy K into
   a test and call that a loop calculation. Test every symmetric transverse
   source basis vector for q along the fourth Euclidean axis, then continue
   Q_E^2=-q^2.
4. Exact null-momentum conserved-source residue: q=(1,0,0,1), six-dimensional
   conserved source space. Prove positivity and rank two; compare with the
   off-shell timelike tensor without using that as proof of a physical ghost
   at the exact massless point.
5. Critical resummation, sign and tadpole tuning; two-species control
   (N1,m1)=(120,1), (N2,m2)=(80,2). No independent species couplings fitted.
6. Diagnostic constants only: N=200, m=1, epsilon=-1/100. Small Euclidean
   Q^2/m^2=(1/100,1/1000,1/10000). These are software fixtures, not measured
   apparatus parameters. Detuning by +/-1/100 is only a sensitivity diagnostic;
   it is not claimed to be a consistent flat-vacuum alternative of the frozen
   clockless theory.
7. Independent heat-kernel coefficient check and complete conventional-EFT
   comparator. Published cubic-graviton results may be cited but must not be
   represented as independently recomputed here.
8. Following the owner's uploaded note during analytic design, derive the
   reparametrization constraints and their Poisson brackets. Distinguish a
   worldvolume deformation algebra from the ADM dynamics of an independent
   spatial metric. Test the sign and the degrees of freedom. Audit physical
   observables, gravitational dressing and the boundary-Hamiltonian premise of
   Marolf's theorem. Do not assume that a pole or Noether identity supplies
   the required quantum physical observable algebra. This adds a hostile
   test to the same action, not a new candidate.

Freeze this action and contract before executing the scientific script. Code
implementation corrections may be logged; they cannot change the physical law
or select a different surviving mechanism.
