# CP-06 — what a connection-first mechanism actually generates

Author-level result: **classical correspondence established; fundamental
selection objective fails**. This is not a new GRUT law. One candidate was
specified and frozen in LAW.md before executable checks. No replacement is
constructed after the obstruction.

## Metric and dynamics derived

For L=sqrt(|det K|)/V, variation with respect to symmetric Ricci gives

    p^{ab}=dL/dR_ab=(M^2/2V)sqrt(|det K|)(K^-1)^{ab}.

Define the nondegenerate conjugate metric by

    p^{ab}=(M^2/2)sqrt(|g|)g^{ab}.

Taking determinants in four dimensions gives sqrt(|det K|)=V^2 sqrt(|g|).
Substitution then yields the exact algebraic relation

    K_ab=V g_ab,
    M^2 R_ab=partial_a phi partial_b phi+V g_ab.

The metric is derived as an invertible conjugate variable on this branch;
its physical signature is inherited from K, not selected by this definition.

Varying the torsion-free connection and integrating the Ricci variation by
parts yields the metric-density compatibility equation. Contracting its free
indices gives nabla_c p^{ab}=0, hence Gamma is Levi-Civita(g). The scalar
variation gives Box_g phi-V'(phi)=0. These statements require appropriate
variational boundary conditions and a regular branch throughout.

The Legendre Hamiltonian density is

    H=p^{ab}R_ab-L
     =sqrt(|g|)[(1/2)g^{ab}partial_a phi partial_b phi+V].

Therefore the first-order dual action is exactly

    S_dual=integral sqrt(|g|)[M^2 g^{ab}R_ab(Gamma)/2
                            -(partial phi)^2/2-V(phi)].

Eliminating Gamma recovers Einstein-scalar gravity, including its scalar
stress tensor and gravitational self-coupling. Boundary terms must be treated
consistently; no equality of arbitrary quantum path integrals is asserted.

On a Lorentzian branch the regular classical constraint and polarization
content is inherited from this invertible equivalence: two tensor helicities
plus one scalar, with conventional gravitational constraints. This is not an
independent quantization or a global existence theorem. Around a curved vacuum,
the flat-space graviton description is only a local short-curvature-radius
approximation; an exact Minkowski vacuum is not present for V>0.

## Two exact hostile tests

**Time is not selected.** At a constant scalar minimum phi=0 with V(0)=V0>0,
the equation is R_ab=(V0/M^2)g_ab. A Euclidean four-sphere and Lorentzian de
Sitter spacetime both satisfy it. Both give an admissible nondegenerate K.
The action has therefore earned a conjugate metric, but has not earned a
physical distinction between time and space. Restricting its domain to one
signature would be an additional assumption. Degenerate K cannot be used as
a demonstrated signature-transition mechanism: the inverse K and Legendre
map fail there.

**Interactions are not selected.** Consider two permitted positive potentials

    V1=V0+m^2 phi^2/2,
    V2=V0+m^2 phi^2/2+lambda phi^4/24, lambda>0.

They have identical vacuum curvature, scalar mass and linear metric coupling,
but different tree-level four-scalar contact amplitudes. Both obey the same
affine mechanism. Thus removing the metric from the starting variables does
not constrain the independent scalar interaction. There is no exclusion of
the second theory from the stated premises.

The single included scalar is minimally coupled in the dual frame. Universal
coupling to a complete physical matter inventory is NOT derived: spinors,
gauge fields and additional matter would require a complete affine prescription
and corresponding physical checks. No arbitrary matter extension is credited.

## Quantum and empirical boundary

Classical Legendre equivalence does not select a quantum measure, prove
unitarity, protect the action from counterterms or predict Newton's constant.
M, V0 and the potential are supplied. Within its earned correspondence the
model gives conventional Einstein-scalar observables; there is no independent
novel empirical prediction. A different quantum measure might change that
conclusion but is not part of this construction and is not evaluated here.

## Computation and literature

check.py evaluates the original diagonal determinant density numerically and
finite-differences its Ricci derivative, then checks the independently derived
Legendre Hamiltonian in both signatures. It corroborates algebra only. Exact
derivations above cover the general nondegenerate symmetric tensor case.

The underlying affine route is established, not invented by GRUT:
- Kijowski and collaborators, Quantum gravity in the Eddington purely affine
  picture, Phys. Rev. D 29, 2746 (1984), DOI 10.1103/PhysRevD.29.2746.
- Poplawski, Gravitation, electromagnetism and the cosmological constant in
  purely affine gravity, https://arxiv.org/abs/gr-qc/0612193.

The sources establish affine/metric correspondence and relevant matter
extensions, not the success of this packet as fundamental new physics.

## Disposition

This attempt changes the primitive variable and genuinely derives a metric
and Einstein-scalar equations on regular branches. Its exact duality also
exposes why that achievement does not explain nature's choice: the causal
branch and interaction function remain independent supplied information.
Do not rescue this construction by imposing a signature or changing its
potential and calling those choices consequences. CP-06 is closed at this
specific obstruction; the broader research mission remains active. Stage 3
and Card 3 are unchanged. No independent scientific validation is claimed.
