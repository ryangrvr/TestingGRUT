# CP-09 independent review: two sources of frustration and the all-pairs fluctuation scale

2026-10-11. Frozen CP-09 source: [LAW.md](https://github.com/ryangrvr/TestingGRUT/blob/1dadbd60ae5719508ac40ce7633768549052691a/research/constructive_physics/attempt09/LAW.md) and [REPORT.md](https://github.com/ryangrvr/TestingGRUT/blob/1dadbd60ae5719508ac40ce7633768549052691a/research/constructive_physics/attempt09/REPORT.md). Separate review branch forked from CP-09. **No frozen construction, Card 3 or canonical branch changed.**

## Verified structure

CP-09 postulates, not derives, the matter-valued connection Laplacian on the graph of all labelled N-event posets:
\[
H=J\sum_{e:C\to C'} L_e^\dagger L_e+\Omega\sum_{i=0}^{N-1}(I-X_i),\quad
(L_e\psi)=\psi_{C'}-U_{ij}\psi_C,\quad
U_{ij}=e^{-i\gamma Z_iX_j}.
\]
Assume J>0, Omega>=0 and fixed finite N. As each term is positive, E0>=0. The exact total-order trial state proves E0<=J(N-1), but **does not** establish convergence of E0/N, a gapped or bounded-excitation limit, manifoldlike ground-state support, clock independence, physical locality, or an Einstein limit. The original CP-08 negative quadratic *lower-energy* obstruction is genuinely removed.

## New exact result: two inequivalent forms of frustration

**Theorem: if Omega>0 and sin gamma !=0, then E0(N)>0 for every finite N>=2.** This is true even if the geometry configuration graph has no loops.

Proof: If <psi|H|psi>=0, positivity forces *both* each L_e psi=0 and the onsite term to vanish. The unique zero-energy matter vector of sum_i(I-X_i) is |+>_X^tensor N. Therefore every nonzero geometry block must be psi_C=a_C|+...+>. Covariant constancy along each transition would require
\[
a_{C'}|+\rangle^{\otimes N}=a_C e^{-i\gamma Z_iX_j}|+\rangle^{\otimes N}
=a_C\{\cos\gamma |+\rangle^{\otimes N}-i\sin\gamma\,Z_i|+\rangle^{\otimes N}\}.
\]
The vectors |+>^tensor N and Z_i|+>^tensor N are orthogonal. For sin gamma !=0, each edge forces a_C=a_C'=0. The connected full poset-reconfiguration graph contains edges covering every vertex (N>=2), so psi=0, contradiction. Finite-dimensional H has a strictly positive smallest eigenvalue. QED.

This theorem is **not evidence of non-Abelian history holonomy**: it holds already at N=2, where the three-configuration graph is a tree and has no connection loops. There the transported geometry connection is flat; with Omega=0 its Laplacian has a zero mode for *all* gamma. Adding Omega>0 pins the matter state and conflicts with the transport. On larger geometry graphs, loop holonomy could cause additional geometric frustration even without Omega; that separate question is not decided here.

At gamma=0 the exact uniform geometry times all-X+ state has zero energy, agreeing with the frozen CP-09 adversarial control. For gamma=pi/4, Omega=1/4 and J=1, the independent N=2 matrix reproduces E0=0.0947552836977.

**Interpretation:** an E0>0 result and a nontrivial matter-transport holonomy must never be used interchangeably as a claimed dynamical geometry selector. The former can follow simply from onsite preference conflicting with local transport.

## New exact result: empty-order spectral variance scales quadratically

Let |empty> be the empty poset, prepare any normalized matter state |phi>, and let N fixed. Each of the N(N-1) oriented one-relation moves produces an orthogonal geometry block and a unitary U_ij. Thus the geometry-hop contribution to the energy variance is exactly J^2N(N-1). The onsite Hamiltonian K=Omega sum_i(I-X_i) contributes independently:
\[
\boxed{\quad {\rm Var}_{|\emptyset,\phi\rangle}(H)
=J^2N(N-1)+{\rm Var}_{|\phi\rangle}(K).\quad}
\]
For the frozen initial all-Z+ product |0...0>, Var(K)=Omega² N. Consequently,
\[
\boxed{\langle H\rangle_{\emptyset,0}=JN(N-1)+\Omega N,\qquad
{\rm Var}_{\emptyset,0}(H)=J^2N(N-1)+\Omega^2N.}
\]
Both identities are exact at arbitrary N; they follow from block orthogonality, not a perturbative numerical fit.

Survival at short time is \(1-{\rm Var}(H)t^2+O_N(t^3)\) for the full **joint** initial quantum state, while the probability that the *geometry alone* remains empty is \(1-J^2N(N-1)t^2+O_N(t^3)\). These are different observables. The latter independently reproduces the CP-09 author's warning that adding a configuration-degree diagonal does not repair the all-pairs escape scaling. Even subtracting an N-dependent c-number from H leaves these variances unchanged. This does not prove that every large-N physical limit fails; collectively scaled observables, a derived normalization or preparation may behave differently.

## Independent N=2 regression

[check.py](check.py) builds the 12x12 N=2 Hamiltonian independently from the LAW (does **not** import CP-08 or CP-09 code). It compares ground energies and zero-mode status for gamma=0, Omega>0; gamma=pi/4, Omega=0; gamma=pi/4, Omega=1/4. It checks the exact all-X+ overlap/orthogonal transport and the exact mean and variance of the frozen all-Z+ empty-geometry initial state. The reported N=2 positive ground energy is independently reproduced in this method. Run with NumPy.

The supplied CP-09 author tests pass at N=2,3,4; this review is a more narrowly targeted algebraic and independently assembled two-event control, not a verification of all those tests, and not a new proposed law.

## Scientific decision

The CP-09 positivity principle is an effective repair of the precise CP-08 energy instability, but **the law still contains global all-pairs reconfiguration**. The spectral variance supplies an invariant, nonzero quantitative obstruction independent of diagonal energy reference. The low-energy configurations must be distinguished by class-level geometric predictions and physical propagation, not by E0 positivity or arbitrary geometric labelling.

A defensible next test on this *frozen law* is to calculate the low-energy spectral gap, state-dependent detector cross-response and energy variance with N increasing, choosing preparation and detector supports before observing the geometry. In parallel compare ground-state support of manifoldlike versus layered configurations under the same Hamiltonian and state. A bounded E0/N is necessary but insufficient.

**Status:** independent exact finite-N theorem and spectral moment, with standalone N=2 control. No CP-10, new theory of reality or Stage-3 advancement claimed.
