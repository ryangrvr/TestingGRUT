# CP-07: a working finite quantum causal-growth mechanism, with a failed spacetime limit

2026-10-11. Author-level construction and tests; independent review outstanding. One frozen candidate, no parameter search. **Verdict: conditional operational success; not a viable fundamental gravitational theory.** No Stage-3 card opened; no canonical integration or banking.

## Outcome

This construction supplies an actual normalized quantum dynamics instead of assigning a determinant to a response matrix. The same gates generate directed relations, propagate matter-state changes and create reversible physical records. It respects ordinary quantum mixture operations. Its probability measure over complete edge histories is defined, and matter changes that measure.

It nevertheless fails the intended extended-spacetime limit for exact reasons: link-count phases do not change diagonal history probabilities, and its all-pairs growth with a strictly positive attachment probability generates macroscopically linear-height chains. Preferred birth order, microscopic tensor factors, gates and Born probabilities are supplied. No physical clock, energy conservation, continuum field equation, massless spin-2 excitation or universal coupling has been derived.

## Explicit mechanism

The law is in LAW.md, frozen in local commit d6645ac before executable testing. Each pair i<j has a fresh relation qubit e, with matter qubits i and j. Define

    Theta_ij = theta0 + g Z_i Z_j,
    R = exp(-i Theta_ij Y_e),
    P = exp(i phi n_e),
    V = exp(-i gamma n_e Z_i X_j).

Apply V P R in birth order, with theta0=pi/4, g=pi/12, gamma=pi/4, phi=pi/3. All initial registers are |0>. Geometry diagnostics use the transitive closure of primitive directed edges; those edges are not necessarily covering links. Distinct primitive graphs with the same closure are retained as distinct physical states.

For one relation, the two conditional matter operators are

    K0 = cos(Theta),
    K1 = exp(i phi) exp(-i gamma Z_i X_j) sin(Theta).

They satisfy K0^dagger K0 + K1^dagger K1 = I. This makes each observed relation a quantum instrument, and its unobserved matter map completely positive and trace preserving. A product over the specified sequence defines all complete-history amplitudes and probabilities. A finite-N state vector is obtained by a unitary circuit on fixed registers; append-a-register growth is an isometry, not a unitary map between different-dimensional spaces without an extension.

The conditional relation probability is tr[rho sin^2(Theta)]. Earlier V operations change the matter state and thereby change later relation probabilities. Matter is not merely painted onto a fixed order. This is a mathematical toy mechanism, not an experimentally supported microscopic law.

## Exact phases-versus-records result

Every relation register is rotated once. Subsequently it is used as a control and never mixed again. Its phase gate commutes with all subsequent operations. Thus the final state has the form

    |Psi(phi)> = sum_E exp(i phi |E|) |E> |psi_E>.

Consequently

    p(E) = ||psi_E||^2,
    p(C) = sum_{E: closure(E)=C} ||psi_E||^2

are independent of phi. Tracing relation registers also makes the matter marginal independent of phi. Coarse-graining different primitive graphs to one closure does not erase their distinguishing microscopic records and does not justify replacing the sum of probabilities by a coherent scalar sum.

Joint geometry/matter coherences can depend on phi. For the two-event initial control, Theta=theta0+g and the conditional matter overlap is cos(gamma). The relation Ramsey observable obeys

    <X_e> = sin(2 Theta) cos(gamma) cos(phi).

Its values are 0.306186217848 at phi=pi/3 and 0.612372435696 at phi=0. Inverse preparation of the relation gives return probability

    cos^4(Theta)+sin^4(Theta)
    +2 cos^2(Theta) sin^2(Theta) cos(gamma) cos(phi),

which is 0.757582521472 for the declared control. This is genuine measurable quantum interference in the toy circuit, with calibrated gate parameters. It is not a preference for a geometric history. Recombining records changes what is measured; it does not retroactively change the original diagonal distribution.

Copying a relation into an observer ancilla suppresses cross-sector coherence after that ancilla is traced out. Inverting the copy restores the full state exactly. These ordinary unitary operations supply a concrete relation between dynamics, fluctuations and records; they do not select a fundamental outcome or insert consciousness.

## Exact geometry obstruction

The eigenvalues of sin^2(Theta) are 1/4 and 3/4. Therefore, even for arbitrary conditional matter states, each fresh edge probability lies in [1/4,3/4]. Relation measurements can be made immediately after each elementary transition for calculating final diagonal-history probabilities: completed registers are never mixed later. Their nonselective measurement leaves those final probabilities unchanged.

Construct a greedy directed path. Begin at event 0. At each later birth j, inspect the fresh edge from the current path endpoint to j; if it is present, append j, otherwise keep the endpoint. Other preceding gates at that birth may change its conditional matter state, but the bound p>=1/4 still holds. Let Y_j be this success indicator. Then

    H_N >= 1+sum_{j=1}^{N-1} Y_j,
    E[H_N] >= 1+(N-1)/4.

With the conditional history filtration, Y_j-E[Y_j|past] is a bounded martingale difference. Azuma's bound gives

    Pr[ H_N < 1+(N-1)/8 ] <= exp[-(N-1)/128].

Thus height is at least linear in N with probability approaching one. No independence of the edges is assumed.

A fixed-shape, fixed-dimensional 4D Minkowski Alexandrov interval sprinkled at increasing density has longest-chain scale N^(1/4), rather than N. The frozen growth law therefore cannot produce that correspondence regime. This does not exclude every finite embedding or all unrestricted curved universes. It excludes the specific intended finite-density extended continuum comparison.

More generally, a uniform conditional attachment lower bound p_min(N) gives E[H_N]>=1+(N-1)p_min(N). Compatibility with an O(N^(1/4)) mean-height target requires p_min(N)=O(N^(-3/4)); this is only a necessary bound, not a construction of geometry. Choosing such a scaling by hand would be another supplied premise. The missing physics is a restriction on actual available interactions, not just an added phase.

## Executed checks

Run `python check.py` with NumPy. All checks passed. Exact state vectors are compared with separately implemented conditional-branch multiplication; the latter does not use the three-qubit relation gate implementation.

| Control | Result |
|---|---:|
| Elementary unitarity error | 4.03e-16 |
| Maximum independent branch/state-vector discrepancy | 6.40e-16 |
| Maximum inverse-evolution discrepancy | 1.04e-15 |
| Maximum phi-dependent history probability difference | 1.12e-16 |
| Affine mixture discrepancy | 8.72e-17 |
| Reduced-state change after copying a record | 0.612372435696 |
| Full-state discrepancy after inverse copying | 0 |

All 8, 64 and 1024 primitive edge assignments for N=3,4,5 are included through the exact state vector and branch calculation. No Monte Carlo or asymptotic extrapolation supplies the long-chain theorem.

| N | Total qubits | Mean height | Total-variation change vs gamma=0 |
|---|---:|---:|---:|
| 3 | 6 | 2.359375 | 0.234375 |
| 4 | 10 | 2.992513 | 0.320408 |
| 5 | 15 | 3.615140 | 0.440190 |

This variation demonstrates matter-dependent causal-history probabilities. It does not establish a gravitational back-reaction law. Full numerical provenance is in results.json. Small matrix errors establish implementation consistency, not truth in nature.

## What construction now requires

The model supplies O(N^2) relation registers, and therefore lies outside the earlier O(N)-variable Gaussian bound. That evasion is not success: the degrees of freedom themselves carry permanent which-history records, defeating the intended diagonal phase selection. Its normalized Born measure is also a new supplied measure; it cannot be credited as derived from the previous uniform poset ensemble.

This attempt identifies two concrete obligations for future construction: physical restrictions must prevent persistent dense all-pairs attachment, and competing microscopic paths must actually recombine if their phases are to affect a geometry measurement. Neither follows just from saying “quantum growth.” A useful future mechanism must also account for birth-order invariance, physical clocks and energy, then recover relativistic gravitational constraints. No such revised law is formulated here; the frozen one is not repaired after failure.

The broader theory-construction mission remains active. This packet is a substantive checkpoint with a functioning conditional mechanism and exact failures, not a claim that the problem is solved or that theory development should stop.

## Prior art and novelty

- Rideout and Sorkin, A Classical Sequential Growth Dynamics for Causal Sets (2000): https://arxiv.org/abs/gr-qc/9904062 . Growth and covariance constraints are established; this circuit does not earn their label-independent physical interpretation.
- Markopoulou, Quantum causal histories (1999): https://arxiv.org/abs/hep-th/9904009 . Causal structure with quantum evolution is an established research direction.
- Hawkins, Markopoulou and Sahlmann, Evolution in Quantum Causal Histories (2003): https://arxiv.org/abs/hep-th/0302111 . Physical subsystem evolution and completely positive maps are established structures.

The published abstracts and bibliographic records were inspected, not all derivations. The present circuit is a specified control construction in ordinary quantum information theory. Originality of its precise assembly is unestablished; no new law of nature is claimed.
